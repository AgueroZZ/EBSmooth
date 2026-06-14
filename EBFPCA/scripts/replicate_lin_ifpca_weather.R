#!/usr/bin/env Rscript

# Approximate R replication of the Canadian weather analysis from the
# interpretable FPCA supplement by Lin, Wang, and Cao. The original supplement
# uses MATLAB code. This script translates the core greedy L0 basis-selection
# routine to R and overlays the resulting iFPCA components with regularized
# FPCA and PACE.

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2 || !startsWith(arg, "--")) {
      stop("Arguments must use --name=value syntax. Got: ", arg, call. = FALSE)
    }
    out[[pieces[[1]]]] <- pieces[[2]]
  }
  out
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript EBFPCA/scripts/replicate_lin_ifpca_weather.R [options]",
      "",
      "Options:",
      "  --zip=PATH          Supplement zip path.",
      "  --out-dir=PATH      Output directory.",
      "  --nbasis=N          B-spline basis size. Default: 181.",
      "  --nharm=N           Number of components. Default: 3.",
      "  --gamma=X           Roughness penalty. Default: 3000.",
      "  --folds=N           CV folds for kappa selection. Default: 10.",
      "  --seed=N            Random seed for CV folds. Default: 983.",
      "  --coef-lambda=X     Small ridge for stable Data2fd coefficients.",
      "  --kappa-min=N       Minimum active basis size. Default: 1.",
      "  --kappa-max=N       Maximum active basis size. Default: nbasis.",
      "  --kappa-step=N      Kappa grid step. Default: 1.",
      sep = "\n"
    )
  )
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
}

find_repo_root <- function() {
  script_path <- get_script_path()
  candidates <- c(getwd())
  if (!is.na(script_path)) {
    script_dir <- dirname(script_path)
    candidates <- c(candidates, script_dir, file.path(script_dir, ".."),
                    file.path(script_dir, "..", ".."))
  }
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "EBFPCA", "data"))) {
      return(candidate)
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
}

require_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

read_weather_temperature <- function(zipfile) {
  member <- "iFPCAcodes/fdaM/examples/weather/temperature.csv"
  label_member <- "iFPCAcodes/fdaM/examples/weather/daillabs.dat"
  con <- unz(zipfile, member, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- as.matrix(utils::read.csv(con, header = FALSE, check.names = FALSE))
  storage.mode(x) <- "double"

  labels <- paste0("site_", seq_len(ncol(x)))
  label_con <- unz(zipfile, label_member, open = "rt")
  raw_labels <- tryCatch(readLines(label_con, warn = FALSE),
                         error = function(e) character())
  close(label_con)
  raw_labels <- trimws(raw_labels)
  raw_labels <- raw_labels[nzchar(raw_labels)]
  if (length(raw_labels) >= ncol(x)) {
    labels <- raw_labels[seq_len(ncol(x))]
  }
  colnames(x) <- make.unique(labels)
  list(x = x, grid = seq(0.5, 364.5, length.out = nrow(x)), member = member)
}

leading_generalized_eigen <- function(A, B, opts = list()) {
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  n <- nrow(A)
  if (n == 1) {
    return(list(value = A[1, 1] / B[1, 1], vector = matrix(1, ncol = 1)))
  }

  fit <- tryCatch({
    RSpectra::eigs_sym(
      A,
      k = 1,
      B = B,
      which = "LA",
      opts = modifyList(list(tol = 1e-8, maxitr = 1000), opts)
    )
  }, error = function(e) NULL)

  if (!is.null(fit) && length(fit$values) >= 1) {
    return(list(value = fit$values[[1]], vector = as.matrix(fit$vectors[, 1])))
  }

  chol_b <- chol(B)
  chol_b_inv <- solve(chol_b)
  transformed <- t(chol_b_inv) %*% A %*% chol_b_inv
  eig <- eigen((transformed + t(transformed)) / 2, symmetric = TRUE)
  list(
    value = eig$values[[1]],
    vector = as.matrix(chol_b_inv %*% eig$vectors[, 1])
  )
}

top_generalized_eigen <- function(A, B, k) {
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  k <- min(k, nrow(A) - 1)
  fit <- tryCatch({
    RSpectra::eigs_sym(A, k = k, B = B, which = "LA",
                       opts = list(tol = 1e-8, maxitr = 2000))
  }, error = function(e) NULL)
  if (!is.null(fit) && length(fit$values) >= k) {
    ord <- order(fit$values, decreasing = TRUE)
    return(list(values = fit$values[ord], vectors = fit$vectors[, ord,
                                                                drop = FALSE]))
  }

  chol_b <- chol(B)
  chol_b_inv <- solve(chol_b)
  transformed <- t(chol_b_inv) %*% A %*% chol_b_inv
  eig <- eigen((transformed + t(transformed)) / 2, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)[seq_len(k)]
  list(values = eig$values[ord], vectors = chol_b_inv %*% eig$vectors[, ord])
}

myopt_l0 <- function(X, W, D2, gamma, kappa, return_path = FALSE) {
  n_basis <- ncol(X)
  n_curve <- nrow(X)
  A <- W %*% crossprod(X) %*% W
  B <- W + gamma * D2
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  B_full <- B

  Bi <- solve(B)
  BiA <- Bi %*% A
  eig <- leading_generalized_eigen(A, B)
  v <- eig$vector[, 1]
  v <- v / sqrt(sum(v * v))
  lambda <- eig$value

  active <- seq_len(n_basis)
  kk <- n_basis
  n_steps <- n_basis - kappa + 1

  if (return_path) {
    path_x <- matrix(NA_real_, nrow = n_basis, ncol = n_steps)
    path_y <- rep(NA_real_, n_steps)
    active_sizes <- rep(NA_integer_, n_steps)
    removed <- rep(NA_integer_, max(n_steps - 1, 0))
    step <- 1L
  }

  while (kk > kappa) {
    if (return_path) {
      b <- rep(0, n_basis)
      b[active] <- v
      b <- b / sqrt(drop(t(b) %*% B_full %*% b))
      path_x[, step] <- b
      path_y[[step]] <- lambda / n_curve
      active_sizes[[step]] <- kk
    }

    diff <- rep(Inf, kk)
    for (i in seq_len(kk)) {
      denom <- 1 - v[[i]]^2
      if (abs(denom) < 1e-10) {
        next
      }
      y <- Bi[i, i]
      g <- Bi[, i]
      p2 <- BiA[, i]
      cval <- BiA[i, i]
      nt1 <- v[[i]] * (sum(v * g) - v[[i]] * y) / y
      nt2 <- v[[i]] * (sum(v * p2) - v[[i]] * cval)
      diff[[i]] <- (nt1 + nt2) / denom
    }

    remove_idx <- which.min(diff)
    if (return_path) {
      removed[[step]] <- active[[remove_idx]]
      step <- step + 1L
    }

    keep <- rep(TRUE, kk)
    keep[[remove_idx]] <- FALSE
    BiA <- BiA[keep, keep, drop = FALSE] -
      Bi[keep, remove_idx, drop = FALSE] %*%
      BiA[remove_idx, keep, drop = FALSE] / Bi[remove_idx, remove_idx]
    Bi <- Bi[keep, keep, drop = FALSE] -
      Bi[keep, remove_idx, drop = FALSE] %*%
      Bi[remove_idx, keep, drop = FALSE] / Bi[remove_idx, remove_idx]
    A <- A[keep, keep, drop = FALSE]
    B <- B[keep, keep, drop = FALSE]
    active <- active[keep]
    kk <- kk - 1L

    eig <- leading_generalized_eigen(A, B)
    v <- eig$vector[, 1]
    v <- v / sqrt(sum(v * v))
    lambda <- eig$value
  }

  b <- rep(0, n_basis)
  b[active] <- v
  b <- b / sqrt(drop(t(b) %*% B_full %*% b))
  opt_y <- lambda / n_curve

  if (!return_path) {
    return(list(value = opt_y, vector = b))
  }

  path_x[, step] <- b
  path_y[[step]] <- opt_y
  active_sizes[[step]] <- kk
  colnames(path_x) <- paste0("kappa_", active_sizes)
  list(
    value = opt_y,
    vector = b,
    path_x = path_x,
    path_y = path_y,
    active_sizes = active_sizes,
    removed = removed
  )
}

make_folds <- function(n, k, seed) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(k), length.out = n))
  lapply(seq_len(k), function(j) which(fold_id == j))
}

mygcv_l0 <- function(A, W, D2, gamma, kappas, folds) {
  kmin <- min(kappas)
  result <- matrix(NA_real_, nrow = length(kappas), ncol = length(folds))
  rownames(result) <- paste0("kappa_", kappas)

  for (fold_idx in seq_along(folds)) {
    test_idx <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(nrow(A)), test_idx)
    path <- myopt_l0(A[train_idx, , drop = FALSE], W, D2, gamma, kmin,
                     return_path = TRUE)
    for (kp_idx in seq_along(kappas)) {
      kappa <- kappas[[kp_idx]]
      path_col <- match(kappa, path$active_sizes)
      xik <- path$path_x[, path_col]
      cQWc <- drop(t(xik) %*% W %*% xik)
      At <- A[test_idx, , drop = FALSE]
      eta <- At - (At %*% (W %*% xik)) %*% t(xik) / cQWc
      result[kp_idx, fold_idx] <- mean(rowSums((eta %*% W) * eta))
    }
    message("  completed CV fold ", fold_idx, "/", length(folds))
  }
  result
}

choose_kappa <- function(cv_result, kappas) {
  mean_cv <- rowMeans(cv_result)
  min_idx <- which.min(mean_cv)
  se <- stats::sd(cv_result[min_idx, ]) / sqrt(max(ncol(cv_result) - 1, 1))
  cutoff <- mean_cv[[min_idx]] + 3 * se
  eligible <- which(mean_cv < cutoff)
  selected_idx <- min(eligible)
  list(
    kappa = kappas[[selected_idx]],
    min_kappa = kappas[[min_idx]],
    mean_cv = mean_cv,
    cutoff = cutoff,
    se = se
  )
}

align_signs_to_ifpca <- function(ifpca_coef, fpca_coef, W) {
  out <- fpca_coef
  for (j in seq_len(nrow(out))) {
    if (drop(out[j, ] %*% W %*% ifpca_coef[j, ]) < 0) {
      out[j, ] <- -out[j, ]
    }
  }
  out
}

coef_to_grid <- function(coef_mat, basis_mat) {
  basis_mat %*% t(coef_mat)
}

estimate_support <- function(values, grid) {
  rows <- list()
  idx <- 1L
  for (j in seq_len(ncol(values))) {
    threshold <- max(abs(values[, j])) * 1e-4
    active <- abs(values[, j]) > threshold
    runs <- rle(active)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    for (r in seq_along(runs$values)) {
      if (isTRUE(runs$values[[r]])) {
        rows[[idx]] <- data.frame(
          component = paste0("PC", j),
          start_day = grid[starts[[r]]],
          end_day = grid[ends[[r]]],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(component = character(), start_day = numeric(),
                      end_day = numeric()))
  }
  do.call(rbind, rows)
}

fit_pace <- function(x, grid, nharm, target_harmonics) {
  ly <- lapply(seq_len(ncol(x)), function(j) x[, j])
  lt <- lapply(seq_len(ncol(x)), function(j) grid)
  fit <- fdapace::FPCA(
    Ly = ly,
    Lt = lt,
    optns = list(
      dataType = "Dense",
      error = TRUE,
      maxK = nharm,
      methodSelectK = nharm,
      methodXi = "CE"
    )
  )
  k <- min(nharm, ncol(fit$phi))
  phi <- apply(fit$phi[, seq_len(k), drop = FALSE], 2, function(y) {
    stats::approx(fit$workGrid, y, xout = grid, rule = 2)$y
  })
  phi <- as.matrix(phi)
  for (j in seq_len(k)) {
    if (stats::cor(phi[, j], target_harmonics[, j]) < 0) {
      phi[, j] <- -phi[, j]
    }
  }
  list(harmonics = phi, lambda = fit$lambda[seq_len(k)], fit = fit)
}

make_raw_plot <- function(x, grid, out_file) {
  df <- data.frame(
    day = rep(grid, times = ncol(x)),
    curve = rep(colnames(x), each = nrow(x)),
    temperature = as.vector(x),
    stringsAsFactors = FALSE
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(day, temperature, group = curve)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.45, color = "#2563EB") +
    ggplot2::labs(title = "Raw Canadian weather temperature curves",
                  x = "Day", y = "Temperature (C)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(out_file, p, width = 7, height = 5, dpi = 160)
}

make_component_plot <- function(grid, ifpca, fpca, pace, variance_table,
                                out_file) {
  rows <- list()
  idx <- 1L
  add_rows <- function(mat, method) {
    local_rows <- list()
    for (j in seq_len(ncol(mat))) {
      local_rows[[j]] <- data.frame(
        day = grid,
        value = mat[, j],
        component = paste0("FPC #", j),
        method = method,
        stringsAsFactors = FALSE
      )
    }
    local_rows
  }
  for (row in add_rows(ifpca, "iFPCA")) rows[[idx <- idx + 1L]] <- row
  for (row in add_rows(fpca, "Regularized FPCA")) rows[[idx <- idx + 1L]] <- row
  if (!is.null(pace)) {
    for (row in add_rows(pace, "PACE")) rows[[idx <- idx + 1L]] <- row
  }
  df <- do.call(rbind, rows)
  labels <- setNames(
    paste0(variance_table$component, "\n",
           "iFPCA ", round(100 * variance_table$ifpca, 1), "%; ",
           "FPCA ", round(100 * variance_table$regularized_fpca, 1), "%; ",
           "PACE ", round(100 * variance_table$pace, 1), "%"),
    paste0("FPC #", seq_len(nrow(variance_table)))
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(day, value, color = method,
                                        linetype = method)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        color = "#9CA3AF") +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::facet_wrap(~ component, nrow = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_color_manual(values = c(
      "iFPCA" = "#2563EB",
      "Regularized FPCA" = "#DC2626",
      "PACE" = "#059669"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "iFPCA" = "solid",
      "Regularized FPCA" = "dashed",
      "PACE" = "dotdash"
    )) +
    ggplot2::labs(title = "Weather eigenfunctions: iFPCA replication",
                  x = "Day", y = "Eigenfunction value") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_file, p, width = 11, height = 4.8, dpi = 160)
}

make_paper_style_component_plot <- function(grid, ifpca, fpca, variance_table,
                                            out_file) {
  rows <- list()
  idx <- 1L
  for (j in seq_len(ncol(ifpca))) {
    rows[[idx]] <- data.frame(
      day = grid,
      value = ifpca[, j],
      component = paste0("FPC #", j),
      method = "iFPCA",
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      day = grid,
      value = fpca[, j],
      component = paste0("FPC #", j),
      method = "FPCA",
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  df <- do.call(rbind, rows)
  labels <- setNames(
    paste0("FPC #", seq_len(ncol(ifpca)), "\n",
           "iFPCA ", round(100 * variance_table$ifpca, 1), "%; ",
           "FPCA ", round(100 * variance_table$regularized_fpca, 1), "%"),
    paste0("FPC #", seq_len(ncol(ifpca)))
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(day, value, color = method,
                                        linetype = method)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        color = "#9CA3AF") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ component, nrow = 1,
                        labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_color_manual(values = c("iFPCA" = "#2563EB",
                                           "FPCA" = "#DC2626")) +
    ggplot2::scale_linetype_manual(values = c("iFPCA" = "solid",
                                              "FPCA" = "dashed")) +
    ggplot2::labs(title = "Weather estimated PCs, paper-style comparison",
                  x = "Day", y = "Temperature") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_file, p, width = 11, height = 4.6, dpi = 160)
}

make_cumulative_plot <- function(variance_table, out_file) {
  methods <- c("ifpca", "regularized_fpca", "pace")
  rows <- list()
  idx <- 1L
  for (method in methods) {
    vals <- variance_table[[method]]
    for (j in seq_along(vals)) {
      rows[[idx]] <- data.frame(
        method = method,
        component_count = j,
        cumulative = sum(vals[seq_len(j)], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  df <- do.call(rbind, rows)
  df$method <- factor(df$method, levels = methods,
                      labels = c("iFPCA", "Regularized FPCA", "PACE"))
  p <- ggplot2::ggplot(df, ggplot2::aes(component_count, cumulative,
                                        color = method, linetype = method)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = seq_len(max(df$component_count))) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                limits = c(0, min(1.05, max(df$cumulative) * 1.1))) +
    ggplot2::labs(title = "Cumulative variance percentage",
                  x = "Number of principal components",
                  y = "Cumulative variance") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_file, p, width = 7, height = 5, dpi = 160)
}

write_summary <- function(out_file, args, selected, support, variance_table) {
  lines <- c(
    "# Lin iFPCA Weather Replication",
    "",
    "## Notes",
    "",
    "This is an R translation of the MATLAB supplement workflow. It uses the",
    "same weather data, greedy L0 support-selection algorithm, CV rule, and",
    "roughness penalty idea, but the R `fda` basis setup is not bit-for-bit",
    "identical to the original MATLAB FDA toolbox.",
    "",
    "## Settings",
    "",
    paste0("- B-spline nbasis: ", args$nbasis),
    paste0("- Number of components: ", args$nharm),
    paste0("- Gamma: ", args$gamma),
    paste0("- CV folds: ", args$folds),
    paste0("- Kappa grid: ", min(args$kappas), " to ", max(args$kappas),
           " by ", args$kappa_step),
    "",
    "## Selected Kappas",
    "",
    paste0("- PC", seq_along(selected$kappa), ": selected kappa ",
           selected$kappa, " (minimum-CV kappa ", selected$min_kappa, ")"),
    "",
    "## Variance Proportions",
    "",
    paste(capture.output(print(variance_table, row.names = FALSE)), collapse = "\n"),
    "",
    "## iFPCA Active Support Intervals",
    "",
    if (nrow(support) == 0) {
      "No active support intervals detected under the grid threshold."
    } else {
      paste(capture.output(print(support, row.names = FALSE)), collapse = "\n")
    },
    "",
    "## Output Files",
    "",
    "- `lin_ifpca_weather_raw_curves.png`",
    "- `lin_ifpca_weather_eigenfunctions.png`",
    "- `lin_ifpca_weather_eigenfunctions_paper_style.png`",
    "- `lin_ifpca_weather_cumulative_variance.png`",
    "- `lin_ifpca_weather_metrics.csv`",
    "- `lin_ifpca_weather_support.csv`",
    "- `lin_ifpca_weather_fits.rds`"
  )
  writeLines(lines, out_file)
}

main <- function() {
  require_installed("fda")
  require_installed("fdapace")
  require_installed("ggplot2")
  require_installed("RSpectra")
  require_installed("scales")

  repo_root <- find_repo_root()
  cli <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) {
    usage()
    return(invisible(NULL))
  }

  args <- list(
    zip = file.path(repo_root, "EBFPCA", "data",
                    "biom12457-sup-0002-suppdatacode.zip"),
    out_dir = file.path(repo_root, "EBFPCA", "results",
                        "lin_ifpca_weather_replication"),
    nbasis = 181L,
    nharm = 3L,
    gamma = 3000,
    folds = 10L,
    seed = 983L,
    coef_lambda = 1e-8,
    kappa_min = 1L,
    kappa_max = NA_integer_,
    kappa_step = 1L
  )
  if (!is.null(cli$zip)) args$zip <- cli$zip
  if (!is.null(cli[["out-dir"]])) args$out_dir <- cli[["out-dir"]]
  if (!is.null(cli$nbasis)) args$nbasis <- as.integer(cli$nbasis)
  if (!is.null(cli$nharm)) args$nharm <- as.integer(cli$nharm)
  if (!is.null(cli$gamma)) args$gamma <- as.numeric(cli$gamma)
  if (!is.null(cli$folds)) args$folds <- as.integer(cli$folds)
  if (!is.null(cli$seed)) args$seed <- as.integer(cli$seed)
  if (!is.null(cli[["coef-lambda"]])) {
    args$coef_lambda <- as.numeric(cli[["coef-lambda"]])
  }
  if (!is.null(cli[["kappa-min"]])) args$kappa_min <- as.integer(cli[["kappa-min"]])
  if (!is.null(cli[["kappa-max"]])) args$kappa_max <- as.integer(cli[["kappa-max"]])
  if (!is.null(cli[["kappa-step"]])) args$kappa_step <- as.integer(cli[["kappa-step"]])
  if (is.na(args$kappa_max)) args$kappa_max <- args$nbasis
  args$kappas <- seq(args$kappa_min, args$kappa_max, by = args$kappa_step)

  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  weather <- read_weather_temperature(args$zip)
  x <- weather$x
  grid <- weather$grid
  centered <- sweep(x, 1, rowMeans(x), "-")

  basis <- fda::create.bspline.basis(
    rangeval = c(min(grid), max(grid)),
    nbasis = args$nbasis,
    norder = 4
  )
  fd <- fda::Data2fd(
    argvals = grid,
    y = centered,
    basisobj = basis,
    lambda = args$coef_lambda
  )
  Acoef <- t(stats::coef(fd))
  W <- fda::eval.penalty(basis, fda::int2Lfd(0))
  D2 <- fda::eval.penalty(basis, fda::int2Lfd(2))
  basis_mat <- fda::eval.basis(grid, basis)

  n_curve <- nrow(Acoef)
  n_total <- max(6L, args$nharm)
  qn <- W %*% crossprod(Acoef) %*% W
  qd <- W + args$gamma * D2
  top_fit <- top_generalized_eigen(qn, qd, n_total)
  total_var <- sum(top_fit$values) / n_curve

  regular_coef <- matrix(NA_real_, nrow = args$nharm, ncol = args$nbasis)
  for (j in seq_len(args$nharm)) {
    coef_j <- top_fit$vectors[, j]
    coef_j <- coef_j / sqrt(drop(t(coef_j) %*% qd %*% coef_j))
    regular_coef[j, ] <- coef_j
  }
  regular_values <- top_fit$values[seq_len(args$nharm)] / n_curve

  folds <- make_folds(n_curve, args$folds, args$seed)
  ifpca_coef <- matrix(NA_real_, nrow = args$nharm, ncol = args$nbasis)
  ifpca_values <- rep(NA_real_, args$nharm)
  selected <- data.frame(
    component = paste0("PC", seq_len(args$nharm)),
    kappa = NA_integer_,
    min_kappa = NA_integer_,
    stringsAsFactors = FALSE
  )
  cv_store <- vector("list", args$nharm)
  A_work <- Acoef

  for (pc in seq_len(args$nharm)) {
    message("Selecting kappa for iFPCA component ", pc, "/", args$nharm)
    if (pc > 1) {
      cc <- ifpca_coef[pc - 1L, ]
      denom <- drop(t(cc) %*% W %*% cc)
      A_work <- A_work - (A_work %*% W %*% cc) %*% t(cc) / denom
    }
    cv <- mygcv_l0(A_work, W, D2, args$gamma, args$kappas, folds)
    choice <- choose_kappa(cv, args$kappas)
    fit_pc <- myopt_l0(A_work, W, D2, args$gamma, choice$kappa,
                       return_path = FALSE)
    ifpca_coef[pc, ] <- fit_pc$vector
    ifpca_values[[pc]] <- fit_pc$value
    selected$kappa[[pc]] <- choice$kappa
    selected$min_kappa[[pc]] <- choice$min_kappa
    cv_store[[pc]] <- list(cv = cv, choice = choice)
    message("  selected kappa=", choice$kappa,
            "; minimum-CV kappa=", choice$min_kappa)
  }

  ifpca_grid <- coef_to_grid(ifpca_coef, basis_mat)
  for (j in seq_len(ncol(ifpca_grid))) {
    if (abs(max(ifpca_grid[, j])) < abs(min(ifpca_grid[, j]))) {
      ifpca_grid[, j] <- -ifpca_grid[, j]
      ifpca_coef[j, ] <- -ifpca_coef[j, ]
    }
  }

  regular_coef <- align_signs_to_ifpca(ifpca_coef, regular_coef, W)
  regular_grid <- coef_to_grid(regular_coef, basis_mat)

  pace <- fit_pace(x, grid, args$nharm, ifpca_grid)

  pace_var <- rep(NA_real_, args$nharm)
  pace_var[seq_along(pace$lambda)] <- pace$lambda / sum(pace$fit$lambda)
  variance_table <- data.frame(
    component = paste0("PC", seq_len(args$nharm)),
    ifpca = ifpca_values / total_var,
    regularized_fpca = regular_values / total_var,
    pace = pace_var,
    selected_kappa = selected$kappa,
    min_cv_kappa = selected$min_kappa,
    stringsAsFactors = FALSE
  )
  support <- estimate_support(ifpca_grid, grid)

  utils::write.csv(variance_table,
                   file.path(args$out_dir, "lin_ifpca_weather_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(support,
                   file.path(args$out_dir, "lin_ifpca_weather_support.csv"),
                   row.names = FALSE)
  saveRDS(
    list(
      args = args,
      weather = weather,
      centered = centered,
      basis = basis,
      coefficients = list(ifpca = ifpca_coef, regularized_fpca = regular_coef),
      harmonics = list(ifpca = ifpca_grid,
                       regularized_fpca = regular_grid,
                       pace = pace$harmonics),
      values = list(ifpca = ifpca_values,
                    regularized_fpca = regular_values,
                    pace = pace$lambda),
      selected = selected,
      cv = cv_store,
      support = support,
      variance_table = variance_table,
      pace_fit = pace$fit
    ),
    file.path(args$out_dir, "lin_ifpca_weather_fits.rds")
  )

  make_raw_plot(x, grid,
                file.path(args$out_dir, "lin_ifpca_weather_raw_curves.png"))
  make_component_plot(
    grid,
    ifpca_grid,
    regular_grid,
    pace$harmonics,
    variance_table,
    file.path(args$out_dir, "lin_ifpca_weather_eigenfunctions.png")
  )
  make_paper_style_component_plot(
    grid,
    ifpca_grid,
    regular_grid,
    variance_table,
    file.path(args$out_dir, "lin_ifpca_weather_eigenfunctions_paper_style.png")
  )
  make_cumulative_plot(
    variance_table,
    file.path(args$out_dir, "lin_ifpca_weather_cumulative_variance.png")
  )
  write_summary(file.path(args$out_dir, "run_summary.md"),
                args, selected, support, variance_table)
  utils::capture.output(utils::sessionInfo(),
                        file = file.path(args$out_dir, "session_info.txt"))

  print(variance_table)
  cat("\nWrote Lin iFPCA weather replication outputs to:\n",
      args$out_dir, "\n", sep = "")
}

main()
