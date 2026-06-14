#!/usr/bin/env Rscript

# Explore several FPCA baselines on the Canadian weather data bundled in the
# Lin et al. iFPCA supplement. The script reads directly from the supplement
# zip file and writes metrics, plots, fitted objects, and session metadata.

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!startsWith(arg, "--")) {
      stop("Arguments must use --name=value syntax. Got: ", arg, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2) {
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
      "  Rscript EBFPCA/scripts/explore_fpca_weather.R [options]",
      "",
      "Options:",
      "  --zip=PATH              Supplement zip path.",
      "  --out-dir=PATH          Output directory.",
      "  --nharm=N               Number of harmonics/components to report.",
      "  --bspline-nbasis=N      Number of B-spline basis functions.",
      "  --fourier-nbasis=N      Number of Fourier basis functions.",
      "  --lambda=X              Smoothing lambda for fda::Data2fd.",
      "",
      "Optional methods:",
      "  fdapace::FPCA and refund::fpca.sc are run only if installed.",
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

has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

open_zip_text <- function(zipfile, member) {
  con <- unz(zipfile, member, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

read_weather_data <- function(zipfile) {
  if (!file.exists(zipfile)) {
    stop("Supplement zip was not found: ", zipfile, call. = FALSE)
  }

  zip_members <- unzip(zipfile, list = TRUE)$Name
  csv_member <- "iFPCAcodes/fdaM/examples/weather/temperature.csv"
  dat_member <- "iFPCAcodes/fdaM/examples/weather/dailtemp.dat"
  label_member <- "iFPCAcodes/fdaM/examples/weather/daillabs.dat"

  if (csv_member %in% zip_members) {
    con <- unz(zipfile, csv_member, open = "rt")
    on.exit(close(con), add = TRUE)
    x <- as.matrix(utils::read.csv(con, header = FALSE, check.names = FALSE))
  } else if (dat_member %in% zip_members) {
    con <- unz(zipfile, dat_member, open = "rt")
    on.exit(close(con), add = TRUE)
    x <- matrix(scan(con, quiet = TRUE), nrow = 365, ncol = 35)
  } else {
    stop("Could not find weather temperature data in supplement zip.",
         call. = FALSE)
  }

  labels <- paste0("site_", seq_len(ncol(x)))
  if (label_member %in% zip_members) {
    raw_labels <- trimws(open_zip_text(zipfile, label_member))
    raw_labels <- raw_labels[nzchar(raw_labels)]
    if (length(raw_labels) >= ncol(x)) {
      labels <- raw_labels[seq_len(ncol(x))]
    }
  }

  storage.mode(x) <- "double"
  colnames(x) <- make.unique(labels)
  rownames(x) <- paste0("day_", seq_len(nrow(x)))

  list(
    x = x,
    grid = seq(0.5, 364.5, length.out = nrow(x)),
    curve_names = colnames(x),
    source = if (csv_member %in% zip_members) csv_member else dat_member
  )
}

apply_sign_convention <- function(harmonics, scores = NULL) {
  harmonics <- as.matrix(harmonics)
  signs <- rep(1, ncol(harmonics))
  for (j in seq_len(ncol(harmonics))) {
    if (sum(harmonics[, j], na.rm = TRUE) < 0) {
      signs[[j]] <- -1
    }
  }
  harmonics <- sweep(harmonics, 2, signs, "*")
  if (!is.null(scores)) {
    scores <- sweep(as.matrix(scores), 2, signs, "*")
  }
  list(harmonics = harmonics, scores = scores)
}

fit_raw_grid_pca <- function(x, grid, nharm) {
  mean_grid <- rowMeans(x)
  centered <- sweep(x, 1, mean_grid, "-")
  pc <- stats::prcomp(t(centered), center = FALSE, scale. = FALSE)
  k <- min(nharm, ncol(pc$rotation), ncol(pc$x))
  scores <- pc$x[, seq_len(k), drop = FALSE]
  signed <- apply_sign_convention(pc$rotation[, seq_len(k), drop = FALSE],
                                  scores)
  harmonics <- signed$harmonics
  scores <- signed$scores
  recon <- t(scores %*% t(harmonics))
  recon <- sweep(recon, 1, mean_grid, "+")
  values <- pc$sdev^2
  list(
    method = "raw_grid_prcomp",
    status = "ok",
    grid = grid,
    mean = mean_grid,
    harmonics = harmonics,
    scores = scores,
    values = values[seq_len(k)],
    varprop = values[seq_len(k)] / sum(values),
    recon = recon,
    message = "Base R prcomp on centered daily grid values."
  )
}

fit_fda_pca <- function(x, grid, nharm, basis_type, nbasis, lambda) {
  rangeval <- c(0, 365)
  if (identical(basis_type, "bspline")) {
    basis <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis,
                                       norder = 4)
    method <- "fda_pca_fd_bspline"
  } else if (identical(basis_type, "fourier")) {
    if (nbasis %% 2 == 0) {
      nbasis <- nbasis + 1
    }
    basis <- fda::create.fourier.basis(rangeval = rangeval, nbasis = nbasis,
                                       period = 365)
    method <- "fda_pca_fd_fourier"
  } else {
    stop("Unknown basis type: ", basis_type, call. = FALSE)
  }

  fdobj <- fda::Data2fd(argvals = grid, y = x, basisobj = basis,
                        lambda = lambda)
  fit <- fda::pca.fd(fdobj, nharm = nharm, centerfns = TRUE)
  harmonics <- as.matrix(fda::eval.fd(grid, fit$harmonics))
  mean_grid <- as.numeric(fda::eval.fd(grid, fit$meanfd))
  scores <- as.matrix(fit$scores[, seq_len(nharm), drop = FALSE])
  signed <- apply_sign_convention(harmonics[, seq_len(nharm), drop = FALSE],
                                  scores)
  harmonics <- signed$harmonics
  scores <- signed$scores
  recon <- harmonics %*% t(scores)
  recon <- sweep(recon, 1, mean_grid, "+")

  list(
    method = method,
    status = "ok",
    grid = grid,
    mean = mean_grid,
    harmonics = harmonics,
    scores = scores,
    values = fit$values[seq_len(nharm)],
    varprop = fit$varprop[seq_len(nharm)],
    recon = recon,
    fdobj = fdobj,
    fda_fit = fit,
    message = paste0("fda::pca.fd with ", basis_type,
                     " basis; nbasis=", nbasis, ", lambda=", lambda, ".")
  )
}

fit_fdapace_optional <- function(x, grid, nharm) {
  method <- "fdapace_FPCA"
  if (!has_pkg("fdapace")) {
    return(list(
      method = method,
      status = "skipped",
      message = "Package fdapace is not installed."
    ))
  }

  out <- tryCatch({
    ly <- lapply(seq_len(ncol(x)), function(j) x[, j])
    lt <- lapply(seq_len(ncol(x)), function(j) grid)
    fit <- fdapace::FPCA(
      Ly = ly,
      Lt = lt,
      optns = list(
        dataType = "Dense",
        error = TRUE,
        maxK = nharm,
        FVEthreshold = 0.99,
        methodXi = "CE"
      )
    )
    work_grid <- fit$workGrid
    k <- min(nharm, ncol(fit$phi), ncol(fit$xiEst))
    harmonics <- apply(fit$phi[, seq_len(k), drop = FALSE], 2, function(y) {
      stats::approx(work_grid, y, xout = grid, rule = 2)$y
    })
    mean_grid <- stats::approx(work_grid, fit$mu, xout = grid, rule = 2)$y
    scores <- fit$xiEst[, seq_len(k), drop = FALSE]
    signed <- apply_sign_convention(harmonics, scores)
    harmonics <- signed$harmonics
    scores <- signed$scores
    recon <- harmonics %*% t(scores)
    recon <- sweep(recon, 1, mean_grid, "+")
    lambda <- fit$lambda[seq_len(k)]
    list(
      method = method,
      status = "ok",
      grid = grid,
      mean = mean_grid,
      harmonics = harmonics,
      scores = scores,
      values = lambda,
      varprop = lambda / sum(fit$lambda),
      recon = recon,
      fdapace_fit = fit,
      message = "fdapace::FPCA on dense list-form trajectories."
    )
  }, error = function(e) {
    list(method = method, status = "error", message = conditionMessage(e))
  })
  out
}

fit_refund_optional <- function(x, grid, nharm) {
  method <- "refund_fpca_sc"
  if (!has_pkg("refund")) {
    return(list(
      method = method,
      status = "skipped",
      message = "Package refund is not installed."
    ))
  }

  out <- tryCatch({
    fit <- refund::fpca.sc(Y = t(x), argvals = grid, npc = nharm, pve = 0.99)
    k <- min(nharm, ncol(fit$efunctions), ncol(fit$scores))
    mean_grid <- as.numeric(fit$mu)
    scores <- as.matrix(fit$scores[, seq_len(k), drop = FALSE])
    signed <- apply_sign_convention(fit$efunctions[, seq_len(k), drop = FALSE],
                                    scores)
    harmonics <- signed$harmonics
    scores <- signed$scores
    recon <- harmonics %*% t(scores)
    recon <- sweep(recon, 1, mean_grid, "+")
    values <- fit$evalues[seq_len(k)]
    list(
      method = method,
      status = "ok",
      grid = grid,
      mean = mean_grid,
      harmonics = harmonics,
      scores = scores,
      values = values,
      varprop = values / sum(fit$evalues),
      recon = recon,
      refund_fit = fit,
      message = "refund::fpca.sc on dense matrix trajectories."
    )
  }, error = function(e) {
    list(method = method, status = "error", message = conditionMessage(e))
  })
  out
}

component_correlation <- function(harmonics, reference, k) {
  if (is.null(harmonics) || is.null(reference)) {
    return(NA_real_)
  }
  if (ncol(harmonics) < k || ncol(reference) < k) {
    return(NA_real_)
  }
  abs(stats::cor(harmonics[, k], reference[, k], use = "pairwise.complete.obs"))
}

compile_metrics <- function(fits, x, nharm) {
  reference <- fits$raw_grid_prcomp$harmonics
  rows <- lapply(fits, function(fit) {
    varprops <- rep(NA_real_, nharm)
    align <- rep(NA_real_, nharm)
    if (identical(fit$status, "ok")) {
      kk <- min(nharm, length(fit$varprop))
      varprops[seq_len(kk)] <- fit$varprop[seq_len(kk)]
      align[seq_len(kk)] <- vapply(seq_len(kk), function(k) {
        component_correlation(fit$harmonics, reference, k)
      }, numeric(1))
    }
    mse <- NA_real_
    if (identical(fit$status, "ok") && !is.null(fit$recon)) {
      mse <- mean((x - fit$recon)^2, na.rm = TRUE)
    }
    data.frame(
      method = fit$method,
      status = fit$status,
      reconstruction_mse = mse,
      varprop_1 = varprops[1],
      varprop_2 = varprops[2],
      varprop_3 = varprops[3],
      abs_cor_with_raw_pc1 = align[1],
      abs_cor_with_raw_pc2 = align[2],
      abs_cor_with_raw_pc3 = align[3],
      message = fit$message,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_curve_plot <- function(x, grid, out_file) {
  curve_id <- rep(colnames(x), each = nrow(x))
  df <- data.frame(
    day = rep(grid, times = ncol(x)),
    curve = curve_id,
    temperature = as.vector(x),
    stringsAsFactors = FALSE
  )
  mean_df <- data.frame(day = grid, temperature = rowMeans(x))
  p <- ggplot2::ggplot(df, ggplot2::aes(day, temperature, group = curve)) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.35, color = "#2B6CB0") +
    ggplot2::geom_line(data = mean_df, ggplot2::aes(day, temperature),
                       inherit.aes = FALSE, linewidth = 1.1,
                       color = "#111827") +
    ggplot2::labs(
      title = "Canadian weather temperature curves",
      x = "Day of year",
      y = "Temperature (C)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(out_file, p, width = 8, height = 5, dpi = 160)
}

make_varprop_plot <- function(fits, nharm, out_file) {
  rows <- list()
  i <- 1
  for (fit in fits) {
    if (!identical(fit$status, "ok")) {
      next
    }
    k <- min(nharm, length(fit$varprop))
    for (j in seq_len(k)) {
      rows[[i]] <- data.frame(
        method = fit$method,
        component = paste0("PC", j),
        varprop = fit$varprop[[j]],
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  df <- do.call(rbind, rows)
  p <- ggplot2::ggplot(df, ggplot2::aes(component, varprop, fill = method)) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      title = "Variance explained by FPCA method",
      x = NULL,
      y = "Variance proportion"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_file, p, width = 8, height = 5, dpi = 160)
}

make_harmonics_plot <- function(fits, nharm, out_file) {
  rows <- list()
  i <- 1
  for (fit in fits) {
    if (!identical(fit$status, "ok")) {
      next
    }
    k <- min(nharm, ncol(fit$harmonics))
    for (j in seq_len(k)) {
      rows[[i]] <- data.frame(
        method = fit$method,
        day = fit$grid,
        component = paste0("PC", j),
        value = fit$harmonics[, j],
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  df <- do.call(rbind, rows)
  p <- ggplot2::ggplot(df, ggplot2::aes(day, value, color = method)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        color = "#9CA3AF") +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::facet_wrap(~ component, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "Estimated eigenfunctions",
      x = "Day of year",
      y = "Eigenfunction value"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_file, p, width = 8, height = 7, dpi = 160)
}

make_scores_plot <- function(fits, out_file) {
  rows <- list()
  i <- 1
  for (fit in fits) {
    if (!identical(fit$status, "ok") || ncol(fit$scores) < 2) {
      next
    }
    curve_names <- rownames(fit$scores)
    if (is.null(curve_names)) {
      curve_names <- paste0("curve_", seq_len(nrow(fit$scores)))
    }
    rows[[i]] <- data.frame(
      method = fit$method,
      curve = curve_names,
      PC1 = fit$scores[, 1],
      PC2 = fit$scores[, 2],
      stringsAsFactors = FALSE
    )
    i <- i + 1
  }
  if (length(rows) == 0) {
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25,
                        color = "#9CA3AF") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25,
                        color = "#9CA3AF") +
    ggplot2::geom_point(color = "#C2410C", size = 1.8, alpha = 0.85) +
    ggplot2::facet_wrap(~ method, scales = "free") +
    ggplot2::labs(
      title = "FPCA score scatter",
      x = "PC1 score",
      y = "PC2 score"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(out_file, p, width = 8, height = 6, dpi = 160)
}

write_run_summary <- function(out_file, data_info, metrics, args) {
  method_lines <- apply(metrics, 1, function(row) {
    paste0("- `", row[["method"]], "`: ", row[["status"]], "; ",
           row[["message"]])
  })
  top_lines <- c(
    "# FPCA Weather Exploration",
    "",
    "## Input",
    "",
    paste0("- Supplement zip: `", args$zip, "`"),
    paste0("- Data member: `", data_info$source, "`"),
    paste0("- Grid points: ", nrow(data_info$x)),
    paste0("- Curves: ", ncol(data_info$x)),
    "",
    "## Methods",
    "",
    method_lines,
    "",
    "## Output Files",
    "",
    "- `metrics.csv`: summary metrics for fitted and skipped methods.",
    "- `fits.rds`: fitted objects and reconstructed curves.",
    "- `weather_curves.png`: raw weather curves plus mean curve.",
    "- `variance_explained.png`: variance proportions by method.",
    "- `harmonics.png`: estimated eigenfunctions.",
    "- `scores_pc1_pc2.png`: PC1/PC2 score scatter, where available.",
    "- `session_info.txt`: R session metadata."
  )
  writeLines(top_lines, out_file)
}

main <- function() {
  require_installed("fda")
  require_installed("ggplot2")
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
                        "fpca_weather_exploration"),
    nharm = 3L,
    bspline_nbasis = 35L,
    fourier_nbasis = 35L,
    lambda = 1e-2
  )
  if (!is.null(cli$zip)) args$zip <- cli$zip
  if (!is.null(cli[["out-dir"]])) args$out_dir <- cli[["out-dir"]]
  if (!is.null(cli$nharm)) args$nharm <- as.integer(cli$nharm)
  if (!is.null(cli[["bspline-nbasis"]])) {
    args$bspline_nbasis <- as.integer(cli[["bspline-nbasis"]])
  }
  if (!is.null(cli[["fourier-nbasis"]])) {
    args$fourier_nbasis <- as.integer(cli[["fourier-nbasis"]])
  }
  if (!is.null(cli$lambda)) args$lambda <- as.numeric(cli$lambda)

  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  data_info <- read_weather_data(args$zip)
  x <- data_info$x
  grid <- data_info$grid

  fits <- list()
  fits$raw_grid_prcomp <- fit_raw_grid_pca(x, grid, args$nharm)
  fits$fda_pca_fd_bspline <- fit_fda_pca(
    x, grid, args$nharm, "bspline", args$bspline_nbasis, args$lambda
  )
  fits$fda_pca_fd_fourier <- fit_fda_pca(
    x, grid, args$nharm, "fourier", args$fourier_nbasis, args$lambda
  )
  fits$fdapace_FPCA <- fit_fdapace_optional(x, grid, args$nharm)
  fits$refund_fpca_sc <- fit_refund_optional(x, grid, args$nharm)

  metrics <- compile_metrics(fits, x, args$nharm)

  utils::write.csv(metrics, file.path(args$out_dir, "metrics.csv"),
                   row.names = FALSE)
  saveRDS(list(args = args, data = data_info, fits = fits, metrics = metrics),
          file.path(args$out_dir, "fits.rds"))

  make_curve_plot(x, grid, file.path(args$out_dir, "weather_curves.png"))
  make_varprop_plot(fits, args$nharm,
                    file.path(args$out_dir, "variance_explained.png"))
  make_harmonics_plot(fits, args$nharm,
                      file.path(args$out_dir, "harmonics.png"))
  make_scores_plot(fits, file.path(args$out_dir, "scores_pc1_pc2.png"))

  session_file <- file.path(args$out_dir, "session_info.txt")
  utils::capture.output(utils::sessionInfo(), file = session_file)

  write_run_summary(file.path(args$out_dir, "run_summary.md"),
                    data_info, metrics, args)

  print(metrics)
  cat("\nWrote FPCA exploration outputs to:\n", args$out_dir, "\n", sep = "")
}

main()
