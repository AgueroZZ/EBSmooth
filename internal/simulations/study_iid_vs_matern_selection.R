args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (!length(file_arg)) {
  stop("Could not determine the script path from commandArgs().")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."))
results_dir <- file.path(script_dir, "results")
fig_dir <- file.path(results_dir, "figures")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

load_ebsmoothr <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(
      file.path(project_root, "EBSmoothr"),
      quiet = TRUE,
      export_all = TRUE
    )
    return(invisible(TRUE))
  }
  stop("This internal selection study requires pkgload.")
}

load_ebsmoothr()

read_key_value_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)

  for (arg in args) {
    arg_clean <- sub("^--", "", arg)
    if (!grepl("=", arg_clean, fixed = TRUE)) {
      out[[arg_clean]] <- TRUE
      next
    }
    pieces <- strsplit(arg_clean, "=", fixed = TRUE)[[1]]
    key <- pieces[1]
    value <- paste(pieces[-1], collapse = "=")
    out[[key]] <- value
  }

  out
}

parse_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))

  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Could not parse logical flag: ", x)
}

parse_int <- function(x, default) {
  if (is.null(x)) return(as.integer(default))
  out <- suppressWarnings(as.integer(x))
  if (!is.finite(out) || is.na(out)) stop("Could not parse integer argument: ", x)
  out
}

parse_num <- function(x, default) {
  if (is.null(x)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(x))
  if (!is.finite(out) || is.na(out)) stop("Could not parse numeric argument: ", x)
  out
}

parse_cli_settings <- function(args) {
  kv <- read_key_value_args(args)
  smoke_flag <- parse_flag(kv$smoke, default = FALSE)
  output_prefix_supplied <- !is.null(kv$output_prefix)

  settings <- list(
    smoke = smoke_flag,
    seed = parse_int(kv$seed, 20260413L),
    pilot_reps = parse_int(kv$pilot_reps, 200L),
    confirmatory_reps = parse_int(kv$confirmatory_reps, 1000L),
    run_confirmatory = parse_flag(kv$run_confirmatory, default = TRUE),
    mesh_sample_size = parse_int(kv$mesh_sample_size, 50L),
    mu0 = parse_num(kv$mu0, 0.3),
    noise_sd = parse_num(kv$noise_sd, 0.2),
    max_edge = parse_num(kv$max_edge, 0.05),
    refined_max_edge = parse_num(kv$refined_max_edge, 0.025),
    tol = parse_num(kv$tol, 1e-6),
    points_per_unit = parse_num(kv$points_per_unit, 80),
    output_prefix = if (output_prefix_supplied) kv$output_prefix else "iid_vs_matern_selection",
    baseline_regime = if (is.null(kv$baseline_regime)) "fixed_domain" else as.character(kv$baseline_regime),
    baseline_tau0 = parse_num(kv$baseline_tau0, 0.5),
    baseline_n = parse_int(kv$baseline_n, 80L)
  )

  if (settings$pilot_reps < 1L) stop("pilot_reps must be >= 1.")
  if (settings$confirmatory_reps < settings$pilot_reps) {
    stop("confirmatory_reps must be >= pilot_reps.")
  }
  if (settings$mesh_sample_size < 0L) stop("mesh_sample_size must be >= 0.")
  if (settings$noise_sd <= 0) stop("noise_sd must be > 0.")
  if (settings$max_edge <= 0 || settings$refined_max_edge <= 0) {
    stop("Both max_edge and refined_max_edge must be > 0.")
  }
  if (settings$points_per_unit <= 0) stop("points_per_unit must be > 0.")
  if (settings$tol < 0) stop("tol must be >= 0.")

  if (settings$smoke) {
    settings$pilot_reps <- 10L
    settings$confirmatory_reps <- 10L
    settings$run_confirmatory <- FALSE
    settings$mesh_sample_size <- min(settings$mesh_sample_size, 5L)
    settings$baseline_regime <- "fixed_domain"
    settings$baseline_tau0 <- 0.5
    settings$baseline_n <- 40L
    if (!output_prefix_supplied) {
      settings$output_prefix <- "iid_vs_matern_selection_smoke"
    }
  }

  settings
}

settings <- parse_cli_settings(commandArgs(trailingOnly = TRUE))

make_cell_id <- function(regime, tau0, n) {
  paste0(regime, "__tau0_", format(tau0, trim = TRUE), "__n_", n)
}

make_locations <- function(regime, n, points_per_unit) {
  if (regime == "fixed_domain") {
    return(seq(0, 1, length.out = n))
  }
  if (regime == "increasing_domain") {
    domain_length <- n / points_per_unit
    return(seq(0, domain_length, length.out = n))
  }
  stop("Unknown regime: ", regime)
}

build_design <- function(settings) {
  if (settings$smoke) {
    out <- data.frame(
      regime = "fixed_domain",
      tau0 = 0.5,
      n = 40L,
      stringsAsFactors = FALSE
    )
  } else {
    out <- rbind(
      expand.grid(
        regime = "fixed_domain",
        tau0 = c(0.2, 0.5),
        n = c(40L, 80L, 160L, 320L),
        stringsAsFactors = FALSE
      ),
      expand.grid(
        regime = "increasing_domain",
        tau0 = c(0.2, 0.5),
        n = c(80L, 160L, 320L, 640L),
        stringsAsFactors = FALSE
      )
    )
  }

  out$cell_id <- mapply(make_cell_id, out$regime, out$tau0, out$n, USE.NAMES = FALSE)
  out$cell_index <- seq_len(nrow(out))
  out$regime_label <- ifelse(out$regime == "fixed_domain", "fixed-domain", "increasing-domain")
  out$domain_length <- vapply(
    seq_len(nrow(out)),
    function(i) {
      loc <- make_locations(out$regime[i], out$n[i], settings$points_per_unit)
      diff(range(loc))
    },
    numeric(1)
  )

  out
}

design <- build_design(settings)

baseline_cell_id <- with(
  subset(
    design,
    regime == settings$baseline_regime &
      tau0 == settings$baseline_tau0 &
      n == settings$baseline_n
  ),
  cell_id
)

if (!length(baseline_cell_id)) {
  baseline_cell_id <- character(0)
}

selection_from_delta <- function(delta, tol) {
  if (!is.finite(delta)) return(NA_character_)
  if (delta > tol) return("matern")
  if (delta < -tol) return("iid")
  "tie"
}

wilson_interval <- function(successes, total, conf = 0.95) {
  if (!is.finite(total) || total <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- stats::qnorm(1 - (1 - conf) / 2)
  p_hat <- successes / total
  denom <- 1 + z^2 / total
  center <- (p_hat + z^2 / (2 * total)) / denom
  half_width <- (z / denom) * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * total)) / total)

  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

safe_quantiles <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    out <- rep(NA_real_, length(probs))
    names(out) <- paste0("q", round(100 * probs))
    return(out)
  }

  stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
}

format_md_table <- function(df, digits = 4) {
  df_fmt <- df
  for (j in seq_along(df_fmt)) {
    if (is.numeric(df_fmt[[j]])) {
      df_fmt[[j]] <- format(round(df_fmt[[j]], digits), nsmall = digits, trim = TRUE)
    }
  }

  header <- paste(names(df_fmt), collapse = " | ")
  sep <- paste(rep("---", ncol(df_fmt)), collapse = " | ")
  body <- apply(df_fmt, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, sep, body), collapse = "\n")
}

simulate_iid_observations <- function(n, mu0, tau0, noise_sd, seed) {
  set.seed(seed)
  theta <- stats::rnorm(n, mean = mu0, sd = tau0)
  s <- rep(noise_sd, n)
  x <- theta + stats::rnorm(n, sd = s)

  list(theta = theta, x = x, s = s)
}

prepare_cell <- function(cell, settings, max_edge = settings$max_edge) {
  locations <- make_locations(cell$regime, cell$n, settings$points_per_unit)
  fit_fun <- ebnm_Matern_generator(
    locations = locations,
    max.edge = max_edge,
    link = "identity"
  )

  list(
    locations = locations,
    fit_fun = fit_fun
  )
}

fit_one_replicate <- function(cell, prepared_cell, rep_id, stage, settings) {
  seed <- as.integer(settings$seed + cell$cell_index * 100000L + rep_id)
  sim <- simulate_iid_observations(
    n = cell$n,
    mu0 = settings$mu0,
    tau0 = cell$tau0,
    noise_sd = settings$noise_sd,
    seed = seed
  )

  matern_fit <- tryCatch(
    prepared_cell$fit_fun(sim$x, sim$s),
    error = function(e) e
  )
  iid_fit <- tryCatch(
    ebnm::ebnm_normal(
      sim$x,
      sim$s,
      mode = "estimate",
      scale = "estimate"
    ),
    error = function(e) e
  )

  if (inherits(matern_fit, "error") || inherits(iid_fit, "error")) {
    status <- if (inherits(matern_fit, "error")) "matern_error" else "iid_error"
    error_message <- if (inherits(matern_fit, "error")) {
      conditionMessage(matern_fit)
    } else {
      conditionMessage(iid_fit)
    }

    return(data.frame(
      cell_id = cell$cell_id,
      cell_index = cell$cell_index,
      regime = cell$regime,
      regime_label = cell$regime_label,
      tau0 = cell$tau0,
      n = cell$n,
      domain_length = cell$domain_length,
      stage = stage,
      rep = rep_id,
      seed = seed,
      status = status,
      error_message = error_message,
      loglik_matern = NA_real_,
      loglik_iid = NA_real_,
      delta = NA_real_,
      selection = NA_character_,
      fitted_range = NA_real_,
      fitted_sigma = NA_real_,
      fitted_beta = NA_real_,
      range_over_mesh = NA_real_,
      iid_mean = NA_real_,
      iid_sigma = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  loglik_matern <- as.numeric(matern_fit$log_likelihood)
  loglik_iid <- as.numeric(iid_fit$log_likelihood)
  delta <- loglik_matern - loglik_iid

  data.frame(
    cell_id = cell$cell_id,
    cell_index = cell$cell_index,
    regime = cell$regime,
    regime_label = cell$regime_label,
    tau0 = cell$tau0,
    n = cell$n,
    domain_length = cell$domain_length,
    stage = stage,
    rep = rep_id,
    seed = seed,
    status = "ok",
    error_message = NA_character_,
    loglik_matern = loglik_matern,
    loglik_iid = loglik_iid,
    delta = delta,
    selection = selection_from_delta(delta, settings$tol),
    fitted_range = exp(matern_fit$fitted_g$theta),
    fitted_sigma = matern_fit$fitted_g$sigma,
    fitted_beta = matern_fit$fitted_beta,
    range_over_mesh = exp(matern_fit$fitted_g$theta) / settings$max_edge,
    iid_mean = iid_fit$fitted_g$mean,
    iid_sigma = iid_fit$fitted_g$sd,
    stringsAsFactors = FALSE
  )
}

run_stage_for_cells <- function(design_subset, rep_ids, stage, settings) {
  if (!length(rep_ids) || !nrow(design_subset)) {
    return(data.frame())
  }

  rows <- vector("list", length = nrow(design_subset) * length(rep_ids))
  row_counter <- 1L

  for (i in seq_len(nrow(design_subset))) {
    cell <- design_subset[i, , drop = FALSE]
    message(
      sprintf(
        "[%s] %s (%d reps)",
        stage,
        cell$cell_id,
        length(rep_ids)
      )
    )

    prepared_cell <- prepare_cell(cell, settings, max_edge = settings$max_edge)
    for (rep_id in rep_ids) {
      rows[[row_counter]] <- fit_one_replicate(
        cell = cell,
        prepared_cell = prepared_cell,
        rep_id = rep_id,
        stage = stage,
        settings = settings
      )
      row_counter <- row_counter + 1L
    }
  }

  do.call(rbind, rows)
}

summarize_draws <- function(draws, settings) {
  if (!nrow(draws)) {
    return(data.frame())
  }

  pieces <- lapply(split(draws, draws$cell_id), function(df_cell) {
    df_valid <- df_cell[df_cell$status == "ok", , drop = FALSE]
    n_valid <- nrow(df_valid)
    n_error <- sum(df_cell$status != "ok")

    n_iid <- sum(df_valid$selection == "iid", na.rm = TRUE)
    n_matern <- sum(df_valid$selection == "matern", na.rm = TRUE)
    n_tie <- sum(df_valid$selection == "tie", na.rm = TRUE)

    p_iid <- if (n_valid > 0) n_iid / n_valid else NA_real_
    p_matern <- if (n_valid > 0) n_matern / n_valid else NA_real_
    p_tie <- if (n_valid > 0) n_tie / n_valid else NA_real_
    ci <- wilson_interval(n_iid, n_valid)

    delta_q <- safe_quantiles(df_valid$delta, c(0.1, 0.5, 0.9))
    range_q <- safe_quantiles(df_valid$fitted_range, c(0.1, 0.5, 0.9))

    df_matern_selected <- df_valid[df_valid$selection == "matern", , drop = FALSE]
    range_over_mesh_selected_q <- safe_quantiles(
      df_matern_selected$range_over_mesh,
      c(0.1, 0.5, 0.9)
    )

    data.frame(
      cell_id = df_cell$cell_id[1],
      regime = df_cell$regime[1],
      regime_label = df_cell$regime_label[1],
      tau0 = df_cell$tau0[1],
      n = df_cell$n[1],
      domain_length = df_cell$domain_length[1],
      n_total = nrow(df_cell),
      n_valid = n_valid,
      n_error = n_error,
      n_iid = n_iid,
      n_matern = n_matern,
      n_tie = n_tie,
      p_iid = p_iid,
      p_matern = p_matern,
      p_tie = p_tie,
      p_iid_ci_low = ci["lower"],
      p_iid_ci_high = ci["upper"],
      roughly_half = isTRUE(
        is.finite(p_iid) &&
          p_iid >= 0.4 &&
          p_iid <= 0.6 &&
          ci["lower"] <= 0.5 &&
          ci["upper"] >= 0.5
      ),
      mean_delta = if (n_valid > 0) mean(df_valid$delta) else NA_real_,
      median_delta = delta_q[2],
      q10_delta = delta_q[1],
      q90_delta = delta_q[3],
      q10_range = range_q[1],
      q50_range = range_q[2],
      q90_range = range_q[3],
      matern_selected_count = nrow(df_matern_selected),
      matern_selected_range_lt_mesh = if (nrow(df_matern_selected) > 0) {
        mean(df_matern_selected$fitted_range < settings$max_edge)
      } else {
        NA_real_
      },
      q10_range_over_mesh_selected = range_over_mesh_selected_q[1],
      q50_range_over_mesh_selected = range_over_mesh_selected_q[2],
      q90_range_over_mesh_selected = range_over_mesh_selected_q[3],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  out[order(out$regime, out$tau0, out$n), , drop = FALSE]
}

select_confirmatory_cells <- function(pilot_summary, baseline_cell_id) {
  if (!nrow(pilot_summary)) return(character(0))

  ambiguous <- pilot_summary$cell_id[
    is.finite(pilot_summary$p_iid) &
      pilot_summary$p_iid >= 0.35 &
      pilot_summary$p_iid <= 0.65
  ]

  unique(c(baseline_cell_id, ambiguous))
}

run_mesh_sensitivity <- function(draws, design, settings) {
  if (!nrow(draws) || settings$mesh_sample_size == 0L) {
    return(data.frame())
  }

  selected <- draws[
    draws$status == "ok" &
      draws$selection == "matern" &
      is.finite(draws$fitted_range) &
      draws$fitted_range < settings$max_edge,
    ,
    drop = FALSE
  ]

  if (!nrow(selected)) {
    return(data.frame())
  }

  out <- list()
  counter <- 1L

  for (cell_id in unique(selected$cell_id)) {
    cell <- design[design$cell_id == cell_id, , drop = FALSE]
    rows_cell <- selected[selected$cell_id == cell_id, , drop = FALSE]
    set.seed(settings$seed + cell$cell_index[1] * 1000L + 77L)
    chosen_idx <- if (nrow(rows_cell) > settings$mesh_sample_size) {
      sample(seq_len(nrow(rows_cell)), size = settings$mesh_sample_size, replace = FALSE)
    } else {
      seq_len(nrow(rows_cell))
    }
    rows_chosen <- rows_cell[chosen_idx, , drop = FALSE]

    message(
      sprintf(
        "[mesh-sensitivity] %s (%d sampled replicates)",
        cell_id,
        nrow(rows_chosen)
      )
    )

    refined_prepared <- prepare_cell(cell, settings, max_edge = settings$refined_max_edge)

    for (k in seq_len(nrow(rows_chosen))) {
      row <- rows_chosen[k, , drop = FALSE]
      sim <- simulate_iid_observations(
        n = cell$n[1],
        mu0 = settings$mu0,
        tau0 = cell$tau0[1],
        noise_sd = settings$noise_sd,
        seed = row$seed[1]
      )

      refined_fit <- tryCatch(
        refined_prepared$fit_fun(sim$x, sim$s),
        error = function(e) e
      )

      if (inherits(refined_fit, "error")) {
        out[[counter]] <- data.frame(
          cell_id = cell_id,
          regime = cell$regime[1],
          regime_label = cell$regime_label[1],
          tau0 = cell$tau0[1],
          n = cell$n[1],
          rep = row$rep[1],
          seed = row$seed[1],
          status = "refined_matern_error",
          error_message = conditionMessage(refined_fit),
          loglik_iid = row$loglik_iid[1],
          loglik_matern_original = row$loglik_matern[1],
          loglik_matern_refined = NA_real_,
          delta_original = row$delta[1],
          delta_refined = NA_real_,
          selection_original = row$selection[1],
          selection_refined = NA_character_,
          fitted_range_original = row$fitted_range[1],
          fitted_range_refined = NA_real_,
          range_ratio_refined = NA_real_,
          flip_to_iid = NA,
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
        next
      }

      loglik_matern_refined <- as.numeric(refined_fit$log_likelihood)
      delta_refined <- loglik_matern_refined - row$loglik_iid[1]
      selection_refined <- selection_from_delta(delta_refined, settings$tol)
      fitted_range_refined <- exp(refined_fit$fitted_g$theta)

      out[[counter]] <- data.frame(
        cell_id = cell_id,
        regime = cell$regime[1],
        regime_label = cell$regime_label[1],
        tau0 = cell$tau0[1],
        n = cell$n[1],
        rep = row$rep[1],
        seed = row$seed[1],
        status = "ok",
        error_message = NA_character_,
        loglik_iid = row$loglik_iid[1],
        loglik_matern_original = row$loglik_matern[1],
        loglik_matern_refined = loglik_matern_refined,
        delta_original = row$delta[1],
        delta_refined = delta_refined,
        selection_original = row$selection[1],
        selection_refined = selection_refined,
        fitted_range_original = row$fitted_range[1],
        fitted_range_refined = fitted_range_refined,
        range_ratio_refined = fitted_range_refined / row$fitted_range[1],
        flip_to_iid = identical(row$selection[1], "matern") && identical(selection_refined, "iid"),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }

  do.call(rbind, out)
}

summarize_mesh_sensitivity <- function(mesh_draws) {
  if (!nrow(mesh_draws)) {
    return(data.frame())
  }

  pieces <- lapply(split(mesh_draws, mesh_draws$cell_id), function(df_cell) {
    df_valid <- df_cell[df_cell$status == "ok", , drop = FALSE]

    data.frame(
      cell_id = df_cell$cell_id[1],
      regime = df_cell$regime[1],
      regime_label = df_cell$regime_label[1],
      tau0 = df_cell$tau0[1],
      n = df_cell$n[1],
      n_sampled = nrow(df_cell),
      n_valid = nrow(df_valid),
      n_flip_to_iid = sum(df_valid$flip_to_iid, na.rm = TRUE),
      p_flip_to_iid = if (nrow(df_valid) > 0) {
        mean(df_valid$flip_to_iid, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_range_ratio_refined = if (nrow(df_valid) > 0) {
        stats::median(df_valid$range_ratio_refined, na.rm = TRUE)
      } else {
        NA_real_
      },
      median_delta_change = if (nrow(df_valid) > 0) {
        stats::median(df_valid$delta_refined - df_valid$delta_original, na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  out[order(out$regime, out$tau0, out$n), , drop = FALSE]
}

draw_selection_probability_plot <- function(summary_df, settings) {
  if (!nrow(summary_df)) return(NULL)

  file_path <- file.path(fig_dir, paste0(settings$output_prefix, "_selection_probability.png"))
  tau_levels <- sort(unique(summary_df$tau0))
  regime_levels <- c("fixed_domain", "increasing_domain")

  grDevices::png(file_path, width = 1400, height = 900, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(length(tau_levels), length(regime_levels)), mar = c(4, 4, 3, 1))

  for (tau0 in tau_levels) {
    for (regime in regime_levels) {
      df_plot <- summary_df[summary_df$tau0 == tau0 & summary_df$regime == regime, , drop = FALSE]
      df_plot <- df_plot[order(df_plot$n), , drop = FALSE]

      if (!nrow(df_plot)) {
        graphics::plot.new()
        next
      }

      graphics::plot(
        df_plot$n,
        df_plot$p_iid,
        type = "b",
        ylim = c(0, 1),
        xlab = "n",
        ylab = "P(select iid)",
        main = sprintf("%s, tau0 = %.1f", df_plot$regime_label[1], tau0),
        pch = 19,
        col = "navy"
      )
      graphics::arrows(
        x0 = df_plot$n,
        y0 = df_plot$p_iid_ci_low,
        x1 = df_plot$n,
        y1 = df_plot$p_iid_ci_high,
        code = 3,
        angle = 90,
        length = 0.04,
        col = "gray40"
      )
      graphics::abline(h = 0.5, col = "firebrick", lty = 2, lwd = 1.5)
    }
  }

  invisible(file_path)
}

draw_delta_vs_range_plot <- function(draws, settings) {
  draws <- draws[draws$status == "ok", , drop = FALSE]
  if (!nrow(draws)) return(NULL)

  file_path <- file.path(fig_dir, paste0(settings$output_prefix, "_delta_vs_range.png"))
  col_map <- c(
    "0.2" = grDevices::rgb(0.05, 0.35, 0.75, 0.35),
    "0.5" = grDevices::rgb(0.80, 0.25, 0.10, 0.35)
  )

  grDevices::png(file_path, width = 1400, height = 700, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  for (regime in c("fixed_domain", "increasing_domain")) {
    df_plot <- draws[draws$regime == regime, , drop = FALSE]
    if (!nrow(df_plot)) {
      graphics::plot.new()
      graphics::title(main = if (regime == "fixed_domain") "fixed-domain" else "increasing-domain")
      next
    }

    x_vals <- log10(pmax(df_plot$fitted_range, .Machine$double.xmin))
    colors <- col_map[format(df_plot$tau0, trim = TRUE)]

    graphics::plot(
      x_vals,
      df_plot$delta,
      pch = 16,
      cex = 0.6,
      col = colors,
      xlab = "log10(fitted range)",
      ylab = "delta = loglik_Matern - loglik_iid",
      main = if (regime == "fixed_domain") "fixed-domain" else "increasing-domain"
    )
    graphics::abline(h = 0, col = "firebrick", lty = 2, lwd = 1.5)
    graphics::abline(v = log10(settings$max_edge), col = "gray40", lty = 3, lwd = 1.2)
    graphics::legend(
      "topright",
      legend = c("tau0 = 0.2", "tau0 = 0.5"),
      col = unname(col_map[c("0.2", "0.5")]),
      pch = 16,
      bty = "n"
    )
  }

  invisible(file_path)
}

write_markdown_summary <- function(summary_df,
                                   pilot_summary_df,
                                   mesh_summary_df,
                                   confirmatory_cells,
                                   settings) {
  file_path <- file.path(results_dir, paste0(settings$output_prefix, "_summary.md"))

  summary_table <- summary_df[, c(
    "regime_label",
    "tau0",
    "n",
    "n_valid",
    "p_iid",
    "p_iid_ci_low",
    "p_iid_ci_high",
    "p_matern",
    "p_tie",
    "roughly_half",
    "q50_range_over_mesh_selected"
  )]
  names(summary_table) <- c(
    "regime",
    "tau0",
    "n",
    "n_valid",
    "p_iid",
    "ci_low",
    "ci_high",
    "p_matern",
    "p_tie",
    "roughly_half",
    "median_range_over_mesh_if_matern_selected"
  )

  lines <- c(
    "# IID vs Matern Selection Study",
    "",
    "## Run Configuration",
    "",
    paste0("- mode: ", if (settings$smoke) "smoke" else "full"),
    paste0("- seed base: ", settings$seed),
    paste0("- pilot reps per cell: ", settings$pilot_reps),
    paste0("- confirmatory reps target: ", settings$confirmatory_reps),
    paste0("- confirmatory enabled: ", settings$run_confirmatory),
    paste0("- mu0: ", settings$mu0),
    paste0("- noise sd: ", settings$noise_sd),
    paste0("- max.edge: ", settings$max_edge),
    paste0("- refined max.edge: ", settings$refined_max_edge),
    paste0("- selection tolerance: ", settings$tol),
    "",
    "## Final Summary",
    "",
    format_md_table(summary_table, digits = 4),
    ""
  )

  if (nrow(pilot_summary_df) && !settings$smoke) {
    pilot_table <- pilot_summary_df[, c("regime_label", "tau0", "n", "p_iid", "p_iid_ci_low", "p_iid_ci_high")]
    names(pilot_table) <- c("regime", "tau0", "n", "pilot_p_iid", "pilot_ci_low", "pilot_ci_high")
    lines <- c(
      lines,
      "## Pilot Summary",
      "",
      format_md_table(pilot_table, digits = 4),
      ""
    )
  }

  if (length(confirmatory_cells)) {
    lines <- c(
      lines,
      "## Confirmatory Cells",
      "",
      paste0("- ", confirmatory_cells),
      ""
    )
  }

  if (nrow(mesh_summary_df)) {
    mesh_table <- mesh_summary_df[, c(
      "regime_label",
      "tau0",
      "n",
      "n_valid",
      "p_flip_to_iid",
      "median_range_ratio_refined",
      "median_delta_change"
    )]
    names(mesh_table)[1] <- "regime"
    lines <- c(
      lines,
      "## Mesh Sensitivity",
      "",
      format_md_table(mesh_table, digits = 4),
      ""
    )
  }

  lines <- c(
    lines,
    "## Figures",
    "",
    paste0("![Selection probability](figures/", settings$output_prefix, "_selection_probability.png)"),
    "",
    paste0("![Delta vs range](figures/", settings$output_prefix, "_delta_vs_range.png)")
  )

  writeLines(lines, con = file_path)
  invisible(file_path)
}

save_csv <- function(df, filename) {
  utils::write.csv(df, file = filename, row.names = FALSE)
}

pilot_rep_ids <- seq_len(settings$pilot_reps)
pilot_draws <- run_stage_for_cells(design, pilot_rep_ids, "pilot", settings)
pilot_summary <- summarize_draws(pilot_draws, settings)

save_csv(
  pilot_draws,
  file.path(results_dir, paste0(settings$output_prefix, "_pilot_draws.csv"))
)
save_csv(
  pilot_summary,
  file.path(results_dir, paste0(settings$output_prefix, "_pilot_summary.csv"))
)

confirmatory_cells <- if (settings$run_confirmatory) {
  select_confirmatory_cells(pilot_summary, baseline_cell_id)
} else {
  character(0)
}

confirmatory_draws <- data.frame()
if (length(confirmatory_cells) && settings$confirmatory_reps > settings$pilot_reps) {
  confirmatory_design <- design[design$cell_id %in% confirmatory_cells, , drop = FALSE]
  confirmatory_rep_ids <- seq.int(settings$pilot_reps + 1L, settings$confirmatory_reps)
  confirmatory_draws <- run_stage_for_cells(
    confirmatory_design,
    confirmatory_rep_ids,
    "confirmatory",
    settings
  )
}

draws <- if (nrow(confirmatory_draws)) {
  rbind(pilot_draws, confirmatory_draws)
} else {
  pilot_draws
}

summary_df <- summarize_draws(draws, settings)
summary_df$confirmatory_selected <- summary_df$cell_id %in% confirmatory_cells

mesh_draws <- run_mesh_sensitivity(draws, design, settings)
mesh_summary <- summarize_mesh_sensitivity(mesh_draws)

save_csv(
  draws,
  file.path(results_dir, paste0(settings$output_prefix, "_draws.csv"))
)
save_csv(
  summary_df,
  file.path(results_dir, paste0(settings$output_prefix, "_summary.csv"))
)

if (nrow(mesh_draws)) {
  save_csv(
    mesh_draws,
    file.path(results_dir, paste0(settings$output_prefix, "_mesh_sensitivity.csv"))
  )
  save_csv(
    mesh_summary,
    file.path(results_dir, paste0(settings$output_prefix, "_mesh_sensitivity_summary.csv"))
  )
}

draw_selection_probability_plot(summary_df, settings)
draw_delta_vs_range_plot(draws, settings)
write_markdown_summary(summary_df, pilot_summary, mesh_summary, confirmatory_cells, settings)

message("Completed iid-vs-Matern selection study.")
