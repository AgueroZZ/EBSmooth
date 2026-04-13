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
  stop("This internal 2D comparison study requires pkgload.")
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

parse_num_pair <- function(x, default) {
  if (is.null(x)) return(as.numeric(default))
  pieces <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  out <- suppressWarnings(as.numeric(pieces))
  if (length(out) != 2L || any(!is.finite(out))) {
    stop("Expected a comma-separated numeric pair, got: ", x)
  }
  out
}

parse_cli_settings <- function(args) {
  kv <- read_key_value_args(args)
  smoke_flag <- parse_flag(kv$smoke, default = FALSE)
  output_prefix_supplied <- !is.null(kv$output_prefix)

  settings <- list(
    smoke = smoke_flag,
    seed = parse_int(kv$seed, 20260413L),
    n_rep = parse_int(kv$n_rep, 6L),
    true_range = parse_num(kv$true_range, 0.22),
    true_sigma = parse_num(kv$true_sigma, 0.35),
    true_beta0 = parse_num(kv$true_beta0, 0.2),
    iid_tau0 = parse_num(kv$iid_tau0, 0.35),
    noise_sd = parse_num(kv$noise_sd, 0.15),
    alpha = parse_num(kv$alpha, 2),
    max_edge = parse_num_pair(kv$max_edge, c(0.08, 0.12)),
    tol = parse_num(kv$tol, 1e-6),
    use_pc_prior = parse_flag(kv$use_pc_prior, default = FALSE),
    pc_range_anchor = parse_num(kv$pc_range_anchor, 0.22),
    pc_range_alpha = parse_num(kv$pc_range_alpha, 0.1),
    pc_sigma_anchor = parse_num(kv$pc_sigma_anchor, 0.35),
    pc_sigma_alpha = parse_num(kv$pc_sigma_alpha, 0.5),
    output_prefix = if (output_prefix_supplied) kv$output_prefix else "2d_matern_vs_iid_comparison"
  )

  if (settings$n_rep < 1L) stop("n_rep must be >= 1.")
  if (settings$true_range <= 0 || settings$true_sigma <= 0 || settings$iid_tau0 <= 0) {
    stop("Range and sigma settings must be > 0.")
  }
  if (settings$noise_sd <= 0) stop("noise_sd must be > 0.")
  if (any(settings$max_edge <= 0)) stop("max_edge entries must be > 0.")
  if (settings$pc_range_anchor <= 0 || settings$pc_sigma_anchor <= 0) {
    stop("PC-prior anchors must be > 0.")
  }
  if (settings$pc_range_alpha <= 0 || settings$pc_range_alpha >= 1 ||
      settings$pc_sigma_alpha <= 0 || settings$pc_sigma_alpha >= 1) {
    stop("PC-prior alphas must satisfy 0 < alpha < 1.")
  }

  if (settings$smoke) {
    settings$n_rep <- 2L
    if (!output_prefix_supplied) {
      settings$output_prefix <- if (settings$use_pc_prior) {
        "2d_matern_vs_iid_comparison_pc_smoke"
      } else {
        "2d_matern_vs_iid_comparison_smoke"
      }
    }
  }
  if (!output_prefix_supplied && settings$use_pc_prior && !settings$smoke) {
    settings$output_prefix <- "2d_matern_vs_iid_comparison_pc"
  }

  settings
}

settings <- parse_cli_settings(commandArgs(trailingOnly = TRUE))

build_grid_sizes <- function(smoke = FALSE) {
  if (smoke) {
    return(list(c(12L, 10L)))
  }
  list(c(18L, 14L), c(24L, 18L))
}

grid_sizes <- build_grid_sizes(settings$smoke)

make_grid_locations <- function(grid_dim) {
  grid <- expand.grid(
    x = seq(0, 1, length.out = grid_dim[1]),
    y = seq(0, 1, length.out = grid_dim[2])
  )
  list(
    grid = grid,
    locations = as.matrix(grid),
    nx = grid_dim[1],
    ny = grid_dim[2],
    n = nrow(grid)
  )
}

simulate_exact_matern <- function(locations,
                                  range,
                                  sigma,
                                  beta0,
                                  noise_sd,
                                  alpha = 2,
                                  max.edge = NULL) {
  loc_info <- .normalize_locations(locations)
  meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)
  spde_template <- INLA::inla.spde2.matern(meshA$mesh, alpha = alpha)

  state <- .exact_matern_state(
    x = rep(beta0, nrow(loc_info$loc)),
    s = rep(noise_sd, nrow(loc_info$loc)),
    A = meshA$A,
    spde_template = spde_template,
    alpha = alpha,
    d = loc_info$d,
    log_range = log(range),
    log_sigma = log(sigma),
    beta0 = beta0
  )

  w <- as.numeric(
    LaplacesDemon::rmvnp(
      1,
      mu = rep(0, nrow(state$Q)),
      Omega = as.matrix(state$Q)
    )
  )

  mean_surface <- as.numeric(beta0 + meshA$A %*% w)
  x <- mean_surface + stats::rnorm(length(mean_surface), sd = noise_sd)

  list(
    x = x,
    s = rep(noise_sd, length(mean_surface)),
    mean_surface = mean_surface
  )
}

simulate_iid_surface <- function(locations,
                                 beta0,
                                 tau0,
                                 noise_sd) {
  n <- nrow(locations)
  mean_surface <- stats::rnorm(n, mean = beta0, sd = tau0)
  x <- mean_surface + stats::rnorm(n, sd = noise_sd)

  list(
    x = x,
    s = rep(noise_sd, n),
    mean_surface = mean_surface
  )
}

selection_from_delta <- function(delta, tol) {
  if (!is.finite(delta)) return(NA_character_)
  if (delta > tol) return("matern")
  if (delta < -tol) return("iid")
  "tie"
}

safe_quantiles <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
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

fit_model_pair <- function(locations, x, s, settings) {
  g_init <- Matern(
    theta = log(settings$true_range * 0.9),
    sigma = max(settings$true_sigma, settings$iid_tau0) * 1.05
  )
  pc_penalty <- if (settings$use_pc_prior) {
    list(
      range = c(settings$pc_range_anchor, settings$pc_range_alpha),
      sigma = c(settings$pc_sigma_anchor, settings$pc_sigma_alpha)
    )
  } else {
    NULL
  }
  fit_fun <- ebnm_Matern_generator(
    locations = locations,
    max.edge = settings$max_edge,
    link = "identity",
    pc.penalty = pc_penalty
  )

  matern_fit <- tryCatch(
    fit_fun(x, s, g_init = g_init),
    error = function(e) e
  )
  iid_fit <- tryCatch(
    ebnm::ebnm_normal(x, s, mode = "estimate", scale = "estimate"),
    error = function(e) e
  )

  list(matern_fit = matern_fit, iid_fit = iid_fit)
}

run_experiment <- function(settings, grid_sizes) {
  rows <- list()
  counter <- 1L
  setting_ids <- c("smooth_truth", "iid_truth")

  for (setting_id in setting_ids) {
    for (grid_dim in grid_sizes) {
      grid_info <- make_grid_locations(grid_dim)
      grid_label <- paste0(grid_info$nx, "x", grid_info$ny)

      message(sprintf("[%s] grid %s, %d reps", setting_id, grid_label, settings$n_rep))

      for (rep_id in seq_len(settings$n_rep)) {
        set.seed(settings$seed + counter * 1000L + rep_id)

        sim <- if (setting_id == "smooth_truth") {
          simulate_exact_matern(
            locations = grid_info$locations,
            range = settings$true_range,
            sigma = settings$true_sigma,
            beta0 = settings$true_beta0,
            noise_sd = settings$noise_sd,
            alpha = settings$alpha,
            max.edge = settings$max_edge
          )
        } else {
          simulate_iid_surface(
            locations = grid_info$locations,
            beta0 = settings$true_beta0,
            tau0 = settings$iid_tau0,
            noise_sd = settings$noise_sd
          )
        }

        fits <- fit_model_pair(grid_info$locations, sim$x, sim$s, settings)
        if (inherits(fits$matern_fit, "error") || inherits(fits$iid_fit, "error")) {
          status <- if (inherits(fits$matern_fit, "error")) "matern_error" else "iid_error"
          error_message <- if (inherits(fits$matern_fit, "error")) {
            conditionMessage(fits$matern_fit)
          } else {
            conditionMessage(fits$iid_fit)
          }

          rows[[counter]] <- data.frame(
            setting = setting_id,
            grid_label = grid_label,
            nx = grid_info$nx,
            ny = grid_info$ny,
            n = grid_info$n,
            rep = rep_id,
            status = status,
            error_message = error_message,
            loglik_matern = NA_real_,
            loglik_iid = NA_real_,
            delta = NA_real_,
            selection = NA_character_,
            est_range = NA_real_,
            est_sigma = NA_real_,
            est_beta = NA_real_,
            range_over_mesh = NA_real_,
            iid_mean = NA_real_,
            iid_sigma = NA_real_,
            surface_rmse = NA_real_,
            surface_corr = NA_real_,
            stringsAsFactors = FALSE
          )
          counter <- counter + 1L
          next
        }

        loglik_matern <- as.numeric(fits$matern_fit$log_likelihood)
        loglik_iid <- as.numeric(fits$iid_fit$log_likelihood)
        est_range <- exp(fits$matern_fit$fitted_g$theta)
        est_sigma <- fits$matern_fit$fitted_g$sigma

        rows[[counter]] <- data.frame(
          setting = setting_id,
          grid_label = grid_label,
          nx = grid_info$nx,
          ny = grid_info$ny,
          n = grid_info$n,
          rep = rep_id,
          status = "ok",
          error_message = NA_character_,
          loglik_matern = loglik_matern,
          loglik_iid = loglik_iid,
          delta = loglik_matern - loglik_iid,
          selection = selection_from_delta(loglik_matern - loglik_iid, settings$tol),
          est_range = est_range,
          est_sigma = est_sigma,
          est_beta = fits$matern_fit$fitted_beta,
          range_over_mesh = est_range / min(settings$max_edge),
          iid_mean = fits$iid_fit$fitted_g$mean,
          iid_sigma = fits$iid_fit$fitted_g$sd,
          surface_rmse = sqrt(mean((fits$matern_fit$posterior$mean - sim$mean_surface)^2)),
          surface_corr = stats::cor(fits$matern_fit$posterior$mean, sim$mean_surface),
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }

  do.call(rbind, rows)
}

summarize_selection <- function(draws) {
  pieces <- lapply(split(draws, interaction(draws$setting, draws$grid_label, drop = TRUE)), function(df) {
    df_valid <- df[df$status == "ok", , drop = FALSE]
    n_valid <- nrow(df_valid)

    data.frame(
      setting = df$setting[1],
      grid_label = df$grid_label[1],
      nx = df$nx[1],
      ny = df$ny[1],
      n = df$n[1],
      n_valid = n_valid,
      p_matern = if (n_valid > 0) mean(df_valid$selection == "matern") else NA_real_,
      p_iid = if (n_valid > 0) mean(df_valid$selection == "iid") else NA_real_,
      p_tie = if (n_valid > 0) mean(df_valid$selection == "tie") else NA_real_,
      mean_delta = if (n_valid > 0) mean(df_valid$delta) else NA_real_,
      median_delta = if (n_valid > 0) stats::median(df_valid$delta) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  out[order(out$setting, out$n), , drop = FALSE]
}

summarize_smooth_truth <- function(draws, settings) {
  df <- draws[draws$setting == "smooth_truth" & draws$status == "ok", , drop = FALSE]
  pieces <- lapply(split(df, df$grid_label), function(dd) {
    data.frame(
      grid_label = dd$grid_label[1],
      nx = dd$nx[1],
      ny = dd$ny[1],
      n = dd$n[1],
      mean_est_range = mean(dd$est_range),
      median_est_range = stats::median(dd$est_range),
      mean_est_sigma = mean(dd$est_sigma),
      median_est_sigma = stats::median(dd$est_sigma),
      mean_abs_err_range = mean(abs(dd$est_range - settings$true_range)),
      mean_abs_err_sigma = mean(abs(dd$est_sigma - settings$true_sigma)),
      mean_surface_rmse = mean(dd$surface_rmse),
      mean_surface_corr = mean(dd$surface_corr),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  out[order(out$n), , drop = FALSE]
}

summarize_iid_truth <- function(draws) {
  df <- draws[draws$setting == "iid_truth" & draws$status == "ok", , drop = FALSE]
  pieces <- lapply(split(df, df$grid_label), function(dd) {
    q_range <- safe_quantiles(dd$est_range, c(0.1, 0.5, 0.9))
    q_range_mesh <- safe_quantiles(dd$range_over_mesh, c(0.1, 0.5, 0.9))

    data.frame(
      grid_label = dd$grid_label[1],
      nx = dd$nx[1],
      ny = dd$ny[1],
      n = dd$n[1],
      mean_est_range = mean(dd$est_range),
      q10_est_range = q_range[1],
      q50_est_range = q_range[2],
      q90_est_range = q_range[3],
      q10_range_over_mesh = q_range_mesh[1],
      q50_range_over_mesh = q_range_mesh[2],
      q90_range_over_mesh = q_range_mesh[3],
      mean_surface_corr = mean(dd$surface_corr),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  out[order(out$n), , drop = FALSE]
}

draw_selection_plot <- function(selection_summary, settings) {
  file_path <- file.path(fig_dir, paste0(settings$output_prefix, "_selection.png"))

  grDevices::png(file_path, width = 1300, height = 700, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))

  for (setting_id in c("smooth_truth", "iid_truth")) {
    df <- selection_summary[selection_summary$setting == setting_id, , drop = FALSE]
    labels <- paste0(df$grid_label, "\n(n=", df$n, ")")
    heights <- rbind(df$p_matern, df$p_iid)

    bp <- graphics::barplot(
      heights,
      beside = TRUE,
      names.arg = labels,
      ylim = c(0, 1),
      col = c("navy", "firebrick"),
      las = 2,
      ylab = "Selection probability",
      main = if (setting_id == "smooth_truth") "Smooth truth" else "IID truth"
    )
    graphics::legend(
      "topright",
      legend = c("Select Matern", "Select iid"),
      fill = c("navy", "firebrick"),
      bty = "n"
    )
    graphics::abline(h = 0.5, lty = 2, col = "gray40")
  }

  invisible(file_path)
}

draw_hyperparameter_plot <- function(draws, settings) {
  file_path <- file.path(fig_dir, paste0(settings$output_prefix, "_hyperparameters.png"))

  grDevices::png(file_path, width = 1400, height = 900, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))

  smooth_df <- draws[draws$setting == "smooth_truth" & draws$status == "ok", , drop = FALSE]
  iid_df <- draws[draws$setting == "iid_truth" & draws$status == "ok", , drop = FALSE]

  if (nrow(smooth_df)) {
    graphics::boxplot(
      est_range ~ grid_label,
      data = smooth_df,
      xlab = "grid",
      ylab = "estimated range",
      main = "Smooth truth: range recovery",
      col = "lightblue"
    )
    graphics::abline(h = settings$true_range, col = "firebrick", lty = 2, lwd = 2)

    graphics::boxplot(
      est_sigma ~ grid_label,
      data = smooth_df,
      xlab = "grid",
      ylab = "estimated sigma",
      main = "Smooth truth: sigma recovery",
      col = "lightblue"
    )
    graphics::abline(h = settings$true_sigma, col = "firebrick", lty = 2, lwd = 2)
  } else {
    graphics::plot.new()
    graphics::plot.new()
  }

  if (nrow(iid_df)) {
    graphics::boxplot(
      log10(est_range) ~ grid_label,
      data = iid_df,
      xlab = "grid",
      ylab = "log10(estimated range)",
      main = "IID truth: fitted range",
      col = "khaki"
    )

    graphics::boxplot(
      range_over_mesh ~ grid_label,
      data = iid_df,
      xlab = "grid",
      ylab = "range / min(max.edge)",
      main = "IID truth: range over mesh",
      col = "khaki"
    )
    graphics::abline(h = 1, col = "firebrick", lty = 2, lwd = 2)
  } else {
    graphics::plot.new()
    graphics::plot.new()
  }

  invisible(file_path)
}

write_markdown_summary <- function(selection_summary,
                                   smooth_summary,
                                   iid_summary,
                                   settings) {
  file_path <- file.path(results_dir, paste0(settings$output_prefix, "_summary.md"))

  selection_table <- selection_summary[, c(
    "setting", "grid_label", "n", "n_valid", "p_matern", "p_iid", "mean_delta", "median_delta"
  )]

  smooth_table <- smooth_summary[, c(
    "grid_label", "n", "mean_est_range", "median_est_range", "mean_est_sigma",
    "median_est_sigma", "mean_abs_err_range", "mean_abs_err_sigma", "mean_surface_corr"
  )]

  iid_table <- iid_summary[, c(
    "grid_label", "n", "mean_est_range", "q10_est_range", "q50_est_range", "q90_est_range",
    "q10_range_over_mesh", "q50_range_over_mesh", "q90_range_over_mesh", "mean_surface_corr"
  )]

  lines <- c(
    "# 2D Matern vs IID Comparison",
    "",
    "## Configuration",
    "",
    paste0("- reps per cell: ", settings$n_rep),
    paste0("- true Matern range: ", settings$true_range),
    paste0("- true Matern sigma: ", settings$true_sigma),
    paste0("- true intercept: ", settings$true_beta0),
    paste0("- iid tau0: ", settings$iid_tau0),
    paste0("- noise sd: ", settings$noise_sd),
    paste0("- max.edge: ", paste(settings$max_edge, collapse = ", ")),
    paste0("- use pc prior: ", settings$use_pc_prior),
    if (settings$use_pc_prior) {
      paste0(
        "- pc.penalty: range=(",
        settings$pc_range_anchor, ", ", settings$pc_range_alpha,
        "), sigma=(",
        settings$pc_sigma_anchor, ", ", settings$pc_sigma_alpha, ")"
      )
    } else {
      "- pc.penalty: none"
    },
    "",
    "## Selection Summary",
    "",
    format_md_table(selection_table, digits = 4),
    "",
    "## Smooth Truth",
    "",
    format_md_table(smooth_table, digits = 4),
    "",
    "## IID Truth",
    "",
    format_md_table(iid_table, digits = 4),
    "",
    "## Figures",
    "",
    paste0("![Selection probabilities](figures/", settings$output_prefix, "_selection.png)"),
    "",
    paste0("![Hyperparameter summaries](figures/", settings$output_prefix, "_hyperparameters.png)")
  )

  writeLines(lines, con = file_path)
  invisible(file_path)
}

draws <- run_experiment(settings, grid_sizes)
selection_summary <- summarize_selection(draws)
smooth_summary <- summarize_smooth_truth(draws, settings)
iid_summary <- summarize_iid_truth(draws)

utils::write.csv(
  draws,
  file = file.path(results_dir, paste0(settings$output_prefix, "_draws.csv")),
  row.names = FALSE
)
utils::write.csv(
  selection_summary,
  file = file.path(results_dir, paste0(settings$output_prefix, "_selection_summary.csv")),
  row.names = FALSE
)
utils::write.csv(
  smooth_summary,
  file = file.path(results_dir, paste0(settings$output_prefix, "_smooth_truth_summary.csv")),
  row.names = FALSE
)
utils::write.csv(
  iid_summary,
  file = file.path(results_dir, paste0(settings$output_prefix, "_iid_truth_summary.csv")),
  row.names = FALSE
)

draw_selection_plot(selection_summary, settings)
draw_hyperparameter_plot(draws, settings)
write_markdown_summary(selection_summary, smooth_summary, iid_summary, settings)

message("Completed 2D Matern-vs-iid comparison study.")
