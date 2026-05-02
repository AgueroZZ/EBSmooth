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
  stop("This internal null study requires pkgload.")
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

parse_num_vector <- function(x, default) {
  if (is.null(x)) return(as.numeric(default))
  pieces <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  out <- suppressWarnings(as.numeric(pieces))
  if (!length(out) || any(!is.finite(out))) {
    stop("Could not parse numeric vector argument: ", x)
  }
  out
}

parse_cli_settings <- function(args) {
  kv <- read_key_value_args(args)
  smoke_flag <- parse_flag(kv$smoke, default = FALSE)
  output_prefix_supplied <- !is.null(kv$output_prefix)

  settings <- list(
    smoke = smoke_flag,
    seed = parse_int(kv$seed, 20260419L),
    pilot_reps = parse_int(kv$pilot_reps, 20L),
    confirmatory_reps = parse_int(kv$confirmatory_reps, 100L),
    nx = parse_int(kv$nx, 18L),
    ny = parse_int(kv$ny, 14L),
    beta0 = parse_num(kv$beta0, 0.2),
    tol = parse_num(kv$tol, 1e-6),
    pilot_noise_sd_grid = parse_num_vector(kv$pilot_noise_sd_grid, c(0.05, 0.1, 0.15, 0.2, 0.3)),
    max_edge_x = parse_num(kv$max_edge_x, 0.08),
    max_edge_y = parse_num(kv$max_edge_y, 0.12),
    pc_range_anchor = parse_num(kv$pc_range_anchor, 0.22),
    pc_range_alpha = parse_num(kv$pc_range_alpha, 0.1),
    pc_sigma_anchor = parse_num(kv$pc_sigma_anchor, 0.2),
    pc_sigma_alpha = parse_num(kv$pc_sigma_alpha, 0.5),
    pc_noise_anchor = parse_num(kv$pc_noise_anchor, 0.2),
    pc_noise_alpha = parse_num(kv$pc_noise_alpha, 0.5),
    output_prefix = if (output_prefix_supplied) kv$output_prefix else "nonspatial_vs_matern_null"
  )

  if (settings$smoke) {
    settings$pilot_reps <- 4L
    settings$confirmatory_reps <- 6L
    settings$nx <- 12L
    settings$ny <- 10L
    if (!output_prefix_supplied) {
      settings$output_prefix <- "nonspatial_vs_matern_null_smoke"
    }
  }

  settings
}

settings <- parse_cli_settings(commandArgs(trailingOnly = TRUE))

make_grid_locations <- function(nx, ny) {
  grid <- expand.grid(
    x = seq(0, 1, length.out = nx),
    y = seq(0, 1, length.out = ny)
  )
  list(
    grid = grid,
    locations = as.matrix(grid),
    n = nrow(grid)
  )
}

grid_info <- make_grid_locations(settings$nx, settings$ny)
matern_setup <- Matern_setup(
  locations = grid_info$locations,
  max.edge = c(settings$max_edge_x, settings$max_edge_y)
)

selection_from_delta <- function(delta, tol) {
  if (!is.finite(delta)) return(NA_character_)
  if (delta > tol) return("matern")
  if (delta < -tol) return("nonspatial")
  "tie"
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

fit_model_pair <- function(x, use_pc_prior, noise_sd) {
  pc_penalty <- if (use_pc_prior) {
    list(
      range = c(settings$pc_range_anchor, settings$pc_range_alpha),
      sigma = c(settings$pc_sigma_anchor, settings$pc_sigma_alpha),
      noise = c(max(noise_sd, settings$pc_noise_anchor), settings$pc_noise_alpha)
    )
  } else {
    NULL
  }

  matern_fit <- tryCatch(
    eb_smoother(
      x = x,
      s = NULL,
      family = "matern",
      setup = matern_setup,
      backend = if (use_pc_prior) "inla_pc" else "exact",
      pc.penalty = pc_penalty
    ),
    error = function(e) e
  )

  nonspatial_fit <- tryCatch(
    eb_smoother(
      x = x,
      s = NULL,
      family = "constant"
    ),
    error = function(e) e
  )

  list(matern_fit = matern_fit, nonspatial_fit = nonspatial_fit)
}

run_replicates <- function(noise_sd, n_rep, use_pc_prior) {
  rows <- vector("list", length = n_rep)

  for (rep_id in seq_len(n_rep)) {
    seed <- settings$seed + rep_id
    set.seed(seed)
    x <- settings$beta0 + stats::rnorm(grid_info$n, sd = noise_sd)

    fits <- fit_model_pair(x = x, use_pc_prior = use_pc_prior, noise_sd = noise_sd)
    if (inherits(fits$matern_fit, "error") || inherits(fits$nonspatial_fit, "error")) {
      rows[[rep_id]] <- data.frame(
        rep = rep_id,
        seed = seed,
        prior_mode = if (use_pc_prior) "pc" else "none",
        noise_sd = noise_sd,
        status = if (inherits(fits$matern_fit, "error")) "matern_error" else "nonspatial_error",
        error_message = if (inherits(fits$matern_fit, "error")) {
          conditionMessage(fits$matern_fit)
        } else {
          conditionMessage(fits$nonspatial_fit)
        },
        loglik_matern = NA_real_,
        loglik_nonspatial = NA_real_,
        delta = NA_real_,
        selection = NA_character_,
        fitted_range = NA_real_,
        fitted_sigma = NA_real_,
        fitted_beta_matern = NA_real_,
        fitted_noise_sd_matern = NA_real_,
        fitted_beta_nonspatial = NA_real_,
        fitted_noise_sd_nonspatial = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    loglik_matern <- as.numeric(fits$matern_fit$log_likelihood)
    loglik_nonspatial <- as.numeric(fits$nonspatial_fit$log_likelihood)

    rows[[rep_id]] <- data.frame(
      rep = rep_id,
      seed = seed,
      prior_mode = if (use_pc_prior) "pc" else "none",
      noise_sd = noise_sd,
      status = "ok",
      error_message = NA_character_,
      loglik_matern = loglik_matern,
      loglik_nonspatial = loglik_nonspatial,
      delta = loglik_matern - loglik_nonspatial,
      selection = selection_from_delta(loglik_matern - loglik_nonspatial, settings$tol),
      fitted_range = exp(fits$matern_fit$fitted_g$theta),
      fitted_sigma = fits$matern_fit$fitted_g$sigma,
      fitted_beta_matern = fits$matern_fit$fitted_beta,
      fitted_noise_sd_matern = fits$matern_fit$fitted_noise_sd,
      fitted_beta_nonspatial = fits$nonspatial_fit$fitted_beta,
      fitted_noise_sd_nonspatial = fits$nonspatial_fit$fitted_noise_sd,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

summarize_draws <- function(draws) {
  pieces <- lapply(split(draws, draws$prior_mode), function(df) {
    df_ok <- df[df$status == "ok", , drop = FALSE]
    data.frame(
      prior_mode = df$prior_mode[1],
      noise_sd = df$noise_sd[1],
      n_total = nrow(df),
      n_valid = nrow(df_ok),
      p_matern = if (nrow(df_ok)) mean(df_ok$selection == "matern") else NA_real_,
      p_nonspatial = if (nrow(df_ok)) mean(df_ok$selection == "nonspatial") else NA_real_,
      p_tie = if (nrow(df_ok)) mean(df_ok$selection == "tie") else NA_real_,
      mean_delta = if (nrow(df_ok)) mean(df_ok$delta) else NA_real_,
      median_delta = if (nrow(df_ok)) stats::median(df_ok$delta) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

plot_selection_summary <- function(summary_df, file_path) {
  png(filename = file_path, width = 900, height = 500, res = 120)
  par(mar = c(4, 5, 4, 1) + 0.1)
  mat <- rbind(summary_df$p_nonspatial, summary_df$p_matern)
  barplot(
    mat,
    beside = FALSE,
    col = c("#4daf4a", "#377eb8"),
    ylim = c(0, 1),
    names.arg = summary_df$prior_mode,
    ylab = "Selection probability",
    main = "Nonspatial vs Matern under null truth"
  )
  legend(
    "topright",
    legend = c("Select nonspatial", "Select matern"),
    fill = c("#4daf4a", "#377eb8"),
    bty = "n"
  )
  dev.off()
}

pilot_rows <- lapply(settings$pilot_noise_sd_grid, function(noise_sd) {
  draws <- run_replicates(noise_sd = noise_sd, n_rep = settings$pilot_reps, use_pc_prior = FALSE)
  df_ok <- draws[draws$status == "ok", , drop = FALSE]
  data.frame(
    noise_sd = noise_sd,
    n_valid = nrow(df_ok),
    p_nonspatial = if (nrow(df_ok)) mean(df_ok$selection == "nonspatial") else NA_real_,
    distance_to_half = if (nrow(df_ok)) abs(mean(df_ok$selection == "nonspatial") - 0.5) else Inf,
    stringsAsFactors = FALSE
  )
})
pilot_summary <- do.call(rbind, pilot_rows)
pilot_summary <- pilot_summary[order(pilot_summary$distance_to_half, pilot_summary$noise_sd), , drop = FALSE]
chosen_noise_sd <- pilot_summary$noise_sd[1]

draws_none <- run_replicates(noise_sd = chosen_noise_sd, n_rep = settings$confirmatory_reps, use_pc_prior = FALSE)
draws_pc <- run_replicates(noise_sd = chosen_noise_sd, n_rep = settings$confirmatory_reps, use_pc_prior = TRUE)
draws <- rbind(draws_none, draws_pc)
summary_df <- summarize_draws(draws)

pilot_csv <- file.path(results_dir, paste0(settings$output_prefix, "_pilot_summary.csv"))
draws_csv <- file.path(results_dir, paste0(settings$output_prefix, "_draws.csv"))
summary_csv <- file.path(results_dir, paste0(settings$output_prefix, "_summary.csv"))
summary_md <- file.path(results_dir, paste0(settings$output_prefix, "_summary.md"))
selection_png <- file.path(fig_dir, paste0(settings$output_prefix, "_selection.png"))

utils::write.csv(pilot_summary, pilot_csv, row.names = FALSE)
utils::write.csv(draws, draws_csv, row.names = FALSE)
utils::write.csv(summary_df, summary_csv, row.names = FALSE)

plot_selection_summary(summary_df, selection_png)

summary_lines <- c(
  "# Nonspatial vs Matern Null Study",
  "",
  paste0("- grid: ", settings$nx, " x ", settings$ny),
  paste0("- confirmatory replicates per prior mode: ", settings$confirmatory_reps),
  paste0("- chosen null noise SD from pilot: ", format(chosen_noise_sd, trim = TRUE)),
  paste0("- selection tolerance: ", settings$tol),
  "",
  "## Pilot Summary",
  "",
  format_md_table(pilot_summary, digits = 4),
  "",
  "## Confirmatory Summary",
  "",
  format_md_table(summary_df, digits = 4),
  "",
  paste0("![Selection summary](figures/", basename(selection_png), ")")
)

writeLines(summary_lines, con = summary_md)

message("Completed nonspatial-vs-Matern null study.")
