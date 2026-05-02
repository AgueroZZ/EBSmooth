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
  stop("This internal MF study requires pkgload.")
}

load_flashier <- function() {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("This internal MF study requires the flashier package.")
  }

  suppressPackageStartupMessages(library(flashier))
  invisible(TRUE)
}

load_ebsmoothr()
load_flashier()

simulation_helpers_path <- file.path(project_root, "attempts_in_EBMF", "simulation_functions.R")
source(simulation_helpers_path, local = globalenv())

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
    seed = parse_int(kv$seed, 20260419L),
    n_reps = parse_int(kv$n_reps, 10L),
    grid_side = parse_int(kv$grid_side, 12L),
    p = parse_int(kv$p, 50L),
    k_spatial = parse_int(kv$k_spatial, 2L),
    k_nonspatial = parse_int(kv$k_nonspatial, 2L),
    sigma_E = parse_num(kv$sigma_E, 0.5),
    tol = parse_num(kv$tol, 1e-6),
    flash_backfit_maxiter = parse_int(kv$flash_backfit_maxiter, 50L),
    flash_conv_tol = parse_num(kv$flash_conv_tol, 1e-3),
    selective_backfit_maxiter = parse_int(kv$selective_backfit_maxiter, 5L),
    output_prefix = if (output_prefix_supplied) kv$output_prefix else "mf_selective_smoothing"
  )

  if (settings$smoke) {
    settings$n_reps <- 2L
    settings$grid_side <- 8L
    settings$p <- 24L
    settings$flash_backfit_maxiter <- 10L
    settings$selective_backfit_maxiter <- 2L
    if (!output_prefix_supplied) {
      settings$output_prefix <- "mf_selective_smoothing_smoke"
    }
  }

  settings
}

settings <- parse_cli_settings(commandArgs(trailingOnly = TRUE))

build_regular_ebnm_list <- function() {
  list(
    flashier::flash_ebnm(
      prior_family = "normal",
      mode = "estimate"
    ),
    flashier::flash_ebnm(
      prior_family = "point_exponential"
    )
  )
}

safe_cor <- function(x, y) {
  value <- suppressWarnings(stats::cor(x, y))
  if (!is.finite(value)) {
    return(NA_real_)
  }
  as.numeric(value)
}

format_md_table <- function(df, digits = 4) {
  if (!nrow(df)) {
    return("(no rows)")
  }

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

make_component_metrics <- function(rep_id,
                                   seed,
                                   stage,
                                   truth_norm,
                                   est_norm,
                                   true_types) {
  K <- min(
    ncol(truth_norm$L_norm),
    ncol(est_norm$L_norm),
    length(true_types)
  )

  if (K == 0L) {
    return(data.frame())
  }

  data.frame(
    rep = rep_id,
    seed = seed,
    stage = stage,
    aligned_component = seq_len(K),
    true_component_type = true_types[seq_len(K)],
    loading_cor = vapply(
      seq_len(K),
      function(k) safe_cor(truth_norm$L_norm[, k], est_norm$L_norm[, k]),
      numeric(1)
    ),
    factor_cor = vapply(
      seq_len(K),
      function(k) safe_cor(truth_norm$F_norm[, k], est_norm$F_norm[, k]),
      numeric(1)
    ),
    loading_mse = colMeans((truth_norm$L_norm[, seq_len(K), drop = FALSE] -
      est_norm$L_norm[, seq_len(K), drop = FALSE])^2),
    factor_mse = colMeans((truth_norm$F_norm[, seq_len(K), drop = FALSE] -
      est_norm$F_norm[, seq_len(K), drop = FALSE])^2),
    stringsAsFactors = FALSE
  )
}

align_screening_summary <- function(screening_summary,
                                    alignment,
                                    true_types,
                                    rep_id,
                                    seed) {
  K <- min(length(alignment$permutation), length(true_types))
  if (K == 0L || !nrow(screening_summary)) {
    return(data.frame())
  }

  matched <- screening_summary[
    match(alignment$permutation[seq_len(K)], screening_summary$factor_k),
    ,
    drop = FALSE
  ]

  matched$rep <- rep_id
  matched$seed <- seed
  matched$estimated_factor_k <- matched$factor_k
  matched$aligned_component <- seq_len(K)
  matched$alignment_sign <- alignment$signs[seq_len(K)]
  matched$true_component_type <- true_types[seq_len(K)]
  matched$label_correct <- matched$status == "ok" &
    matched$selection == matched$true_component_type

  matched[, c(
    "rep",
    "seed",
    "aligned_component",
    "estimated_factor_k",
    "alignment_sign",
    "true_component_type",
    "status",
    "selection",
    "label_correct",
    "delta",
    "loglik_matern",
    "loglik_nonspatial",
    "matern_beta",
    "nonspatial_beta",
    "matern_noise_sd",
    "nonspatial_noise_sd",
    "matern_range",
    "matern_sigma",
    "error_message"
  )]
}

mean_metric <- function(df, stage, true_type, column) {
  sel <- df$stage == stage & df$true_component_type == true_type
  if (!any(sel)) {
    return(NA_real_)
  }
  mean(df[[column]][sel], na.rm = TRUE)
}

summarize_label_accuracy <- function(screening_aligned) {
  if (!nrow(screening_aligned)) {
    return(data.frame())
  }

  groups <- c(unique(as.character(screening_aligned$true_component_type)), "overall")
  rows <- lapply(groups, function(group_name) {
    df <- if (identical(group_name, "overall")) {
      screening_aligned
    } else {
      screening_aligned[screening_aligned$true_component_type == group_name, , drop = FALSE]
    }

    ok <- !is.na(df$label_correct)
    data.frame(
      component_group = group_name,
      n_components = nrow(df),
      n_scored = sum(ok),
      p_correct = if (sum(ok)) mean(df$label_correct[ok]) else NA_real_,
      mean_delta = if (nrow(df)) mean(df$delta, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

summarize_recovery <- function(component_metrics) {
  if (!nrow(component_metrics)) {
    return(data.frame())
  }

  keys <- unique(component_metrics[, c("stage", "true_component_type"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    df <- component_metrics[
      component_metrics$stage == key$stage &
        component_metrics$true_component_type == key$true_component_type,
      ,
      drop = FALSE
    ]

    data.frame(
      stage = key$stage,
      true_component_type = key$true_component_type,
      n_components = nrow(df),
      mean_loading_cor = mean(df$loading_cor, na.rm = TRUE),
      mean_factor_cor = mean(df$factor_cor, na.rm = TRUE),
      mean_loading_mse = mean(df$loading_mse, na.rm = TRUE),
      mean_factor_mse = mean(df$factor_mse, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  overall_rows <- lapply(unique(component_metrics$stage), function(stage_name) {
    df <- component_metrics[component_metrics$stage == stage_name, , drop = FALSE]
    data.frame(
      stage = stage_name,
      true_component_type = "overall",
      n_components = nrow(df),
      mean_loading_cor = mean(df$loading_cor, na.rm = TRUE),
      mean_factor_cor = mean(df$factor_cor, na.rm = TRUE),
      mean_loading_mse = mean(df$loading_mse, na.rm = TRUE),
      mean_factor_mse = mean(df$factor_mse, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, c(rows, overall_rows))
}

bind_rows_safe <- function(frames) {
  nonempty <- Filter(
    function(x) is.data.frame(x) && nrow(x) > 0L,
    frames
  )
  if (!length(nonempty)) {
    return(data.frame())
  }
  do.call(rbind, nonempty)
}

plot_loading_mse <- function(component_metrics, file_path) {
  df <- component_metrics[is.finite(component_metrics$loading_mse), , drop = FALSE]
  if (!nrow(df)) {
    return(invisible(FALSE))
  }

  groups <- interaction(df$true_component_type, df$stage, sep = " : ", drop = TRUE)
  grouped_values <- split(df$loading_mse, groups)
  group_names <- names(grouped_values)

  png(filename = file_path, width = 1100, height = 600, res = 120)
  par(mar = c(8, 5, 4, 1) + 0.1)
  boxplot(
    grouped_values,
    names = group_names,
    las = 2,
    col = ifelse(grepl("original$", group_names), "#9ecae1", "#31a354"),
    ylab = "Normalized loading MSE",
    main = "Selective smoothing in spatial MF"
  )
  abline(h = 0, lty = 2, col = "gray70")
  dev.off()

  invisible(TRUE)
}

run_single_rep <- function(rep_id, settings, ebnm_list_regular) {
  seed <- settings$seed + rep_id - 1L
  simulation <- simulate_spatial_flash_data(
    seed = seed,
    grid_side = settings$grid_side,
    p = settings$p,
    k_spatial = settings$k_spatial,
    k_nonspatial = settings$k_nonspatial,
    sigma_E = settings$sigma_E
  )

  fit_res <- tryCatch(
    run_flash_model(
      Y = simulation$Y,
      S = simulation$sigma_E,
      ebnm_list = ebnm_list_regular,
      greedy_Kmax = simulation$K,
      backfit_maxiter = settings$flash_backfit_maxiter,
      conv_tol = settings$flash_conv_tol
    ),
    error = function(e) e
  )

  if (inherits(fit_res, "error")) {
    run_summary <- data.frame(
      rep = rep_id,
      seed = seed,
      status = "flash_error",
      error_message = conditionMessage(fit_res),
      n_estimated_components = NA_integer_,
      n_aligned_components = NA_integer_,
      n_screened_spatial = NA_integer_,
      label_accuracy = NA_real_,
      original_spatial_loading_mse = NA_real_,
      selective_spatial_loading_mse = NA_real_,
      original_nonspatial_loading_mse = NA_real_,
      selective_nonspatial_loading_mse = NA_real_,
      max_abs_nonspatial_loading_change = NA_real_,
      flash_elapsed_sec = NA_real_,
      selective_elapsed_sec = NA_real_,
      stringsAsFactors = FALSE
    )

    return(list(
      run_summary = run_summary,
      screening = data.frame(),
      screening_aligned = data.frame(),
      component_metrics = data.frame()
    ))
  }

  fit <- fit_res$fit
  flash_elapsed_sec <- as.numeric(fit_res$elapsed[["elapsed"]])

  screening <- screen_flash_loadings_with_eb_smoother(
    fit = fit,
    locations = simulation$locations,
    tol = settings$tol
  )
  screening_raw <- screening$screening_summary
  screening_raw$rep <- rep_id
  screening_raw$seed <- seed

  selective_time <- system.time({
    selective <- tryCatch(
      selective_smooth_flash_loadings(
        fit = fit,
        Y = simulation$Y,
        S = simulation$sigma_E,
        locations = simulation$locations,
        ebnm_fn = ebnm_list_regular,
        tol = settings$tol,
        refit_factors = TRUE,
        backfit_maxiter = settings$selective_backfit_maxiter
      ),
      error = function(e) e
    )
  })
  selective_elapsed_sec <- as.numeric(selective_time[["elapsed"]])

  truth_norm <- normalize_component_matrices(
    L_mat = simulation$L_true,
    F_mat = simulation$F_true
  )

  original_alignment <- align_flash_components_to_truth(
    L_est = fit$L_pm,
    F_est = fit$F_pm,
    L_true = simulation$L_true
  )
  original_norm <- normalize_component_matrices(
    L_mat = original_alignment$L,
    F_mat = original_alignment$F
  )

  screening_aligned <- align_screening_summary(
    screening_summary = screening$screening_summary,
    alignment = original_alignment,
    true_types = simulation$component_type,
    rep_id = rep_id,
    seed = seed
  )

  component_metrics <- make_component_metrics(
    rep_id = rep_id,
    seed = seed,
    stage = "original",
    truth_norm = truth_norm,
    est_norm = original_norm,
    true_types = simulation$component_type
  )

  status <- "ok"
  error_message <- NA_character_
  max_abs_nonspatial_loading_change <- NA_real_

  if (inherits(selective, "error")) {
    status <- "selective_error"
    error_message <- conditionMessage(selective)
  } else {
    selective_alignment <- align_flash_components_to_truth(
      L_est = selective$L,
      F_est = selective$F,
      L_true = simulation$L_true
    )
    selective_norm <- normalize_component_matrices(
      L_mat = selective_alignment$L,
      F_mat = selective_alignment$F
    )

    component_metrics <- rbind(
      component_metrics,
      make_component_metrics(
        rep_id = rep_id,
        seed = seed,
        stage = "selective",
        truth_norm = truth_norm,
        est_norm = selective_norm,
        true_types = simulation$component_type
      )
    )

    unchanged_idx <- screening$screening_summary$factor_k[
      screening$screening_summary$selection != "spatial" |
        is.na(screening$screening_summary$selection)
    ]
    if (length(unchanged_idx)) {
      max_abs_nonspatial_loading_change <- max(
        abs(selective$L[, unchanged_idx, drop = FALSE] - fit$L_pm[, unchanged_idx, drop = FALSE])
      )
    } else {
      max_abs_nonspatial_loading_change <- 0
    }
  }

  run_summary <- data.frame(
    rep = rep_id,
    seed = seed,
    status = status,
    error_message = error_message,
    n_estimated_components = ncol(fit$L_pm),
    n_aligned_components = min(ncol(fit$L_pm), simulation$K),
    n_screened_spatial = sum(screening$screening_summary$selection == "spatial", na.rm = TRUE),
    label_accuracy = if (nrow(screening_aligned)) mean(screening_aligned$label_correct, na.rm = TRUE) else NA_real_,
    original_spatial_loading_mse = mean_metric(component_metrics, "original", "spatial", "loading_mse"),
    selective_spatial_loading_mse = mean_metric(component_metrics, "selective", "spatial", "loading_mse"),
    original_nonspatial_loading_mse = mean_metric(component_metrics, "original", "nonspatial", "loading_mse"),
    selective_nonspatial_loading_mse = mean_metric(component_metrics, "selective", "nonspatial", "loading_mse"),
    max_abs_nonspatial_loading_change = max_abs_nonspatial_loading_change,
    flash_elapsed_sec = flash_elapsed_sec,
    selective_elapsed_sec = selective_elapsed_sec,
    stringsAsFactors = FALSE
  )

  list(
    run_summary = run_summary,
    screening = screening_raw,
    screening_aligned = screening_aligned,
    component_metrics = component_metrics
  )
}

ebnm_list_regular <- build_regular_ebnm_list()
rep_results <- lapply(seq_len(settings$n_reps), run_single_rep, settings = settings, ebnm_list_regular = ebnm_list_regular)

run_summary <- do.call(rbind, lapply(rep_results, `[[`, "run_summary"))
screening_raw <- bind_rows_safe(lapply(rep_results, `[[`, "screening"))
screening_aligned <- bind_rows_safe(lapply(rep_results, `[[`, "screening_aligned"))
component_metrics <- bind_rows_safe(lapply(rep_results, `[[`, "component_metrics"))

label_summary <- summarize_label_accuracy(
  if (nrow(screening_aligned) && "status" %in% names(screening_aligned)) {
    screening_aligned[screening_aligned$status == "ok", , drop = FALSE]
  } else {
    data.frame()
  }
)
recovery_summary <- summarize_recovery(component_metrics)

run_summary_csv <- file.path(results_dir, paste0(settings$output_prefix, "_run_summary.csv"))
screening_csv <- file.path(results_dir, paste0(settings$output_prefix, "_screening.csv"))
screening_aligned_csv <- file.path(results_dir, paste0(settings$output_prefix, "_screening_aligned.csv"))
component_metrics_csv <- file.path(results_dir, paste0(settings$output_prefix, "_component_metrics.csv"))
label_summary_csv <- file.path(results_dir, paste0(settings$output_prefix, "_label_summary.csv"))
recovery_summary_csv <- file.path(results_dir, paste0(settings$output_prefix, "_recovery_summary.csv"))
summary_md <- file.path(results_dir, paste0(settings$output_prefix, "_summary.md"))
loading_mse_png <- file.path(fig_dir, paste0(settings$output_prefix, "_loading_mse.png"))

utils::write.csv(run_summary, run_summary_csv, row.names = FALSE)
utils::write.csv(screening_raw, screening_csv, row.names = FALSE)
utils::write.csv(screening_aligned, screening_aligned_csv, row.names = FALSE)
utils::write.csv(component_metrics, component_metrics_csv, row.names = FALSE)
utils::write.csv(label_summary, label_summary_csv, row.names = FALSE)
utils::write.csv(recovery_summary, recovery_summary_csv, row.names = FALSE)

has_plot <- plot_loading_mse(component_metrics, loading_mse_png)

summary_lines <- c(
  "# Spatial Matrix-Factorization Selective Smoothing Study",
  "",
  paste0("- replicates: ", settings$n_reps),
  paste0("- grid side: ", settings$grid_side),
  paste0("- features (p): ", settings$p),
  paste0("- spatial components: ", settings$k_spatial),
  paste0("- nonspatial components: ", settings$k_nonspatial),
  paste0("- observation noise SD: ", format(settings$sigma_E, trim = TRUE)),
  paste0("- screening tolerance: ", format(settings$tol, trim = TRUE)),
  "",
  "## Run Summary",
  "",
  format_md_table(run_summary[, c(
    "rep",
    "status",
    "n_estimated_components",
    "n_screened_spatial",
    "label_accuracy",
    "original_spatial_loading_mse",
    "selective_spatial_loading_mse",
    "original_nonspatial_loading_mse",
    "selective_nonspatial_loading_mse",
    "max_abs_nonspatial_loading_change"
  )], digits = 4),
  "",
  "## Label Accuracy",
  "",
  format_md_table(label_summary, digits = 4),
  "",
  "## Recovery Summary",
  "",
  format_md_table(recovery_summary, digits = 4)
)

if (isTRUE(has_plot)) {
  summary_lines <- c(
    summary_lines,
    "",
    paste0("![Loading MSE comparison](figures/", basename(loading_mse_png), ")")
  )
}

writeLines(summary_lines, con = summary_md)

message("Completed spatial MF selective smoothing study.")
