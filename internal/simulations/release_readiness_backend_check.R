#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EBSmoothr)
})

results_dir <- file.path("internal", "simulations", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
fast_mode <- identical(Sys.getenv("EBSMOOTH_RELEASE_BENCH_FAST"), "1")

time_fit <- function(label, expr) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch(
    {
      fit <- force(expr)
      list(status = "ok", fit = fit, error = NA_character_)
    },
    error = function(e) {
      list(status = "error", fit = NULL, error = conditionMessage(e))
    }
  )
  out$elapsed_sec <- proc.time()[["elapsed"]] - start
  out$label <- label
  out
}

fit_backend <- function(x, s, family, locations = NULL, setup = NULL, link,
                        backend, beta_fixed = NULL, beta_prec = NULL,
                        pc.penalty = NULL) {
  args <- list(
    x = x,
    s = s,
    family = family,
    locations = locations,
    setup = setup,
    link = link,
    backend = backend,
    beta_fixed = beta_fixed,
    beta_prec = beta_prec,
    pc.penalty = pc.penalty
  )
  args <- args[!vapply(args, is.null, logical(1))]
  do.call(eb_smoother, args)
}

metric_row <- function(scenario, backend, result, reference = NULL) {
  fit <- result$fit
  raw <- if (!is.null(fit) && !is.null(fit$raw_fit)) fit$raw_fit else fit
  get_num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)[1]
  row <- data.frame(
    scenario = scenario,
    backend = backend,
    status = result$status,
    elapsed_sec = result$elapsed_sec,
    actual_backend = if (is.null(fit)) NA_character_ else as.character(fit$backend),
    implementation = if (is.null(raw) || is.null(raw$laplace_implementation)) NA_character_ else raw$laplace_implementation,
    validation_status = if (is.null(raw) || is.null(raw$inla_validation)) NA_character_ else raw$inla_validation$status,
    log_likelihood = get_num(if (is.null(fit)) NULL else fit$log_likelihood),
    fitted_beta = get_num(if (is.null(fit)) NULL else fit$fitted_beta),
    fitted_range = get_num(if (is.null(fit) || is.null(fit$fitted_g$theta)) NULL else exp(fit$fitted_g$theta)),
    fitted_sigma = get_num(if (is.null(fit) || is.null(fit$fitted_g$sigma)) NULL else fit$fitted_g$sigma),
    fitted_noise_sd = get_num(if (is.null(fit) || is.null(fit$fitted_noise_sd)) NULL else fit$fitted_noise_sd),
    error = if (is.na(result$error)) NA_character_ else result$error,
    stringsAsFactors = FALSE
  )

  if (!is.null(reference) && result$status == "ok" && reference$status == "ok") {
    row$delta_log_likelihood <- abs(row$log_likelihood - as.numeric(reference$fit$log_likelihood)[1])
    row$delta_beta <- abs(row$fitted_beta - as.numeric(reference$fit$fitted_beta)[1])
    row$delta_noise_sd <- abs(row$fitted_noise_sd - get_num(reference$fit$fitted_noise_sd))
    row$posterior_mean_max_abs <- max(abs(fit$posterior$mean - reference$fit$posterior$mean), na.rm = TRUE)
  } else {
    row$delta_log_likelihood <- NA_real_
    row$delta_beta <- NA_real_
    row$delta_noise_sd <- NA_real_
    row$posterior_mean_max_abs <- NA_real_
  }
  row
}

make_1d_log_data <- function(n, noise = 0.1) {
  loc <- seq(0, 1, length.out = n)
  eta <- 0.1 + 0.25 * sin(2 * pi * loc) + 0.08 * cos(6 * pi * loc)
  list(loc = loc, x = exp(eta) + rnorm(n, sd = noise), s = rep(noise, n))
}

make_1d_identity_data <- function(n, noise = 0.1) {
  loc <- seq(0, 1, length.out = n)
  mean <- 0.15 + 0.3 * sin(2 * pi * loc) + 0.1 * cos(4 * pi * loc)
  list(loc = loc, x = mean + rnorm(n, sd = noise), s = rep(noise, n))
}

make_2d_log_data <- function(k, noise = 0.1) {
  loc <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = k),
    y = seq(0, 1, length.out = k)
  ))
  eta <- 0.08 + 0.18 * sin(2 * pi * loc[, 1]) + 0.10 * cos(2 * pi * loc[, 2])
  list(loc = loc, x = exp(eta) + rnorm(nrow(loc), sd = noise), s = rep(noise, nrow(loc)))
}

set.seed(20260430)

accuracy_specs <- list()

d1_id <- make_1d_identity_data(40)
accuracy_specs[["identity_known_1d_eb"]] <- list(
  reference = "exact",
  fits = list(
    exact = function() fit_backend(d1_id$x, d1_id$s, "matern", locations = d1_id$loc, link = "identity", backend = "exact"),
    laplace_tmb = function() fit_backend(d1_id$x, d1_id$s, "matern", locations = d1_id$loc, link = "identity", backend = "laplace_tmb"),
    inla = function() fit_backend(d1_id$x, d1_id$s, "matern", locations = d1_id$loc, link = "identity", backend = "inla")
  )
)

d1_log <- make_1d_log_data(24)
accuracy_specs[["log_known_1d_eb"]] <- list(
  reference = "laplace_fisher",
  fits = list(
    auto = function() fit_backend(d1_log$x, d1_log$s, "matern", locations = d1_log$loc, link = "log", backend = "auto"),
    laplace_fisher = function() fit_backend(d1_log$x, d1_log$s, "matern", locations = d1_log$loc, link = "log", backend = "laplace_fisher"),
    laplace_tmb = function() fit_backend(d1_log$x, d1_log$s, "matern", locations = d1_log$loc, link = "log", backend = "laplace_tmb"),
    laplace_r = function() fit_backend(d1_log$x, d1_log$s, "matern", locations = d1_log$loc, link = "log", backend = "laplace_r"),
    inla = function() fit_backend(d1_log$x, d1_log$s, "matern", locations = d1_log$loc, link = "log", backend = "inla")
  )
)

d1_learn <- make_1d_log_data(24)
pc_learn <- list(range = c(0.3, 0.5), sigma = c(0.25, 0.5), noise = c(0.1, 0.5))
accuracy_specs[["log_learned_1d_eb"]] <- list(
  reference = "laplace_fisher",
  fits = list(
    auto = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "auto"),
    laplace_fisher = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_fisher"),
    laplace_tmb = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_tmb"),
    laplace_r = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_r"),
    inla = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "inla")
  )
)

accuracy_specs[["log_learned_1d_eb_pc"]] <- list(
  reference = "laplace_fisher_pc",
  fits = list(
    auto_pc = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "auto", pc.penalty = pc_learn),
    laplace_fisher_pc = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_fisher", pc.penalty = pc_learn),
    laplace_tmb_pc = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_tmb", pc.penalty = pc_learn),
    laplace_r_pc = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "laplace_r", pc.penalty = pc_learn),
    inla_pc = function() fit_backend(d1_learn$x, NULL, "matern", locations = d1_learn$loc, link = "log", backend = "inla_pc", pc.penalty = pc_learn)
  )
)

d2_log <- make_2d_log_data(4)
d2_setup <- Matern_setup(d2_log$loc, max.edge = 0.6)
accuracy_specs[["log_known_2d_fixed_beta"]] <- list(
  reference = "laplace_fisher",
  fits = list(
    auto = function() fit_backend(d2_log$x, d2_log$s, "matern", setup = d2_setup, link = "log", backend = "auto", beta_fixed = 0),
    laplace_fisher = function() fit_backend(d2_log$x, d2_log$s, "matern", setup = d2_setup, link = "log", backend = "laplace_fisher", beta_fixed = 0),
    laplace_tmb = function() fit_backend(d2_log$x, d2_log$s, "matern", setup = d2_setup, link = "log", backend = "laplace_tmb", beta_fixed = 0),
    inla = function() fit_backend(d2_log$x, d2_log$s, "matern", setup = d2_setup, link = "log", backend = "inla", beta_fixed = 0)
  )
)

accuracy_rows <- list()
for (scenario in names(accuracy_specs)) {
  spec <- accuracy_specs[[scenario]]
  results <- lapply(names(spec$fits), function(backend) {
    time_fit(paste(scenario, backend, sep = "::"), spec$fits[[backend]]())
  })
  names(results) <- names(spec$fits)
  reference <- results[[spec$reference]]
  for (backend in names(results)) {
    accuracy_rows[[length(accuracy_rows) + 1L]] <- metric_row(scenario, backend, results[[backend]], reference = reference)
  }
}
accuracy <- do.call(rbind, accuracy_rows)

runtime_specs <- list()
for (n in if (fast_mode) c(100, 300) else c(100, 300, 1000)) {
  dat_n <- make_1d_log_data(n)
  runtime_specs[[paste0("1d_log_learned_eb_n", n, "_auto")]] <- local({
    dat <- dat_n
    function() fit_backend(dat$x, NULL, "matern", locations = dat$loc, link = "log", backend = "auto")
  })
  runtime_specs[[paste0("1d_log_learned_eb_n", n, "_laplace_tmb")]] <- local({
    dat <- dat_n
    function() fit_backend(dat$x, NULL, "matern", locations = dat$loc, link = "log", backend = "laplace_tmb")
  })
  if (n == 100) {
    runtime_specs[[paste0("1d_log_learned_eb_n", n, "_laplace_r")]] <- local({
      dat <- dat_n
      function() fit_backend(dat$x, NULL, "matern", locations = dat$loc, link = "log", backend = "laplace_r")
    })
  }
}

for (k in if (fast_mode) c(16, 24) else c(16, 24, 32)) {
  dat_k <- make_2d_log_data(k)
  setup_k <- Matern_setup(dat_k$loc, max.edge = 0.08)
  runtime_specs[[paste0("2d_log_known_fixed_k", k, "_auto")]] <- local({
    dat <- dat_k
    setup <- setup_k
    function() fit_backend(dat$x, dat$s, "matern", setup = setup, link = "log", backend = "auto", beta_fixed = 0)
  })
  runtime_specs[[paste0("2d_log_known_fixed_k", k, "_laplace_tmb")]] <- local({
    dat <- dat_k
    setup <- setup_k
    function() fit_backend(dat$x, dat$s, "matern", setup = setup, link = "log", backend = "laplace_tmb", beta_fixed = 0)
  })
  runtime_specs[[paste0("2d_log_known_fixed_k", k, "_inla")]] <- local({
    dat <- dat_k
    setup <- setup_k
    function() fit_backend(dat$x, dat$s, "matern", setup = setup, link = "log", backend = "inla", beta_fixed = 0)
  })
}

runtime_rows <- list()
for (label in names(runtime_specs)) {
  result <- time_fit(label, runtime_specs[[label]]())
  runtime_rows[[length(runtime_rows) + 1L]] <- metric_row(label, sub("^.*_", "", label), result)
}
runtime <- do.call(rbind, runtime_rows)

accuracy_path <- file.path(results_dir, "backend_release_accuracy.csv")
runtime_path <- file.path(results_dir, "backend_release_runtime.csv")
summary_path <- file.path(results_dir, "backend_release_summary.md")
write.csv(accuracy, accuracy_path, row.names = FALSE)
write.csv(runtime, runtime_path, row.names = FALSE)

ok_accuracy <- subset(accuracy, status == "ok")
failed_accuracy <- subset(accuracy, status != "ok")
runtime_ok <- subset(runtime, status == "ok")

summary_lines <- c(
  "# Backend Release-Readiness Summary",
  "",
  paste("Generated:", timestamp),
  "",
  "## Accuracy",
  "",
  paste("- Successful fits:", nrow(ok_accuracy)),
  paste("- Failed or guarded fits:", nrow(failed_accuracy)),
  paste("- Maximum successful objective delta:", if (nrow(ok_accuracy)) signif(max(ok_accuracy$delta_log_likelihood, na.rm = TRUE), 4) else "NA"),
  paste("- Maximum successful posterior-mean delta:", if (nrow(ok_accuracy)) signif(max(ok_accuracy$posterior_mean_max_abs, na.rm = TRUE), 4) else "NA"),
  "",
  "Failed or guarded accuracy rows are expected only when explicit INLA cannot be validated against the package reference.",
  "",
  "## Runtime",
  "",
  paste("- Successful runtime rows:", nrow(runtime_ok)),
  paste("- Fastest successful row:", if (nrow(runtime_ok)) runtime_ok$scenario[which.min(runtime_ok$elapsed_sec)] else "NA"),
  paste("- Slowest successful row:", if (nrow(runtime_ok)) runtime_ok$scenario[which.max(runtime_ok$elapsed_sec)] else "NA"),
  "",
  "## Auto Policy",
  "",
  "- Identity-link Matern: exact.",
  "- Nonspatial: exact.",
  "- Log-link Matern learned-noise: Laplace/TMB.",
  "- Log-link Matern known-noise 2D fixed beta without PC prior: INLA when validation evidence remains consistent.",
  "- R Laplace: reference/debug implementation, not the auto default when TMB is supported.",
  "",
  paste("Accuracy CSV:", accuracy_path),
  paste("Runtime CSV:", runtime_path)
)
writeLines(summary_lines, summary_path)

cat(paste(summary_lines, collapse = "\n"))
cat("\n")
