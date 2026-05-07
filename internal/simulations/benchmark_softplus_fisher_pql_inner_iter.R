#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EBSmoothr)
})

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all("EBSmoothr", export_all = TRUE, quiet = TRUE)
}

timed_fit <- function(expr) {
  gc()
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = as.numeric(elapsed))
}

rmse <- function(x) sqrt(mean(as.numeric(x)^2))

softplus <- function(x) log1p(exp(-abs(x))) + pmax(x, 0)

extract_matern_summary <- function(fit) {
  raw <- if (!is.null(fit$raw_fit)) fit$raw_fit else fit
  if (is.null(raw$posterior_eta)) {
    stop("The fit object does not contain `posterior_eta`.")
  }
  list(
    log_range = as.numeric(raw$fitted_g$theta),
    log_sigma = log(as.numeric(raw$fitted_g$sigma)),
    beta = as.numeric(raw$fitted_beta),
    noise_sd = if (is.null(raw$fitted_noise_sd)) NA_real_ else as.numeric(raw$fitted_noise_sd),
    log_likelihood = as.numeric(raw$log_likelihood),
    eta_mean = as.numeric(raw$posterior_eta$mean),
    eta_sd = as.numeric(raw$posterior_eta$sd),
    response_mean = as.numeric(raw$posterior$mean),
    response_sd = sqrt(pmax(as.numeric(raw$posterior$var), 0))
  )
}

compare_to_reference <- function(case, pql_inner_iter, pql_fit, pql_time, ref_fit, ref_time, n) {
  pql <- extract_matern_summary(pql_fit)
  ref <- extract_matern_summary(ref_fit)
  data.frame(
    case = case,
    pql_inner_iter = as.integer(pql_inner_iter),
    pql_time_sec = pql_time,
    reference_time_sec = ref_time,
    speedup = ref_time / pql_time,
    log_range_pql = pql$log_range,
    log_range_ref = ref$log_range,
    log_range_diff = pql$log_range - ref$log_range,
    log_sigma_pql = pql$log_sigma,
    log_sigma_ref = ref$log_sigma,
    log_sigma_diff = pql$log_sigma - ref$log_sigma,
    beta_pql = pql$beta,
    beta_ref = ref$beta,
    beta_diff = pql$beta - ref$beta,
    noise_sd_pql = pql$noise_sd,
    noise_sd_ref = ref$noise_sd,
    noise_sd_diff = pql$noise_sd - ref$noise_sd,
    log_noise_diff = if (is.finite(pql$noise_sd) && is.finite(ref$noise_sd)) {
      log(pql$noise_sd) - log(ref$noise_sd)
    } else {
      NA_real_
    },
    eta_mean_max_abs = max(abs(pql$eta_mean - ref$eta_mean)),
    eta_mean_rmse = rmse(pql$eta_mean - ref$eta_mean),
    eta_sd_max_abs = max(abs(pql$eta_sd - ref$eta_sd)),
    eta_sd_rmse = rmse(pql$eta_sd - ref$eta_sd),
    response_mean_max_abs = max(abs(pql$response_mean - ref$response_mean)),
    response_mean_rmse = rmse(pql$response_mean - ref$response_mean),
    log_likelihood_diff_per_obs = (pql$log_likelihood - ref$log_likelihood) / n,
    stringsAsFactors = FALSE
  )
}

fit_fisher_pql_known <- function(x, s, setup, pql_inner_iter) {
  resolved <- EBSmoothr:::.resolve_matern_g_init(
    x = x,
    s = s,
    g_init = NULL,
    beta_fixed = NULL,
    beta_prec = NULL,
    penalty_range0 = setup$penalty_range,
    pc.penalty = NULL,
    allow_noise = FALSE,
    link = "softplus"
  )
  EBSmoothr:::.fit_matern_fisher_pql_known_noise(
    x = x,
    s = s,
    A = setup$A,
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    theta_init = resolved$theta_init,
    sigma_init = resolved$sigma_init,
    beta_init = resolved$beta_init,
    beta_mode = "empirical_bayes",
    pc_penalty = resolved$pc_penalty,
    link = "softplus",
    suppress_warnings = TRUE,
    pql_max_iter = pql_inner_iter
  )
}

fit_fisher_pql_unknown <- function(x, setup, pql_inner_iter) {
  resolved <- EBSmoothr:::.resolve_matern_g_init(
    x = x,
    s = NULL,
    g_init = NULL,
    beta_fixed = NULL,
    beta_prec = NULL,
    penalty_range0 = setup$penalty_range,
    pc.penalty = NULL,
    allow_noise = TRUE,
    link = "softplus"
  )
  EBSmoothr:::.fit_matern_fisher_pql_unknown_noise(
    x = x,
    A = setup$A,
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    theta_init = resolved$theta_init,
    sigma_init = resolved$sigma_init,
    noise_sd_init = resolved$noise_sd_init,
    beta_init = resolved$beta_init,
    beta_mode = "empirical_bayes",
    pc_penalty = resolved$pc_penalty,
    link = "softplus",
    suppress_warnings = TRUE,
    pql_max_iter = pql_inner_iter
  )
}

set.seed(20260507)
n <- 1000L
pql_grid <- c(1L, 2L, 3L, 5L, 10L)
noise_sd_true <- 0.15
locations <- cbind(runif(n), runif(n))
eta_true <- 0.45 +
  0.85 * sin(2 * pi * locations[, 1]) +
  0.55 * cos(2 * pi * locations[, 2]) +
  0.35 * sin(2 * pi * (locations[, 1] + locations[, 2]))
x <- softplus(eta_true) + stats::rnorm(n, sd = noise_sd_true)
s_known <- rep(noise_sd_true, n)
setup <- Matern_setup(locations, max.edge = c(0.08, 0.24), suppress_warnings = TRUE)

message("Fitting softplus known-noise Laplace reference...")
known_ref <- timed_fit(
  ebnm_Matern_generator(setup = setup, link = "softplus", backend = "laplace")(x, s_known)
)

message("Fitting softplus learned-noise Laplace reference...")
unknown_ref <- timed_fit(
  eb_smoother(
    x,
    s = NULL,
    family = "matern",
    setup = setup,
    backend = "laplace",
    link = "softplus",
    suppress_warnings = TRUE
  )
)

rows <- list()
row_id <- 1L
for (pql_inner_iter in pql_grid) {
  message("Fitting known-noise Fisher-PQL, inner_iter = ", pql_inner_iter, "...")
  pql_known <- timed_fit(fit_fisher_pql_known(x, s_known, setup, pql_inner_iter))
  rows[[row_id]] <- compare_to_reference(
    case = "known_s_ebnm",
    pql_inner_iter = pql_inner_iter,
    pql_fit = pql_known$value,
    pql_time = pql_known$elapsed,
    ref_fit = known_ref$value,
    ref_time = known_ref$elapsed,
    n = n
  )
  row_id <- row_id + 1L

  message("Fitting learned-noise Fisher-PQL, inner_iter = ", pql_inner_iter, "...")
  pql_unknown <- timed_fit(fit_fisher_pql_unknown(x, setup, pql_inner_iter))
  rows[[row_id]] <- compare_to_reference(
    case = "learned_s_eb_smoother",
    pql_inner_iter = pql_inner_iter,
    pql_fit = pql_unknown$value,
    pql_time = pql_unknown$elapsed,
    ref_fit = unknown_ref$value,
    ref_time = unknown_ref$elapsed,
    n = n
  )
  row_id <- row_id + 1L
}

results <- do.call(rbind, rows)
results_dir <- file.path("internal", "simulations", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(results_dir, "softplus_fisher_pql_inner_iter_n1000.csv")
utils::write.csv(results, out_file, row.names = FALSE)

print(results[, c(
  "case", "pql_inner_iter", "pql_time_sec", "reference_time_sec", "speedup",
  "log_range_diff", "log_sigma_diff", "beta_diff", "noise_sd_diff",
  "eta_mean_rmse", "eta_sd_rmse", "log_likelihood_diff_per_obs"
)], row.names = FALSE)
message("Wrote results to ", out_file)
