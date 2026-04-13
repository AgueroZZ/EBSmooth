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
  stop("The internal comparison script requires pkgload.")
}

load_ebsmoothr()

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
  x <- mean_surface + rnorm(length(mean_surface), sd = noise_sd)

  list(
    x = x,
    s = rep(noise_sd, length(mean_surface)),
    mean_surface = mean_surface,
    mesh = meshA$mesh,
    A = meshA$A,
    spde_template = spde_template
  )
}

legacy_inla_joint_pc_fit <- function(locations,
                                     x,
                                     s,
                                     max.edge = NULL,
                                     alpha = 2,
                                     prior_range = NULL,
                                     prior_sigma = NULL,
                                     suppress_warnings = TRUE) {
  loc_info <- .normalize_locations(locations)
  loc_mat <- loc_info$loc
  d <- loc_info$d
  meshA <- .build_mesh_A(loc_mat, max.edge = max.edge)
  mesh <- meshA$mesh
  A <- meshA$A
  n <- nrow(loc_mat)

  prior_range0 <- if (is.null(prior_range)) .default_penalty_range(loc_mat) else prior_range
  prior_sigma0 <- if (is.null(prior_sigma)) stats::sd(x) else prior_sigma
  if (!is.finite(prior_sigma0) || prior_sigma0 <= 0) prior_sigma0 <- 1

  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = alpha,
    prior.range = c(prior_range0, 0.5),
    prior.sigma = c(prior_sigma0, 0.5)
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

  stackA <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(A, matrix(1, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )

  formulaA <- Y ~ 0 + beta0 + f(spatial.field, model = spde)

  runA <- function() {
    INLA::inla(
      formulaA,
      scale = (1 / s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      control.fixed = INLA::control.fixed(prec = 0),
      data = INLA::inla.stack.data(stackA),
      control.predictor = list(A = INLA::inla.stack.A(stackA), link = 1),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE),
      silent = TRUE
    )
  }

  resA <- if (suppress_warnings) suppressWarnings(runA()) else runA()
  theta_hat <- as.numeric(resA$mode$theta)
  beta0_hat <- as.numeric(resA$summary.fixed$mean)

  stackB <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(A, matrix(beta0_hat, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )

  formulaB <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

  runB <- function() {
    INLA::inla(
      formulaB,
      scale = (1 / s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      control.fixed = INLA::control.fixed(prec = 0),
      data = INLA::inla.stack.data(stackB),
      control.predictor = list(A = INLA::inla.stack.A(stackB), link = 1),
      control.mode = INLA::control.mode(theta = theta_hat, fixed = TRUE),
      control.compute = list(mlik = TRUE),
      silent = TRUE
    )
  }

  resB <- if (suppress_warnings) suppressWarnings(runB()) else runB()
  ii <- INLA::inla.stack.index(stackB, tag = "est")$data

  list(
    posterior = data.frame(
      mean = resB$summary.fitted.values$mean[ii],
      var = resB$summary.fitted.values$sd[ii]^2
    ),
    fitted_range = unname(exp(theta_hat[1])),
    fitted_sigma = unname(exp(theta_hat[2])),
    fitted_beta = beta0_hat,
    raw_mlik = as.numeric(resB$mlik[[1]]),
    summary_hyperpar = resA$summary.hyperpar,
    inla_step_a = resA,
    inla_step_b = resB,
    mesh = mesh,
    A = A,
    d = d
  )
}

profile_exact_given_hyperparameters <- function(x,
                                                s,
                                                A,
                                                spde_template,
                                                alpha,
                                                d,
                                                range_fixed,
                                                sigma_fixed,
                                                beta_init = NULL) {
  if (is.null(beta_init)) {
    beta_init <- stats::weighted.mean(x, 1 / (s^2))
  }

  objective <- function(beta0) {
    st <- tryCatch(
      suppressWarnings(
        .exact_matern_state(
          x = x,
          s = s,
          A = A,
          spde_template = spde_template,
          alpha = alpha,
          d = d,
          log_range = log(range_fixed),
          log_sigma = log(sigma_fixed),
          beta0 = beta0
        )
      ),
      error = function(e) e
    )
    if (inherits(st, "error")) return(Inf)
    -st$log_marginal
  }

  beta_span <- max(4 * stats::sd(x), 1)
  opt <- optimize(
    f = objective,
    interval = c(beta_init - beta_span, beta_init + beta_span)
  )

  st <- suppressWarnings(
    .exact_matern_state(
      x = x,
      s = s,
      A = A,
      spde_template = spde_template,
      alpha = alpha,
      d = d,
      log_range = log(range_fixed),
      log_sigma = log(sigma_fixed),
      beta0 = opt$minimum
    )
  )

  list(
    state = st,
    fitted_sigma = sigma_fixed,
    fitted_beta = opt$minimum
  )
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

set.seed(20260412)

grid <- expand.grid(
  x = seq(0, 1, length.out = 20),
  y = seq(0, 1, length.out = 16)
)
locations <- as.matrix(grid)
alpha <- 2
true_range <- 0.30
true_sigma <- 0.35
true_beta0 <- 0.2
noise_sd <- 0.22
max.edge <- c(0.08, 0.12)

sim <- simulate_exact_matern(
  locations = locations,
  range = true_range,
  sigma = true_sigma,
  beta0 = true_beta0,
  noise_sd = noise_sd,
  alpha = alpha,
  max.edge = max.edge
)

exact_fit_fun <- ebnm_Matern_generator(locations = locations, max.edge = max.edge, alpha = alpha)
exact_time <- system.time(
  exact_fit <- exact_fit_fun(
    sim$x,
    sim$s,
    g_init = Matern(theta = log(true_range * 0.9), sigma = true_sigma * 1.1)
  )
)[["elapsed"]]

legacy_time <- system.time(
  legacy_fit <- legacy_inla_joint_pc_fit(
    locations = locations,
    x = sim$x,
    s = sim$s,
    max.edge = max.edge,
    alpha = alpha,
    prior_range = true_range,
    prior_sigma = true_sigma
  )
)[["elapsed"]]

loc_info <- .normalize_locations(locations)
meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)
spde_template <- INLA::inla.spde2.matern(meshA$mesh, alpha = alpha)
exact_profile_legacy_hyper <- profile_exact_given_hyperparameters(
  x = sim$x,
  s = sim$s,
  A = meshA$A,
  spde_template = spde_template,
  alpha = alpha,
  d = loc_info$d,
  range_fixed = legacy_fit$fitted_range,
  sigma_fixed = legacy_fit$fitted_sigma,
  beta_init = legacy_fit$fitted_beta
)

posterior_compare <- data.frame(
  metric = c(
    "corr_exact_vs_legacy",
    "rmse_exact_vs_legacy",
    "rmse_exact_vs_truth",
    "rmse_legacy_vs_truth",
    "corr_exact_vs_truth",
    "corr_legacy_vs_truth"
  ),
  value = c(
    stats::cor(exact_fit$posterior$mean, legacy_fit$posterior$mean),
    sqrt(mean((exact_fit$posterior$mean - legacy_fit$posterior$mean)^2)),
    sqrt(mean((exact_fit$posterior$mean - sim$mean_surface)^2)),
    sqrt(mean((legacy_fit$posterior$mean - sim$mean_surface)^2)),
    stats::cor(exact_fit$posterior$mean, sim$mean_surface),
    stats::cor(legacy_fit$posterior$mean, sim$mean_surface)
  )
)

hyper_compare <- data.frame(
  method = c("Exact EB", "Legacy INLA + Joint PC prior"),
  fitted_range = c(exp(exact_fit$fitted_g$theta), legacy_fit$fitted_range),
  fitted_sigma = c(exact_fit$fitted_g$sigma, legacy_fit$fitted_sigma),
  fitted_beta = c(exact_fit$fitted_beta, legacy_fit$fitted_beta)
)

runtime_compare <- data.frame(
  method = c("Exact EB", "Legacy INLA + Joint PC prior"),
  elapsed_sec = c(exact_time, legacy_time)
)

objective_compare <- data.frame(
  quantity = c(
    "Exact log marginal likelihood at exact optimum",
    "Exact log marginal likelihood at legacy fitted range/sigma (beta profiled)",
    "Gap: exact optimum minus profiled legacy-hyperparameter objective",
    "Legacy INLA raw mlik (not directly comparable)"
  ),
  value = c(
    as.numeric(exact_fit$log_likelihood),
    exact_profile_legacy_hyper$state$log_marginal,
    as.numeric(exact_fit$log_likelihood) - exact_profile_legacy_hyper$state$log_marginal,
    legacy_fit$raw_mlik
  )
)

x_vals <- sort(unique(grid$x))
y_vals <- sort(unique(grid$y))

png(file.path(fig_dir, "matern_exact_vs_legacy_inla.png"), width = 1400, height = 1100, res = 150)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 5))
image(
  x = x_vals,
  y = y_vals,
  z = matrix(sim$x, nrow = length(x_vals), ncol = length(y_vals)),
  main = "Observed noisy surface",
  xlab = "x",
  ylab = "y",
  col = hcl.colors(20, "BluYl")
)
image(
  x = x_vals,
  y = y_vals,
  z = matrix(sim$mean_surface, nrow = length(x_vals), ncol = length(y_vals)),
  main = "True latent mean",
  xlab = "x",
  ylab = "y",
  col = hcl.colors(20, "BluYl")
)
image(
  x = x_vals,
  y = y_vals,
  z = matrix(exact_fit$posterior$mean, nrow = length(x_vals), ncol = length(y_vals)),
  main = "Exact EB posterior mean",
  xlab = "x",
  ylab = "y",
  col = hcl.colors(20, "BluYl")
)
image(
  x = x_vals,
  y = y_vals,
  z = matrix(legacy_fit$posterior$mean, nrow = length(x_vals), ncol = length(y_vals)),
  main = "Legacy INLA posterior mean",
  xlab = "x",
  ylab = "y",
  col = hcl.colors(20, "BluYl")
)
par(op)
dev.off()

write.csv(hyper_compare, file.path(results_dir, "matern_exact_vs_legacy_hyperparameters.csv"), row.names = FALSE)
write.csv(runtime_compare, file.path(results_dir, "matern_exact_vs_legacy_runtime.csv"), row.names = FALSE)
write.csv(posterior_compare, file.path(results_dir, "matern_exact_vs_legacy_posterior_metrics.csv"), row.names = FALSE)
write.csv(objective_compare, file.path(results_dir, "matern_exact_vs_legacy_objectives.csv"), row.names = FALSE)

summary_lines <- c(
  "---",
  'title: "Exact Matern vs Legacy INLA PC-Prior Sanity Check"',
  "---",
  "",
  "# Exact Matern vs Legacy INLA PC-Prior Sanity Check",
  "",
  paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Setup",
  "",
  "- Single 2D simulation on a 20 x 16 grid (n = 320).",
  "- Truth: range = 0.30, sigma = 0.35, beta0 = 0.2, noise sd = 0.22, alpha = 2.",
  "- Exact method: current package implementation with exact Gaussian marginal likelihood.",
  "- Legacy method: an INLA route using `inla.spde2.pcmatern()` and INLA's internal EB optimization with joint PC priors on both range and sigma.",
  "",
  "## Important Caveat",
  "",
  "- The exact method optimizes an unpenalized exact Gaussian marginal likelihood, while the legacy INLA route uses a penalized objective induced by the PC priors and INLA's internal approximation machinery.",
  "- Because of that objective mismatch, the legacy `mlik` reported by INLA is still not directly comparable to the exact method's raw marginal log-likelihood, even though both methods now fit range and sigma.",
  "",
  "## Hyperparameters",
  "",
  format_md_table(hyper_compare, digits = 4),
  "",
  "## Runtime",
  "",
  format_md_table(runtime_compare, digits = 4),
  "",
  "## Posterior Surface Comparison",
  "",
  format_md_table(posterior_compare, digits = 6),
  "",
  "Figure:",
  "",
  "![Exact vs legacy Matern fit](figures/matern_exact_vs_legacy_inla.png)",
  "",
  "## Objective Comparison",
  "",
  format_md_table(objective_compare, digits = 6),
  "",
  "Interpretation:",
  "",
  paste0(
    "- The posterior correlation between the two fitted surfaces is ",
    format(round(posterior_compare$value[posterior_compare$metric == "corr_exact_vs_legacy"], 6), nsmall = 6),
    "."
  ),
  paste0(
    "- The exact objective gap between the exact optimum and the best exact-model fit constrained to the legacy fitted range/sigma pair is ",
    format(round(objective_compare$value[3], 6), nsmall = 6),
    "."
  ),
  "- This exact-objective gap is the more meaningful same-scale comparison, because the raw INLA `mlik` includes a different model/objective configuration."
)

writeLines(summary_lines, con = file.path(results_dir, "matern_exact_vs_legacy_summary.md"))

message("Saved exact-vs-legacy Matern comparison to: ", results_dir)
