args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (!length(file_arg)) {
  stop("Could not determine the script path from commandArgs().")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."))
results_dir <- file.path(script_dir, "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

load_ebsmoothr <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(
      file.path(project_root, "EBSmoothr"),
      quiet = TRUE,
      export_all = TRUE
    )
    return(invisible(TRUE))
  }
  stop("The internal objective-inspection script requires pkgload.")
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

run_stepA_stepB <- function(locations,
                            x,
                            s,
                            prior_range,
                            prior_sigma,
                            alpha = 2,
                            max.edge = c(0.08, 0.12)) {
  loc_info <- .normalize_locations(locations)
  meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)
  spde_template <- INLA::inla.spde2.matern(meshA$mesh, alpha = alpha)

  spde <- INLA::inla.spde2.pcmatern(
    mesh = meshA$mesh,
    alpha = alpha,
    prior.range = prior_range,
    prior.sigma = prior_sigma
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)
  n <- nrow(loc_info$loc)

  stackA <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(meshA$A, matrix(1, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )
  formulaA <- Y ~ 0 + beta0 + f(spatial.field, model = spde)

  resA <- suppressWarnings(
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
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  )

  theta_hat <- as.numeric(resA$mode$theta)
  beta0_hat <- as.numeric(resA$summary.fixed$mean)

  exact_at_A <- .exact_matern_state(
    x = x,
    s = s,
    A = meshA$A,
    spde_template = spde_template,
    alpha = alpha,
    d = loc_info$d,
    log_range = theta_hat[1],
    log_sigma = theta_hat[2],
    beta0 = beta0_hat
  )

  stackB <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(meshA$A, matrix(beta0_hat, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )
  formulaB <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

  resB <- suppressWarnings(
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
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  )

  list(
    resA = resA,
    resB = resB,
    theta_hat = theta_hat,
    beta0_hat = beta0_hat,
    exact_at_A = exact_at_A,
    spde_template = spde_template,
    A = meshA$A,
    d = loc_info$d
  )
}

run_fixed_theta_stepB <- function(locations,
                                  x,
                                  s,
                                  prior_range,
                                  prior_sigma,
                                  theta_fixed,
                                  beta_fixed,
                                  alpha = 2,
                                  max.edge = c(0.08, 0.12)) {
  loc_info <- .normalize_locations(locations)
  meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)

  spde <- INLA::inla.spde2.pcmatern(
    mesh = meshA$mesh,
    alpha = alpha,
    prior.range = prior_range,
    prior.sigma = prior_sigma
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde$n.spde)
  n <- nrow(loc_info$loc)

  stackB <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(meshA$A, matrix(beta_fixed, nrow = n, ncol = 1)),
    effects = list(
      spatial.field = idx$spatial.field,
      beta0 = 1
    ),
    tag = "est"
  )
  formulaB <- Y ~ 0 + offset(beta0) + f(spatial.field, model = spde)

  resB <- suppressWarnings(
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
      control.mode = INLA::control.mode(theta = theta_fixed, fixed = TRUE),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  )

  resB
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

scenario_df <- data.frame(
  scenario = c("balanced", "smooth_lowvar", "rough_highvar"),
  range_prob = c(0.5, 0.1, 0.9),
  sigma_prob = c(0.5, 0.1, 0.9),
  stringsAsFactors = FALSE
)

step_table <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_df)), function(i) {
    fit <- run_stepA_stepB(
      locations = locations,
      x = sim$x,
      s = sim$s,
      prior_range = c(true_range, scenario_df$range_prob[i]),
      prior_sigma = c(true_sigma, scenario_df$sigma_prob[i]),
      alpha = alpha,
      max.edge = max.edge
    )

    data.frame(
      scenario = scenario_df$scenario[i],
      range_prob = scenario_df$range_prob[i],
      sigma_prob = scenario_df$sigma_prob[i],
      fitted_range = exp(fit$theta_hat[1]),
      fitted_sigma = exp(fit$theta_hat[2]),
      fitted_beta = fit$beta0_hat,
      stepA_mlik_integration = as.numeric(fit$resA$mlik["log marginal-likelihood (integration)", 1]),
      stepA_mlik_gaussian = as.numeric(fit$resA$mlik["log marginal-likelihood (Gaussian)", 1]),
      stepA_log_posterior_mode = as.numeric(fit$resA$misc$log.posterior.mode),
      stepA_max_log_posterior = as.numeric(fit$resA$misc$configs$max.log.posterior),
      stepA_joint_log_posterior = as.numeric(fit$resA$joint.hyper[, "Log posterior density"]),
      exact_loglik_at_stepA_mode = as.numeric(fit$exact_at_A$log_marginal),
      stepB_mlik_integration = as.numeric(fit$resB$mlik["log marginal-likelihood (integration)", 1]),
      stepB_mlik_gaussian = as.numeric(fit$resB$mlik["log marginal-likelihood (Gaussian)", 1])
    )
  })
)

theta_fixed <- c(log(step_table$fitted_range[1]), log(step_table$fitted_sigma[1]))
beta_fixed <- step_table$fitted_beta[1]

fixed_theta_table <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_df)), function(i) {
    resB <- run_fixed_theta_stepB(
      locations = locations,
      x = sim$x,
      s = sim$s,
      prior_range = c(true_range, scenario_df$range_prob[i]),
      prior_sigma = c(true_sigma, scenario_df$sigma_prob[i]),
      theta_fixed = theta_fixed,
      beta_fixed = beta_fixed,
      alpha = alpha,
      max.edge = max.edge
    )

    data.frame(
      scenario = scenario_df$scenario[i],
      range_prob = scenario_df$range_prob[i],
      sigma_prob = scenario_df$sigma_prob[i],
      stepB_fixed_theta_mlik_integration = as.numeric(resB$mlik["log marginal-likelihood (integration)", 1]),
      stepB_fixed_theta_mlik_gaussian = as.numeric(resB$mlik["log marginal-likelihood (Gaussian)", 1]),
      stepB_fixed_theta_log_posterior_mode = as.numeric(resB$misc$log.posterior.mode),
      stepB_fixed_theta_max_log_posterior = as.numeric(resB$misc$configs$max.log.posterior)
    )
  })
)

write.csv(step_table, file.path(results_dir, "matern_inla_stepA_stepB_objectives.csv"), row.names = FALSE)
write.csv(fixed_theta_table, file.path(results_dir, "matern_inla_stepB_fixed_theta_objectives.csv"), row.names = FALSE)

summary_lines <- c(
  "---",
  'title: "Inspect INLA Step-A and Step-B Objective Quantities"',
  "---",
  "",
  "# Inspect INLA Step-A and Step-B Objective Quantities",
  "",
  paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Setup",
  "",
  "- Fixed single 2D dataset on a 20 x 16 grid (n = 320).",
  "- Truth: range = 0.30, sigma = 0.35, beta0 = 0.20, noise sd = 0.22, alpha = 2.",
  "- Three legacy INLA scenarios are compared by changing the PC-prior tail probabilities only.",
  "",
  "## Step A vs Step B",
  "",
  format_md_table(step_table, digits = 6),
  "",
  "## Step B With Fixed Theta/Beta Across Different Priors",
  "",
  format_md_table(fixed_theta_table, digits = 12),
  "",
  "## Interpretation",
  "",
  "- `stepB_mlik_*` matches the exact Gaussian marginal likelihood evaluated at the fitted hyperparameters.",
  "- `stepB_mlik_*` and `stepB_fixed_theta_*` are invariant to the PC-prior choice once theta is fixed.",
  "- The Step A quantities that change with the prior are `stepA_mlik_*`, `stepA_log_posterior_mode`, `stepA_max_log_posterior`, and `stepA_joint_log_posterior`.",
  "- Empirically, `stepA_log_posterior_mode` and `stepA_max_log_posterior` behave like penalized optimization targets, whereas Step B does not retain a prior-dependent criterion once theta is fixed."
)

writeLines(summary_lines, con = file.path(results_dir, "matern_inla_objectives_summary.md"))

message("Saved INLA objective inspection results to: ", results_dir)
