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
  stop("The internal sensitivity script requires pkgload.")
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
                                     prior_range,
                                     prior_sigma,
                                     suppress_warnings = TRUE) {
  loc_info <- .normalize_locations(locations)
  loc_mat <- loc_info$loc
  meshA <- .build_mesh_A(loc_mat, max.edge = max.edge)
  mesh <- meshA$mesh
  A <- meshA$A
  n <- nrow(loc_mat)

  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = alpha,
    prior.range = prior_range,
    prior.sigma = prior_sigma
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
    inla_step_a = resA,
    inla_step_b = resB
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

loc_info <- .normalize_locations(locations)
meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)
spde_template <- INLA::inla.spde2.matern(meshA$mesh, alpha = alpha)

scenario_df <- data.frame(
  scenario = c(
    "weak_balanced",
    "strong_smooth_lowvar",
    "strong_rough_highvar"
  ),
  range_prob = c(0.5, 0.1, 0.9),
  sigma_prob = c(0.5, 0.1, 0.9),
  stringsAsFactors = FALSE
)

scenario_results <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_df)), function(i) {
    prior_range <- c(true_range, scenario_df$range_prob[i])
    prior_sigma <- c(true_sigma, scenario_df$sigma_prob[i])

    elapsed <- system.time(
      legacy_fit <- legacy_inla_joint_pc_fit(
        locations = locations,
        x = sim$x,
        s = sim$s,
        max.edge = max.edge,
        alpha = alpha,
        prior_range = prior_range,
        prior_sigma = prior_sigma
      )
    )[["elapsed"]]

    exact_at_legacy <- profile_exact_given_hyperparameters(
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

    data.frame(
      scenario = scenario_df$scenario[i],
      range_prob = scenario_df$range_prob[i],
      sigma_prob = scenario_df$sigma_prob[i],
      fitted_range = legacy_fit$fitted_range,
      fitted_sigma = legacy_fit$fitted_sigma,
      fitted_beta = legacy_fit$fitted_beta,
      runtime_sec = elapsed,
      raw_mlik = legacy_fit$raw_mlik,
      exact_loglik_at_legacy = exact_at_legacy$state$log_marginal,
      exact_gap_to_optimum = as.numeric(exact_fit$log_likelihood) - exact_at_legacy$state$log_marginal,
      posterior_corr_vs_exact = stats::cor(legacy_fit$posterior$mean, exact_fit$posterior$mean),
      posterior_rmse_vs_exact = sqrt(mean((legacy_fit$posterior$mean - exact_fit$posterior$mean)^2)),
      posterior_corr_vs_truth = stats::cor(legacy_fit$posterior$mean, sim$mean_surface),
      posterior_rmse_vs_truth = sqrt(mean((legacy_fit$posterior$mean - sim$mean_surface)^2))
    )
  })
)

exact_summary <- data.frame(
  method = "Exact EB",
  fitted_range = exp(exact_fit$fitted_g$theta),
  fitted_sigma = exact_fit$fitted_g$sigma,
  fitted_beta = exact_fit$fitted_beta,
  runtime_sec = exact_time,
  exact_loglik = as.numeric(exact_fit$log_likelihood),
  posterior_corr_vs_truth = stats::cor(exact_fit$posterior$mean, sim$mean_surface),
  posterior_rmse_vs_truth = sqrt(mean((exact_fit$posterior$mean - sim$mean_surface)^2))
)

plot_df <- rbind(
  data.frame(
    label = "Truth",
    type = "reference",
    range = true_range,
    sigma = true_sigma,
    stringsAsFactors = FALSE
  ),
  data.frame(
    label = "Exact EB",
    type = "exact",
    range = exp(exact_fit$fitted_g$theta),
    sigma = exact_fit$fitted_g$sigma,
    stringsAsFactors = FALSE
  ),
  data.frame(
    label = scenario_results$scenario,
    type = "legacy",
    range = scenario_results$fitted_range,
    sigma = scenario_results$fitted_sigma,
    stringsAsFactors = FALSE
  )
)

png(file.path(fig_dir, "matern_pc_prior_sensitivity.png"), width = 1500, height = 550, res = 150)
op <- par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

plot(
  plot_df$range,
  plot_df$sigma,
  type = "n",
  xlab = "Fitted range",
  ylab = "Fitted sigma",
  main = "Hyperparameter sensitivity"
)
points(
  plot_df$range[plot_df$type == "reference"],
  plot_df$sigma[plot_df$type == "reference"],
  pch = 4,
  cex = 1.5,
  lwd = 2,
  col = "firebrick"
)
points(
  plot_df$range[plot_df$type == "exact"],
  plot_df$sigma[plot_df$type == "exact"],
  pch = 19,
  cex = 1.2,
  col = "navy"
)
points(
  plot_df$range[plot_df$type == "legacy"],
  plot_df$sigma[plot_df$type == "legacy"],
  pch = 17,
  cex = 1.2,
  col = "darkgreen"
)
text(plot_df$range, plot_df$sigma, labels = plot_df$label, pos = 3, cex = 0.8)

barplot(
  height = scenario_results$exact_gap_to_optimum,
  names.arg = scenario_results$scenario,
  las = 2,
  col = "grey75",
  border = "grey25",
  ylab = "Exact objective gap",
  main = "How far legacy fits move from exact EB"
)
abline(h = 0, lwd = 1.2)

barplot(
  height = scenario_results$posterior_rmse_vs_exact,
  names.arg = scenario_results$scenario,
  las = 2,
  col = "grey75",
  border = "grey25",
  ylab = "Posterior RMSE vs exact",
  main = "Surface difference from exact EB"
)
abline(h = 0, lwd = 1.2)

par(op)
dev.off()

write.csv(exact_summary, file.path(results_dir, "matern_pc_prior_sensitivity_exact.csv"), row.names = FALSE)
write.csv(scenario_results, file.path(results_dir, "matern_pc_prior_sensitivity_legacy.csv"), row.names = FALSE)

summary_lines <- c(
  "---",
  'title: "Matern PC Prior Sensitivity Check"',
  "---",
  "",
  "# Matern PC Prior Sensitivity Check",
  "",
  paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Setup",
  "",
  "- Fixed single 2D dataset on a 20 x 16 grid (n = 320).",
  "- Truth: range = 0.30, sigma = 0.35, beta0 = 0.20, noise sd = 0.22, alpha = 2.",
  "- Exact baseline: current exact Gaussian EB Matern implementation.",
  "- Legacy variants: INLA + joint PC priors on both range and sigma.",
  "",
  "## Prior Scenario Interpretation",
  "",
  "- For the range prior, INLA uses `P(range < rho0) = alpha_r`.",
  "- For the sigma prior, INLA uses `P(sigma > sigma0) = alpha_sigma`.",
  "- Small `(alpha_r, alpha_sigma)` pushes toward smoother, lower-variance fields.",
  "- Large `(alpha_r, alpha_sigma)` pushes toward rougher, higher-variance fields.",
  "",
  "## Exact Baseline",
  "",
  format_md_table(exact_summary, digits = 6),
  "",
  "## Legacy INLA Scenarios",
  "",
  format_md_table(scenario_results, digits = 6),
  "",
  "Figure:",
  "",
  "![Matern PC prior sensitivity](figures/matern_pc_prior_sensitivity.png)",
  "",
  "## Interpretation",
  "",
  paste0(
    "- The weakest balanced prior (`0.5`, `0.5`) gives fitted `(range, sigma) = (",
    format(round(scenario_results$fitted_range[scenario_results$scenario == "weak_balanced"], 4), nsmall = 4),
    ", ",
    format(round(scenario_results$fitted_sigma[scenario_results$scenario == "weak_balanced"], 4), nsmall = 4),
    ")`."
  ),
  paste0(
    "- The smooth/low-variance prior and the rough/high-variance prior move the legacy fit by exact-objective gaps of ",
    format(round(scenario_results$exact_gap_to_optimum[scenario_results$scenario == "strong_smooth_lowvar"], 4), nsmall = 4),
    " and ",
    format(round(scenario_results$exact_gap_to_optimum[scenario_results$scenario == "strong_rough_highvar"], 4), nsmall = 4),
    ", respectively."
  ),
  "- If these gaps and posterior RMSE values stay small, then the exact EB fit and the legacy INLA fit are practically hard to distinguish on this dataset.",
  "- If they become large under strong priors, that is direct evidence that the prior is driving the legacy fit away from the exact EB objective."
)

writeLines(summary_lines, con = file.path(results_dir, "matern_pc_prior_sensitivity_summary.md"))

message("Saved Matern PC-prior sensitivity study to: ", results_dir)
