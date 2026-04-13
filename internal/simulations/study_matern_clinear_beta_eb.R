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
  stop("This internal spike requires pkgload.")
}

load_ebsmoothr()

fit_matern_clinear_pc <- function(locations,
                                  x,
                                  s,
                                  pc_penalty,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE) {
  loc_info <- .normalize_locations(locations)
  meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)
  mesh <- meshA$mesh
  A <- meshA$A
  d <- loc_info$d

  spde_exact <- INLA::inla.spde2.matern(mesh = mesh, alpha = alpha)
  spde_pc <- .build_matern_pc_spde(
    mesh = mesh,
    alpha = alpha,
    pc_penalty = pc_penalty,
    suppress_warnings = suppress_warnings
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde_pc$n.spde)
  n <- length(x)

  stack <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(A, 1),
    effects = list(
      spatial.field = idx$spatial.field,
      data.frame(beta_cov = rep(1, n))
    ),
    tag = "est"
  )

  formula <- Y ~ -1 +
    f(
      beta_cov,
      model = "clinear",
      hyper = list(
        beta = list(
          initial = stats::weighted.mean(x, 1 / (s^2)),
          fixed = FALSE,
          prior = "flat",
          param = numeric(0)
        )
      )
    ) +
    f(spatial.field, model = spde_pc)

  run_fit <- function() {
    INLA::inla(
      formula,
      scale = (1 / s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      data = INLA::inla.stack.data(stack),
      control.predictor = list(A = INLA::inla.stack.A(stack), link = 1, compute = TRUE),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  }

  res <- if (suppress_warnings) suppressWarnings(run_fit()) else run_fit()

  theta <- res$mode$theta
  beta_mode <- as.numeric(theta[grep("Beta_intern", names(theta))[1]])
  log_range_mode <- as.numeric(theta[grep("log\\(Range\\)", names(theta))[1]])
  log_sigma_mode <- as.numeric(theta[grep("log\\(Stdev\\)", names(theta))[1]])

  exact_objective <- function(par) {
    state <- .exact_matern_state(
      x = as.numeric(x),
      s = as.numeric(s),
      A = A,
      spde_template = spde_exact,
      alpha = alpha,
      d = d,
      log_range = par[1],
      log_sigma = par[2],
      beta0 = par[3]
    )

    state$log_marginal +
      .log_pc_prior_matern_internal(
        log_range = par[1],
        log_sigma = par[2],
        range_spec = pc_penalty$range,
        sigma_spec = pc_penalty$sigma,
        d = d
      )
  }

  opt <- stats::optim(
    par = c(log_range_mode, log_sigma_mode, beta_mode),
    fn = function(par) {
      out <- tryCatch(exact_objective(par), error = function(e) NA_real_)
      if (!is.finite(out)) return(Inf)
      -out
    },
    method = "BFGS"
  )

  exact_at_mode <- exact_objective(c(log_range_mode, log_sigma_mode, beta_mode))
  exact_opt <- -opt$value

  list(
    result = res,
    d = d,
    exact_at_mode = exact_at_mode,
    exact_opt = exact_opt,
    beta_mode = beta_mode,
    beta_exact = opt$par[3],
    range_mode = exp(log_range_mode),
    range_exact = exp(opt$par[1]),
    sigma_mode = exp(log_sigma_mode),
    sigma_exact = exp(opt$par[2]),
    stepA_penalized = as.numeric(res$misc$log.posterior.mode),
    stepA_mlik_integration = as.numeric(res$mlik["log marginal-likelihood (integration)", 1]),
    stepA_mlik_gaussian = as.numeric(res$mlik["log marginal-likelihood (Gaussian)", 1])
  )
}

format_md_table <- function(df, digits = 6) {
  df_fmt <- df
  for (j in seq_along(df_fmt)) {
    if (is.numeric(df_fmt[[j]])) {
      df_fmt[[j]] <- format(round(df_fmt[[j]], digits), nsmall = digits, trim = TRUE)
    }
  }

  header <- paste(names(df_fmt), collapse = " | ")
  sep <- paste(rep("---", ncol(df_fmt)), collapse = " | ")
  rows <- apply(df_fmt, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, sep, rows), collapse = "\n")
}

set.seed(9101)
locations_1d <- seq(0, 10, length.out = 60)
s_1d <- rep(0.22, length(locations_1d))
x_1d <- 0.3 + sin(locations_1d / 1.8) + stats::rnorm(length(locations_1d), sd = s_1d)
pc_1d <- list(
  range = c(anchor = 1.0, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
fit_1d <- fit_matern_clinear_pc(
  locations = locations_1d,
  x = x_1d,
  s = s_1d,
  pc_penalty = pc_1d
)

set.seed(9102)
grid_2d <- expand.grid(
  x = seq(0, 1, length.out = 7),
  y = seq(0, 1, length.out = 6)
)
locations_2d <- as.matrix(grid_2d)
s_2d <- rep(0.2, nrow(locations_2d))
x_2d <- with(grid_2d, sin(2 * pi * x) * cos(2 * pi * y)) +
  stats::rnorm(nrow(locations_2d), sd = s_2d)
pc_2d <- list(
  range = c(anchor = 0.25, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
fit_2d <- fit_matern_clinear_pc(
  locations = locations_2d,
  x = x_2d,
  s = s_2d,
  pc_penalty = pc_2d,
  max.edge = c(0.25, 0.4)
)

summary_df <- data.frame(
  case = c("1d", "2d"),
  stepA_penalized = c(fit_1d$stepA_penalized, fit_2d$stepA_penalized),
  exact_profile_plus_prior_at_mode = c(fit_1d$exact_at_mode, fit_2d$exact_at_mode),
  exact_profile_plus_prior_optimum = c(fit_1d$exact_opt, fit_2d$exact_opt),
  gap_mode = c(
    fit_1d$stepA_penalized - fit_1d$exact_at_mode,
    fit_2d$stepA_penalized - fit_2d$exact_at_mode
  ),
  gap_optimum = c(
    fit_1d$stepA_penalized - fit_1d$exact_opt,
    fit_2d$stepA_penalized - fit_2d$exact_opt
  ),
  beta_clinear = c(fit_1d$beta_mode, fit_2d$beta_mode),
  beta_exact = c(fit_1d$beta_exact, fit_2d$beta_exact),
  range_clinear = c(fit_1d$range_mode, fit_2d$range_mode),
  range_exact = c(fit_1d$range_exact, fit_2d$range_exact),
  sigma_clinear = c(fit_1d$sigma_mode, fit_2d$sigma_mode),
  sigma_exact = c(fit_1d$sigma_exact, fit_2d$sigma_exact),
  stepA_mlik_integration = c(fit_1d$stepA_mlik_integration, fit_2d$stepA_mlik_integration),
  stepA_mlik_gaussian = c(fit_1d$stepA_mlik_gaussian, fit_2d$stepA_mlik_gaussian)
)

csv_path <- file.path(results_dir, "matern_clinear_beta_eb_summary.csv")
utils::write.csv(summary_df, csv_path, row.names = FALSE)

md_path <- file.path(results_dir, "matern_clinear_beta_eb_summary.md")
md_lines <- c(
  "# Matern `clinear` Beta-EB Spike",
  "",
  "This internal experiment checks whether an INLA `pcmatern` fit with",
  "`clinear` used for the intercept can reproduce the exact profiled objective",
  "used by the current exact Gaussian Matern implementation.",
  "",
  "## Setup",
  "",
  "- `range` and `sigma` use PC priors via `inla.spde2.pcmatern()`.",
  "- The intercept is represented through `f(beta_cov, model = \"clinear\")`",
  "  with `beta_cov = 1` for every observation.",
  "- The `clinear` hyperparameter uses a flat prior, so this spike targets the",
  "  `betaprec < 0` interpretation in which `beta0` is optimized rather than",
  "  integrated out.",
  "",
  "## Results",
  "",
  format_md_table(summary_df, digits = 6),
  "",
  "## Takeaway",
  "",
  "With the correct `clinear` specification, the INLA Step A penalized",
  "objective is numerically almost identical to the exact profiled objective",
  "plus the PC prior contribution. The fitted `(beta0, range, sigma)` values",
  "are also nearly identical to the exact optimum in both the 1D and 2D",
  "examples.",
  "",
  "This supports the idea that `betaprec < 0` may be implementable within the",
  "fast INLA framework by representing the intercept as a `clinear`",
  "hyperparameter instead of a fixed effect."
)
writeLines(md_lines, md_path)

message("Saved summary CSV to: ", csv_path)
message("Saved summary markdown to: ", md_path)
