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
  stop("This internal benchmark requires pkgload.")
}

load_ebsmoothr()

build_benchmark_setup <- function(locations,
                                  pc_penalty,
                                  max.edge = NULL,
                                  alpha = 2,
                                  suppress_warnings = TRUE) {
  loc_info <- .normalize_locations(locations)
  meshA <- .build_mesh_A(loc_info$loc, max.edge = max.edge)

  spde_exact <- if (suppress_warnings) {
    suppressWarnings(INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha))
  } else {
    INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha)
  }

  spde_pc <- .build_matern_pc_spde(
    mesh = meshA$mesh,
    alpha = alpha,
    pc_penalty = pc_penalty,
    suppress_warnings = suppress_warnings
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde_pc$n.spde)

  list(
    d = loc_info$d,
    A = meshA$A,
    mesh = meshA$mesh,
    spde_exact = spde_exact,
    spde_pc = spde_pc,
    idx = idx,
    pc_penalty = pc_penalty,
    alpha = alpha,
    max.edge = max.edge,
    n = nrow(loc_info$loc)
  )
}

fit_exact_profile_pc <- function(setup, x, s) {
  d <- setup$d
  eval_count <- 0L

  objective <- function(par) {
    eval_count <<- eval_count + 1L
    state <- tryCatch(
      .exact_matern_state(
        x = as.numeric(x),
        s = as.numeric(s),
        A = setup$A,
        spde_template = setup$spde_exact,
        alpha = setup$alpha,
        d = d,
        log_range = par[1],
        log_sigma = par[2],
        beta0 = par[3]
      ),
      error = function(e) e
    )
    if (inherits(state, "error")) return(Inf)

    -(
      state$log_marginal +
        .log_pc_prior_matern_internal(
          log_range = par[1],
          log_sigma = par[2],
          range_spec = setup$pc_penalty$range,
          sigma_spec = setup$pc_penalty$sigma,
          d = d
        )
    )
  }

  par0 <- c(
    log(as.numeric(setup$pc_penalty$range["anchor"])),
    log(as.numeric(setup$pc_penalty$sigma["anchor"])),
    stats::weighted.mean(x, 1 / (s^2))
  )

  opt <- stats::optim(
    par = par0,
    fn = objective,
    method = "BFGS"
  )

  list(
    log_likelihood = -opt$value,
    beta = opt$par[3],
    range = exp(opt$par[1]),
    sigma = exp(opt$par[2]),
    eval_count = eval_count
  )
}

fit_clinear_pc <- function(setup, x, s, suppress_warnings = TRUE) {
  n <- setup$n

  stack <- INLA::inla.stack(
    data = list(Y = as.numeric(x)),
    A = list(setup$A, 1),
    effects = list(
      spatial.field = setup$idx$spatial.field,
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
    f(spatial.field, model = setup$spde_pc)

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
  list(
    log_likelihood = as.numeric(res$misc$log.posterior.mode),
    beta = as.numeric(theta[grep("Beta_intern", names(theta))[1]]),
    range = exp(as.numeric(theta[grep("log\\(Range\\)", names(theta))[1]])),
    sigma = exp(as.numeric(theta[grep("log\\(Stdev\\)", names(theta))[1]]))
  )
}

time_one <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(value = value, elapsed = elapsed)
}

benchmark_case <- function(case_name,
                           locations,
                           x,
                           s,
                           pc_penalty,
                           max.edge = NULL,
                           reps = 3L) {
  setup <- build_benchmark_setup(
    locations = locations,
    pc_penalty = pc_penalty,
    max.edge = max.edge
  )

  invisible(fit_exact_profile_pc(setup, x, s))
  invisible(fit_clinear_pc(setup, x, s))

  exact_runs <- vector("list", reps)
  clinear_runs <- vector("list", reps)

  for (i in seq_len(reps)) {
    exact_runs[[i]] <- time_one(fit_exact_profile_pc(setup, x, s))
    clinear_runs[[i]] <- time_one(fit_clinear_pc(setup, x, s))
  }

  exact_elapsed <- vapply(exact_runs, function(z) z$elapsed, numeric(1))
  clinear_elapsed <- vapply(clinear_runs, function(z) z$elapsed, numeric(1))

  exact_last <- exact_runs[[reps]]$value
  clinear_last <- clinear_runs[[reps]]$value

  list(
    summary = data.frame(
      case = case_name,
      n = length(x),
      reps = reps,
      exact_mean_seconds = mean(exact_elapsed),
      exact_sd_seconds = stats::sd(exact_elapsed),
      clinear_mean_seconds = mean(clinear_elapsed),
      clinear_sd_seconds = stats::sd(clinear_elapsed),
      speedup_exact_over_clinear = mean(exact_elapsed) / mean(clinear_elapsed),
      exact_eval_count = exact_last$eval_count,
      objective_gap = clinear_last$log_likelihood - exact_last$log_likelihood,
      beta_gap = clinear_last$beta - exact_last$beta,
      range_gap = clinear_last$range - exact_last$range,
      sigma_gap = clinear_last$sigma - exact_last$sigma
    ),
    timings = data.frame(
      case = case_name,
      method = rep(c("exact_profile", "clinear_inla"), each = reps),
      rep = rep(seq_len(reps), times = 2L),
      elapsed_seconds = c(exact_elapsed, clinear_elapsed)
    )
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
  rows <- apply(df_fmt, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, sep, rows), collapse = "\n")
}

set.seed(9201)
locations_1d <- seq(0, 12, length.out = 120)
s_1d <- rep(0.22, length(locations_1d))
x_1d <- 0.3 + sin(locations_1d / 1.8) + stats::rnorm(length(locations_1d), sd = s_1d)
pc_1d <- list(
  range = c(anchor = 1.0, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
bench_1d <- benchmark_case(
  case_name = "1d_n120",
  locations = locations_1d,
  x = x_1d,
  s = s_1d,
  pc_penalty = pc_1d,
  reps = 3L
)

set.seed(9202)
grid_2d <- expand.grid(
  x = seq(0, 1, length.out = 18),
  y = seq(0, 1, length.out = 14)
)
locations_2d <- as.matrix(grid_2d)
s_2d <- rep(0.2, nrow(locations_2d))
x_2d <- with(grid_2d, sin(2 * pi * x) * cos(2 * pi * y)) +
  stats::rnorm(nrow(locations_2d), sd = s_2d)
pc_2d <- list(
  range = c(anchor = 0.25, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
bench_2d <- benchmark_case(
  case_name = "2d_n252",
  locations = locations_2d,
  x = x_2d,
  s = s_2d,
  pc_penalty = pc_2d,
  max.edge = c(0.12, 0.18),
  reps = 2L
)

set.seed(9203)
grid_2d_big <- expand.grid(
  x = seq(0, 1, length.out = 24),
  y = seq(0, 1, length.out = 18)
)
locations_2d_big <- as.matrix(grid_2d_big)
s_2d_big <- rep(0.2, nrow(locations_2d_big))
x_2d_big <- with(grid_2d_big, sin(2 * pi * x) * cos(2 * pi * y)) +
  stats::rnorm(nrow(locations_2d_big), sd = s_2d_big)
pc_2d_big <- list(
  range = c(anchor = 0.25, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
bench_2d_big <- benchmark_case(
  case_name = "2d_n432",
  locations = locations_2d_big,
  x = x_2d_big,
  s = s_2d_big,
  pc_penalty = pc_2d_big,
  max.edge = c(0.10, 0.15),
  reps = 1L
)

summary_df <- rbind(bench_1d$summary, bench_2d$summary, bench_2d_big$summary)
timings_df <- rbind(bench_1d$timings, bench_2d$timings, bench_2d_big$timings)

summary_csv <- file.path(results_dir, "matern_clinear_vs_exact_betaeb_benchmark_summary.csv")
timings_csv <- file.path(results_dir, "matern_clinear_vs_exact_betaeb_benchmark_timings.csv")

utils::write.csv(summary_df, summary_csv, row.names = FALSE)
utils::write.csv(timings_df, timings_csv, row.names = FALSE)

summary_md <- file.path(results_dir, "matern_clinear_vs_exact_betaeb_benchmark_summary.md")
md_lines <- c(
  "# Matern `clinear` vs Exact Beta-EB Benchmark",
  "",
  "This benchmark compares two implementations of the `betaprec < 0`",
  "interpretation for the Matern smoother:",
  "",
  "- `exact_profile`: exact Gaussian profile objective plus the PC prior;",
  "- `clinear_inla`: INLA Step A with `clinear` used to represent the",
  "  intercept as a hyperparameter.",
  "",
  "The timings below are fit-phase timings after the spatial setup has already",
  "been built. This mirrors the package pattern in which the generator is built",
  "once for a fixed set of locations and then reused for fitting.",
  "",
  "## Summary",
  "",
  format_md_table(summary_df, digits = 6),
  "",
  "## Interpretation",
  "",
  "- `speedup_exact_over_clinear > 1` means the exact/profile method is slower",
  "  than the `clinear`-based INLA fit by that factor.",
  "- The parameter and objective gaps are included to confirm that both methods",
  "  still converge to essentially the same fitted solution."
)
writeLines(md_lines, summary_md)

message("Saved summary CSV to: ", summary_csv)
message("Saved timings CSV to: ", timings_csv)
message("Saved summary markdown to: ", summary_md)
