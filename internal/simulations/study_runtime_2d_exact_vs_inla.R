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
  stop("This internal runtime study requires pkgload.")
}

load_ebsmoothr()

make_benchmark_data_2d <- function(nx = 60L,
                                   ny = 50L,
                                   noise_sd = 0.18) {
  grid <- expand.grid(
    x = seq(0, 1, length.out = nx),
    y = seq(0, 1, length.out = ny)
  )

  truth <- with(
    grid,
    0.2 +
      sin(2 * pi * x) * cos(1.5 * pi * y) +
      0.15 * cos(pi * (x + y))
  )
  s <- rep(noise_sd, nrow(grid))
  x <- truth + stats::rnorm(nrow(grid), sd = s)

  list(
    locations = as.matrix(grid),
    truth = truth,
    x = x,
    s = s,
    n = nrow(grid),
    d = 2L
  )
}

time_one <- function(expr) {
  gc(verbose = FALSE)
  t0 <- proc.time()[["elapsed"]]

  result <- tryCatch(
    list(status = "completed", value = force(expr), detail = NA_character_),
    error = function(e) list(status = "error", value = NULL, detail = conditionMessage(e))
  )

  list(
    elapsed = proc.time()[["elapsed"]] - t0,
    status = result$status,
    detail = result$detail,
    value = result$value
  )
}

safe_number <- function(x) {
  if (length(x) == 0L || !is.finite(x)) NA_real_ else as.numeric(x)
}

set.seed(4172026)

dat <- make_benchmark_data_2d()
setup <- Matern_setup(locations = dat$locations, max.edge = c(0.1, 0.2))
n_spde <- setup$spde_template$n.spde
sigma_anchor <- stats::sd(dat$x)
if (!is.finite(sigma_anchor) || sigma_anchor <= 0) sigma_anchor <- 1
pc_penalty <- list(
  range = c(0.2, 0.5),
  sigma = c(sigma_anchor, 0.5),
  noise = c(sigma_anchor, 0.5)
)

exact_known <- time_one(
  eb_smoother(
    x = dat$x,
    s = dat$s,
    family = "matern",
    setup = setup,
    pc.penalty = NULL
  )
)

exact_learned <- time_one(
  eb_smoother(
    x = dat$x,
    s = NULL,
    family = "matern",
    setup = setup,
    pc.penalty = NULL
  )
)

pc_learned <- time_one(
  eb_smoother(
    x = dat$x,
    s = NULL,
    family = "matern",
    setup = setup,
    backend = "inla_pc",
    pc.penalty = pc_penalty
  )
)

exact_objective <- time_one(
  .exact_matern_profile_objective_unknown_noise(
    x = dat$x,
    A = setup$A,
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    log_range = log(0.2),
    log_sigma = log(sigma_anchor),
    log_noise_sd = log(median(dat$s))
  )
)

timings <- data.frame(
  case = c(
    "matern_exact_known_2d",
    "matern_exact_learned_2d",
    "matern_pc_learned_2d",
    "matern_exact_objective_eval_2d"
  ),
  elapsed_seconds = c(
    exact_known$elapsed,
    exact_learned$elapsed,
    pc_learned$elapsed,
    exact_objective$elapsed
  ),
  status = c(
    exact_known$status,
    exact_learned$status,
    pc_learned$status,
    exact_objective$status
  ),
  detail = c(
    exact_known$detail,
    exact_learned$detail,
    pc_learned$detail,
    exact_objective$detail
  ),
  n = dat$n,
  d = dat$d,
  n_spde = n_spde,
  max_edge = I(list(c(0.1, 0.2), c(0.1, 0.2), c(0.1, 0.2), c(0.1, 0.2))),
  stringsAsFactors = FALSE
)

timings$max_edge <- vapply(
  timings$max_edge,
  function(z) paste(z, collapse = ", "),
  character(1)
)

csv_path <- file.path(results_dir, "runtime_2d_exact_vs_inla.csv")
utils::write.csv(timings, csv_path, row.names = FALSE)

summary_lines <- c(
  "# Runtime Benchmark: 2D Exact Matern vs INLA-PC",
  "",
  sprintf("- Date: %s", as.character(Sys.Date())),
  sprintf("- Data geometry: %d x %d regular grid (%d observations)", 60L, 50L, dat$n),
  sprintf("- Matern setup max.edge: %s", paste(c(0.1, 0.2), collapse = ", ")),
  sprintf("- SPDE nodes: %d", n_spde),
  "",
  "## Results",
  "",
  sprintf("- Exact known-noise fit: %.3f s (%s)", safe_number(exact_known$elapsed), exact_known$status),
  sprintf("- Exact learned-noise fit: %.3f s (%s)", safe_number(exact_learned$elapsed), exact_learned$status),
  sprintf("- INLA-PC learned-noise fit: %.3f s (%s)", safe_number(pc_learned$elapsed), pc_learned$status),
  sprintf("- Single exact objective evaluation: %.3f s (%s)", safe_number(exact_objective$elapsed), exact_objective$status),
  "",
  "## Notes",
  "",
  "- The exact objective evaluation isolates the shared sparse-factorization path used by the 2D exact optimizer.",
  "- Compare this file against earlier slow-path timings to catch regressions in permutation-aware sparse linear algebra."
)

summary_path <- file.path(results_dir, "runtime_2d_exact_vs_inla_summary.md")
writeLines(summary_lines, con = summary_path)

message("Wrote benchmark outputs to:")
message("  - ", csv_path)
message("  - ", summary_path)
