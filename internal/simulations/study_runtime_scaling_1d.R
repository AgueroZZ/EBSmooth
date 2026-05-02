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
  stop("This internal scaling study requires pkgload.")
}

load_ebsmoothr()

make_benchmark_data_1d <- function(n,
                                   points_per_unit = 100,
                                   noise_sd = 0.22) {
  domain_length <- n / points_per_unit
  t <- seq(0, domain_length, length.out = n)
  truth <- 0.3 + sin(t / 1.8) + 0.15 * cos(t / 3.5)
  s <- rep(noise_sd, n)
  x <- truth + stats::rnorm(n, sd = s)

  list(
    t = t,
    x = x,
    s = s,
    truth = truth,
    domain_length = domain_length
  )
}

choose_num_knots <- function(domain_length,
                             knot_edge = 0.2,
                             min_knots = 30L,
                             max_knots = 500L) {
  as.integer(min(max(min_knots, ceiling(domain_length / knot_edge) + 1L), max_knots))
}

time_one <- function(expr, timeout_seconds = Inf) {
  gc(verbose = FALSE)
  t0 <- proc.time()[["elapsed"]]

  result <- tryCatch(
    {
      if (is.finite(timeout_seconds)) {
        value <- R.utils::withTimeout(
          expr = force(expr),
          timeout = timeout_seconds,
          onTimeout = "error"
        )
      } else {
        value <- force(expr)
      }
      list(status = "completed", value = value, detail = NA_character_)
    },
    TimeoutException = function(e) {
      list(status = "timed_out", value = NULL, detail = conditionMessage(e))
    },
    error = function(e) {
      list(status = "error", value = NULL, detail = conditionMessage(e))
    }
  )

  elapsed <- proc.time()[["elapsed"]] - t0
  if (identical(result$status, "timed_out") && is.finite(timeout_seconds)) {
    elapsed <- timeout_seconds
  }

  list(
    value = result$value,
    elapsed = elapsed,
    status = result$status,
    detail = result$detail
  )
}

build_exact_matern_setup <- function(locations,
                                     max.edge = 0.2,
                                     alpha = 2,
                                     suppress_warnings = TRUE) {
  loc_info <- .normalize_locations(locations)
  meshA <- if (suppress_warnings) {
    suppressWarnings(.build_mesh_A(loc_info$loc, max.edge = max.edge))
  } else {
    .build_mesh_A(loc_info$loc, max.edge = max.edge)
  }

  spde_template <- if (suppress_warnings) {
    suppressWarnings(INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha))
  } else {
    INLA::inla.spde2.matern(mesh = meshA$mesh, alpha = alpha)
  }

  list(
    d = loc_info$d,
    A = meshA$A,
    mesh = meshA$mesh,
    spde_template = spde_template,
    alpha = alpha,
    max.edge = max.edge
  )
}

build_pc_matern_setup <- function(locations,
                                  pc_penalty,
                                  max.edge = 0.2,
                                  alpha = 2,
                                  suppress_warnings = TRUE) {
  exact_setup <- build_exact_matern_setup(
    locations = locations,
    max.edge = max.edge,
    alpha = alpha,
    suppress_warnings = suppress_warnings
  )

  spde_pc <- .build_matern_pc_spde(
    mesh = exact_setup$mesh,
    alpha = alpha,
    pc_penalty = pc_penalty,
    suppress_warnings = suppress_warnings
  )

  idx <- INLA::inla.spde.make.index("spatial.field", n.spde = spde_pc$n.spde)

  c(
    exact_setup,
    list(
      spde_pc = spde_pc,
      idx = idx,
      pc_penalty = pc_penalty,
      n = length(locations)
    )
  )
}

exact_matern_profile_objective <- function(setup, x, s, par, pc_penalty = NULL) {
  state <- .exact_matern_state(
    x = as.numeric(x),
    s = as.numeric(s),
    A = setup$A,
    spde_template = setup$spde_template,
    alpha = setup$alpha,
    d = setup$d,
    log_range = par[1],
    log_sigma = par[2],
    beta0 = par[3]
  )

  out <- state$log_marginal
  if (!is.null(pc_penalty)) {
    out <- out + .log_pc_prior_matern_internal(
      log_range = par[1],
      log_sigma = par[2],
      range_spec = pc_penalty$range,
      sigma_spec = pc_penalty$sigma,
      d = setup$d
    )
  }
  out
}

exact_matern_integrated_flat_objective <- function(setup, x, s, par, pc_penalty = NULL) {
  log_range <- par[1]
  log_sigma <- par[2]

  ts <- .matern_tau_from_range_sigma(
    range = exp(log_range),
    sigma = exp(log_sigma),
    alpha = setup$alpha,
    d = setup$d
  )
  theta_spde <- c(log(ts$tau), log(ts$kappa))

  Q <- INLA::inla.spde2.precision(setup$spde_template, theta = theta_spde)
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))

  w_prec_diag <- 1 / (s^2)
  W <- Matrix::Diagonal(x = w_prec_diag)
  AtW <- Matrix::t(setup$A) %*% W

  Q_post <- Q + AtW %*% setup$A
  Q_post <- Matrix::forceSymmetric(Matrix::Matrix(Q_post, sparse = TRUE))
  chol_post <- Matrix::Cholesky(Q_post, LDL = FALSE, perm = FALSE)

  u <- rep(1, length(x))
  c_u <- as.numeric(AtW %*% u)
  c_x <- as.numeric(AtW %*% x)

  solve_post <- function(rhs) as.numeric(Matrix::solve(chol_post, rhs, system = "A"))

  quad_u <- sum(w_prec_diag * u^2) - sum(c_u * solve_post(c_u))
  if (!is.finite(quad_u) || quad_u <= 0) {
    stop("Integrated-flat beta term is not positive.")
  }

  quad_x <- sum(w_prec_diag * x^2) - sum(c_x * solve_post(c_x))
  cross_ux <- sum(w_prec_diag * x * u) - sum(c_u * solve_post(c_x))
  quad_resid <- quad_x - cross_ux^2 / quad_u

  logdet_D <- sum(log(s^2))
  logdet_Q <- .compute_logdet_spd(Q)
  logdet_Q_post <- .compute_logdet_spd(Q_post)
  logdet_Sigma <- logdet_D - logdet_Q + logdet_Q_post

  out <- -0.5 * (
    (length(x) - 1) * log(2 * pi) +
      logdet_Sigma +
      log(quad_u) +
      quad_resid
  )

  if (!is.null(pc_penalty)) {
    out <- out + .log_pc_prior_matern_internal(
      log_range = log_range,
      log_sigma = log_sigma,
      range_spec = pc_penalty$range,
      sigma_spec = pc_penalty$sigma,
      d = setup$d
    )
  }

  out
}

fit_matern_exact_profile <- function(dat,
                                     g_init,
                                     max.edge = 0.2) {
  fit_fun <- ebnm_Matern_generator(
    locations = dat$t,
    max.edge = max.edge,
    link = "identity"
  )
  fit_fun(dat$x, dat$s, g_init = g_init)
}

fit_matern_pc_package_b0 <- function(dat,
                                     g_init,
                                     pc_penalty,
                                     max.edge = 0.2) {
  fit_fun <- ebnm_Matern_generator(
    locations = dat$t,
    max.edge = max.edge,
    link = "identity",
    pc.penalty = pc_penalty
  )
  fit_fun(dat$x, dat$s, g_init = g_init)
}

fit_matern_exact_integrated_flat <- function(dat,
                                             range_init = 1,
                                             sigma_init = 1,
                                             max.edge = 0.2) {
  setup <- build_exact_matern_setup(dat$t, max.edge = max.edge)
  opt <- stats::optim(
    par = c(log(range_init), log(sigma_init)),
    fn = function(par) {
      val <- tryCatch(
        exact_matern_integrated_flat_objective(setup, dat$x, dat$s, par, pc_penalty = NULL),
        error = function(e) NA_real_
      )
      if (!is.finite(val)) return(Inf)
      -val
    },
    method = "BFGS"
  )
  list(log_likelihood = -opt$value)
}

fit_matern_pc_exact_profile <- function(dat,
                                        pc_penalty,
                                        range_init = 1,
                                        sigma_init = 1,
                                        max.edge = 0.2) {
  setup <- build_exact_matern_setup(dat$t, max.edge = max.edge)
  par0 <- c(log(range_init), log(sigma_init), stats::weighted.mean(dat$x, 1 / (dat$s^2)))
  opt <- stats::optim(
    par = par0,
    fn = function(par) {
      val <- tryCatch(
        exact_matern_profile_objective(setup, dat$x, dat$s, par, pc_penalty = pc_penalty),
        error = function(e) NA_real_
      )
      if (!is.finite(val)) return(Inf)
      -val
    },
    method = "BFGS"
  )
  list(log_likelihood = -opt$value)
}

fit_matern_pc_clinear <- function(dat,
                                  pc_penalty,
                                  range_init = 1,
                                  sigma_init = 1,
                                  max.edge = 0.2) {
  setup <- build_pc_matern_setup(
    locations = dat$t,
    pc_penalty = pc_penalty,
    max.edge = max.edge
  )

  stack <- INLA::inla.stack(
    data = list(Y = as.numeric(dat$x)),
    A = list(setup$A, 1),
    effects = list(
      spatial.field = setup$idx$spatial.field,
      data.frame(beta_cov = rep(1, setup$n))
    ),
    tag = "est"
  )

  formula <- Y ~ -1 +
    f(
      beta_cov,
      model = "clinear",
      hyper = list(
        beta = list(
          initial = stats::weighted.mean(dat$x, 1 / (dat$s^2)),
          fixed = FALSE,
          prior = "flat",
          param = numeric(0)
        )
      )
    ) +
    f(spatial.field, model = setup$spde_pc)

  res <- suppressWarnings(
    INLA::inla(
      formula,
      scale = 1 / (dat$s^2),
      control.inla = list(int.strategy = "eb", strategy = "gaussian"),
      control.family = list(
        control.link = list(model = "identity"),
        hyper = list(prec = list(fixed = TRUE, initial = 0))
      ),
      control.mode = INLA::control.mode(theta = c(
        stats::weighted.mean(dat$x, 1 / (dat$s^2)),
        log(range_init),
        log(sigma_init)
      ), restart = TRUE),
      data = INLA::inla.stack.data(stack),
      control.predictor = list(A = INLA::inla.stack.A(stack), link = 1, compute = TRUE),
      control.compute = list(mlik = TRUE, hyperpar = TRUE, return.marginals = FALSE, config = TRUE),
      silent = TRUE
    )
  )

  list(log_likelihood = as.numeric(res$misc$log.posterior.mode))
}

fit_lgp_package <- function(dat,
                            betaprec) {
  domain_length <- max(dat$t) - min(dat$t)
  num_knots <- choose_num_knots(domain_length)
  setup <- LGP_setup(
    t = dat$t,
    p = 2,
    num_knots = num_knots,
    betaprec = betaprec,
    link = "identity"
  )
  fit_fun <- ebnm_LGP_generator(setup, link = "identity")
  fit_fun(dat$x, dat$s, g_init = LGP(0))
}

fit_public_lgp_known <- function(dat) {
  domain_length <- max(dat$t) - min(dat$t)
  num_knots <- choose_num_knots(domain_length)
  setup <- LGP_setup(
    t = dat$t,
    p = 2,
    num_knots = num_knots,
    betaprec = 0,
    link = "identity"
  )
  eb_smoother(dat$x, s = dat$s, family = "lgp", setup = setup)
}

fit_public_lgp_learned <- function(dat) {
  domain_length <- max(dat$t) - min(dat$t)
  num_knots <- choose_num_knots(domain_length)
  setup <- LGP_setup(
    t = dat$t,
    p = 2,
    num_knots = num_knots,
    betaprec = 0,
    link = "identity"
  )
  eb_smoother(dat$x, s = NULL, family = "lgp", setup = setup)
}

fit_public_matern_exact_known <- function(dat,
                                          max.edge = 0.2) {
  eb_smoother(
    dat$x,
    s = dat$s,
    family = "matern",
    locations = dat$t,
    backend = "exact",
    max.edge = max.edge
  )
}

fit_public_matern_exact_learned <- function(dat,
                                            max.edge = 0.2) {
  eb_smoother(
    dat$x,
    s = NULL,
    family = "matern",
    locations = dat$t,
    backend = "exact",
    max.edge = max.edge
  )
}

fit_public_matern_pc_known <- function(dat,
                                       pc_penalty,
                                       max.edge = 0.2) {
  eb_smoother(
    dat$x,
    s = dat$s,
    family = "matern",
    locations = dat$t,
    backend = "inla_pc",
    pc.penalty = pc_penalty,
    max.edge = max.edge
  )
}

fit_public_matern_pc_learned <- function(dat,
                                         pc_penalty,
                                         max.edge = 0.2) {
  eb_smoother(
    dat$x,
    s = NULL,
    family = "matern",
    locations = dat$t,
    backend = "inla_pc",
    pc.penalty = pc_penalty,
    max.edge = max.edge
  )
}

benchmark_public_api_n3000 <- function(pc_penalty_known,
                                       pc_penalty_learned,
                                       max.edge = 0.2,
                                       timeout_seconds = Inf) {
  dat <- make_benchmark_data_1d(n = 3000L)
  specs <- list(
    list(id = "lgp_known", fit = function() fit_public_lgp_known(dat)),
    list(id = "lgp_learned", fit = function() fit_public_lgp_learned(dat)),
    list(id = "matern_exact_known", fit = function() fit_public_matern_exact_known(dat, max.edge = max.edge)),
    list(id = "matern_exact_learned", fit = function() fit_public_matern_exact_learned(dat, max.edge = max.edge)),
    list(id = "matern_pc_known", fit = function() fit_public_matern_pc_known(dat, pc_penalty = pc_penalty_known, max.edge = max.edge)),
    list(id = "matern_pc_learned", fit = function() fit_public_matern_pc_learned(dat, pc_penalty = pc_penalty_learned, max.edge = max.edge))
  )

  do.call(
    rbind,
    lapply(specs, function(spec) {
      tm <- time_one(spec$fit(), timeout_seconds = timeout_seconds)
      data.frame(
        n = 3000L,
        benchmark = "public_api_single_fit",
        method_id = spec$id,
        elapsed_seconds = tm$elapsed,
        status = tm$status,
        detail = if (is.na(tm$detail)) "" else tm$detail,
        stringsAsFactors = FALSE
      )
    })
  )
}

benchmark_repeated_matern_setup_reuse <- function(pc_penalty,
                                                  max.edge = 0.2,
                                                  timeout_seconds = Inf,
                                                  n = 3000L,
                                                  n_components = 4L) {
  dat <- make_benchmark_data_1d(n = n)
  component_matrix <- sapply(seq_len(n_components), function(k) {
    dat$truth + 0.1 * cos((k + 1) * dat$t / 5) + stats::rnorm(n, sd = dat$s)
  })
  if (is.null(dim(component_matrix))) {
    component_matrix <- matrix(component_matrix, ncol = 1)
  }

  fit_loop_with_locations <- function() {
    lapply(seq_len(ncol(component_matrix)), function(k) {
      eb_smoother(
        component_matrix[, k],
        s = NULL,
        family = "matern",
        locations = dat$t,
        backend = "inla_pc",
        pc.penalty = pc_penalty,
        max.edge = max.edge
      )
    })
  }

  fit_loop_with_setup <- function() {
    setup <- Matern_setup(locations = dat$t, max.edge = max.edge)
    lapply(seq_len(ncol(component_matrix)), function(k) {
      eb_smoother(
        component_matrix[, k],
        s = NULL,
        family = "matern",
        setup = setup,
        backend = "inla_pc",
        pc.penalty = pc_penalty
      )
    })
  }

  results <- list(
    rebuild_each_call = time_one(fit_loop_with_locations(), timeout_seconds = timeout_seconds),
    reuse_matern_setup = time_one(fit_loop_with_setup(), timeout_seconds = timeout_seconds)
  )

  do.call(
    rbind,
    lapply(names(results), function(id) {
      tm <- results[[id]]
      data.frame(
        n = n,
        benchmark = "repeated_matern_pc_learned",
        method_id = id,
        n_components = n_components,
        elapsed_seconds = tm$elapsed,
        status = tm$status,
        detail = if (is.na(tm$detail)) "" else tm$detail,
        stringsAsFactors = FALSE
      )
    })
  )
}

build_method_specs <- function(dat,
                               pc_penalty,
                               max.edge = 0.2,
                               range_init = 1,
                               sigma_init = 1) {
  list(
    list(
      family = "LGP",
      variant = "betaprec = 0",
      implementation = "package",
      id = "lgp_b0",
      fit = function() fit_lgp_package(dat, betaprec = 0)
    ),
    list(
      family = "LGP",
      variant = "betaprec < 0",
      implementation = "package",
      id = "lgp_blt0",
      fit = function() fit_lgp_package(dat, betaprec = -1)
    ),
    list(
      family = "Matern",
      variant = "betaprec = 0",
      implementation = "exact",
      id = "matern_b0_exact",
      fit = function() fit_matern_exact_integrated_flat(
        dat,
        range_init = range_init,
        sigma_init = sigma_init,
        max.edge = max.edge
      )
    ),
    list(
      family = "Matern",
      variant = "betaprec < 0",
      implementation = "exact",
      id = "matern_blt0_exact",
      fit = function() fit_matern_exact_profile(
        dat,
        g_init = Matern(theta = log(range_init), sigma = sigma_init),
        max.edge = max.edge
      )
    ),
    list(
      family = "Matern + PC prior",
      variant = "betaprec = 0",
      implementation = "INLA",
      id = "matern_pc_b0_inla",
      fit = function() fit_matern_pc_package_b0(
        dat,
        g_init = Matern(theta = log(range_init), sigma = sigma_init),
        pc_penalty = pc_penalty,
        max.edge = max.edge
      )
    ),
    list(
      family = "Matern + PC prior",
      variant = "betaprec < 0",
      implementation = "exact",
      id = "matern_pc_blt0_exact",
      fit = function() fit_matern_pc_exact_profile(
        dat,
        pc_penalty = pc_penalty,
        range_init = range_init,
        sigma_init = sigma_init,
        max.edge = max.edge
      )
    ),
    list(
      family = "Matern + PC prior",
      variant = "betaprec < 0",
      implementation = "clinear",
      id = "matern_pc_blt0_clinear",
      fit = function() fit_matern_pc_clinear(
        dat,
        pc_penalty = pc_penalty,
        range_init = range_init,
        sigma_init = sigma_init,
        max.edge = max.edge
      )
    )
  )
}

benchmark_one_method_for_n <- function(n,
                                       method_id,
                                       pc_penalty,
                                       max.edge = 0.2,
                                       range_init = 1,
                                       sigma_init = 1,
                                       timeout_seconds = Inf) {
  dat <- make_benchmark_data_1d(n = n)
  methods <- build_method_specs(
    dat = dat,
    pc_penalty = pc_penalty,
    max.edge = max.edge,
    range_init = range_init,
    sigma_init = sigma_init
  )
  method_ids <- vapply(methods, function(m) m$id, character(1))
  idx <- match(method_id, method_ids)
  if (is.na(idx)) {
    stop("Unknown method_id: ", method_id)
  }

  m <- methods[[idx]]
  tm <- time_one(m$fit(), timeout_seconds = timeout_seconds)

  data.frame(
    n = n,
    domain_length = dat$domain_length,
    family = m$family,
    variant = m$variant,
    implementation = m$implementation,
    method_id = m$id,
    elapsed_seconds = tm$elapsed,
    status = tm$status,
    detail = if (is.na(tm$detail)) "" else tm$detail,
    stringsAsFactors = FALSE
  )
}

benchmark_methods_for_n <- function(n,
                                    pc_penalty,
                                    max.edge = 0.2,
                                    range_init = 1,
                                    sigma_init = 1,
                                    timeout_seconds = Inf) {
  dat <- make_benchmark_data_1d(n = n)

  methods <- build_method_specs(
    dat = dat,
    pc_penalty = pc_penalty,
    max.edge = max.edge,
    range_init = range_init,
    sigma_init = sigma_init
  )

  do.call(
    rbind,
    lapply(methods, function(m) {
      message("  - ", m$id)
      tm <- time_one(m$fit(), timeout_seconds = timeout_seconds)
      data.frame(
        n = n,
        domain_length = dat$domain_length,
        family = m$family,
        variant = m$variant,
        implementation = m$implementation,
        method_id = m$id,
        elapsed_seconds = tm$elapsed,
        status = tm$status,
        detail = if (is.na(tm$detail)) "" else tm$detail,
        stringsAsFactors = FALSE
      )
    })
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

save_outputs <- function(results_df,
                         results_dir,
                         timeout_seconds,
                         public_api_df = NULL,
                         repeated_setup_df = NULL) {
  csv_path <- file.path(results_dir, "runtime_scaling_1d_methods.csv")
  utils::write.csv(results_df, csv_path, row.names = FALSE)

  public_api_csv_path <- if (is.null(public_api_df)) {
    NA_character_
  } else {
    path <- file.path(results_dir, "runtime_public_api_n3000.csv")
    utils::write.csv(public_api_df, path, row.names = FALSE)
    path
  }

  repeated_setup_csv_path <- if (is.null(repeated_setup_df)) {
    NA_character_
  } else {
    path <- file.path(results_dir, "runtime_repeated_matern_n3000.csv")
    utils::write.csv(repeated_setup_df, path, row.names = FALSE)
    path
  }

  results_df$series <- paste0(results_df$variant, " [", results_df$implementation, "]")
  results_df$elapsed_plot_seconds <- ifelse(
    results_df$status == "timed_out",
    timeout_seconds,
    results_df$elapsed_seconds
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_path <- file.path(results_dir, "runtime_scaling_1d_methods.png")

    p <- ggplot2::ggplot(
      results_df,
      ggplot2::aes(
        x = n,
        y = elapsed_plot_seconds,
        color = series,
        shape = status,
        group = series
      )
    ) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::facet_wrap(~ family, scales = "free_y") +
      ggplot2::labs(
        title = "1D Runtime Scaling Across LGP and Matern Variants",
        subtitle = "End-to-end runtime for one fit at each sample size",
        x = "Sample size n (log scale)",
        y = "Elapsed seconds (log scale)",
        color = NULL,
        shape = "Fit status",
        caption = paste0(
          "Timed-out fits are plotted at the timeout threshold (",
          timeout_seconds,
          " seconds)."
        )
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.minor = ggplot2::element_blank()
      )

    ggplot2::ggsave(
      filename = plot_path,
      plot = p,
      width = 12,
      height = 7,
      dpi = 160
    )
  } else {
    plot_path <- NA_character_
  }

  summary_path <- file.path(results_dir, "runtime_scaling_1d_methods_summary.md")

  wide_summary <- reshape(
    results_df[, c("n", "method_id", "elapsed_seconds", "status")],
    idvar = "n",
    timevar = "method_id",
    direction = "wide"
  )
  names(wide_summary) <- sub("^elapsed_seconds\\.", "", names(wide_summary))
  names(wide_summary) <- sub("^status\\.", "status_", names(wide_summary))
  wide_summary <- wide_summary[order(wide_summary$n), , drop = FALSE]

  status_counts <- as.data.frame(table(results_df$status), stringsAsFactors = FALSE)
  names(status_counts) <- c("status", "count")

  md_lines <- c(
    "# 1D Runtime Scaling Study",
    "",
    "This internal benchmark studies how runtime scales with sample size `n`",
    "for the main 1D algorithmic variants currently available or internally",
    "prototyped in this repository.",
    "",
    "## Design",
    "",
    "- Sample sizes: `100, 300, 1000, 3000, 10000`.",
    "- Data are generated on a 1D increasing domain with approximately 100",
    "  observation points per unit length.",
    "- Noise level is fixed at `0.22`.",
    "- The runtime includes the full end-to-end fit for each method, including",
    "  setup/generator construction.",
    "- For Matern methods, the mesh edge is fixed at `0.2` to let latent-field",
    "  complexity grow with the domain length.",
    "- For L-GP, the knot count is chosen from the domain length with a minimum",
    "  of 30 and a maximum of 500.",
    paste0("- Each individual fit is capped at `", timeout_seconds, "` seconds."),
    "",
    "## Methods in the plot",
    "",
    "- `LGP`, `betaprec = 0`",
    "- `LGP`, `betaprec < 0`",
    "- `Matern`, `betaprec = 0`, exact integrated-flat beta objective",
    "- `Matern`, `betaprec < 0`, exact profiled beta objective",
    "- `Matern + PC prior`, `betaprec = 0`, INLA Step A with beta integrated out",
    "- `Matern + PC prior`, `betaprec < 0`, exact profiled objective + PC prior",
    "- `Matern + PC prior`, `betaprec < 0`, INLA `clinear` prototype",
    "",
    "## Status counts",
    "",
    format_md_table(status_counts, digits = 0),
    "",
    if (is.na(plot_path)) {
      "Plot was not generated because `ggplot2` was not available."
    } else {
      paste0("Plot: `", plot_path, "`")
    },
    "",
    "## Timing Table",
    "",
    format_md_table(wide_summary, digits = 4),
    "",
    "## Public API Benchmark at `n = 3000`",
    "",
    if (is.null(public_api_df)) {
      "Not yet rendered in this partial save."
    } else {
      format_md_table(public_api_df[, c("method_id", "elapsed_seconds", "status")], digits = 4)
    },
    "",
    "## Repeated Matern Smoothing Benchmark at `n = 3000`",
    "",
    if (is.null(repeated_setup_df)) {
      "Not yet rendered in this partial save."
    } else {
      format_md_table(repeated_setup_df[, c("method_id", "n_components", "elapsed_seconds", "status")], digits = 4)
    }
  )
  writeLines(md_lines, summary_path)

  list(
    csv_path = csv_path,
    plot_path = plot_path,
    summary_path = summary_path,
    public_api_csv_path = public_api_csv_path,
    repeated_setup_csv_path = repeated_setup_csv_path
  )
}

get_flag_value <- function(args, name) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) {
    return(NULL)
  }
  sub(paste0("^", prefix), "", hit[[1]])
}

ns <- c(100L, 300L, 1000L, 3000L, 10000L)
pc_penalty <- list(
  range = c(anchor = 1.0, alpha = 0.5),
  sigma = c(anchor = 0.7, alpha = 0.5)
)
mesh_edge <- 0.2
range_init <- 1
sigma_init <- 1
timeout_seconds <- 120

single_n <- get_flag_value(commandArgs(trailingOnly = TRUE), "single-n")
single_method <- get_flag_value(commandArgs(trailingOnly = TRUE), "single-method")
single_timeout <- get_flag_value(commandArgs(trailingOnly = TRUE), "timeout-seconds")
render_from_csv <- get_flag_value(commandArgs(trailingOnly = TRUE), "render-from-csv")

pc_penalty_public_known <- pc_penalty
pc_penalty_public_learned <- c(pc_penalty, list(noise = c(anchor = 0.22, alpha = 0.5)))

if (!is.null(single_timeout)) {
  timeout_seconds <- as.numeric(single_timeout)
}

if (!is.null(render_from_csv)) {
  set.seed(999)
  public_api_df <- benchmark_public_api_n3000(
    pc_penalty_known = pc_penalty_public_known,
    pc_penalty_learned = pc_penalty_public_learned,
    max.edge = mesh_edge,
    timeout_seconds = timeout_seconds
  )
  set.seed(1001)
  repeated_setup_df <- benchmark_repeated_matern_setup_reuse(
    pc_penalty = pc_penalty_public_learned,
    max.edge = mesh_edge,
    timeout_seconds = timeout_seconds
  )
  rendered_paths <- save_outputs(
    results_df = utils::read.csv(render_from_csv, stringsAsFactors = FALSE),
    results_dir = results_dir,
    timeout_seconds = timeout_seconds,
    public_api_df = public_api_df,
    repeated_setup_df = repeated_setup_df
  )
  message("Saved timing CSV to: ", rendered_paths$csv_path)
  message("Saved summary markdown to: ", rendered_paths$summary_path)
  if (!is.na(rendered_paths$plot_path)) {
    message("Saved plot to: ", rendered_paths$plot_path)
  }
  quit(save = "no", status = 0)
}

if (!is.null(single_n) || !is.null(single_method)) {
  if (is.null(single_n) || is.null(single_method)) {
    stop("Both --single-n and --single-method must be supplied together.")
  }

  single_df <- benchmark_one_method_for_n(
    n = as.integer(single_n),
    method_id = single_method,
    pc_penalty = pc_penalty,
    max.edge = mesh_edge,
    range_init = range_init,
    sigma_init = sigma_init,
    timeout_seconds = timeout_seconds
  )
  utils::write.csv(single_df, stdout(), row.names = FALSE)
  quit(save = "no", status = 0)
}

warmup_n <- 80L
set.seed(101)
invisible(benchmark_methods_for_n(
  n = warmup_n,
  pc_penalty = pc_penalty,
  max.edge = mesh_edge,
  range_init = range_init,
  sigma_init = sigma_init,
  timeout_seconds = timeout_seconds
))

results_df <- data.frame()
for (i in seq_along(ns)) {
  set.seed(200 + i)
  message("Benchmarking n = ", ns[i], " ...")
  batch_df <- benchmark_methods_for_n(
    n = ns[i],
    pc_penalty = pc_penalty,
    max.edge = mesh_edge,
    range_init = range_init,
    sigma_init = sigma_init,
    timeout_seconds = timeout_seconds
  )
  results_df <- rbind(results_df, batch_df)
  save_outputs(
    results_df = results_df,
    results_dir = results_dir,
    timeout_seconds = timeout_seconds
  )
}
set.seed(999)
public_api_df <- benchmark_public_api_n3000(
  pc_penalty_known = pc_penalty_public_known,
  pc_penalty_learned = pc_penalty_public_learned,
  max.edge = mesh_edge,
  timeout_seconds = timeout_seconds
)
set.seed(1001)
repeated_setup_df <- benchmark_repeated_matern_setup_reuse(
  pc_penalty = pc_penalty_public_learned,
  max.edge = mesh_edge,
  timeout_seconds = timeout_seconds
)
paths <- save_outputs(
  results_df = results_df,
  results_dir = results_dir,
  timeout_seconds = timeout_seconds,
  public_api_df = public_api_df,
  repeated_setup_df = repeated_setup_df
)
message("Saved timing CSV to: ", paths$csv_path)
message("Saved summary markdown to: ", paths$summary_path)
if (!is.na(paths$public_api_csv_path)) {
  message("Saved public API timing CSV to: ", paths$public_api_csv_path)
}
if (!is.na(paths$repeated_setup_csv_path)) {
  message("Saved repeated-setup timing CSV to: ", paths$repeated_setup_csv_path)
}
if (!is.na(paths$plot_path)) {
  message("Saved plot to: ", paths$plot_path)
}
