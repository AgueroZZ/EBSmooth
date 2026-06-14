load("EBFPCA/data/weather_data.rda")

library(flashier)
library(EBSmoothr)


# Fit EBMF model, unconstrained
# pc.penalty anchors are in the units of `locations` (days here):
# P(range < 60 days) = 0.05, P(sigma > 1) = 0.5
# ebnm_factor <- EBSmoothr::ebnm_Matern_generator(
#   locations = seq(1, 365),
#   link = "softplus",
#   backend = "laplace_fisher"
# )
ebnm_factor <- EBSmoothr::ebnm_LGP_generator(
  LGP_setup = LGP_setup(t = seq(1, 365))
)

ebnm_loading <- flash_ebnm(
  prior_family = "point_exponential",
  mode = 0,
  scale = "estimate"
)

Y_centered <- Y - rowMeans(Y)
# Y_pos <- Y - min(Y) + 1 # shift to positive for EBSmoothr

# # Try EBSmoothr on a single location to see if it captures the seasonal pattern
# location_idx <- 1
# location_data <- Y[location_idx, ]
# location_dates <- as.Date("2023-01-01") + 0:364
# # ebnm_fit <- EBSmoothr::eb_smoother(x = location_data, setup = LGP_setup(t = location_dates), s = NULL, family = "lgp", link = "identity")
# ebnm_fit <- EBSmoothr::eb_smoother(
#   x = location_data,
#   locations = location_dates,
#   s = NULL,
#   family = "matern",
#   link = "identity",
#   pc.penalty = list(range = c(60, 0.05), sigma = 1)
# )
# plot(
#   ebnm_fit$posterior$mean ~ location_dates,
#   type = "l",
#   col = "blue",
#   main = "EBSmoothr Fit for One Location",
#   xlab = "Date",
#   ylab = "Temperature",
#   ylim = range(c(location_data, ebnm_fit$posterior$mean))
# )
# points(location_data ~ location_dates, col = "red", pch = 16, cex = 0.5)

# Step-by-step equivalent of
flashier_fit <- flash_init(Y_centered, S = NULL, var_type = 0) |>
  flash_set_conv_crit(flash_conv_crit_max_chg, tol = 1e-5) |>
  flash_set_verbose(verbose = 3) |>
  flash_greedy(
    Kmax = 5,
    ebnm_fn = c(ebnm_loading, ebnm_factor),
    maxiter = 100
  ) |>
  flash_backfit(maxiter = 50) |>
  flash_nullcheck()


K <- ncol(flashier_fit$F_pm)
factor_cols <- hcl.colors(K, "Dark 3")
dates <- as.Date("2023-01-01") + 0:364
month_starts <- which(format(dates, "%d") == "01")

matplot(
  flashier_fit$F_pm,
  type = "l",
  lty = 1,
  col = factor_cols,
  ylim = range(flashier_fit$F_pm),
  main = "Factors from EBMF model",
  xlab = "Date",
  ylab = "Factor",
  xaxt = "n"
)
axis(1, at = month_starts, labels = format(dates[month_starts], "%d/%b"))
legend(
  "topright",
  legend = paste("Factor", 1:K),
  col = factor_cols,
  lty = 1,
  cex = 0.8
)


plot_stations <- function(stations, Y, dates, cols = NULL) {
  idx <- match(stations, rownames(Y))
  if (anyNA(idx)) {
    stop("Unknown station(s): ", paste(stations[is.na(idx)], collapse = ", "))
  }
  if (is.null(cols)) {
    cols <- hcl.colors(length(stations), "Dark 3")
  }
  plot(
    dates,
    rep(NA_real_, length(dates)),
    ylim = range(Y),
    main = "Temperature",
    xlab = "Date",
    ylab = "Temperature"
  )
  for (i in seq_along(idx)) {
    points(dates, Y[idx[i], ], col = cols[i], pch = 16, cex = 0.4)
  }
  legend(
    "topright",
    legend = stations,
    col = cols,
    pch = 16,
    cex = 0.6,
    pt.cex = 0.8,
    bty = "n"
  )
}
plot_stations(
  c("sydney", "stjohns", "victoria", "princeru", "inuvik", "resolute"),
  Y = Y_centered,
  dates = dates
)


# Heatmap of the loadings: stations x factors, diverging palette centered at 0
L <- flashier_fit$L_pm
rownames(L) <- rownames(Y)
colnames(L) <- paste("Factor", 1:K)
loading_lim <- max(abs(L))
heatmap(
  L,
  Colv = NA,
  scale = "none",
  col = hcl.colors(101, "Blue-Red 3"),
  breaks = seq(-loading_lim, loading_lim, length.out = 102),
  margins = c(8, 6),
  main = "Loadings (L_pm)"
)


##### FPCA via fdapace (PACE), compared with the EBMF factors
library(fdapace)
day <- 1:365

# PACE with kernel-smoothed mean and covariance
# (default for dense data is unsmoothed cross-sectional estimates)
pace_fit <- FPCA(
  Ly = asplit(Y, 1),
  Lt = rep(list(day), nrow(Y)),
  optns = list(
    dataType = "Dense",
    methodMuCovEst = "smooth",
    # keep as many PCs as EBMF factors (default selects by 99% FVE)
    methodSelectK = K
  )
)
K_cmp <- min(K, ncol(pace_fit$phi))
pace_pcs <- pace_fit$phi[, 1:K_cmp, drop = FALSE]
pace_fve <- c(pace_fit$cumFVE[1], diff(pace_fit$cumFVE[1:K_cmp]))

# EBMF factors have arbitrary scale and sign; normalize to unit norm and
# align signs with the PACE eigenfunctions
ebmf_pcs <- apply(
  flashier_fit$F_pm[, 1:K_cmp, drop = FALSE],
  2,
  function(f) f / sqrt(sum(f^2))
)
for (k in 1:K_cmp) {
  if (sum(ebmf_pcs[, k] * pace_pcs[, k]) < 0) {
    ebmf_pcs[, k] <- -ebmf_pcs[, k]
  }
}

ticks <- month_starts[seq(1, 12, by = 2)]
yl <- range(ebmf_pcs, pace_pcs)
par(mfrow = c(ceiling(K_cmp / 3), min(K_cmp, 3)))
for (k in 1:K_cmp) {
  plot(
    day,
    ebmf_pcs[, k],
    type = "l",
    col = "red",
    ylim = yl,
    main = sprintf("Component %d (PACE FVE %.1f%%)", k, 100 * pace_fve[k]),
    xlab = "Date",
    ylab = "Eigenfunction",
    xaxt = "n"
  )
  axis(1, at = ticks, labels = format(dates[ticks], "%d/%b"))
  lines(pace_fit$workGrid, pace_pcs[, k], col = "blue")
}
par(mfrow = c(1, 1))


# Heatmap of the PACE scores: stations x PCs, diverging palette centered at 0
pace_scores <- pace_fit$xiEst[, 1:K_cmp, drop = FALSE]
rownames(pace_scores) <- rownames(Y)
colnames(pace_scores) <- paste("PC", 1:K_cmp)
score_lim <- max(abs(pace_scores))
heatmap(
  pace_scores,
  Colv = NA,
  scale = "none",
  col = hcl.colors(101, "Blue-Red 3"),
  breaks = seq(-score_lim, score_lim, length.out = 102),
  margins = c(8, 6),
  main = "PACE FPC scores"
)
