library(EBSmoothr)
library(smashr)
library(splines)

# simulate a function
f <- function(x) {
  # bspline basis with local support
  basis <- bs(x, df = 10, degree = 3)
  # coefficients for the basis functions
  coefs <- c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0)
  # linear combination of the basis functions
  return(basis %*% coefs)
}

# simulate data
set.seed(123)
n <- 100
loc <- seq(0, 1, length.out = n)
x <- f(loc) + rnorm(n, sd = 0.2)

# fit the EBSmoothr model
ebsmoothr_nn <- ebnm_Matern_generator(
  locations = loc,
  link = "log", # softplus
  backend = "fisher_pql", # laplace_fisher
  pql_inner_iter = 10
)(x = x, s = rep(0.2, n))

# plot the result
plot(
  loc,
  x,
  main = "EBSmoothr Fit",
  cex = 0.5,
  xlab = "Location",
  ylab = "Value"
)
lines(loc, ebsmoothr_nn$posterior$mean, col = "blue", lwd = 2)
lines(loc, f(loc), col = "red", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("EBSmoothr Fit", "True Function"),
  col = c("blue", "red"),
  lwd = 2,
  lty = c(1, 2)
)


# What is the ML?
vec <- c()
for (iter in 1:10) {
  ebsmoothr_nn <- ebnm_Matern_generator(
    locations = loc,
    link = "softplus", # softplus
    backend = "fisher_pql", # laplace_fisher
    pql_inner_iter = iter
  )(x = x, s = rep(0.2, n))
  vec <- c(vec, ebsmoothr_nn$log_likelihood)
}

ebsmoothr_nn_2 <- ebnm_Matern_generator(
  locations = loc,
  link = "softplus", # softplus
  backend = "laplace_fisher" # laplace_fisher
)(x = x, s = rep(0.2, n))
ebsmoothr_nn_2$log_likelihood

ebsmoothr_nn_3 <- ebnm_Matern_generator(
  locations = loc,
  link = "softplus", # softplus
  backend = "laplace" # laplace_fisher
)(x = x, s = rep(0.2, n))
ebsmoothr_nn_3$log_likelihood

plot(
  vec,
  log = "y",
  type = "b",
  xlab = "PQL Inner Iteration",
  ylab = "Log-Likelihood"
)
abline(h = ebsmoothr_nn_2$log_likelihood, col = "red", lty = 2)
abline(h = ebsmoothr_nn_3$log_likelihood, col = "blue", lty = 2)
