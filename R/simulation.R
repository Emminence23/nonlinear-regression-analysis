# Simulation checks: generate fake data sets from a fitted model (same t,
# same noise level), refit with exactly the pipeline used on the real data,
# and see whether the methods behave as they should.

# Fit the candidate models to one data set. Model 1's two-exponential fit is
# included only to see whether it breaks down (rates merging).
fit_candidates <- function(t, y, two_exp = TRUE) {
  fits <- list(
    M1c = lm_fit(MODELS$M1c, t, y, vp_grid_M1c(t, y, rate_grid(150))$theta),
    M1s = lm_fit(MODELS$M1s, t, y, vp_grid_M1s(t, y, rate_grid(150))$theta),
    M2  = lm_fit(MODELS$M2,  t, y, vp_grid_M2(t, y)$theta),
    M3  = ols_fit(MODELS$M3, t, y)
  )
  if (two_exp) {
    f1 <- lm_fit(MODELS$M1, t, y, vp_grid_M1(t, y, rate_grid(40))$theta)
    # The two-exponential fit "breaks down" when it cannot beat the confluent
    # model (its boundary limit) at a genuinely separate pair of rates.
    gap <- abs(f1$theta[["b2"]] - f1$theta[["b1"]])
    fits$breakdown <- !(f1$converged && gap > 1e-3 && f1$rss < fits$M1c$rss - 1e-9)
  }
  fits
}

# The rule used in the report: pick the family by AIC, then inside Model 1
# keep one exponential unless the F test says the extra term is needed.
choose_model <- function(fits) {
  aic <- vapply(fits[c("M1c", "M1s", "M2", "M3")], function(f) info_criteria(f)$AIC, 0)
  best <- names(which.min(aic))
  if (best %in% c("M1c", "M1s"))
    best <- if (nested_f_test(fits$M1s, fits$M1c)$p_value > 0.05) "M1s" else "M1c"
  best
}

# One simulated data set from `truth` (a fitted model). Returns the chosen
# model, whether Model 1 broke down, and whether each 95% Wald interval of
# the true model covers the true parameter value.
simulate_once <- function(truth, name, t) {
  sigma <- sqrt(sigma2_hat(truth))
  y <- truth$fitted + rnorm(length(t), 0, sigma)
  fits <- fit_candidates(t, y)
  refit <- fits[[name]]
  ci <- wald_ci(refit)
  list(chosen = choose_model(fits), breakdown = fits$breakdown,
       covered = setNames(ci$lower <= truth$theta & truth$theta <= ci$upper, names(truth$theta)),
       converged = refit$converged)
}

# Is the wave real? Under the simple model (one exponential), how often does
# the wave model improve the RSS by at least the observed amount?
wave_once <- function(simple, t) {
  y <- simple$fitted + rnorm(length(t), 0, sqrt(sigma2_hat(simple)))
  f0 <- lm_fit(MODELS$M1s, t, y, vp_grid_M1s(t, y, rate_grid(150))$theta)
  f1 <- lm_fit(MODELS$M1w, t, y, vp_grid_M1w(t, y)$theta)
  f0$rss - f1$rss
}
