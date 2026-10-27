# Inference for a fitted model (result of lm_fit / ols_fit).
#
# Under y_t = mu(t; theta) + e_t, e_t iid N(0, sigma^2):
#   sigma^2-hat      = RSS / (n - p)                        (unbiased; Q4)
#   Fisher info      I(theta) = J'J / sigma^2, J = d mu / d theta at theta-hat
#   Wald 95% CI      theta_j -/+ t_{n-p, 0.975} * sqrt([sigma^2-hat (J'J)^{-1}]_jj)   (Q5)
# Wald intervals rely on a locally linear model; profile-likelihood and
# bootstrap intervals are reported alongside as checks on that approximation.

sigma2_hat <- function(fit) fit$rss / (length(fit$resid) - length(fit$theta))

# Gaussian log-likelihood at the MLE (sigma^2 = RSS/n). k counts sigma too.
info_criteria <- function(fit) {
  n <- length(fit$resid); p <- length(fit$theta); k <- p + 1
  ll <- -n / 2 * (log(2 * pi) + log(fit$rss / n) + 1)
  data.frame(RSS = fit$rss, p = p, sigma2 = sigma2_hat(fit), logLik = ll,
             AIC = -2 * ll + 2 * k,
             AICc = -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1),
             BIC = -2 * ll + log(n) * k)
}

# (J'J)^{-1} via the column-scaled Jacobian (raw powers of t up to t^4 and
# near-collinear exponentials make the unscaled matrix badly conditioned).
# Returns all-Inf when J'J is numerically singular.
jtj_inverse <- function(J) {
  s <- sqrt(colSums(J^2))
  inv <- tryCatch(solve(crossprod(sweep(J, 2, s, "/")), tol = 1e-20), error = function(e) NULL)
  if (is.null(inv)) return(matrix(Inf, ncol(J), ncol(J), dimnames = list(colnames(J), colnames(J))))
  inv / outer(s, s)
}

wald_ci <- function(fit, level = 0.95) {
  n <- length(fit$resid); p <- length(fit$theta)
  cov <- sigma2_hat(fit) * jtj_inverse(fit$J)
  se <- sqrt(diag(cov))
  q <- qt(1 - (1 - level) / 2, n - p)
  data.frame(parameter = names(fit$theta), estimate = unname(fit$theta), se = unname(se),
             lower = unname(fit$theta - q * se), upper = unname(fit$theta + q * se))
}

# Identifiability: condition number of the column-scaled J'J (scale-free) and
# the largest absolute correlation between parameter estimates.
identifiability <- function(fit) {
  J <- fit$J
  Js <- sweep(J, 2, sqrt(colSums(J^2)), "/")
  ev <- eigen(crossprod(Js), symmetric = TRUE, only.values = TRUE)$values
  inv <- jtj_inverse(J)
  maxcor <- if (all(is.finite(inv))) { cor <- cov2cor(inv); diag(cor) <- 0; max(abs(cor)) } else 1
  data.frame(condition_number = max(ev) / max(min(ev), .Machine$double.eps * max(ev)),
             max_abs_correlation = maxcor)
}

# Profile-likelihood CI for each parameter: {c : S(c) <= S-hat (1 + F_{1,n-p}/(n-p))},
# where S(c) is the RSS minimised over the other parameters with theta_j = c.
profile_ci <- function(fit, t, y, level = 0.95, max_se = 30) {
  n <- length(y); p <- length(fit$theta)
  cutoff <- fit$rss * (1 + qf(level, 1, n - p) / (n - p))
  se <- wald_ci(fit)$se
  m <- fit$model
  profile_rss <- function(j, c, start) {
    sub <- list(pars = m$pars[-j],
                mu = function(th, t) { full <- append(th, c, after = j - 1); m$mu(full, t) },
                jac = function(th, t) { full <- append(th, c, after = j - 1); m$jac(full, t)[, -j, drop = FALSE] })
    f <- lm_fit(sub, t, y, start[-j], max_iter = 200)
    list(rss = f$rss, theta = append(f$theta, c, after = j - 1))
  }
  bound <- function(j, dir) {
    # Walk outward in steps of se_j until RSS crosses the cutoff, then bisect.
    prev <- list(c = fit$theta[j], theta = fit$theta)
    for (k in seq_len(max_se)) {
      c <- fit$theta[j] + dir * k * se[j]
      pr <- profile_rss(j, c, prev$theta)
      if (!is.finite(pr$rss)) return(NA_real_)
      if (pr$rss > cutoff) {
        g <- function(cc) profile_rss(j, cc, prev$theta)$rss - cutoff
        return(uniroot(g, sort(c(prev$c, c)), tol = 1e-8 * se[j])$root)
      }
      prev <- list(c = c, theta = pr$theta)
    }
    NA_real_                                          # not bounded within max_se SEs
  }
  data.frame(parameter = m$pars,
             lower = vapply(seq_len(p), bound, 0, dir = -1),
             upper = vapply(seq_len(p), bound, 0, dir = +1))
}

# Residual bootstrap: resample centred, leverage-free rescaled residuals,
# refit from theta-hat, take percentile intervals.
bootstrap_ci <- function(fit, t, y, B = 2000, level = 0.95, seed = 686) {
  set.seed(seed)
  n <- length(y); p <- length(fit$theta)
  e <- (fit$resid - mean(fit$resid)) * sqrt(n / (n - p))
  draws <- t(replicate(B, {
    ystar <- fit$fitted + sample(e, n, replace = TRUE)
    f <- if (fit$iterations == 0) ols_fit(fit$model, t, ystar) else lm_fit(fit$model, t, ystar, fit$theta)
    if (f$converged) f$theta else rep(NA_real_, p)
  }))
  ok <- stats::complete.cases(draws)
  a <- (1 - level) / 2
  data.frame(parameter = names(fit$theta),
             lower = apply(draws[ok, , drop = FALSE], 2, quantile, a),
             upper = apply(draws[ok, , drop = FALSE], 2, quantile, 1 - a),
             converged_share = mean(ok), row.names = NULL)
}

# Extra-sum-of-squares F test for nested models (small model `f0` inside `f1`).
nested_f_test <- function(f0, f1) {
  n <- length(f1$resid); p0 <- length(f0$theta); p1 <- length(f1$theta)
  F <- ((f0$rss - f1$rss) / (p1 - p0)) / (f1$rss / (n - p1))
  data.frame(df1 = p1 - p0, df2 = n - p1, F = F,
             p_value = pf(F, p1 - p0, n - p1, lower.tail = FALSE),
             LR = n * log(f0$rss / f1$rss),
             p_value_LR = pchisq(n * log(f0$rss / f1$rss), p1 - p0, lower.tail = FALSE))
}
