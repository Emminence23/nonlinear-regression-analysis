# Least-squares fitting.
#
# Every nonlinear model here is *partially linear*: once the nonlinear
# parameters (the rates b, or the denominator slope) are fixed, the mean is
# linear in the remaining amplitudes, which then have a closed-form OLS
# solution. Minimising the RSS over the nonlinear parameters alone, with the
# linear ones profiled out ("variable projection", Golub & Pereyra 1973),
# shrinks a 5-D non-convex search for Model 1 to a 2-D one that can be
# scanned exhaustively. The best grid point then seeds a full
# Levenberg-Marquardt refinement of all parameters jointly.

# OLS of y on X via QR; returns coefficients and RSS.
ols <- function(X, y) {
  q <- qr(X)
  if (q$rank < ncol(X)) return(list(coef = rep(NA_real_, ncol(X)), rss = Inf))
  list(coef = qr.coef(q, y), rss = sum(qr.resid(q, y)^2))
}

# ---- Variable-projection grid searches -------------------------------------

# Rate grid for the exponential models: negative rates, log-spaced magnitudes.
rate_grid <- function(n = 80, lo = 1e-3, hi = 3) -exp(seq(log(hi), log(lo), length.out = n))

vp_grid_M1 <- function(t, y, grid = rate_grid()) {
  best <- list(rss = Inf)
  surf <- expand.grid(b1 = grid, b2 = grid)
  surf <- surf[surf$b1 < surf$b2, ]               # order the rates: b1 < b2
  surf$rss <- NA_real_
  for (i in seq_len(nrow(surf))) {
    b <- c(surf$b1[i], surf$b2[i])
    f <- ols(cbind(1, exp(b[1] * t), exp(b[2] * t)), y)
    surf$rss[i] <- f$rss
    if (f$rss < best$rss) best <- list(theta = c(f$coef, b), rss = f$rss)
  }
  best$surface <- surf
  best
}

vp_grid_M1s <- function(t, y, grid = rate_grid(400)) {
  rss <- vapply(grid, function(b) ols(cbind(1, exp(b * t)), y)$rss, 0)
  b <- grid[which.min(rss)]
  f <- ols(cbind(1, exp(b * t)), y)
  list(theta = c(f$coef, b), rss = f$rss)
}

vp_grid_M1c <- function(t, y, grid = rate_grid(400)) {
  rss <- vapply(grid, function(b) ols(cbind(1, exp(b * t), t * exp(b * t)), y)$rss, 0)
  b <- grid[which.min(rss)]
  f <- ols(cbind(1, exp(b * t), t * exp(b * t)), y)
  list(theta = c(f$coef, b), rss = f$rss)
}

# Given the decay rho and frequency w, the wave is linear in (a0, c, d).
vp_grid_M1w <- function(t, y, rho = rate_grid(60), w = seq(0.02, pi, length.out = 60)) {
  best <- list(rss = Inf)
  for (r in rho) for (om in w) {
    e <- exp(r * t)
    f <- ols(cbind(1, e * cos(om * t), e * sin(om * t)), y)
    if (f$rss < best$rss) best <- list(theta = c(f$coef, r, om), rss = f$rss)
  }
  best
}

# Model 2 with b0 = 1: the denominator 1 + b1 t must stay positive on the data
# range, i.e. b1 > -1 / max(t).
vp_grid_M2 <- function(t, y, grid = NULL) {
  if (is.null(grid)) {
    lo <- -1 / max(t)
    grid <- c(lo * (1 - exp(seq(log(0.999), log(1e-3), length.out = 50))),
              exp(seq(log(1e-4), log(100), length.out = 400)))
  }
  rss <- vapply(grid, function(b) { d <- 1 + b * t; ols(cbind(1 / d, t / d), y)$rss }, 0)
  b <- grid[which.min(rss)]
  d <- 1 + b * t
  f <- ols(cbind(1 / d, t / d), y)
  list(theta = c(f$coef, b), rss = f$rss, grid = data.frame(b1 = grid, rss = rss))
}

# ---- Levenberg-Marquardt ----------------------------------------------------

# Minimise ||y - mu(theta)||^2 with Marquardt's diagonal scaling.
# Returns the estimate plus convergence diagnostics.
lm_fit <- function(model, t, y, theta0, max_iter = 500, tol = 1e-12) {
  th <- unname(theta0)
  r <- y - model$mu(th, t)
  rss <- sum(r^2)
  lambda <- 1e-3
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    J <- model$jac(th, t)
    A <- crossprod(J); g <- crossprod(J, r)
    repeat {
      step <- tryCatch(solve(A + lambda * diag(diag(A)), g), error = function(e) NULL)
      if (!is.null(step)) {
        th_new <- th + drop(step)
        r_new <- y - model$mu(th_new, t)
        rss_new <- sum(r_new^2)
        if (is.finite(rss_new) && rss_new <= rss) break
      }
      lambda <- lambda * 10
      if (lambda > 1e16) break
    }
    if (lambda > 1e16) { converged <- TRUE; break }   # no downhill step left: at a minimum
    done <- (rss - rss_new) <= tol * (rss + tol)
    th <- th_new; r <- r_new; rss <- rss_new
    lambda <- max(lambda / 10, 1e-12)
    if (done) { converged <- TRUE; break }
  }
  J <- model$jac(th, t)
  names(th) <- model$pars
  list(model = model, theta = th, theta0 = setNames(unname(theta0), model$pars),
       rss = rss, resid = r, fitted = y - r, J = J,
       iterations = iter, converged = converged,
       grad_norm = sqrt(sum(crossprod(J, r)^2)))
}

# Linear model: closed form, wrapped in the same result structure.
ols_fit <- function(model, t, y) {
  X <- model$jac(NULL, t)
  f <- ols(X, y)
  th <- setNames(f$coef, model$pars)
  r <- y - drop(X %*% th)
  list(model = model, theta = th, theta0 = th, rss = sum(r^2), resid = r,
       fitted = y - r, J = X, iterations = 0L, converged = TRUE,
       grad_norm = sqrt(sum(crossprod(X, r)^2)))
}
