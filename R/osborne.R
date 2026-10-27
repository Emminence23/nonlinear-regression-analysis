# Osborne's modified Prony algorithm for sums of exponentials
# (Osborne 1975; Osborne & Smyth 1991; MTH 686 course notes, not included here).
#
# For equally spaced t = 1..n, any mu_t = sum_k a_k x_k^t (x_k = e^{b_k})
# satisfies the linear recurrence  g_0 mu_t + g_1 mu_{t+1} + ... + g_p mu_{t+p} = 0,
# where the x_k are the roots of g(x) = sum_j g_j x^j. Writing that recurrence
# as G' mu = 0 (G is n x (n-p), banded), the profiled least-squares criterion
# becomes
#     R(g) = Y' P_G Y,      P_G = G (G'G)^{-1} G',
# a function of the p+1 recurrence coefficients only. Setting dR/dg = 0 gives
# B(g) g = 0 with B symmetric; Osborne iterates g <- eigenvector of B(g) with
# the smallest |eigenvalue| until it converges.
#
# Model 1 has a constant term a0 = a0 * 1^t, i.e. a root fixed at x = 1. We
# impose it exactly by writing g(x) = (x - 1) h(x), so g = C h is linear in the
# free coefficients h, and the iteration runs on h.
#
# Note: differentiating R(g) gives
#   b_jk = Y'G_j W G_k'Y + Y'G_k W G_j'Y - v'(G_j'G_k + G_k'G_j) v,
#   W = (G'G)^{-1}, v = W G'Y,
# i.e. the middle term carries a minus sign (the handout prints a plus).

# Banded (n-p) x n matrix whose rows apply the recurrence with coefficients g.
recurrence_rows <- function(g, n) {
  p <- length(g) - 1
  Gt <- matrix(0, n - p, n)
  for (j in 0:p) Gt[cbind(seq_len(n - p), seq_len(n - p) + j)] <- g[j + 1]
  Gt
}

# Map h -> g for g(x) = (x - 1) h(x); q = length(h), coefficients in increasing powers.
root_at_one <- function(q) {
  C <- matrix(0, q + 1, q)
  for (k in seq_len(q)) { C[k, k] <- -1; C[k + 1, k] <- 1 }
  C
}

osborne_B <- function(h, y, C) {
  n <- length(y)
  g <- drop(C %*% h)
  Gt <- recurrence_rows(g, n)
  # dG'/dh_k = sum_j C[j, k] * (shift matrix j)
  Ht <- lapply(seq_len(ncol(C)), function(k) recurrence_rows(C[, k], n))
  W <- solve(tcrossprod(Gt))
  v <- W %*% (Gt %*% y)
  HY <- vapply(Ht, function(H) drop(H %*% y), numeric(nrow(Gt)))   # columns: H_k' Y
  q <- ncol(C)
  B <- matrix(0, q, q)
  for (j in seq_len(q)) for (k in seq_len(q)) {
    HjHk <- Ht[[j]] %*% t(Ht[[k]])
    B[j, k] <- 2 * drop(crossprod(HY[, j], W %*% HY[, k])) -
      drop(crossprod(v, (HjHk + t(HjHk)) %*% v))
  }
  B
}

# Classical Prony start: minimise ||M C h|| over |h| = 1, where M has rows
# (y_t, ..., y_{t+p}). This is the unweighted (biased) version of the problem
# Osborne's iteration solves exactly.
prony_start <- function(y, C) {
  p <- nrow(C) - 1; n <- length(y)
  M <- sapply(0:p, function(j) y[(1 + j):(n - p + j)])
  sv <- svd(M %*% C)
  sv$v[, ncol(C)]
}

# Profiled criterion R(h) = Y' P_G Y = RSS of the best fit within the
# recurrence class defined by h.
osborne_R <- function(h, y, C) {
  Gt <- recurrence_rows(drop(C %*% h), length(y))
  z <- Gt %*% y
  drop(crossprod(z, solve(tcrossprod(Gt), z)))
}

# Converged when the unit eigenvector moves by less than `tol`. Near a root
# close to 1, G'G is poorly conditioned and rounding alone moves h by ~1e-9
# per step, so a much smaller tolerance is never met.
osborne_iterate <- function(h, y, C, max_iter = 500, tol = 1e-8) {
  h <- h / sqrt(sum(h^2))
  for (iter in seq_len(max_iter)) {
    e <- eigen(osborne_B(h, y, C), symmetric = TRUE)
    h_new <- e$vectors[, which.min(abs(e$values))]
    h_new <- h_new * sign(sum(h_new * h))  # resolve the eigenvector sign
    delta <- sqrt(sum((h_new - h)^2))
    h <- h_new
    if (delta < tol) break
  }
  list(h = h, R = osborne_R(h, y, C), roots = polyroot(h),
       iterations = iter, converged = delta < tol)
}

# Check B(h) h against a central-difference gradient of R(h) at the points h
# (columns of H). B(h) h equals dR/dh, so the two should agree to rounding.
# Returns the largest relative difference.
osborne_gradient_check <- function(H, y, C, eps = 1e-6) {
  max(apply(H, 2, function(h) {
    h <- h / sqrt(sum(h^2))
    num <- vapply(seq_along(h), function(k) {
      e <- replace(numeric(length(h)), k, eps)
      (osborne_R(h + e, y, C) - osborne_R(h - e, y, C)) / (2 * eps)
    }, 0)
    an <- drop(osborne_B(h, y, C) %*% h)
    max(abs(an - num)) / max(abs(num))
  }))
}

# Coefficients of h(x) = (x - e^{b1}) (x - e^{b2}).
h_from_rates <- function(b) c(exp(sum(b)), -sum(exp(b)), 1)

# The iteration finds *a* stationary point of R(h), not necessarily the global
# minimum, so it is run from the Prony start plus a spread of rate pairs and
# the best stationary point is kept.
osborne_M1 <- function(t, y, starts = NULL) {
  stopifnot(all(diff(t) == 1))
  C <- root_at_one(3)                     # h quadratic: two free exponentials + the constant
  if (is.null(starts)) {
    rates <- c(-0.05, -0.2, -0.5, -1, -2)
    pairs <- t(combn(rates, 2))
    starts <- c(list(prony = prony_start(y, C)),
                setNames(lapply(seq_len(nrow(pairs)), function(i) h_from_rates(pairs[i, ])),
                         sprintf("b=(%g,%g)", pairs[, 1], pairs[, 2])))
  }
  runs <- lapply(starts, osborne_iterate, y = y, C = C)
  summary <- data.frame(
    start = names(starts),
    R = vapply(runs, `[[`, 0, "R"),
    iterations = vapply(runs, `[[`, 0L, "iterations"),
    converged = vapply(runs, `[[`, NA, "converged"),
    roots = vapply(runs, function(r) paste(format(r$roots, digits = 4), collapse = ", "), ""),
    real_positive = vapply(runs, function(r) all(abs(Im(r$roots)) < 1e-8) && all(Re(r$roots) > 0), NA)
  )
  best <- runs[[which.min(summary$R)]]
  x <- best$roots
  best$real_positive <- all(abs(Im(x)) < 1e-8) && all(Re(x) > 0)
  # Complex pair x = r e^{+-i w}  <=>  mu_t = a0 + r^t (c cos(w t) + d sin(w t)).
  best$modulus <- Mod(x[1]); best$angle <- abs(Arg(x[1]))
  list(best = best, runs = summary,
       gradient_check = osborne_gradient_check(do.call(cbind, starts), y, C))
}
