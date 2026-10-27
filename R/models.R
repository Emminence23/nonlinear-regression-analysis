# Model definitions: mean functions and analytic Jacobians d mu / d theta.
#
#   M1   bi-exponential   mu = a0 + a1 exp(b1 t) + a2 exp(b2 t)          theta = (a0, a1, a2, b1, b2)
#   M1c  confluent        mu = a0 + (a1 + g t) exp(b1 t)                  theta = (a0, a1, g, b1)
#   M1w  decaying wave    mu = a0 + e^{rho t} (c cos(w t) + d sin(w t))  theta = (a0, c, d, rho, w)
#   M1s  mono-exponential mu = a0 + a1 exp(b1 t)                         theta = (a0, a1, b1)
#   M2   rational         mu = (a0 + a1 t) / (b0 + b1 t),  b0 := 1        theta = (a0, a1, b1)
#   M3   quartic          mu = sum_k b_k t^k, k = 0..4                    theta = (b0, ..., b4)
#
# Model 2 as stated in the brief is not identifiable: (a0, a1, b0, b1) and
# c * (a0, a1, b0, b1) give the same curve for any c != 0, so the Fisher
# information matrix is singular. We fix b0 = 1 (valid because b0 != 0 for
# any curve that is finite at t = 0) and estimate the remaining three
# parameters. Estimates under any other normalisation follow by rescaling.

MODELS <- list(
  M1 = list(
    label = "Model 1: bi-exponential",
    pars  = c("a0", "a1", "a2", "b1", "b2"),
    mu = function(th, t) th[1] + th[2] * exp(th[4] * t) + th[3] * exp(th[5] * t),
    jac = function(th, t) {
      e1 <- exp(th[4] * t); e2 <- exp(th[5] * t)
      cbind(a0 = 1, a1 = e1, a2 = e2, b1 = th[2] * t * e1, b2 = th[3] * t * e2)
    }
  ),
  M1s = list(
    label = "Model 1 (reduced): mono-exponential",
    pars  = c("a0", "a1", "b1"),
    mu = function(th, t) th[1] + th[2] * exp(th[3] * t),
    jac = function(th, t) {
      e1 <- exp(th[3] * t)
      cbind(a0 = 1, a1 = e1, b1 = th[2] * t * e1)
    }
  ),
  # Limit of M1 as b2 -> b1 with a1 = -a2 -> infinity (a1 (b2 - b1) -> g):
  # the closure of the M1 family, where the M1 infimum of the RSS is attained.
  M1c = list(
    label = "Model 1 (boundary): confluent exponential",
    pars  = c("a0", "a1", "g", "b1"),
    mu = function(th, t) th[1] + (th[2] + th[3] * t) * exp(th[4] * t),
    jac = function(th, t) {
      e <- exp(th[4] * t)
      cbind(a0 = 1, a1 = e, g = t * e, b1 = (th[2] + th[3] * t) * t * e)
    }
  ),
  # Decaying wave: Model 1 with a complex-conjugate pair of rates rho +- i w.
  # Osborne's algorithm finds the unconstrained optimum of Model 1 here.
  M1w = list(
    label = "Decaying wave",
    pars  = c("a0", "c", "d", "rho", "w"),
    mu = function(th, t) th[1] + exp(th[4] * t) * (th[2] * cos(th[5] * t) + th[3] * sin(th[5] * t)),
    jac = function(th, t) {
      e <- exp(th[4] * t); co <- cos(th[5] * t); si <- sin(th[5] * t)
      wave <- th[2] * co + th[3] * si
      cbind(a0 = 1, c = e * co, d = e * si, rho = t * e * wave,
            w = t * e * (th[3] * co - th[2] * si))
    }
  ),
  M2 = list(
    label = "Model 2: rational",
    pars  = c("a0", "a1", "b1"),
    mu = function(th, t) (th[1] + th[2] * t) / (1 + th[3] * t),
    jac = function(th, t) {
      d <- 1 + th[3] * t
      cbind(a0 = 1 / d, a1 = t / d, b1 = -t * (th[1] + th[2] * t) / d^2)
    }
  ),
  M3 = list(
    label = "Model 3: quartic polynomial",
    pars  = paste0("b", 0:4),
    mu = function(th, t) drop(outer(t, 0:4, `^`) %*% th),
    jac = function(th, t) {
      X <- outer(t, 0:4, `^`); colnames(X) <- paste0("b", 0:4); X
    }
  )
)
