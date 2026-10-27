#!/usr/bin/env Rscript
# Full analysis: fits, inference, diagnostics, figures.
# Writes results/*.csv, figures/*.png and report/_generated/*.md (the tables the
# report includes), so every number in the report comes from this run.
#
#   Rscript run_analysis.R

suppressPackageStartupMessages({
  for (f in list.files("R", full.names = TRUE)) source(f)
})
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("report/_generated", showWarnings = FALSE, recursive = TRUE)

d <- read.table("data/set-94.dat", col.names = c("t", "y"))
t <- d$t; y <- d$y; n <- length(y)
cat(sprintf("n = %d observations, t = %g..%g\n\n", n, min(t), max(t)))

# ---- 1. Fits -----------------------------------------------------------------

g1  <- vp_grid_M1(t, y)
g1c <- vp_grid_M1c(t, y)
g1s <- vp_grid_M1s(t, y)
g1w <- vp_grid_M1w(t, y)
g2  <- vp_grid_M2(t, y)

fit1  <- lm_fit(MODELS$M1,  t, y, g1$theta)        # does not converge: see below
fit1c <- lm_fit(MODELS$M1c, t, y, g1c$theta)
fit1s <- lm_fit(MODELS$M1s, t, y, g1s$theta)
fit1w <- lm_fit(MODELS$M1w, t, y, g1w$theta)
fit2  <- lm_fit(MODELS$M2,  t, y, g2$theta)
fit3  <- ols_fit(MODELS$M3, t, y)

# Model 1 drifts towards the boundary b1 = b2: follow LM for more iterations.
m1_path <- do.call(rbind, lapply(c(10, 50, 200, 1000, 5000), function(k) {
  f <- lm_fit(MODELS$M1, t, y, g1$theta, max_iter = k, tol = 0)
  data.frame(iterations = k, a1 = f$theta[["a1"]], a2 = f$theta[["a2"]],
             b1 = f$theta[["b1"]], b2 = f$theta[["b2"]],
             rate_gap = f$theta[["b2"]] - f$theta[["b1"]], RSS = f$rss)
}))

# Cross-check: R's own Golub-Pereyra variable projection (algorithm = "plinear").
plinear <- tryCatch({
  f <- nls(y ~ cbind(1, exp(b1 * t), exp(b2 * t)), start = list(b1 = g1$theta[4], b2 = g1$theta[5]),
           algorithm = "plinear", control = nls.control(maxiter = 1000))
  list(converged = TRUE, rss = deviance(f), message = "converged")
}, error = function(e) list(converged = FALSE, rss = NA, message = conditionMessage(e)))

osb <- osborne_M1(t, y)

fits <- list(M1c = fit1c, M1s = fit1s, M1w = fit1w, M2 = fit2, M3 = fit3)

# ---- 2. Tables ---------------------------------------------------------------

estimates <- do.call(rbind, lapply(c(list(M1 = fit1), fits), function(f)
  data.frame(parameter = names(f$theta), start = unname(f$theta0), estimate = unname(f$theta))))
estimates$model <- sub("\\..*", "", rownames(estimates)); rownames(estimates) <- NULL

convergence <- do.call(rbind, lapply(names(c(list(M1 = fit1), fits)), function(k) {
  f <- c(list(M1 = fit1), fits)[[k]]
  data.frame(model = k, iterations = f$iterations, converged = f$converged,
             grad_norm = f$grad_norm, RSS = f$rss)
}))

selection <- do.call(rbind, lapply(names(fits), function(k) cbind(model = k, info_criteria(fits[[k]]))))
selection <- selection[order(selection$AIC), ]

ident <- do.call(rbind, lapply(names(c(list(M1 = fit1), fits)), function(k)
  cbind(model = k, identifiability(c(list(M1 = fit1), fits)[[k]]))))

cat("Computing confidence intervals (Wald, profile, bootstrap) ...\n")
cis <- do.call(rbind, lapply(names(fits), function(k) {
  f <- fits[[k]]
  w <- wald_ci(f); pr <- profile_ci(f, t, y); bs <- bootstrap_ci(f, t, y)
  data.frame(model = k, parameter = w$parameter, estimate = w$estimate, se = w$se,
             wald_lower = w$lower, wald_upper = w$upper,
             profile_lower = pr$lower, profile_upper = pr$upper,
             boot_lower = bs$lower, boot_upper = bs$upper)
}))
wald_m1 <- wald_ci(fit1)

tests <- residual_tests(fits)
nested <- rbind(cbind(test = "M1s (one exponential) vs M1c (confluent)", nested_f_test(fit1s, fit1c)))

write.csv(estimates,   "results/estimates.csv",   row.names = FALSE)
write.csv(convergence, "results/convergence.csv", row.names = FALSE)
write.csv(selection,   "results/model_selection.csv", row.names = FALSE)
write.csv(ident,       "results/identifiability.csv", row.names = FALSE)
write.csv(cis,         "results/confidence_intervals.csv", row.names = FALSE)
write.csv(tests,       "results/residual_tests.csv", row.names = FALSE)
write.csv(nested,      "results/nested_test.csv", row.names = FALSE)
write.csv(m1_path,     "results/model1_divergence.csv", row.names = FALSE)
write.csv(osb$runs,    "results/osborne_runs.csv", row.names = FALSE)
write.csv(data.frame(max_relative_difference = osb$gradient_check),
          "results/osborne_gradient_check.csv", row.names = FALSE)

# ---- 3. Figures --------------------------------------------------------------

ggsave("figures/fits.png", plot_fits(t, y, fits), width = 7, height = 4.6, dpi = 200, bg = "white")
ggsave("figures/residuals.png", plot_residuals(t, fits), width = 7.5, height = 10.5, dpi = 200, bg = "white")
ggsave("figures/model1_rss_surface.png", plot_m1_surface(g1$surface, fit1c$rss, fit1c$theta[["b1"]]),
       width = 5.6, height = 4.4, dpi = 200, bg = "white")

# ---- 4. Report fragments -----------------------------------------------------

# Numbers are formatted to `digits` significant figures (fixed rounding would
# print the quartic's 1.29e-06 as 0.000001).
md <- function(x, file, digits = 4, ...) {
  num <- vapply(x, is.numeric, NA)
  x[num] <- lapply(x[num], function(v) trimws(formatC(v, digits = digits, format = "g")))
  writeLines(knitr::kable(x, format = "pipe", row.names = FALSE,
                          align = ifelse(num, "r", "l"), ...),
             file.path("report/_generated", file))
}
fmt <- function(x, d = 4) trimws(formatC(x, digits = d, format = "g"))

md(estimates[, c("model", "parameter", "start", "estimate")], "estimates.md", digits = 6)
md(convergence, "convergence.md", digits = 6)
md(selection[, c("model", "p", "RSS", "sigma2", "AIC", "AICc", "BIC")], "selection.md", digits = 4)
md(ident, "identifiability.md", digits = 3)
md(m1_path, "model1_divergence.md", digits = 4)
md(wald_m1, "wald_m1.md", digits = 4)
md(cis[, c("model", "parameter", "estimate", "se", "wald_lower", "wald_upper",
           "profile_lower", "profile_upper", "boot_lower", "boot_upper")],
   "confidence_intervals.md", digits = 4,
   col.names = c("model", "param", "estimate", "SE", "Wald lo", "Wald hi",
                 "profile lo", "profile hi", "boot lo", "boot hi"))
md(tests[, c("model", "shapiro_wilk_W", "shapiro_wilk_p", "anderson_darling_p", "jarque_bera_p",
             "runs", "runs_expected", "runs_p")], "residual_tests.md", digits = 3,
   col.names = c("model", "SW W", "SW p", "AD p", "JB p", "runs", "E[runs]", "runs p"))
md(nested, "nested_test.md", digits = 4)
md(osb$runs[, c("start", "R", "iterations", "converged", "roots", "real_positive")],
   "osborne_runs.md", digits = 6)

# Simulation results (written by simulate.R; skipped if it has not been run).
sim <- list()
if (file.exists("results/simulation/selection.csv")) {
  sim_sel <- read.csv("results/simulation/selection.csv")
  sim_cov <- read.csv("results/simulation/coverage.csv")
  sim_brk <- read.csv("results/simulation/breakdown.csv")
  sim_wave <- read.csv("results/simulation/wave_test.csv")
  wide <- reshape(sim_sel, idvar = "truth", timevar = "chosen", direction = "wide")
  names(wide) <- sub("count.", "chose ", names(wide))
  wide$truth <- paste("true", wide$truth)
  md(wide, "sim_selection.md", digits = 4)
  cov_wide <- do.call(rbind, lapply(split(sim_cov, sim_cov$truth), function(g)
    data.frame(truth = paste("true", g$truth[1]), parameters = paste(g$parameter, collapse = ", "),
               coverage = paste(sprintf("%.1f%%", 100 * g$coverage), collapse = ", "),
               undefined = max(g$undefined))))
  md(cov_wide, "sim_coverage.md")
  cov_range <- function(x, d = 0) sprintf("%.*f–%.*f%%", d, 100 * min(x), d, 100 * max(x))
  lin_cov <- read.csv("results/simulation/coverage_check_M3.csv")
  md(data.frame(truth = paste("true", sim_brk$truth), `breaks down` = sprintf("%.0f%%", 100 * sim_brk$breakdown_share),
                `fake data sets` = sim_brk$reps, check.names = FALSE), "sim_breakdown.md")
  ggsave("figures/wave_test.png", plot_wave_null(read.csv("results/simulation/wave_null.csv")$improvement,
                                                 sim_wave$observed),
         width = 6, height = 3.2, dpi = 200, bg = "white")
  pct_right <- function(k) {
    g <- sim_sel[sim_sel$truth == k, ]
    sprintf("%.0f%%", 100 * g$count[g$chosen == k] / sum(g$count))
  }
  sim <- list(wave_p = sprintf("%.0f%%", 100 * sim_wave$p_value), wave_obs = fmt(sim_wave$observed, 3),
              sim_reps = sim_wave$reps,
              right_M1s = pct_right("M1s"), right_M2 = pct_right("M2"), right_M3 = pct_right("M3"),
              brk_M1s = sprintf("%.0f%%", 100 * sim_brk$breakdown_share[sim_brk$truth == "M1s"]),
              cov_M3 = cov_range(sim_cov$coverage[sim_cov$truth == "M3"]),
              cov_M3_check = cov_range(lin_cov$coverage, 1), cov_M3_check_reps = format(lin_cov$reps[1], big.mark = ","))
}

# Scalars quoted in the text.
# Selected model: the F test cannot reject M1s inside M1c, and BIC prefers M1s.
best <- if (nested$p_value > 0.05) "M1s" else "M1c"
fmt3 <- function(x) sprintf("%.3f", x)
sci_tex <- function(x) { e <- floor(log10(abs(x))); sprintf("$%.1f \\times 10^{%d}$", x / 10^e, e) }
aic <- setNames(selection$AIC, selection$model)
best_m1_aic <- min(aic[c("M1c", "M1s")])
kappa <- setNames(ident$condition_number, ident$model)
conv <- setNames(convergence$iterations, convergence$model)
grad <- setNames(convergence$grad_norm, convergence$model)
fs <- fits[[best]]$theta
sel_eq <- sprintf("\\hat y(t) = %s + %s\\, e^{%s\\, t}", fmt3(fs[["a0"]]), fmt3(fs[["a1"]]), fmt3(fs[["b1"]]))
sel_s2 <- sigma2_hat(fits[[best]])
writeLines(c("$$", sprintf("%s, \\qquad \\hat\\sigma^2 = %s .", sel_eq, fmt3(sel_s2)), "$$"),
           "report/_generated/selected_model.md")
vals <- list(
  n = n,
  m1_grid_rss = fmt(g1$rss, 7), m1c_rss = fmt(fit1c$rss, 7), m1_lm_rss = fmt(fit1$rss, 7),
  m1c_b = fmt(fit1c$theta[["b1"]]), m1s_b = fmt(fit1s$theta[["b1"]]),
  plinear = if (plinear$converged) sprintf("converged (RSS %s)", fmt(plinear$rss, 7)) else
    sprintf("fails (%s)", gsub("['\"]", "", plinear$message)),
  osb_R = fmt(osb$best$R, 7), osb_mod = fmt(osb$best$modulus), osb_ang = fmt(osb$best$angle),
  osb_roots = paste(format(osb$best$roots, digits = 4), collapse = ", "),
  osb_iter = osb$best$iterations,
  best = best, best_aic = selection$model[1],
  nested_p = fmt(nested$p_value, 3), nested_F = fmt(nested$F, 3),
  wave_rss = fmt(fit1w$rss, 7), wave_decay = fmt(exp(fit1w$theta[["rho"]]), 3),
  wave_period = fmt(2 * pi / fit1w$theta[["w"]], 3),
  wave_gain = fmt(fit1s$rss - fit1w$rss, 3),
  sel_eq = paste0("$", sel_eq, "$"), sel_sigma2 = fmt3(sel_s2), sel_sigma = fmt3(sqrt(sel_s2)),
  daic_m1 = sprintf("%.1f", abs(aic[["M1s"]] - aic[["M1c"]])),
  daic_m2 = sprintf("%.1f", aic[["M2"]] - best_m1_aic), daic_m3 = sprintf("%.0f", aic[["M3"]] - best_m1_aic),
  kappa_best = fmt(kappa[[best]], 2), kappa_m1 = sci_tex(kappa[["M1"]]),
  m1_iter = fit1$iterations,
  m1_se_amp = sci_tex(max(wald_m1$se[wald_m1$parameter %in% c("a1", "a2")])),
  m1_se_rate = fmt(max(wald_m1$se[wald_m1$parameter %in% c("b1", "b2")]), 2),
  sw_p = fmt(tests$shapiro_wilk_p[tests$model == best], 2), runs_p = fmt(tests$runs_p[tests$model == best], 1),
  iter_max = max(conv[c("M1c", "M1s", "M2")]), iter_m1w = conv[["M1w"]],
  grad_bound = sprintf("$10^{%d}$", ceiling(log10(max(grad[c("M1c", "M1s", "M2", "M1w")])))),
  osb_fd = sci_tex(osb$gradient_check), osb_starts = nrow(osb$runs),
  osb_n_conv = sum(osb$runs$converged)
)
vals <- c(vals, sim)
writeLines(sprintf("%s: %s", names(vals), unlist(vals)), "results/key_values.txt")
yaml_str <- function(v) gsub('"', '\\"', gsub("\\", "\\\\", v, fixed = TRUE), fixed = TRUE)
writeLines(paste0(names(vals), ": \"", yaml_str(unlist(vals)), "\""), "report/_variables.yml")

saveRDS(list(fits = c(list(M1 = fit1), fits), osborne = osb, grids = list(g1 = g1, g1c = g1c, g1s = g1s, g1w = g1w, g2 = g2)),
        "results/fits.rds")

# ---- 5. Console summary ------------------------------------------------------

cat("\nConvergence\n"); print(convergence, digits = 6, row.names = FALSE)
cat("\nModel 1 drift towards b1 = b2\n"); print(m1_path, digits = 5, row.names = FALSE)
cat("\nplinear cross-check:", vals$plinear, "\n")
cat("\nOsborne best stationary point: R =", vals$osb_R, " roots", vals$osb_roots, "\n")
cat("\nModel selection\n"); print(selection, digits = 5, row.names = FALSE)
cat("\nIdentifiability\n"); print(ident, digits = 3, row.names = FALSE)
cat("\nNested test\n"); print(nested, digits = 4, row.names = FALSE)
cat("\nResidual tests\n"); print(tests, digits = 3, row.names = FALSE)
cat("\nConfidence intervals\n"); print(cis, digits = 4, row.names = FALSE)

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
