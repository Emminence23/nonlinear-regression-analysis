#!/usr/bin/env Rscript
# Simulation checks (slower than run_analysis.R; results are cached as CSV).
#
#   Rscript simulate.R            # 1000 data sets per check, all cores
#
# Writes results/simulation/*.csv, which run_analysis.R turns into report
# tables and figures. Reproducible: fixed seed and one random-number stream
# per data set, so any number of cores gives the same results
# (CORES=2 Rscript simulate.R to limit the cores).

suppressPackageStartupMessages({
  for (f in list.files("R", full.names = TRUE)) source(f)
  library(parallel)
})
REPS  <- 1000
CORES <- as.integer(Sys.getenv("CORES", max(1, detectCores())))
RNGkind("L'Ecuyer-CMRG")
dir.create("results/simulation", showWarnings = FALSE, recursive = TRUE)

d <- read.table("data/set-94.dat", col.names = c("t", "y"))
t <- d$t; y <- d$y
real <- fit_candidates(t, y, two_exp = FALSE)

# Data set i always uses random-number stream i, so the results do not depend
# on how many cores share the work.
run <- function(f, reps = REPS, seed = 686) {
  set.seed(seed)
  seeds <- vector("list", reps)
  seeds[[1]] <- .Random.seed
  for (i in seq_len(reps - 1)) seeds[[i + 1]] <- nextRNGStream(seeds[[i]])
  mclapply(seq_len(reps), function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    f()
  }, mc.cores = CORES)
}

# ---- 1. Is the wave real? -------------------------------------------------------

wave <- lm_fit(MODELS$M1w, t, y, vp_grid_M1w(t, y)$theta)
observed <- real$M1s$rss - wave$rss
cat(sprintf("Wave test: %d data sets ...\n", REPS))
null <- unlist(run(function() wave_once(real$M1s, t)))
write.csv(data.frame(improvement = null), "results/simulation/wave_null.csv", row.names = FALSE)
write.csv(data.frame(observed = observed, p_value = mean(null >= observed), reps = REPS),
          "results/simulation/wave_test.csv", row.names = FALSE)
cat(sprintf("  observed improvement %.4f, p = %.3f\n", observed, mean(null >= observed)))

# ---- 2. Selection, breakdown, interval coverage -------------------------------

selection <- list(); coverage <- list(); breakdown <- list()
for (name in c("M1s", "M2", "M3")) {
  cat(sprintf("Truth %s: %d data sets ...\n", name, REPS))
  out <- run(function() simulate_once(real[[name]], name, t))
  ok <- !vapply(out, inherits, NA, "try-error")
  out <- out[ok]
  chosen <- factor(vapply(out, `[[`, "", "chosen"), levels = c("M1s", "M1c", "M2", "M3"))
  selection[[name]] <- data.frame(truth = name, chosen = levels(chosen), count = as.vector(table(chosen)))
  breakdown[[name]] <- data.frame(truth = name, breakdown_share = mean(vapply(out, `[[`, NA, "breakdown")),
                                  reps = length(out))
  cov <- do.call(rbind, lapply(out, `[[`, "covered"))
  # A few fits can give an undefined interval (NaN standard error); count
  # those separately rather than as misses.
  coverage[[name]] <- data.frame(truth = name, parameter = colnames(cov),
                                 coverage = colMeans(cov, na.rm = TRUE),
                                 undefined = colSums(is.na(cov)), reps = nrow(cov), row.names = NULL)
}
write.csv(do.call(rbind, selection), "results/simulation/selection.csv", row.names = FALSE)
write.csv(do.call(rbind, coverage),  "results/simulation/coverage.csv",  row.names = FALSE)
write.csv(do.call(rbind, breakdown), "results/simulation/breakdown.csv", row.names = FALSE)

# ---- 3. Exact check for the linear model ---------------------------------------

# The quartic is linear, so its 95% Wald intervals cover exactly 95%. With
# 1000 data sets the estimated coverage still varies by about +-1.4 points,
# so check it once more on 10,000 fresh data sets (cheap: no nonlinear fits).
cat("Quartic coverage check: 10000 data sets ...\n")
lin <- run(function() {
  ystar <- real$M3$fitted + rnorm(length(t), 0, sqrt(sigma2_hat(real$M3)))
  ci <- wald_ci(ols_fit(MODELS$M3, t, ystar))
  ci$lower <= real$M3$theta & real$M3$theta <= ci$upper
}, reps = 10000, seed = 687)
write.csv(data.frame(parameter = names(real$M3$theta), coverage = colMeans(do.call(rbind, lin)), reps = 10000),
          "results/simulation/coverage_check_M3.csv", row.names = FALSE)

print(do.call(rbind, selection), row.names = FALSE)
print(do.call(rbind, coverage), row.names = FALSE, digits = 3)
print(do.call(rbind, breakdown), row.names = FALSE, digits = 3)
