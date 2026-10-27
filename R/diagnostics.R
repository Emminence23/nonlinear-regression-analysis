# Residual diagnostics (Q6, Q7) and figures (Q8).

library(ggplot2)
library(patchwork)

# Fixed categorical order (validated reference palette, slots 1-5) plus a
# distinct linetype per model so identity never rests on colour alone.
MODEL_COLOURS   <- c(M1c = "#2a78d6", M1s = "#eb6834", M2 = "#1baf7a", M3 = "#eda100", M1w = "#e87ba4")
MODEL_LINETYPES <- c(M1c = "solid",   M1s = "longdash", M2 = "dotdash", M3 = "dotted",  M1w = "twodash")
MODEL_SHORT     <- c(M1c = "Model 1 (confluent limit)", M1s = "Model 1 (one exponential)",
                     M2 = "Model 2 (rational)", M3 = "Model 3 (quartic)", M1w = "Decaying wave")

theme_report <- function() {
  theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom", legend.title = element_blank())
}

# Jarque-Bera: n/6 (S^2 + (K-3)^2 / 4) ~ chi^2_2 asymptotically.
jarque_bera <- function(e) {
  n <- length(e); z <- e - mean(e); m2 <- mean(z^2)
  S <- mean(z^3) / m2^1.5; K <- mean(z^4) / m2^2
  stat <- n / 6 * (S^2 + (K - 3)^2 / 4)
  list(statistic = stat, p.value = pchisq(stat, 2, lower.tail = FALSE))
}

# Wald-Wolfowitz runs test on residual signs: too few runs means systematic
# lack of fit (residuals clustered above/below zero along t).
runs_test <- function(e) {
  s <- e[e != 0] > 0
  n1 <- sum(s); n2 <- sum(!s); n <- n1 + n2
  runs <- 1 + sum(s[-1] != s[-n])
  mu <- 2 * n1 * n2 / n + 1
  v <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
  z <- (runs - mu) / sqrt(v)
  list(runs = runs, expected = mu, z = z, p.value = 2 * pnorm(-abs(z)))
}

residual_tests <- function(fits) {
  do.call(rbind, lapply(names(fits), function(k) {
    e <- fits[[k]]$resid
    data.frame(model = k,
               shapiro_wilk_W = unname(shapiro.test(e)$statistic),
               shapiro_wilk_p = shapiro.test(e)$p.value,
               anderson_darling_p = nortest::ad.test(e)$p.value,
               jarque_bera_p = jarque_bera(e)$p.value,
               runs = runs_test(e)$runs, runs_expected = runs_test(e)$expected,
               runs_p = runs_test(e)$p.value)
  }))
}

plot_fits <- function(t, y, fits) {
  grid <- seq(min(t), max(t), length.out = 400)
  curves <- do.call(rbind, lapply(names(fits), function(k)
    data.frame(model = k, t = grid, mu = fits[[k]]$model$mu(fits[[k]]$theta, grid))))
  curves$model <- factor(curves$model, levels = names(fits))
  ggplot() +
    geom_point(aes(t, y), data = data.frame(t, y), colour = "grey25", size = 1.4) +
    geom_line(aes(t, mu, colour = model, linetype = model), data = curves, linewidth = 0.7) +
    scale_colour_manual(values = MODEL_COLOURS, labels = MODEL_SHORT) +
    scale_linetype_manual(values = MODEL_LINETYPES, labels = MODEL_SHORT) +
    labs(x = "t", y = "y(t)", title = "Observed data and least-squares fits") +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    theme_report()
}

plot_residuals <- function(t, fits) {
  panels <- lapply(names(fits), function(k) {
    f <- fits[[k]]
    d <- data.frame(t = t, fitted = f$fitted, r = f$resid / sqrt(sigma2_hat(f)))
    col <- MODEL_COLOURS[[k]]
    p1 <- ggplot(d, aes(t, r)) + geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
      geom_point(colour = col, size = 1.2) + labs(x = "t", y = "std. residual", title = MODEL_SHORT[[k]]) +
      theme_report()
    p2 <- ggplot(d, aes(fitted, r)) + geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
      geom_point(colour = col, size = 1.2) + labs(x = "fitted", y = NULL, title = " ") + theme_report()
    p3 <- ggplot(d, aes(sample = r)) + stat_qq_line(colour = "grey50", linewidth = 0.3) +
      stat_qq(colour = col, size = 1.2) + labs(x = "normal quantile", y = NULL, title = " ") + theme_report()
    p1 | p2 | p3
  })
  wrap_plots(panels, ncol = 1)
}

# Excess profiled RSS over the two rates of Model 1, relative to the
# infimum rss_inf (log10 scale). The valley floor runs into the diagonal
# b1 = b2, where Model 1 degenerates into the confluent model.
plot_m1_surface <- function(surface, rss_inf, limit_rate) {
  s <- surface
  s$excess <- log10(pmax(s$rss - rss_inf, 1e-4))
  s$x <- log10(-s$b1); s$y <- log10(-s$b2)
  ggplot(s, aes(x, y)) +
    geom_raster(aes(fill = excess), interpolate = TRUE) +
    geom_contour(aes(z = excess), breaks = c(-3, -2, -1, 0), colour = "white", linewidth = 0.25, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.4, linetype = "dashed") +
    annotate("point", x = log10(-limit_rate), y = log10(-limit_rate), colour = "white", fill = "#eb6834",
             shape = 21, size = 2.8, stroke = 0.8) +
    annotate("text", x = log10(-limit_rate) - 0.08, y = log10(-limit_rate) + 0.12, hjust = 1, size = 3,
             label = "infimum: b1 = b2") +
    scale_x_continuous(breaks = -3:0, labels = c("0.001", "0.01", "0.1", "1")) +
    scale_y_continuous(breaks = -3:0, labels = c("0.001", "0.01", "0.1", "1")) +
    scale_fill_gradient(low = "#0d366b", high = "#cde2fb", name = expression(log[10]~(RSS - RSS[inf])),
                        breaks = -4:0, labels = c("\u2264 -4", -3:0)) +
    labs(x = expression(-b[1]~"(faster rate)"), y = expression(-b[2]~"(slower rate)"),
         title = "Model 1: excess RSS profiled over the two rates") +
    coord_equal() +
    theme_report() + theme(legend.position = "right", legend.title = element_text(size = 8))
}

# How much the wave model improves the RSS on data simulated from the simple
# model, versus the improvement seen on the real data.
plot_wave_null <- function(null, observed) {
  ggplot(data.frame(x = null), aes(x)) +
    geom_histogram(bins = 40, fill = "#86b6ef", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = observed, colour = "#eb6834", linewidth = 0.8) +
    annotate("text", x = observed + 0.02, y = Inf, vjust = 1.5, hjust = 0, size = 3,
             label = sprintf("real data: %.3f\n%.0f%% of fake data sets do better",
                             observed, 100 * mean(null >= observed))) +
    labs(x = "RSS improvement of the wave over one exponential",
         y = "fake data sets", title = "Is the wave real? Improvements on data with no wave") +
    theme_report()
}
