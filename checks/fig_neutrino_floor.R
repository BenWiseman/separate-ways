#!/usr/bin/env Rscript
# fig_neutrino_floor.R -- Section 7: what the matter rule does to the neutrino spectrum.
#
# When the rule is exact the sterile partner cannot decay, its Yukawa column is forbidden, and
# the light mass matrix has rank at most two. The two measured splittings then leave one light
# neutrino EXACTLY massless, so for normal ordering the sum has no free parameter:
#     Sigma m_nu = sqrt(Dm21) + sqrt(Dm31) = 58.78 +- 0.32 meV.
# The weakly broken version adds a small gauge-invariant Yukawa column and lifts the zero, but
# for the example lifetime of 1e28 s the lift is 2.1e-55 eV, which is fifty-three orders of
# magnitude below the lighter splitting. The prediction of a massless state is not fragile.
# Arithmetic: checks/calc/weak_breaking_neutrino_bound.R

d21 <- 7.53e-5; s21 <- 0.18e-5
d31 <- 2.510e-3; s31 <- 0.030e-3
m2 <- sqrt(d21); m3 <- sqrt(d31)
sum_meV <- (m2 + m3) * 1e3
err_meV <- 1e3 * sqrt((s21/(2*m2))^2 + (s31/(2*m3))^2)
hbar <- 6.582119569e-25; v <- 246.22; M1 <- 4.916e8; tau <- 1e28
m1_broken <- 4*pi*v^2*hbar/(M1^2*tau) * 1e9          # eV
cat("=== the spectrum the matter rule forces ===\n\n")
cat(sprintf("  m_2 = %.4f meV,  m_3 = %.3f meV\n", m2*1e3, m3*1e3))
cat(sprintf("  Sigma m_nu = %.2f +- %.2f meV   (paper: 58.78 +- 0.32)\n", sum_meV, err_meV))
cat(sprintf("  m_1 = 0 exactly when the rule is exact; at most %.3e eV when weakly broken\n", m1_broken))
cat(sprintf("  that is %.0f orders of magnitude below m_2\n", log10(m2/m1_broken)))
stopifnot(abs(sum_meV - 58.78) < 0.05, abs(err_meV - 0.32) < 0.05,
          abs(m1_broken - 2.075e-55) < 5e-58)

ink <- "grey15"; c1 <- "#1f4e79"; mark <- "#a8400f"; grey <- "#8a97a4"
draw <- function() {
  # Two panels, because 8.68 and 50.1 meV are three quarters of a decade apart on an axis
  # fifty-seven decades wide and land on top of each other. Left, the range; right, the zoom.
  par(mfrow = c(1, 2), mar = c(3.9, 1.2, 2.6, 1.2), mgp = c(2.5, 0.7, 0), xpd = FALSE)

  lo <- -57; hi <- 0.5
  plot(NA, xlim = c(lo, hi), ylim = c(0, 1), axes = FALSE, xlab = "neutrino mass  (eV)", ylab = "")
  rect(lo, 0, hi, 1, col = "#fbfbfa", border = NA)
  at <- seq(-55, 0, 10); axis(1, at = at, labels = parse(text = sprintf("10^%d", at)))
  segments(lo, 0.50, hi, 0.50, col = "grey75", lwd = 1)
  rect(log10(m2), 0.455, log10(m3), 0.545, col = c1, border = NA)
  points(log10(m1_broken), 0.50, pch = 19, cex = 1.6, col = mark)
  text(log10(m3) + 0.6, 0.66, "the two", col = c1, cex = 0.82, adj = 1)
  text(log10(m3) + 0.6, 0.59, "measured", col = c1, cex = 0.82, adj = 1)
  text(log10(m1_broken) + 1.5, 0.87, "zero when the rule is exact;", col = mark, cex = 0.82, adj = 0)
  text(log10(m1_broken) + 1.5, 0.80, "at most", col = mark, cex = 0.82, adj = 0)
  text(log10(m1_broken) + 1.5, 0.73, expression(paste("2.1 ", {} %*% {}, " 10"^-55, " eV when")),
       col = mark, cex = 0.82, adj = 0)
  text(log10(m1_broken) + 1.5, 0.66, "weakly broken", col = mark, cex = 0.82, adj = 0)
  arrows(log10(m1_broken) + 1.0, 0.30, log10(m2) - 1.0, 0.30, code = 3, length = 0.06,
         col = grey, lwd = 1.2)
  text((log10(m1_broken) + log10(m2))/2, 0.375, "fifty-three orders of magnitude",
       col = grey, cex = 0.84, adj = 0.5)
  text((log10(m1_broken) + log10(m2))/2, 0.205, "so it stays massless for anything",
       col = grey, cex = 0.82, adj = 0.5)
  text((log10(m1_broken) + log10(m2))/2, 0.135, "an experiment can reach",
       col = grey, cex = 0.82, adj = 0.5)
  mtext("one state is left massless", side = 3, line = 0.8, cex = 0.88, col = ink)

  plot(NA, xlim = c(-3, 63), ylim = c(0, 1), axes = FALSE, xlab = "neutrino mass  (meV)", ylab = "")
  rect(-3, 0, 63, 1, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 60, 20))
  segments(-3, 0.50, 63, 0.50, col = "grey75", lwd = 1)
  for (m in c(0, m2*1e3, m3*1e3)) {
    col <- if (m == 0) mark else c1
    points(m, 0.50, pch = 19, cex = 1.6, col = col); segments(m, 0.50, m, 0.63, col = col, lty = 3) }
  text(0,      0.695, expression(italic(m)[1]), col = mark, cex = 0.95, adj = 0.5)
  text(m2*1e3, 0.695, expression(italic(m)[2]), col = c1, cex = 0.95, adj = 0.5)
  text(m3*1e3, 0.695, expression(italic(m)[3]), col = c1, cex = 0.95, adj = 0.5)
  text(0,      0.785, "0", col = mark, cex = 0.82, adj = 0.5)
  text(m2*1e3, 0.785, "8.68", col = c1, cex = 0.82, adj = 0.5)
  text(m3*1e3, 0.785, "50.1", col = c1, cex = 0.82, adj = 0.5)
  # a stacked bar, not an arrow: the sum is m_2 + m_3 laid end to end, and an arrow spanning
  # 0 to 58.78 would read as a range on the axis instead.
  rect(0, 0.265, m2*1e3, 0.335, col = "#4a78a8", border = NA)
  rect(m2*1e3, 0.265, sum_meV, 0.335, col = c1, border = NA)
  segments(sum_meV, 0.245, sum_meV, 0.355, col = ink, lwd = 1.4)
  text(sum_meV/2, 0.375, expression(paste(Sigma, italic(m)[nu], " = 58.78 ", {} %+-% {}, " 0.32 meV")),
       col = ink, cex = 0.86, adj = 0.5)
  text(sum_meV/2, 0.195, expression(paste(italic(m)[2], " and ", italic(m)[3], " laid end to end: the two")),
       col = ink, cex = 0.82, adj = 0.5)
  text(sum_meV/2, 0.120, "splittings fix the sum once the third mass is zero",
       col = ink, cex = 0.82, adj = 0.5)
  mtext("and the sum then has nothing left to choose", side = 3, line = 0.8, cex = 0.88, col = ink)
}
for (f in c("paper/fig_companion_neutrino.pdf", "paper/fig_companion_neutrino.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 9.0, height = 3.7) else png(f, width = 1390, height = 575, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
