#!/usr/bin/env Rscript
# fig_two_sided_growth.R -- Section 9: what feeding from both sheets actually changes.
#
# A hole with a past region takes matter from our side and antimatter from the mirror. If the
# mirror sheet supplies accretion at T times our own rate, with T the fraction of its infall the
#     seam passes, M(t) = M_seed exp[ (1+T)(t - t_i) / t_S ],   t_S = 450 Myr * eps/(1-eps),
# so the mirror sheet does not add to the MASS, it multiplies the EXPONENT. Over the window
# the overmassive early holes occupy, z = 20 to z = 7, that is the difference between 11.6
# e-folds and 23.3, and a factor of 1.1e5 in final mass.
#
# The right panel is the caveat and belongs in the figure, not only in the text: T is not
# derived anywhere in this paper. T = 1 is a seam fully open to matter. What does not depend on
# its value is the direction, since any T > 0 raises the exponent.
# Arithmetic and assertions: checks/calc/two_sided_growth.R

H0 <- 67.36; Om <- 0.3153; OL <- 0.6847
invH0_Gyr <- 9.77792 / (H0/100)
age <- function(z) (2/3) * invH0_Gyr / sqrt(OL) * asinh(sqrt(OL/Om) * (1+z)^(-1.5))
eps <- 0.1; tS <- 0.450 * eps/(1-eps)              # Gyr
z_i <- 20; z_f <- 7; seed <- 100
t_i <- age(z_i); t_f <- age(z_f); dt <- t_f - t_i
efold <- function(k) k * dt / tS
mass  <- function(z, k) seed * exp(k * (age(z) - t_i) / tS)

stopifnot(abs(t_i - 0.179) < 2e-3, abs(t_f - 0.761) < 2e-3, abs(tS*1e3 - 50.1) < 0.3,
          abs(efold(1) - 11.63) < 0.05, abs(efold(2) - 23.26) < 0.1)
cat("=== growth over z = 20 to 7, seed 100 Msun, eps = 0.1 ===\n\n")
cat(sprintf("  window %.4f Gyr, Salpeter time %.1f Myr\n", dt, tS*1e3))
cat(sprintf("  one-sided: %.2f e-folds, final %.3e Msun\n", efold(1), mass(z_f, 1)))
cat(sprintf("  two-sided: %.2f e-folds, final %.3e Msun, ratio %.2e\n",
            efold(2), mass(z_f, 2), mass(z_f, 2)/mass(z_f, 1)))
kneed <- function(x) 1 + log(x)/efold(1)
for (x in c(10, 100, 1000))
  cat(sprintf("  a factor %5.0f over the same window needs T = %.3f\n", x, kneed(x) - 1))
stopifnot(abs(kneed(100) - 1 - 0.396) < 2e-3)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"
draw <- function() {
  par(mfrow = c(1, 2), mar = c(3.9, 4.4, 2.1, 1.0), mgp = c(2.5, 0.7, 0), xpd = FALSE)

  zs <- seq(z_i, z_f, length.out = 400)
  plot(NA, xlim = c(z_i, z_f), ylim = c(2, 12.6), axes = FALSE,
       xlab = "redshift", ylab = expression(paste("black-hole mass  (", italic(M)[sun], ")")))
  rect(z_i, 2, z_f, 12.6, col = "#fbfbfa", border = NA)
  axis(1, at = c(20, 16, 12, 10, 8, 7)); axis(2, at = seq(2, 12, 2),
       labels = parse(text = sprintf("10^%d", seq(2, 12, 2))), las = 1)
  rect(z_i, 7, z_f, 9, col = "#eceff2", border = NA)
  text(19.8, 7.35, "where JWST's overmassive holes sit", col = grey, cex = 0.80, adj = 0)

  lines(zs, log10(mass(zs, 1)), col = c1, lwd = 2.6)
  lines(zs, log10(mass(zs, 2)), col = c2, lwd = 2.6)
  points(z_i, 2, pch = 19, cex = 1.1, col = ink)

  segments(19.85, 12.25, 19.15, 12.25, col = c2, lwd = 2.6)
  text(18.95, 12.25, "both sheets: 23.3 e-folds", col = c2, cex = 0.82, adj = 0)
  segments(19.85, 11.35, 19.15, 11.35, col = c1, lwd = 2.6)
  text(20.10, 11.35, "one sheet: 11.6 e-folds", col = c1, cex = 0.82, adj = 0)
  mtext("the exponent doubles, not the rate", side = 3, line = 0.6, cex = 0.88, col = ink)

  ps <- seq(0, 1, length.out = 400)
  plot(NA, xlim = c(0, 1), ylim = c(0, 2.35), axes = FALSE,
       xlab = expression(paste("fraction ", italic(p), " of each stream that crosses the seam")),
       ylab = expression(paste("what our hole receives, in units of  ", italic(F))))
  rect(0, 0, 1, 2.35, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1, 0.25))
  axis(2, at = seq(0, 2, 0.5), las = 1)
  lines(ps, 1 + ps, col = grey, lwd = 2.2, lty = 2)
  lines(ps, (1 - ps) + ps, col = c2, lwd = 2.8)
  lines(ps, 1 - ps, col = c1, lwd = 1.6, lty = 3)
  lines(ps, ps, col = c1, lwd = 1.6, lty = 4)
  text(0.05, 1.93, "ours kept and theirs added,", col = grey, cex = 0.78, adj = 0)
  text(0.05, 1.79, "which creates 2pF from nothing", col = grey, cex = 0.78, adj = 0)
  text(0.31, 0.86, "swap: (1-p)F + pF = F", col = c2, cex = 0.82, adj = 0)
  text(0.13, 0.60, "(1-p)F ours", col = c1, cex = 0.80, adj = 0)
  # two diagonals cross this panel and the label meets one or the other wherever it goes
  source("checks/fig_label.R")
  lab_on(0.78, 0.60, "pF theirs", bg = "#fbfbfa", col = c1, cex = 0.80, adj = 1)
  points(c(0, 1), c(1, 1), pch = 19, cex = 0.9, col = c2)
  mtext("and the sum that takes it away again", side = 3, line = 0.6,
        cex = 0.88, col = ink)
}
for (f in c("paper/fig_companion_growth.pdf", "paper/fig_companion_growth.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 8.6, height = 4.2) else png(f, width = 1330, height = 650, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
