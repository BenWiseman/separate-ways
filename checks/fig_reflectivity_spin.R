#!/usr/bin/env Rscript
# fig_reflectivity_spin.R -- Section 4's correction, drawn.
#
# The generalized Boltzmann reflectivity of the black-mirror horizon carries the co-rotation
# term, |R| = exp(-|omega - m Omega_H| / 2 T_H). Dropping m Omega_H does not shift the answer,
# it reverses it: without the term the reflectivity FALLS with spin to 1e-6 and below; with it
# the reflectivity RISES to nearly a half. Everything here uses the same Berti-Cardoso-Will
# fit for the fundamental l = m = 2 mode that Section 4 uses.
# Arithmetic and assertions: checks/calc/horizon_reflectivity_spin.R

rp   <- function(a) 1 + sqrt(1 - a^2)
OmH  <- function(a) a / (2 * rp(a))
TH   <- function(a) sqrt(1 - a^2) / (4 * pi * rp(a))
Mw   <- function(a) 1.5251 - 1.1568 * (1 - a)^0.1292
Rful <- function(a) exp(-abs(Mw(a) - 2 * OmH(a)) / (2 * TH(a)))
Rnai <- function(a) exp(-abs(Mw(a))              / (2 * TH(a)))

need <- sqrt(1 - 1/1.14)            # |R| that delta-tau-hat_220 = 0.14 asks for
stopifnot(abs(need - 0.3504) < 1e-3, abs(Rful(0.9) - 0.4444) < 1e-3,
          abs(Rnai(0.9) - 1.03e-6) < 1e-7)
cat(sprintf("  |R| needed for the damping excess %.4f; full form at a=0.9 %.4f; naive %.2e\n",
            need, Rful(0.9), Rnai(0.9)))

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#8c8c8c"; mark <- "#a8400f"
draw <- function() {
  par(mar = c(3.9, 4.6, 1.6, 1.4), mgp = c(2.5, 0.7, 0), xpd = FALSE)
  a <- seq(0, 0.995, by = 0.0005)
  ylo <- -8.4; yhi <- 0.25
  plot(NA, xlim = c(0, 1), ylim = c(ylo, yhi), axes = FALSE,
       xlab = expression(paste("spin  ", italic(a), "/", italic(M))),
       ylab = expression(paste("amplitude reflectivity  |", italic(R), "|")))
  rect(0, ylo, 1, yhi, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1, 0.25))
  axis(2, at = seq(-8, 0, 2), labels = parse(text = sprintf("10^%d", seq(-8, 0, 2))), las = 1)

  rect(0, log10(need), 1, yhi, col = "#f6ece4", border = NA)
  segments(0, log10(need), 1, log10(need), col = mark, lwd = 1.6)
  text(0.03, log10(need) + 0.35, expression(paste("what a damping excess ", delta, hat(tau)[220],
       " = 0.14 would need")), col = mark, cex = 0.80, adj = 0)

  lines(a, log10(Rnai(a)), col = c2, lwd = 2.4, lty = 2)
  lines(a, log10(Rful(a)), col = c1, lwd = 2.6)

  # Label positions are chosen against computed curve values, not by eye: the naive
  # curve sits above y = -4.2 for every a < 0.795, so the legend block below it is clear.
  segments(0.05, -4.30, 0.11, -4.30, col = c1, lwd = 2.6)
  text(0.13, -4.30, "with the co-rotation term, as the formula has it", col = c1, cex = 0.82, adj = 0)
  segments(0.05, -4.95, 0.11, -4.95, col = c2, lwd = 2.4, lty = 2)
  text(0.13, -4.95, "without it, as the bound assumed", col = c2, cex = 0.82, adj = 0)

  points(0.9, log10(Rful(0.9)), pch = 19, cex = 1.15, col = c1)
  text(0.875, -1.30, "0.444 at a/M = 0.9, a fifth", col = c1, cex = 0.80, adj = 1)
  text(0.875, -1.72, "of the energy returned", col = c1, cex = 0.80, adj = 1)
  points(0.9, log10(Rnai(0.9)), pch = 19, cex = 1.15, col = c2)
  text(0.865, -6.60, expression(paste("10"^-6, " at the same spin")), col = c2, cex = 0.80, adj = 1)
  points(0, log10(Rful(0)), pch = 19, cex = 1.0, col = ink)
  text(0.04, -7.55, "the two agree only at zero spin, where the horizon is not rotating;",
       col = ink, cex = 0.82, adj = 0)
  text(0.04, -7.97, "away from it the term reverses the answer rather than shifting it.",
       col = ink, cex = 0.82, adj = 0)
}
for (f in c("papers/2_over_the_horizon/fig_companion_reflectivity.pdf", "papers/2_over_the_horizon/fig_companion_reflectivity.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 7.6, height = 4.4) else png(f, width = 1180, height = 690, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
