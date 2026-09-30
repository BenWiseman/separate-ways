# The two ends of the construction, and the one frequency where they meet.
#
# Left: what the crossing at the bang selects. Theta-invariance forces half-and-half occupancy, so
# the pair-block floor is n_*(P) = (1 - sqrt(1-P))/2, and since a squeeze carries <N> = sinh^2 r,
# the selected occupation is a selected squeeze r = arcsinh(sqrt(n_*)). It tops out at P = 1.
#
# Right: what equilibrium at a horizon fixes. The fold's map is the half-period thermal shift, so
# the cross-sheet correlator is the direct one shifted by i beta/2, the state is the thermofield
# double, and tanh r = exp(-beta omega/2).
#
# The two are the same number at beta omega = ln 3, because sinh r = 1/sqrt2 gives cosh r =
# sqrt(3/2) and hence tanh r = 1/sqrt3. Nothing was arranged to make that happen.

n_star <- function(P) (1 - sqrt(1 - P)) / 2
r_bang <- function(P) asinh(sqrt(n_star(P)))
r_hor  <- function(bw) atanh(exp(-bw / 2))

ceiling_r <- asinh(1 / sqrt(2))
cross_bw  <- log(3)

cat("=== the two ends ===\n\n")
cat(sprintf("  bang ceiling  arcsinh(1/sqrt2)      = %.12f\n", ceiling_r))
cat(sprintf("  horizon value at beta*omega = ln 3  = %.12f\n", r_hor(cross_bw)))
cat(sprintf("  difference                          = %.3e\n", abs(ceiling_r - r_hor(cross_bw))))
stopifnot(abs(ceiling_r - r_hor(cross_bw)) < 1e-12,
          abs(r_bang(1) - ceiling_r) < 1e-12,
          abs(sinh(ceiling_r)^2 - 0.5) < 1e-12)
cat("  planted: a horizon value at beta*omega = 1 instead gives %.6f, which misses it by %.4f\n" |>
      sprintf(r_hor(1), abs(r_hor(1) - ceiling_r)))

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.0, 4.4, 2.2, 1.0), mgp = c(2.5, 0.7, 0), xpd = FALSE)

  Ps <- seq(0.001, 1, length.out = 400)
  plot(NA, xlim = c(0, 1), ylim = c(0, 0.75), axes = FALSE,
       xlab = expression(paste("crossing probability  ", italic(P))),
       ylab = expression(paste("selected squeeze  ", italic(r))))
  rect(0, 0, 1, 0.75, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1, 0.25)); axis(2, at = seq(0, 0.75, 0.25), las = 1)
  abline(h = ceiling_r, col = grey, lty = 2, lwd = 1.8)
  lines(Ps, r_bang(Ps), col = c1, lwd = 2.8)
  points(1, ceiling_r, pch = 19, cex = 1.0, col = c1)
  text(0.40, ceiling_r + 0.045, expression(paste("arcsinh(1/", sqrt(2), ") = 0.6585")),
       col = grey, cex = 0.86, adj = 0)
  text(0.20, 0.09, expression(paste(italic(r), " = arcsinh ", sqrt(italic(n)["*"]))),
       col = c1, cex = 0.90, adj = 0)
  mtext("what the bang selects", side = 3, line = 0.6, cex = 0.92, col = ink)

  bws <- seq(0.35, 6, length.out = 400)
  plot(NA, xlim = c(0, 6), ylim = c(0, 1.6), axes = FALSE,
       xlab = expression(paste("mode frequency  ", beta * omega)),
       ylab = expression(paste("horizon squeeze  ", italic(r))))
  rect(0, 0, 6, 1.6, col = "#fbfbfa", border = NA)
  axis(1, at = 0:6); axis(2, at = seq(0, 1.5, 0.5), las = 1)
  abline(h = ceiling_r, col = grey, lty = 2, lwd = 1.8)
  segments(cross_bw, 0, cross_bw, ceiling_r, col = grey, lty = 3, lwd = 1.5)
  lines(bws, r_hor(bws), col = c2, lwd = 2.8)
  points(cross_bw, ceiling_r, pch = 19, cex = 1.1, col = ink)
  text(cross_bw + 0.18, ceiling_r + 0.20, expression(paste(beta * omega, " = ln 3")),
       col = ink, cex = 0.90, adj = 0)
  text(2.5, 1.18, expression(paste("tanh ", italic(r), " = ", e^{-beta * omega / 2})),
       col = c2, cex = 0.92, adj = 0)
  mtext("what a horizon fixes, and where they agree", side = 3, line = 0.6, cex = 0.92, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("papers/2_over_the_horizon/fig_companion_two_ends.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.6, height = 4.3, pointsize = 12)
  else png(out, width = 9.6, height = 4.3, units = "in", res = 150, pointsize = 12)
  draw(); dev.off(); cat(sprintf("  wrote %s\n", out))
}
