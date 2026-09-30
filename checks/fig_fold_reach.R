# How far a point is from its own fold image, and how it reaches its floor.
#
# Section 5.4's field, drawn. Left: the separation s(r) across the exterior, against the floor
# pi r_h that the angular half-turn costs. Right: the approach to that floor, which is the local
# clock rate squared, flattening onto the closed-form coefficient.
#
# The integrals are redone here in R rather than read from the Python that computes them for the
# text, so the figure is an independent check of its own numbers as well as a picture of them.
# On the t=0 slice the bridge is dl^2 + r(l)^2 dOmega^2 with dl = dr/sqrt(1 - rh/r); the
# substitution r = rh/(1-u^2) removes the endpoint singularity exactly, giving
# dl = 2 rh du/(1-u^2)^2. Clairaut conserves L = r sin(psi), so
#     ds = dl/sqrt(1 - L^2/r^2),  dtheta = L dl/(r^2 sqrt(1 - L^2/r^2)).

rh <- 1

# Near the throat the Clairaut factor sqrt(1 - L^2/r^2) collapses like sqrt(delta + u^2) with
# delta = 1 - L/rh, so the integrand has a boundary layer of width sqrt(delta). Substituting
# u = A sinh(v) with A = sqrt(delta) absorbs it exactly and leaves a bounded integrand: near
# v = 0 the ratio du/dv over the Clairaut factor tends to 1/sqrt(2).
halves <- function(r0, L) {
  A <- sqrt(max(1 - L / rh, 1e-14))
  u0 <- sqrt(1 - rh / r0)
  vmax <- asinh(u0 / A)
  parts <- function(v) {
    u <- A * sinh(v)
    r <- rh / (1 - u^2)
    w <- sqrt(pmax(1 - L^2 / r^2, 0))
    dl <- 2 * rh / (1 - u^2)^2
    dudv <- A * cosh(v)
    list(ang = L * dl * dudv / (r^2 * w), len = dl * dudv / w)
  }
  # Simpson on a smooth bounded integrand, so there is no adaptive routine to fail and the
  # accuracy is checked against the Python the text quotes rather than asserted.
  n <- 4000
  v <- seq(0, vmax, length.out = n + 1)
  wsimp <- c(1, rep(c(4, 2), length.out = n - 1), 1) * (vmax / n) / 3
  pv <- parts(v)
  ang <- sum(wsimp * pv$ang)
  len <- sum(wsimp * pv$len)
  c(ang, len)
}

separation <- function(r0) {
  if (r0 <= rh * (1 + 1e-12)) return(pi * rh)
  f <- function(L) 2 * halves(r0, L)[1] - pi
  L <- uniroot(f, c(1e-9, rh * (1 - 1e-9)), tol = 1e-12)$root
  2 * halves(r0, L)[2]
}

rs <- exp(seq(log(1.0005), log(10), length.out = 90))
ss <- vapply(rs, separation, numeric(1))

cat("=== separation to one's own fold image, units 2M = 1 ===\n\n")
for (r0 in c(2, 3, 5, 10))
  cat(sprintf("  s(%-4g) = %.4f\n", r0, separation(r0)))
cat(sprintf("  floor pi r_h = %.6f; smallest sampled s = %.6f\n", pi * rh, min(ss)))
stopifnot(abs(separation(2) - 6.0539) < 2e-3, abs(separation(3) - 8.5194) < 2e-3,
          abs(separation(10) - 23.7996) < 5e-3, min(ss) >= pi * rh - 1e-6)
cat("  matches the Python that the text quotes, to the tolerances above.\n")

C <- 2 * sqrt(2) / tanh(pi / (2 * sqrt(2)))
near <- exp(seq(log(1e-4), log(0.1), length.out = 40))
ratio <- vapply(near, function(e) (separation(rh * (1 + e)) - pi * rh) / (1 - rh / (rh * (1 + e))),
                numeric(1))
cat(sprintf("\n  closed form 2 sqrt2 coth(pi/(2 sqrt2)) = %.9f\n", C))
cat(sprintf("  ratio at eps = 1e-4 is %.9f, differing by %.2e\n", ratio[1], abs(ratio[1] - C)))
stopifnot(abs(ratio[1] - C) < 3e-4)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.0, 4.5, 2.2, 1.0), mgp = c(2.6, 0.7, 0), xpd = FALSE)

  plot(NA, xlim = c(1, 10), ylim = c(0, 25), axes = FALSE,
       xlab = expression(paste("radius  ", italic(r), "   (units ", 2*italic(M), " = 1)")),
       ylab = "separation to one's own image")
  rect(1, 0, 10, 25, col = "#fbfbfa", border = NA)
  axis(1, at = c(1, 2, 4, 6, 8, 10)); axis(2, at = seq(0, 25, 5), las = 1)
  abline(h = pi * rh, col = grey, lty = 2, lwd = 1.8)
  text(6.6, pi * rh - 1.5, expression(paste("floor  ", pi, italic(r)[h],
       "  (half the horizon circumference)")), col = grey, cex = 0.80, adj = 0.5)
  lines(rs, ss, col = c1, lwd = 2.8)
  points(1, pi * rh, pch = 19, cex = 1.0, col = c1)
  segments(7, 2 * (7 - 1) + pi - 2.4, 10, 2 * (10 - 1) + pi - 2.4, col = c2, lty = 3, lwd = 1.8)
  text(8.3, 2 * (8.3 - 1) + pi + 0.6, expression(paste("slope ", 2)), col = c2, cex = 0.86)
  mtext("nothing is nearer its image than half a turn", side = 3, line = 0.6,
        cex = 0.90, col = ink)

  plot(NA, xlim = c(-4, -1), ylim = c(3.50, 3.60), axes = FALSE,
       xlab = expression(paste(log[10], "  of  ", italic(r) / italic(r)[h] - 1)),
       ylab = expression(paste("(", italic(s) - pi * italic(r)[h], ") / ", italic(f))))
  rect(-4, 3.50, -1, 3.60, col = "#fbfbfa", border = NA)
  axis(1, at = -4:-1); axis(2, at = seq(3.50, 3.60, 0.02), las = 1)
  abline(h = C, col = c2, lty = 2, lwd = 1.8)
  lines(log10(near), ratio, col = c1, lwd = 2.8)
  text(-2.4, C - 0.012, expression(paste(2 * sqrt(2), " coth(", pi / (2 * sqrt(2)), ") = 3.5166")),
       col = c2, cex = 0.86, adj = 0.5)
  mtext("and it arrives as the square of the clock rate", side = 3, line = 0.6,
        cex = 0.90, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("papers/2_over_the_horizon/fig_companion_fold_reach.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.6, height = 4.3, pointsize = 12)
  else png(out, width = 9.6, height = 4.3, units = "in", res = 150, pointsize = 12)
  draw(); dev.off(); cat(sprintf("  wrote %s\n", out))
}
