#!/usr/bin/env Rscript
# fig_desitter_silence.R -- Section 1: why the two sheets are silent, drawn.
#
# de Sitter is conformal to a piece of the Einstein static cylinder:
#     ds^2 ~ -dT^2 + d chi^2 + sin^2(chi) d Omega_2^2,   T in (0, pi),
# with the spatial slices round three-spheres. The antipodal map of the embedding,
# X -> -X, is the antipodal map of S^3 together with T -> pi - T. Two points antipodal
# on S^3 are at geodesic distance exactly pi, whatever their position.
#
# So plot conformal time T against distance d ALONG S^3 from the point p, not against the
# polar angle chi: the conformal square hides the angular part of the separation, and that
# part is the whole argument. Radial null rays are 45 degree lines in these variables, so
# the causal region of p is the pair of cones |T - T_p| >= d, and p's image sits at
#     Theta p = (d = pi, T = pi - T_p).
# Reaching d = pi from p costs conformal time pi, and only pi - 2 T_p is available, which is
# less than pi for every T_p > 0. The image is outside both cones at every T_p, touching only
# in the limit T_p -> 0, which is conformal infinity and not a point of the spacetime.
#
# Equivalently and with no picture: Z = X.Y/l^2 is 1 at coincidence and greater than 1 for
# timelike separation, and Z(X, -X) = -1. The commutator vanishes identically.

gap <- function(Tp) pi - abs(pi - 2*Tp)      # how much short of the image the cone falls
stopifnot(all(gap(seq(0.01, pi - 0.01, length.out = 999)) > 0))
cat("=== the antipodal image is spacelike at every conformal time ===\n\n")
for (Tp in c(0.05, pi/4, pi/2, 3*pi/4, pi - 0.05))
  cat(sprintf("  T_p = %.4f   image at (d, T) = (pi, %.4f)   cone falls short by %.4f\n",
              Tp, pi - Tp, gap(Tp)))
cat(sprintf("\n  The shortfall is LARGEST at the midpoint, %.4f at T_p = pi/2, and smallest near\n", gap(pi/2)))
cat("  the ends: it tends to zero as T_p -> 0 or pi, which is conformal infinity and not a\n")
cat("  point of the spacetime. So the tight case is a point created or destroyed at the\n")
cat("  very beginning or end, and even there the separation is spacelike at every T_p > 0.\n")

# The check must be able to fail: halve the S^3 distance the map moves through and it does not.
bad <- function(Tp) pi/2 - abs(pi - 2*Tp)
cat(sprintf("\n  control, a map moving only pi/2 around the sphere: at T_p = 0.10 the shortfall is\n"))
cat(sprintf("  %.4f, negative, so THAT map puts a point in causal contact with its image. It is\n", bad(0.10)))
cat("  the full half-turn of the sphere that makes the separation spacelike everywhere, so\n")
cat("  the check can fail and the antipodal map is what stops it failing.\n")
stopifnot(bad(0.10) < 0, bad(pi/2) > 0)

ink <- "grey15"; c1 <- "#1f4e79"; mark <- "#a8400f"; cone <- "#dde5ee"
draw <- function() {
  par(mar = c(3.9, 4.4, 2.2, 1.2), mgp = c(2.5, 0.7, 0), xpd = FALSE)
  plot(NA, xlim = c(0, pi*1.06), ylim = c(0, pi), axes = FALSE,
       xlab = expression(paste("geodesic distance ", italic(d), " from ", italic(p), " along the spatial ",
                               italic(S)^3)),
       ylab = expression(paste("conformal time  ", italic(T))))
  rect(0, 0, pi*1.06, pi, col = "#fbfbfa", border = NA)
  axis(1, at = c(0, pi/4, pi/2, 3*pi/4, pi),
       labels = c("0", expression(pi/4), expression(pi/2), expression(3*pi/4), expression(pi)))
  axis(2, at = c(0, pi/2, pi), labels = c("0", expression(pi/2), expression(pi)), las = 1)

  Tp <- pi/2
  # Within the spacetime the cones are TRIANGLES, not half-planes: the future one runs out of
  # conformal time at T = pi having covered only d = pi - T_p, and the past one at T = 0 having
  # covered T_p. Drawing them as half-planes wraps them off the top of the panel and hides the
  # entire point, which is how the first version of this figure was wrong.
  polygon(c(0, pi - Tp, 0), c(Tp, pi, pi), col = cone, border = NA)
  polygon(c(0, Tp, 0),      c(Tp, 0, 0),   col = cone, border = NA)
  segments(0, Tp, pi - Tp, pi, col = "#9fb3c8", lwd = 1.3)
  segments(0, Tp, Tp, 0,      col = "#9fb3c8", lwd = 1.3)

  segments(pi, 0, pi, pi, col = mark, lwd = 2.6)
  text(pi - 0.05, 0.30, expression(paste("the image ", Theta, italic(p), " lies here, at")),
       col = mark, cex = 0.82, adj = 1)
  text(pi - 0.05, 0.17, expression(paste(italic(d), " = ", pi, " exactly, for every ", italic(T)[p])),
       col = mark, cex = 0.82, adj = 1)

  points(0, Tp, pch = 19, cex = 1.4, col = c1)
  text(0.08, Tp + 0.11, expression(italic(p)), col = c1, cex = 1.05, adj = 0)
  points(pi, pi - Tp, pch = 19, cex = 1.4, col = mark)
  text(pi - 0.05, pi - Tp + 0.15, expression(paste(Theta, italic(p))), col = mark, cex = 1.05, adj = 1)

  text(0.10, pi - 0.17, "everything", col = "#5b6b7d", cex = 0.82, adj = 0)
  text(0.10, pi - 0.32, "p can reach", col = "#5b6b7d", cex = 0.82, adj = 0)
  text(0.10, 0.30, "and everything", col = "#5b6b7d", cex = 0.82, adj = 0)
  text(0.10, 0.15, "that can reach p", col = "#5b6b7d", cex = 0.82, adj = 0)

  arrows(Tp + 0.06, pi - Tp, pi - 0.06, pi - Tp, code = 2, length = 0.07, col = ink, lwd = 1.1)
  text(pi*0.50, pi - Tp + 0.56, "the cones run out of conformal time", col = ink, cex = 0.82, adj = 0.5)
  text(pi*0.50, pi - Tp + 0.41, "having covered half the distance to the", col = ink, cex = 0.82, adj = 0.5)
  text(pi*0.50, pi - Tp + 0.26, "antipode. The rest of the sphere is", col = ink, cex = 0.82, adj = 0.5)
  text(pi*0.50, pi - Tp + 0.11, "spacelike, so the commutator vanishes.", col = ink, cex = 0.82, adj = 0.5)
  text(pi*0.50, 0.72, expression(paste("shortfall ", pi, " - |", pi, " - 2", italic(T)[p],
       "|, zero only at ", italic(T)[p], " = 0 or ", pi, ",")), col = ink, cex = 0.82, adj = 0.5)
  # the forty-five degree cone crosses this wherever it sits, because the label is wider
  # than the gap between the cone and the panel edge at every height. Masked.
  source("checks/fig_label.R")
  lab_on(pi*0.50, 0.44, "which is conformal infinity and not a point of the spacetime",
         bg = "#fbfbfa",
       col = ink, cex = 0.82, adj = 0.5)
  mtext("a point and its fold image are spacelike separated, always", side = 3, line = 0.7,
        cex = 0.88, col = ink)
}
for (f in c("paper/fig_companion_silence.pdf", "paper/fig_companion_silence.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 7.0, height = 4.6) else png(f, width = 1090, height = 715, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}
