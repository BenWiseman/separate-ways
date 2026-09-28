# Can the Group F machinery see the contact divergence at all? A decisive diagnostic.
#
# A.18 establishes that the image two-point function diverges at the contact configuration,
# as s^-3/2. That is not in doubt: it is computed there from a mode sum in a model geometry,
# with a control returning the known generic power and a plant that removes the divergence.
#
# Every calculation in Group F builds the image correlator from Schwarzschild-coordinate
# mode functions evaluated at the same (t, r) for both points, because Theta preserves both.
# If that machinery is right it must reproduce the divergence at r = M, where A.15 puts the
# contact boundary. If it returns something finite there, it cannot see the contact at all
# and every number Group F has produced is void, whatever its convergence looked like.
#
# The cheapest version of the test is <phi^2>_img rather than the stress: same mode sum, no
# derivatives, same parity, and A.18 says it must blow up.

M <- 1
f <- function(r) 1 - 2 * M / r
V <- function(r, L) f(r) * (L * (L + 1) / r^2 + 2 * M / r^3)
p <- function(r, L, w) sqrt(w^2 - V(r, L))

# WKB |psi|^2 ~ 1/p, so <phi^2>_img ~ sum_l c_l (-1)^l / p_l
phi2 <- function(r, w, eps) {
  l <- 0:ceiling(45 / eps)
  sum(((2 * l + 1) / (4 * pi)) * (-1)^l / p(r, l, w) * exp(-l * eps))
}
rich <- function(v) { while (length(v) > 1) v <- 2 * v[-1] - v[-length(v)]; v }

cat("=== 1. <phi^2>_img from this machinery, approaching the contact boundary ===\n")
cat("   A.18 requires a divergence at r = M. Watch whether one appears.\n")
cat("      r/M        eps=0.04      0.02        0.01       0.005     extrapolated\n")
for (r in c(1.20, 1.05, 1.00, 0.95, 0.80)) {
  es <- c(0.04, 0.02, 0.01, 0.005)
  vs <- sapply(es, function(e) phi2(r, 1, e))
  cat(sprintf("   %6.2f  %11.5f %11.5f %11.5f %11.5f   %+12.6f\n",
              r, vs[1], vs[2], vs[3], vs[4], rich(vs)))
}

cat("\n=== 2. squeeze right onto the boundary ===\n")
cat("      |r - M|        extrapolated <phi^2>_img\n")
for (d in c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)) {
  vs <- sapply(c(0.04, 0.02, 0.01, 0.005), function(e) phi2(M - d, 1, e))
  cat(sprintf("   %9.0e      %+16.8f\n", d, rich(vs)))
}

cat("\n=== 3. the verdict ===\n")
v_far  <- rich(sapply(c(0.04, 0.02, 0.01, 0.005), function(e) phi2(M - 1e-1, 1, e)))
v_near <- rich(sapply(c(0.04, 0.02, 0.01, 0.005), function(e) phi2(M - 1e-5, 1, e)))
cat(sprintf("   at |r-M| = 1e-1: %+.8f\n", v_far))
cat(sprintf("   at |r-M| = 1e-5: %+.8f\n", v_near))
cat(sprintf("   ratio over four decades of approach: %.4f\n", v_near / v_far))
if (abs(v_near / v_far) < 10) {
  cat("\n   FINITE. Four decades closer to the contact boundary and the value barely moves.\n")
  cat("   A.18 requires a divergence here. This machinery does not produce one, so it is\n")
  cat("   not computing the image correlator at the contact configuration.\n")
  cat("   The reason is structural: Theta maps Kruskal region II to region IV, the two are\n")
  cat("   distinguished only by the sign of the Kruskal time, and the Schwarzschild chart\n")
  cat("   assigns them the same (t, r). Mode functions written in that chart cannot tell\n")
  cat("   the two points apart, so no separation enters and no divergence can appear.\n")
  cat("   Everything Group F computed with it is void. The fix is Kruskal coordinates.\n")
} else {
  cat("\n   DIVERGES. The machinery does see the contact, and the earlier Group F numbers\n")
  cat("   need re-examining rather than discarding.\n")
}
