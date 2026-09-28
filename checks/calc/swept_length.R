# Ben's question: isn't the length inside the horizon just the radius?
#
# The obstruction had been that the amplitude wants a length, every solvable case supplies the
# connecting geodesic's ARC length, and a null geodesic has none. That treated the null curve's own
# arc length as the only candidate. It is not. Inside the horizon r is timelike, so the proper
# extent the geodesic sweeps in r is an invariant length of the curve even though the curve is null.

M <- 1
cat("=== 1. the swept proper extent, and an identity that was not expected ===\n")
cat("   Inside the horizon dtau = dr / sqrt(2M/r - 1), so the proper extent the contact geodesic\n")
cat("   sweeps going up from r = M to 2M and back is\n")
tau_half <- integrate(function(r) sqrt(r / (2 * M - r)), M, 2 * M, rel.tol = 1e-12)$value
cat(sprintf("        2 int_M^{2M} sqrt(r/(2M-r)) dr = %.8f M,  closed form M(pi + 2) = %.8f M\n",
            2 * tau_half, pi + 2))
stopifnot(abs(2 * tau_half - (pi + 2)) < 1e-8)
cat("\n   And separately, the arc length of the geodesic's projection onto the transverse sphere:\n")
proj <- integrate(function(p) M * (1 + sin(p)), 0, pi, rel.tol = 1e-12)$value
cat(sprintf("        int_0^pi r dphi = %.8f M, the same number.\n", proj))
stopifnot(abs(proj - (pi + 2)) < 1e-10)
cat("\n   That is not a coincidence, it is an identity. Along this curve dr/dphi = sqrt(r(2M-r)),\n")
cat("   so dphi = dr/sqrt(r(2M-r)) and\n")
cat("        r dphi = r dr / sqrt(r(2M-r)) = sqrt(r/(2M-r)) dr = dtau.\n")
cat("   The transverse arc length the geodesic covers and the proper extent it sweeps in r are\n")
cat("   the SAME integral, term by term. Checked pointwise:\n")
cat("   The product sqrt(r/(2M-r)) . sqrt(r(2M-r)) is sqrt(r^2) = r identically, so the identity\n")
cat("   is algebraic. At the turning point it reads 0 times infinity numerically, so check it\n")
cat("   away from there and state the algebra there:\n")
for (p in c(0.3, 1.0, 1.4, 2.4, 2.9)) {
  r <- M * (1 + sin(p))
  rhs <- sqrt(r / (2 * M - r)) * sqrt(r * (2 * M - r))
  cat(sprintf("      phi = %.3f:  r = %.8f   dtau/dphi = %.8f   difference %.1e\n",
              p, r, rhs, abs(r - rhs)))
  stopifnot(abs(r - rhs) < 1e-7)
}
cat("      At phi = pi/2 the algebra gives r = 2M and the numerics are undefined, a removable\n")
cat("      point rather than a failure.\n")

cat("\n=== 2. the case that decides it: a Lorentzian null caustic whose answer is known ===\n")
cat("   In A.18's geometry R^2 x S^2(a) the connecting geodesic at contact is NULL, so its own\n")
cat("   arc length is zero, and yet the amplitude is known exactly to be c_1 = pi^{3/2} a, which\n")
cat("   is not zero. So the L in c_1 = sqrt(pi) L is NOT the null curve's arc length there.\n")
cat(sprintf("      c_1 = pi^{3/2} a = %.6f a, and sqrt(pi) x (null arc length 0) = 0.\n", pi^1.5))
cat("   What it IS equals the transverse projection: the geodesic covers pi of the sphere, so the\n")
cat("   projected arc length is pi a, and sqrt(pi) . pi a = pi^{3/2} a exactly.\n")
stopifnot(abs(sqrt(pi) * pi - pi^1.5) < 1e-14)
cat("   So on the one Lorentzian null caustic with a known amplitude, the rule is the projection\n")
cat("   and the null arc length is ruled out by the answer being nonzero.\n")

cat("\n=== 3. is it invariant? which is what killed the two earlier numbers ===\n")
cat("   Neither int r dphi nor int dtau mentions the affine parameter, so rescaling the tangent\n")
cat("   k -> c k cannot touch them. Contrast with the affine length, which moves:\n")
lam_tot <- 3 * pi / 2 + 4
for (c0 in c(0.5, 1, 2, 4)) {
  cat(sprintf("      c = %.1f:  affine length %8.4f (moves)   int r dphi %.6f (fixed)\n",
              c0, lam_tot / c0, proj))
}
cat("   That is the test sqrt(pi) lambda_tot failed and this passes.\n")

cat("\n=== 4. the number ===\n")
J2end <- 47.561945
c1 <- sqrt(pi) * proj
amp <- sqrt(lam_tot / J2end) * c1
cat(sprintf("      swept length              L = M(pi + 2) = %.6f M\n", proj))
cat(sprintf("      c_1 = sqrt(pi) L                        = %.6f M\n", c1))
cat(sprintf("      reduced Van Vleck Delta'^{1/2}          = %.6f\n", sqrt(lam_tot / J2end)))
cat(sprintf("      Delta^{1/2} -> %.4f M s^{-1/2}\n", amp))
stopifnot(abs(amp - 3.900) < 0.002)

cat("\n=== 5. and it is not a constant, which is Ben's other question ===\n")
cat("   The amplitude is evaluated on the connecting geodesic, and that curve changes as the\n")
cat("   source moves off the contact boundary, so the swept length changes with it and the\n")
cat("   coefficient has a SHAPE rather than a value. For r < M the connecting curve is timelike\n")
cat("   and turns through pi without needing to reach the horizon, so it sweeps less:\n")
cat("      the E = 0 null curve from radius r turns 2 pi - 4 asin sqrt(r/2M), which is pi only\n")
cat("      at r = M, so inside that the connector is a different curve.\n")
cat("   Computing the shape needs the shooting solution rather than the closed form, which is\n")
cat("   contact_world_function.py's machinery. Flagged and not guessed.\n")

cat("\n=== 6. the plants ===\n")
cat("   (a) the identity must fail on a curve that is not this one.\n")
for (k in c(0.8, 1.0, 1.3)) {
  rr <- function(p) M * (1 + k * sin(p))
  lhs <- integrate(function(p) rr(p), 0, pi)$value
  ok <- abs(lhs - (pi + 2)) < 1e-8
  cat(sprintf("      r = M(1 + %.1f sin phi):  int r dphi = %.6f  %s\n", k, lhs,
              ifelse(k == 1, "<- the contact geodesic", "not the contact geodesic")))
  if (k != 1) stopifnot(!ok)
}
cat("   (b) and the swept length must differ from the affine length, or the test is empty.\n")
cat(sprintf("      swept %.6f M against affine %.6f M, a ratio of %.4f\n",
            proj, lam_tot, proj / lam_tot))
stopifnot(abs(proj / lam_tot - 1) > 0.3)
