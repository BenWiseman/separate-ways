# Why the plane wave has no invariant length, and whether a black hole inherits that.
#
# planewave_null_caustic.R concluded "a vacuum null caustic has no invariant length at all" and
# the companion says so. The argument ran on a plane wave. This checks whether the property being
# used is VACUUM or something narrower, because if it is narrower the conclusion does not transfer
# and a route that was closed is open again.

cat("=== 1. the plane wave's curvature invariants ===\n")
cat("   Brinkmann form ds^2 = 2 du dv + A_ij x^i x^j du^2 + dx^2, whose only independent Riemann\n")
cat("   components are R_{u i u j} = -A_ij. Raising an index pair needs g^{u a}, and the inverse\n")
cat("   metric has g^{uv} = 1 with g^{uu} = 0, so every raised u index becomes a v index and\n")
cat("   R^{u i u j} is proportional to R_{v i v j}, which vanishes identically.\n")
cat("   So every polynomial scalar invariant is zero. Plane waves are VSI.\n")
cat("      Kretschmann R_abcd R^abcd = 0\n")
cat("      Ricci scalar R            = 0\n")
cat("      R_ab R^ab                 = 0\n")
cat("   Numerically, with A = diag(-w^2, +w^2) and the only nonzero components being R_{uiuj}:\n")
w <- 1.3
Riem_uiuj <- c(w^2, -w^2)                 # -A_11, -A_22
cat(sprintf("      R_{u1u1} = %+.4f and R_{u2u2} = %+.4f are nonzero,\n", Riem_uiuj[1], Riem_uiuj[2]))
cat("      yet every full contraction pairs a u with a u through g^{uu} = 0, so all vanish.\n")
K_pw <- 0
cat(sprintf("      Kretschmann of the plane wave: %.1f\n", K_pw))

cat("\n=== 2. Schwarzschild's, which do not ===\n")
M <- 1
K_sch <- function(r) 48 * M^2 / r^6
cat("      r/M     Kretschmann 48 M^2/r^6     K^{-1/4}, an invariant LENGTH\n")
for (r in c(1.0, 1.5, 2.0)) {
  cat(sprintf("   %7.2f  %22.6f  %25.6f M\n", r, K_sch(r), K_sch(r)^(-1/4)))
}
stopifnot(K_sch(1) > 1)
cat("   Schwarzschild is Ricci-flat and its Kretschmann is not zero, so it carries invariant\n")
cat("   lengths at every point. K^{-1/4} is one, and it is 0.3799 M at the contact boundary.\n")
cat(sprintf("      K^{-1/4} at r = M: %.6f M, and 48^{-1/4} = %.6f\n",
            K_sch(1)^(-1/4), 48^(-1/4)))

cat("\n=== 3. so the conclusion was drawn too widely ===\n")
cat("   The plane wave has no invariant length because it is VSI, which is a statement about\n")
cat("   plane waves and not about vacuum. Vacuum means R_ab = 0; it does not mean the Weyl\n")
cat("   invariants vanish, and for a black hole they emphatically do not.\n")
cat("   What survives of the argument: the affine parameter and the Jacobi fields still supply no\n")
cat("   invariant length, since both scale under k -> c k. What does NOT survive: the claim that\n")
cat("   the geometry supplies none. At a hole it supplies K^{-1/4} and every other Weyl invariant\n")
cat("   of the right dimension.\n")

cat("\n=== 4. and the sphere family, to see whether the candidate can reduce correctly ===\n")
cat("   For S^N of radius a the Kretschmann is 2 N (N-1)/a^4, so K^{-1/4} is proportional to a\n")
cat("   and therefore to L = pi a. If the enhancement were built from K^{-1/4} it would reduce\n")
cat("   correctly on the family, with an N-dependent coefficient:\n")
cat("      N    K = 2N(N-1)/a^4     L / K^{-1/4} = pi a K^{1/4}\n")
for (N in 2:5) {
  K <- 2 * N * (N - 1); cat(sprintf("   %4d  %17.4f /a^4  %22.6f\n", N, K, pi * K^(1/4)))
}
cat("   The coefficient moves with N, so K^{-1/4} alone is not L. At the case that matters,\n")
cat("   though, only n = 1 is needed, and that is N = 2 on the family and Schwarzschild at a\n")
cat("   hole, so a single coefficient would do. What blocks it is that Schwarzschild's K varies\n")
cat("   along the contact geodesic while the sphere's does not, so the reduction fixes the\n")
cat("   coefficient and not the functional.\n")
cat(sprintf("      K^{-1/4} along the contact geodesic runs from %.4f M at r = M to %.4f M at r = 2M\n",
            K_sch(1)^(-1/4), K_sch(2)^(-1/4)))

cat("\n=== 5. where this leaves fork 2 ===\n")
cat("   The route is REOPENED. It was closed on the ground that a vacuum null caustic has no\n")
cat("   invariant length; that ground was a property of plane waves. The enhancement can be\n")
cat("   built from Weyl invariants along the contact geodesic, and the remaining question is\n")
cat("   which functional of them, which the sphere family constrains but does not determine\n")
cat("   because its K is constant where Schwarzschild's is not.\n")
cat("   NEAREST ROUTE: a compact-symmetry caustic whose curvature VARIES along the geodesic, so\n")
cat("   that the functional is pinned rather than just the coefficient. A warped product\n")
cat("   du^2 + f(u)^2 dOmega^2 with f not a sine has exactly that, and its heat kernel is the\n")
cat("   Sturm-Liouville problem caustic_amplitude_profile.R already solves.\n")

cat("\n=== 6. the plant ===\n")
cat("   The VSI claim must distinguish a plane wave from something that is not one.\n")
for (nm in list(list("plane wave", 0), list("Schwarzschild at r = M", K_sch(1)),
                list("S^2 of radius 1", 4), list("flat space", 0))) {
  cat(sprintf("      %-24s Kretschmann %10.4f   invariant length available: %s\n",
              nm[[1]], nm[[2]], ifelse(nm[[2]] > 0, "yes", "no")))
}
cat("   Flat space and the plane wave both return no length, for different reasons, and the two\n")
cat("   curved cases both return one. The test separates them.\n")

cat("\n=== 7. the decisive test, from machinery already built ===\n")
cat("   If the enhancement were any functional of the curvature along the geodesic, then\n")
cat("   changing the curvature profile at fixed L would change it. caustic_amplitude_profile.R\n")
cat("   did exactly that: surfaces of revolution du^2 + f^2 dphi^2 with f = a[sin th + eps h],\n")
cat("   h odd at both poles, which hold L = pi a and V_fam = 2 pi fixed while moving f.\n")
cat("   Their Gaussian curvature is K = -f''/f, and it moves a lot:\n\n")
a <- 1; L <- pi * a
Kg <- function(th, e) {
  f   <- function(t) sin(t) + e * sin(t)^3
  fpp <- function(t) -sin(t) + e * (6 * sin(t) * cos(t)^2 - 3 * sin(t)^3)
  -fpp(th) / f(th) / a^2
}
cat("      theta     K (eps=0)   K (eps=0.30)   K (eps=0.50)   ratio at 0.50\n")
for (th in c(0.4, 0.9, pi/2, 2.2, 2.7)) {
  cat(sprintf("   %8.3f  %11.5f  %13.5f  %13.5f  %14.4f\n",
              th, Kg(th, 0), Kg(th, 0.30), Kg(th, 0.50), Kg(th, 0.50) / Kg(th, 0)))
}
sp <- max(abs(sapply(seq(0.3, 2.8, length.out = 60), function(t) Kg(t, 0.5) / Kg(t, 0) - 1)))
cat(sprintf("\n   worst fractional change in the curvature at eps = 0.5: %.1f per cent\n", 100 * sp))
stopifnot(sp > 0.3)
cat("   And the amplitude did not move: the kernel ratios to the round case extrapolated to\n")
cat("   0.9991, 0.9978, 0.9956, 0.9982, 0.9959 and 0.9912 across six perturbations, against a\n")
cat("   curvature-dependent amplitude needing between 0.73 and 1.50.\n")

cat("\n=== 8. so the obstruction is structural, and is now nailed down ===\n")
cat("   The enhancement is not a functional of the curvature along the geodesic, because the\n")
cat("   curvature was changed by tens of per cent at fixed L and the answer did not move. It is\n")
cat("   the ARC LENGTH and nothing else. K^{-1/4} is therefore ruled out along with every other\n")
cat("   curvature functional, and reopening the route on the strength of Schwarzschild not being\n")
cat("   VSI does not survive contact with the profile test.\n")
cat("   A null geodesic has no arc length. That is the obstruction, and it is structural rather\n")
cat("   than a calculation not yet done.\n")
cat("   What the correction to section 3 is worth: the REASON in the manuscript was wrong even\n")
cat("   though the conclusion was right. 'A vacuum null caustic has no invariant length' is a\n")
cat("   property of plane waves, which are VSI. The right statement is that the amplitude needs\n")
cat("   an arc length, which a null geodesic lacks whatever invariants its spacetime carries.\n")

cat("\n=== 9. the plant for the decisive test ===\n")
cat("   The curvature-variation measurement must return zero on the unperturbed case.\n")
for (e in c(0, 0.15, 0.5)) {
  s0 <- max(abs(sapply(seq(0.3, 2.8, length.out = 60), function(t) Kg(t, e) / Kg(t, 0) - 1)))
  cat(sprintf("      eps = %.2f:  worst fractional curvature change %.4f%s\n", e, s0,
              ifelse(e == 0, "   <- zero, as it must be", "")))
  if (e == 0) stopifnot(s0 < 1e-12)
}
