# Deriving the caustic replacement rule instead of guessing an inner product.
#
# collective_jacobian_check.R settled the measure and left the collective-coordinate Jacobian
# short by exactly c^2 under the affine rescaling. Rather than guess again, solve for the rule on
# the family where the answer is known at four orders, then translate.

Vol <- function(n) 2 * pi^((n + 1) / 2) / gamma((n + 1) / 2)
a <- 1; L <- pi * a
c_n <- function(n) Vol(n) * (L / (2 * sqrt(pi)))^n     # established, two independent routes

cat("=== 1. solving for the degenerate transverse factor ===\n")
cat("   On S^N the heat kernel at the antipode is (4 pi s)^{-N/2} c_n s^{-n/2} e^{-L^2/4s} with\n")
cat("   n = N - 1. Split off the longitudinal (4 pi s)^{-1/2} and what is left covers the n\n")
cat("   degenerate transverse directions:\n")
cat("        X_n = (4 pi s)^{-n/2} c_n s^{-n/2}.\n")
cat("   Evaluate it and look for the pattern.\n\n")
cat("      n   X_n at s = 1        Vol(S^n) (L/4 pi s)^n     difference\n")
for (n in 1:4) for (s in c(1)) {
  X <- (4 * pi * s)^(-n / 2) * c_n(n) * s^(-n / 2)
  guess <- Vol(n) * (L / (4 * pi * s))^n
  cat(sprintf("   %4d  %18.10f  %22.10f   %.1e\n", n, X, guess, abs(X - guess)))
  stopifnot(abs(X - guess) < 1e-12)
}
cat("\n   and at other s, so it is not a coincidence at s = 1:\n")
for (s in c(0.3, 2.5)) {
  ok <- all(sapply(1:4, function(n)
    abs((4*pi*s)^(-n/2) * c_n(n) * s^(-n/2) - Vol(n) * (L/(4*pi*s))^n) < 1e-12))
  cat(sprintf("      s = %.1f: all four orders agree: %s\n", s, ok))
  stopifnot(ok)
}
cat("\n   So each degenerate direction contributes L/(4 pi s) and the family contributes its\n")
cat("   volume once:   X_n = Vol(S^n) ( L / 4 pi s )^n.\n")

cat("\n=== 2. reading that as a replacement for the vanishing Jacobi field ===\n")
cat("   A NON-degenerate transverse direction contributes (4 pi s)^{-1/2} (lambda_tot/J)^{1/2},\n")
cat("   the ratio to the flat Jacobi field being what makes it dimensionless. Matching:\n")
cat("        (4 pi s)^{-1/2} (lambda_tot/J)^{1/2}  ->  L/(4 pi s)\n")
cat("   which is the same as\n")
cat("        lambda_tot / J   ->   lambda_tot^2 / (4 pi s),    i.e.   J -> 4 pi s / L.\n")
for (s in c(0.5, 1, 2)) {
  lhs <- (4 * pi * s)^(-1/2) * (L / (4 * pi * s / L))^(1/2)
  cat(sprintf("      s = %.1f:  rebuilt %.8f   against L/(4 pi s) = %.8f\n", s, lhs, L / (4 * pi * s)))
  stopifnot(abs(lhs - L / (4 * pi * s)) < 1e-12)
}
cat("   The vanishing Jacobi field is replaced by 4 pi s / L. Nothing about the profile enters,\n")
cat("   which is caustic_amplitude_profile.R's result arrived at from the other side.\n")

cat("\n=== 3. the dictionary to the null case, checked on the direction that is NOT degenerate ===\n")
cat("   In the Schwinger representation the transverse problem has mass p and time DU, and the\n")
cat("   heat-kernel s corresponds to s = i DU / 2p. Feed that into the NON-degenerate factor:\n")
cat("        (4 pi s)^{-1/2} (lambda_tot/J)^{1/2}  with s = i lambda_tot / 2p\n")
lam_tot <- 3 * pi / 2 + 4; J2 <- 47.561945
for (p in c(0.7, 3.0)) {
  s_eq <- 1i * lam_tot / (2 * p)
  lhs <- (4 * pi * s_eq)^(-1/2) * (lam_tot / J2)^(1/2)
  rhs <- sqrt(p / (2i * pi * J2))
  cat(sprintf("      p = %.1f:  dictionary %s   standard sqrt(p/2 pi i J2) %s\n",
              p, format(lhs, digits = 8), format(rhs, digits = 8)))
  stopifnot(abs(lhs - rhs) < 1e-10)
}
cat("   It reproduces the standard sqrt(p / 2 pi i J) exactly, so the dictionary is right.\n")

cat("\n=== 4. and what it gives for the degenerate direction, with the obstruction now a single factor ===\n")
cat("   Writing the degenerate factor as the non-degenerate one times an enhancement,\n")
cat("        L/(4 pi s) = (4 pi s)^{-1/2} . [ L / sqrt(4 pi s) ],\n")
cat("   the bracket is dimensionless and is the whole of the caustic's effect. In null variables\n")
cat("   it is sqrt( p lambda_tot / 2 pi i ).\n")
for (p in c(0.7, 3.0)) {
  s_eq <- 1i * lam_tot / (2 * p)
  cat(sprintf("      p = %.1f:  L/sqrt(4 pi s) = %s   sqrt(p lambda_tot/2 pi i) = %s\n",
              p, format(lam_tot / sqrt(4 * pi * s_eq), digits = 8),
              format(sqrt(p * lam_tot / (2i * pi)), digits = 8)))
  stopifnot(abs(lam_tot / sqrt(4*pi*s_eq) - sqrt(p * lam_tot / (2i*pi))) < 1e-10)
}
cat("\n   Under the rescaling p -> c p and lambda_tot -> c lambda_tot that enhancement scales as\n")
cat("   c, and everything else in the assembly is invariant. So the residual mismatch is one\n")
cat("   factor and it is this one.\n")
for (c0 in c(0.5, 2, 4)) {
  e1 <- sqrt(0.7 * lam_tot); e2 <- sqrt((c0 * 0.7) * (c0 * lam_tot))
  cat(sprintf("      c = %.1f:  enhancement scales by %.4f\n", c0, e2 / e1))
  stopifnot(abs(e2 / e1 - c0) < 1e-10)
}

cat("\n=== 5. where that leaves it, stated exactly ===\n")
cat("   The Riemannian rule is derived rather than guessed and is exact at four orders. The\n")
cat("   dictionary is verified on the direction that is not degenerate. What has no invariant\n")
cat("   null counterpart is a single dimensionless enhancement sqrt(p lambda_tot / 2 pi i),\n")
cat("   because p lambda_tot is not invariant while p/J and p dlambda are.\n")
cat("   The one invariant of the right kind available is p times the WORLD FUNCTION's partner:\n")
cat("   lambda_tot DV is invariant, since lambda_tot scales as c and DV as 1/c, and\n")
cat("   lambda_tot DV is twice sigma. So the enhancement wants to be sqrt(p^2 lambda_tot DV) or\n")
cat("   equivalently p sqrt(2 sigma) rather than sqrt(p lambda_tot), and those differ by\n")
cat("   sqrt(p DV), which is invariant. That is a candidate and it is NOT yet a result: it is\n")
cat("   fixed by requiring invariance, which is one condition, and the earlier lesson is that one\n")
cat("   condition does not determine a functional form.\n")
cat("   NEAREST ROUTE: the same derivation on a Lorentzian case with a compact symmetry zero mode\n")
cat("   and a known propagator, which pins the enhancement instead of constraining it.\n")

cat("\n=== 6. the plant ===\n")
cat("   The pattern hunt in section 1 must reject a wrong pattern.\n")
for (k in c(2, 4, 8)) {
  ok <- all(sapply(1:4, function(n)
    abs((4*pi)^(-n/2) * c_n(n) - Vol(n) * (L/(k*pi))^n) < 1e-12))
  cat(sprintf("      Vol(S^n)(L/%d pi s)^n at s = 1:  %s\n", k,
              ifelse(ok, "matches all four", "rejected, as it must be")))
  if (k != 4) stopifnot(!ok)
}

cat("\n=== 7. testing the one natural invariant that could stand in for L ===\n")
cat("   The enhancement needs a dimensionless invariant. The obvious candidate built from the\n")
cat("   geodesic alone is int sqrt(T) dlambda, with T the focusing tidal eigenvalue: sqrt(T)\n")
cat("   scales as c and dlambda as 1/c, so the integral is invariant and dimensionless. On the\n")
cat("   sphere family T = |k_S|^2/a^2 is constant and the integral is exactly pi, which is also\n")
cat("   the angle the geodesic turns. If the contact geodesic returned pi as well, the\n")
cat("   enhancement could be written with it and the obstruction would be gone.\n")
sph <- pi                                        # int sqrt(T) dlambda on R x S^N
# contact geodesic: T = 3 M L^2 / r^5, dlambda = r^2 dphi / L, r = M(1 + sin phi), M = L = 1
bh <- integrate(function(p) sqrt(3) * (1 + sin(p))^(-1/2), 0, pi, rel.tol = 1e-12)$value
cat(sprintf("\n      sphere family:      int sqrt(T) dlambda = %.8f  (= pi)\n", sph))
cat(sprintf("      contact geodesic:   int sqrt(T) dlambda = %.8f\n", bh))
cat(sprintf("      ratio to pi: %.6f\n", bh / pi))
stopifnot(abs(bh - pi) > 0.5)
cat("   It is not pi. The geodesic turns through pi geometrically, but the tidal integral is a\n")
cat("   different quantity and it is 4.318. So this candidate is ruled OUT, and ruled out by a\n")
cat("   number rather than by taste.\n")
cat("   That is worth having: it was the only invariant available from the geodesic alone, so\n")
cat("   the enhancement cannot be built from the contact curve's own data. It needs something\n")
cat("   the product geometries supply and a hole does not, which is what the plane-wave argument\n")
cat("   said from the other direction.\n")
