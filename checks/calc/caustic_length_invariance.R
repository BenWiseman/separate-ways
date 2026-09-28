# Which length the caustic amplitude uses, and why the number I put in yesterday is not one.
#
# caustic_amplitude_profile.R established c_n = V_fam (L / 2 sqrt(pi))^n, profile-independent,
# and applied it to the contact geodesic with L = lambda_tot = 3 pi/2 + 4, the affine length
# GR57 computes. That step does not survive inspection.
#
# Delta and s are both invariant under rescaling a null geodesic's affine parameter k -> c k,
# so c_1 must be too. But lambda scales as 1/c and so does every Jacobi field, so the affine
# length is NOT invariant and sqrt(pi) lambda_tot is a choice of normalisation rather than a
# number. On a Riemannian caustic L is arc length and the problem does not arise, which is why
# the sphere tests never saw it.

a <- 1

cat("=== 1. the scaling, made explicit on a case with a known answer ===\n")
cat("   Take A.18's geometry R^2 x S^2(a). The caustic lives in the sphere factor and the true\n")
cat("   answer is c_1 = pi^{3/2} a, fixed by the sphere alone. Now look at the NULL geodesic\n")
cat("   from a point to its antipodal image in the 4D spacetime. Its sphere-direction Jacobi\n")
cat("   field obeys J'' = -(|k_S|^2/a^2) J, so J = (a/|k_S|) sin(|k_S| lambda / a) and\n")
cat("        lambda_tot = pi a / |k_S|,\n")
cat("   which depends on how the tangent was normalised. The answer does not.\n\n")
cat("      |k_S|   lambda_tot    sqrt(pi) lambda_tot    true c_1 = pi^{3/2} a\n")
for (kS in c(0.5, 1, 2, 4)) {
  lt <- pi * a / kS
  cat(sprintf("   %7.2f  %11.6f  %20.6f  %22.6f\n", kS, lt, sqrt(pi) * lt, pi^1.5 * a))
}
cat("   The middle column moves by a factor of eight and the right one does not move at all.\n")
cat("   So sqrt(pi) lambda_tot is not the amplitude, and the 6.61 M it gave is WITHDRAWN.\n")

cat("\n=== 2. and arc length does not rescue it either ===\n")
cat("   The connecting geodesic for a point just off contact is spacelike with proper length\n")
cat("   l = sqrt(2 sigma), and sigma -> 0 AT contact, so l -> 0. Parametrised by arc length the\n")
cat("   sphere projection still turns pi, so |k_S| = pi a / l and the Jacobi field is\n")
cat("   (l/pi) sin(pi u / l), vanishing at u = 0 and u = l. Feeding L = l into the formula:\n")
cat("        delta        l = sqrt(2 sigma)     sqrt(pi) l      true c_1\n")
for (d in c(1e-1, 1e-2, 1e-3)) {
  sig <- pi * a^2 * d; l <- sqrt(2 * sig)
  cat(sprintf("   %9.0e  %18.6f  %14.6f  %12.6f\n", d, l, sqrt(pi) * l, pi^1.5 * a))
}
cat("   It goes to zero while the answer stays at 5.568. Neither the affine length nor the\n")
cat("   arc length of the connecting geodesic is the L in the formula.\n")

cat("\n=== 3. what the L actually was, on every case that worked ===\n")
cat("   Round sphere S^N: L = pi a, the arc length of the caustic geodesic IN THE SPHERE.\n")
cat("   Surfaces of revolution: L = pi a, the meridian's arc length, same thing.\n")
cat("   R^2 x S^2: L = pi a again, the length of the geodesic's PROJECTION onto the sphere,\n")
cat("   which is the factor the caustic lives in and is where the Jacobi field is degenerate.\n")
cat("   Every case is the same invariant: the arc length of the projection onto the focusing\n")
cat("   factor. For a product that is unambiguous. Schwarzschild is not a product.\n")

cat("\n=== 4. the candidate that reduces correctly, offered as a candidate ===\n")
cat("   The transverse sphere is still there at a black hole, with radius r(phi), and the\n")
cat("   contact geodesic's projection onto it has arc length int r dphi. On a product r is\n")
cat("   constant and that is pi a exactly, so the candidate reduces correctly on every case\n")
cat("   above. On the contact geodesic, with r(phi) = M(1 + sin phi):\n")
Lproj <- integrate(function(p) 1 + sin(p), 0, pi)$value
cat(sprintf("        int_0^pi r dphi = M (pi + 2) = %.6f M\n", Lproj))
cat(sprintf("        against the affine length lambda_tot = 3pi/2 + 4 = %.6f M\n", 3*pi/2 + 4))
J2end <- 47.561945; lam_tot <- 3 * pi / 2 + 4
for (nm in c("projection", "affine")) {
  Lx <- if (nm == "projection") Lproj else lam_tot
  cat(sprintf("        %-11s  c_1 = %9.4f M   Delta'^{1/2} c_1 = %8.4f M\n",
              nm, sqrt(pi) * Lx, sqrt(lam_tot / J2end) * sqrt(pi) * Lx))
}
cat("   The projection reading gives 3.90 M. It is a CANDIDATE and not a result: it reduces\n")
cat("   correctly on the one family available, which is exactly the position the profile\n")
cat("   question was in before a non-symmetric surface settled it. A second non-product\n")
cat("   caustic would settle this one the same way.\n")

cat("\n=== 5. what survives untouched ===\n")
cat("   The exponent. -(D - 2 + n)/2 is a count of degenerate directions and carries no\n")
cat("   length, so nothing here touches it, and it still rests on five points.\n")
cat("   The profile-independence. Six perturbed surfaces of revolution, ratios to 1.\n")
cat("   The V_fam linearity, tested across four caustic orders.\n")
cat("   The reduced Van Vleck 0.183180, which is a ratio of like quantities and invariant:\n")
for (c0 in c(0.5, 1, 2)) {
  cat(sprintf("      k -> %.1f k:  lambda_tot/J2 = %.6f\n", c0, (lam_tot/c0) / (J2end/c0)))
}
cat("   What is withdrawn is one step, the identification of L, and the number it produced.\n")

cat("\n=== 6. the plant ===\n")
cat("   The invariance test must pass something that IS invariant and fail something that is\n")
cat("   not, or it is not a test.\n")
for (c0 in c(0.5, 1, 2, 4)) {
  cat(sprintf("      k -> %.1f k:  affine length %8.4f (moves)   lambda/J2 %.6f (fixed)   int r dphi %.6f (fixed)\n",
              c0, lam_tot / c0, (lam_tot/c0)/(J2end/c0), Lproj))
}
stopifnot(abs((lam_tot/2)/(J2end/2) - lam_tot/J2end) < 1e-12)
cat("   The affine length moves by a factor of eight across those rows. That is the whole\n")
cat("   argument, and it is arithmetic rather than opinion.\n")
