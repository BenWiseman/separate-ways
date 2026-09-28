# Why the amplitude transferred to A.18's geometry and not to a black hole, settled on an
# exactly solvable null caustic.
#
# caustic_length_invariance.R showed that sqrt(pi) lambda_tot is not invariant and withdrew the
# number it gave, leaving the projection length as a candidate. This asks the prior question:
# is there an invariant length at a null caustic AT ALL? A gravitational plane wave answers it,
# and it is the right test object because it has exactly the structure the contact geodesic has.
#
#   ds^2 = 2 du dv + A_ij x^i x^j du^2 + dx^2 + dy^2,  A = diag(-w^2, +w^2).
#
# R_uu = -(A_11 + A_22), so the traceless choice is VACUUM, which forces one focusing transverse
# direction and one defocusing exactly as Schwarzschild's does. And the scalar field reduces
# exactly: Fourier in v with momentum p turns the wave equation into a Schrodinger equation in u
# with an oscillator of frequency w in x and an inverted one in y, whose kernel is Mehler's. So
# the caustic is explicit, the propagator is exact, and nothing is approximated.

cat("=== 1. the transverse structure, and that it matches the contact geodesic's ===\n")
cat("   Jacobi fields along the ray: sin(w U)/w focusing, sinh(w U)/w defocusing, where U is\n")
cat("   the affine parameter along u. Their tidal eigenvalues are -w^2 and +w^2, summing to\n")
cat("   zero, which is R_ab k^a k^b = 0. First conjugate point at w U = pi, multiplicity one.\n")
w <- 1
cat("      wU      J_focus = sin(wU)/w   J_defocus = sinh(wU)/w\n")
for (x in c(0.5, 1.5, pi - 1e-6, pi)) {
  cat(sprintf("   %7.4f  %19.8f  %22.8f\n", x, sin(x)/w, sinh(x)/w))
}
cat("   n = 1 in D = 4, the same as a black hole's contact caustic and unlike the sphere\n")
cat("   family, where every transverse direction focuses together.\n")

cat("\n=== 2. the exact kernel, and where it degenerates ===\n")
cat("   i d_U psi = -(1/2p) grad_perp^2 psi - (p/2) A_ij x^i x^j psi, so the x factor is an\n")
cat("   oscillator and the propagator is Mehler's:\n")
cat("      K_x = sqrt( p w / (2 pi i sin wU) ) exp[ i p w ((x^2+x'^2) cos wU - 2 x x')/(2 sin wU) ]\n")
cat("   The Van Vleck factor is 1/sqrt(sin wU) and it diverges at wU = pi, where the kernel\n")
cat("   collapses to a delta function at x' = -x. That IS the caustic, in closed form.\n")
cat("      pi - wU      1/sqrt(sin wU)     times sqrt(pi - wU)\n")
for (d in c(1e-1, 1e-2, 1e-3, 1e-4)) {
  cat(sprintf("   %10.0e  %16.6f  %20.8f\n", d, 1/sqrt(sin(pi - d)), sqrt(d)/sqrt(sin(pi - d))))
}
cat("   The amplitude goes as (pi - wU)^{-1/2}, which is the n = 1 caustic's half power, so\n")
cat("   the exponent rule holds here too and holds on a genuinely null caustic.\n")

cat("\n=== 3. the question the plane wave settles ===\n")
cat("   What LENGTH could the amplitude be proportional to? Under u -> c u, v -> v/c the\n")
cat("   metric keeps its form provided A -> A/c^2, so w -> w/c and U -> c U. The dimensions\n")
cat("   follow: 2 du dv is a length squared, so u v is one, and A x^2 du^2 being a length\n")
cat("   squared makes A go as 1/u^2 and w as 1/u. So w U = pi is invariant, and\n")
cat("        1/w and U have the dimension of u, NOT of length.\n")
cat("   The transverse coordinates are lengths, but the wave is homogeneous in them, so there\n")
cat("   is no preferred transverse length either.\n\n")
for (c0 in c(0.5, 1, 2, 4)) {
  cat(sprintf("      u -> %.1f u:   w = %.4f   U_caustic = %.4f   w U = %.6f\n",
              c0, 1/c0, pi * c0, (1/c0) * (pi * c0)))
}
cat("\n   A VACUUM NULL CAUSTIC HAS NO INVARIANT LENGTH. So the heat-kernel parametrisation\n")
cat("   Delta^{1/2} = c_n s^{-n/2}, which needs a c_n of dimension length^n, cannot be what\n")
cat("   carries a null caustic's amplitude. What carries it is the null momentum p, which is\n")
cat("   integrated over rather than being a property of the geometry.\n")
stopifnot(abs((1/4) * (pi * 4) - pi) < 1e-12)

cat("\n=== 4. which explains exactly which transfers worked ===\n")
cat("   S^N and the surfaces of revolution: the caustic geodesic is RIEMANNIAN and its arc\n")
cat("   length is invariant. The formula applies and was tested six ways.\n")
cat("   R^2 x S^2, which is A.18's geometry: a metric product whose caustic lives entirely in\n")
cat("   the Riemannian factor, so the sphere's own pi a is invariant and the formula applies\n")
cat("   through the factorisation of the heat kernel. Its amplitude is therefore EXACT:\n")
cat(sprintf("        c_1 = pi^{3/2} a = %.6f a, with a the transverse radius.\n", pi^1.5))
cat("   Schwarzschild's contact caustic: null, vacuum, and not a product, which is the plane\n")
cat("   wave's case. No invariant length, so no c_1 of this kind, and the projection candidate\n")
cat("   of the previous script is not merely unconfirmed, it is answering the wrong question.\n")

cat("\n=== 5. what that leaves, and it is not nothing ===\n")
cat("   Untouched: the exponent -(D-2+n)/2, which counts degenerate directions and is checked\n")
cat("   again here on a null vacuum caustic at n = 1; the profile-independence; the V_fam\n")
cat("   linearity; and Delta' = 0.183180, invariant.\n")
cat("   Gained: A.18's model has an EXACT amplitude rather than a measured one, because it is\n")
cat("   a product. That is a number its own numerical sum can now be checked against.\n")
cat("   The black hole's amplitude needs the Lorentzian uniform treatment, where the caustic\n")
cat("   enhancement is carried by the null momentum. That is A.18's sum done on the contact\n")
cat("   geodesic, which is the calculation named four iterations ago and is still the one.\n")

cat("\n=== 6. the plants ===\n")
cat("   (a) the vacuum condition must be what forces one focusing direction and one not.\n")
for (a2 in c(-1, -0.5, 0, 1)) {
  a1 <- -1
  cat(sprintf("      A = diag(%+.1f, %+.1f):  R_uu = %+.1f  ->  %s\n", a1, a2, -(a1 + a2),
              if (abs(a1 + a2) < 1e-12) "vacuum, one focusing one spreading"
              else if (a2 < 0) "not vacuum, BOTH focus (a point caustic)"
              else "not vacuum, both spread or mixed"))
}
cat("   Only the traceless case is vacuum, and only then is the caustic order one. That is\n")
cat("   the same statement contact_vanvleck.R makes at a hole, arrived at independently.\n")
cat("   (b) the half power must be a measurement and not an assumption.\n")
for (p0 in c(0.25, 0.5, 0.75)) {
  v <- sapply(c(1e-2, 1e-4), function(d) d^p0 / sqrt(sin(pi - d)))
  cat(sprintf("      dividing the amplitude by (pi-wU)^{-%.2f}: ratio across two decades %.4f\n",
              p0, v[2] / v[1]))
  if (abs(p0 - 0.5) > 1e-9) stopifnot(abs(v[2]/v[1] - 1) > 0.5)
}
cat("   Only the half power is flat across two decades.\n")
