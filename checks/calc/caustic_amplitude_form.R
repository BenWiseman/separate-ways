# What the caustic amplitude is made of, and the one factor four points cannot settle.
#
# caustic_amplitude.R fixed c_n at n = 1,2,3,4 on the sphere family, two independent ways, and
# left it written as L^n pi^{n/2} 4 pi / (n! Vol(S^{n+1})). That form references the volume of
# a sphere one dimension above the geodesic family, which is not something a local caustic can
# know about, so it is the wrong way to write a right answer. This finds a form built only
# from the family and then asks which parts of it the four points actually test.

a <- 1
Vol <- function(N) 2 * pi^((N + 1) / 2) / gamma((N + 1) / 2)
c_ref <- function(n) { L <- pi * a; L^n * pi^(n/2) * 4 * pi / (factorial(n) * Vol(n + 1)) }

cat("=== 1. first simplification: it is a function of the geodesic length alone ===\n")
cat("   Vol(S^{n+1}) = 2 pi^{(n+2)/2}/Gamma((n+2)/2) turns the reference form into\n")
cat("        c_n = L^n . 2 Gamma(n/2 + 1) / n! .\n")
cat("      n     c_n reference        2 Gamma(n/2+1)/n! . L^n      difference\n")
for (n in 1:4) {
  L <- pi * a; v <- L^n * 2 * gamma(n / 2 + 1) / factorial(n)
  cat(sprintf("   %4d  %18.10f   %24.10f   %.1e\n", n, c_ref(n), v, abs(v - c_ref(n))))
  stopifnot(abs(v - c_ref(n)) < 1e-12)
}
cat("   Cleaner, and no stray sphere. But L is still the only length a round sphere has.\n")

cat("\n=== 2. the form a collective-coordinate integral would produce ===\n")
cat("   The zero modes of the Jacobi operator are the directions the family moves in. In the\n")
cat("   heat-kernel path integral the action is (1/4s) int_0^1 |zdot|^2 dtau, so the norm the\n")
cat("   Gaussian measure uses is the tau-parametrised mean square,\n")
cat("        <J^2> = (1/L) int_0^L J(u)^2 du,\n")
cat("   with J the Jacobi field normalised to J'(0) = 1 and u arc length. Replacing each zero\n")
cat("   mode by its collective coordinate should then give a factor per mode and the family's\n")
cat("   own volume once. Test the shape\n")
cat("        c_n = V_fam . ( pi <J^2> / 2 )^{n/2}.\n")
J_sph  <- function(u) a * sin(u / a)
L_sph  <- pi * a
J2_sph <- integrate(function(u) J_sph(u)^2, 0, L_sph)$value / L_sph
cat(sprintf("\n   round sphere: <J^2> = %.10f, and a^2/2 = %.10f\n", J2_sph, a^2 / 2))
cat("      n   V_fam = Vol(S^n)    V_fam (pi <J^2>/2)^{n/2}     c_n reference     difference\n")
for (n in 1:4) {
  V <- Vol(n); v <- V * (pi * J2_sph / 2)^(n / 2)
  cat(sprintf("   %4d  %16.8f  %26.10f  %16.10f   %.1e\n", n, V, v, c_ref(n), abs(v - c_ref(n))))
  stopifnot(abs(v - c_ref(n)) < 1e-10)
}
cat("   All four, with V_fam genuinely different at each n (2pi, 4pi, 2pi^2, 8pi^2/3), so the\n")
cat("   LINEAR dependence on the family's volume is tested and not assumed.\n")

cat("\n=== 3. what the four points do NOT test ===\n")
cat("   Every one of them has the same Jacobi profile, a sin(u/a), so the bracket is only\n")
cat("   ever evaluated at one number. Any functional agreeing with <J^2> on a sine passes.\n")
cat("   The obvious rival is max(J)^2/2, which for a sine is the same thing exactly:\n")
cat(sprintf("      <J^2> = %.10f   max(J)^2/2 = %.10f   on the round sphere\n",
            J2_sph, max(sapply(seq(0, L_sph, length.out = 20001), J_sph))^2 / 2))
cat("   They differ on any other profile, and the contact geodesic is another profile.\n")

cat("\n=== 4. the contact geodesic's own numbers ===\n")
# J1(phi) = M (1 + sin phi) sin phi, with dlambda = r^2 dphi / L_ang and M = L_ang = 1
M <- 1
lam_tot <- 3 * pi / 2 + 4
J1   <- function(p) (1 + sin(p)) * sin(p)
dlam <- function(p) (1 + sin(p))^2
J2m  <- integrate(function(p) J1(p)^2 * dlam(p), 0, pi)$value / lam_tot
Jmax <- max(sapply(seq(0, pi, length.out = 200001), J1))
cat(sprintf("   affine length              L       = %.8f M\n", lam_tot))
cat(sprintf("   path-space mean square     <J^2>   = %.8f M^2\n", J2m))
cat(sprintf("   maximum of the field       max J   = %.8f M,  max^2/2 = %.8f M^2\n", Jmax, Jmax^2/2))
cat("   On a sine those last two agree; here they differ by a factor of\n")
cat(sprintf("        %.4f, so the two readings of the formula genuinely separate.\n", (Jmax^2/2)/J2m))
cat("\n   The reduced Van Vleck for the ONE non-degenerate transverse direction is not 1 here,\n")
cat("   unlike the sphere family where every transverse direction is degenerate. It is\n")
cat("   lambda_tot / J2(lambda_tot), which contact_vanvleck.R computed:\n")
J2end <- 47.561945
cat(sprintf("        Delta' = %.6f / %.6f = %.6f\n", lam_tot, J2end, lam_tot / J2end))
cat("\n      reading                         c_1                Delta'^{1/2} c_1\n")
for (nm in c("<J^2>", "max^2/2")) {
  br <- if (nm == "<J^2>") J2m else Jmax^2 / 2
  c1 <- 2 * pi * (pi * br / 2)^(1 / 2)
  cat(sprintf("      %-28s %12.6f M     %12.6f M\n", nm, c1, sqrt(lam_tot / J2end) * c1))
}
cat("\n   The two readings differ by 9.4 per cent in the bracket and, since c_1 goes as its\n")
cat("   square root, by 4.8 per cent in the amplitude. So the contact caustic's amplitude is\n")
br1 <- J2m; br2 <- Jmax^2 / 2
A1 <- sqrt(lam_tot / J2end) * 2 * pi * sqrt(pi * br1 / 2)
A2 <- sqrt(lam_tot / J2end) * 2 * pi * sqrt(pi * br2 / 2)
cat(sprintf("        Delta^{1/2} -> %.2f M . s^{-1/2},  the two readings giving %.3f and %.3f,\n",
            (A1 + A2) / 2, A1, A2))
cat("   which is a number with a five per cent uncertainty rather than an open question. The\n")
cat("   uncertainty is the functional form of the bracket and nothing else.\n")

cat("\n=== 5. a second example that looked free and is not ===\n")
cat("   S^2/Z_k, the sphere quotiented by a rotation of order k about the polar axis, has a\n")
cat("   heat kernel by images: K_quot(x,y) = sum_j K_{S^2}(x, g_j y). The poles are FIXED by\n")
cat("   every g_j, so at the pole pair all k terms coincide and K_quot = k K_{S^2}. Meanwhile\n")
cat("   the meridians downstairs run over phi in [0, 2pi/k), so V_fam = 2 pi / k.\n")
for (k in 2:4) {
  cat(sprintf("      k = %d:  kernel is %d times the sphere's, but V_fam is %.4f = (1/%d) of it\n",
              k, k, 2 * pi / k, k))
}
cat("   The formula would need the amplitude to fall like 1/k and it rises like k, a factor\n")
cat("   of k^2 out. That is not a refutation of section 2: the quotient puts CONICAL points\n")
cat("   exactly where the geodesics converge, so the caustic's local structure is not the\n")
cat("   smooth one the formula describes. Recorded as a test that does not apply, with the\n")
cat("   reason, rather than as a discrepancy left lying around.\n")

cat("\n=== 6. what would settle it ===\n")
cat("   A compact symmetric space whose Jacobi fields along a geodesic do NOT all share one\n")
cat("   frequency. Complex projective space is the standard example: along a Fubini-Study\n")
cat("   geodesic the Hopf direction has twice the frequency of the rest, so conjugate points\n")
cat("   come at two distances and the profile is not a single sine. Its spectrum is known in\n")
cat("   closed form, so the same two routes would run. That is the calculation, and it is\n")
cat("   the one thing standing between this and an amplitude at a black hole.\n")

cat("\n=== 7. the plants ===\n")
cat("   (a) the section 2 form must fail if the family volume enters any other way.\n")
for (p in c(0.5, 1, 1.5, 2)) {
  err <- max(sapply(1:4, function(n) abs(Vol(n)^p * (pi * J2_sph / 2)^(n / 2) / c_ref(n) - 1)))
  cat(sprintf("      V_fam^%.1f:  worst relative error over n = 1..4 is %.3f%s\n",
              p, err, if (abs(p - 1) < 1e-9) "   <- the only exponent that works" else ""))
  if (abs(p - 1) > 1e-9) stopifnot(err > 0.1)
}
cat("   (b) and if the per-mode bracket carried a different power of n.\n")
for (q in c(0.25, 0.5, 0.75)) {
  err <- max(sapply(1:4, function(n) abs(Vol(n) * (pi * J2_sph / 2)^(n * q) / c_ref(n) - 1)))
  cat(sprintf("      exponent n*%.2f:  worst relative error %.3f%s\n",
              q, err, if (abs(q - 0.5) < 1e-9) "   <- the only one that works" else ""))
  if (abs(q - 0.5) > 1e-9) stopifnot(err > 0.1)
}
