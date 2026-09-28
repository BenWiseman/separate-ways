# Checking the collective-coordinate structure, and finding where the bookkeeping actually breaks.
#
# M3 walked this branch and, before its budget ran out, proposed a structure. It is worth checking
# rather than inheriting, and the check turns up something useful: the structure is right, the
# power it gives agrees with the one derived independently, and the step that is NOT yet right is
# one nobody had named.
#
# The proposal. Extract the rotation zero mode as a collective coordinate. The path-space inner
# product with the mass p explicit is <f,g> = p int f g dlambda, so the Jacobian is
# sqrt(p int xi^2 dlambda). The degenerate direction's remaining Gaussian is the REDUCED
# determinant, which by Gelfand-Yaglom is governed by the Jacobi field's derivative at the
# endpoint rather than by the field itself, so it contributes sqrt(p / 2 pi i xidot(T)) with
# xidot(T) standing in for the vanishing J1(T). The spreading direction contributes
# sqrt(p / 2 pi i J2). The orbit contributes its volume, 2 pi.

M <- 1; lam_tot <- 3 * pi / 2 + 4; J2end <- 47.561945
xi  <- function(p) (1 + sin(p)) * sin(p)          # the Killing displacement per unit rotation
dlam<- function(p) (1 + sin(p))^2
I2  <- integrate(function(p) xi(p)^2 * dlam(p), 0, pi, rel.tol = 1e-12)$value
xidT <- -1                                         # dxi/dlambda at the far end, exactly -L/M

cat("=== 1. the pieces, computed ===\n")
cat(sprintf("   orbit volume                              2 pi        = %.6f\n", 2 * pi))
cat(sprintf("   int xi^2 dlambda along the contact curve              = %.6f M^3\n", I2))
cat(sprintf("   xidot at the far end (Gelfand-Yaglom stand-in)        = %.1f\n", xidT))
cat(sprintf("   spreading partner                          J2(T)      = %.6f\n", J2end))
cat(sprintf("   affine length                              lambda_tot = %.6f\n", lam_tot))

cat("\n=== 2. the power, which agrees with the route derived separately ===\n")
cat("   Jacobian sqrt(p .) times reduced Gaussian sqrt(p/.) times spreading sqrt(p/.) gives\n")
cat("   p^{3/2}. null_momentum_route.R derived p^{1 + n/2} independently, which at n = 1 is the\n")
cat("   same. Two routes, same power, so the structure is not obviously wrong.\n")
for (n in 0:2) cat(sprintf("      n = %d:  collective structure p^{%.1f}   rule p^{1+n/2} = p^{%.1f}\n",
                           n, 1 + n/2, 1 + n/2))

cat("\n=== 3. the invariance test, which is where two earlier numbers died ===\n")
cat("   The residual freedom of Brinkmann coordinates is U -> c U, V -> V/c, under which\n")
cat("   lambda -> lambda/c, every Jacobi field -> J/c, p -> c p, and xi is UNCHANGED because it\n")
cat("   is a Killing displacement per unit rotation angle and knows nothing of the affine scale.\n")
scal <- function(c0) list(
  jac   = c0 * (1 / c0),                  # p int xi^2 dlambda : p up by c, dlambda down by c
  red   = c0 / c0,                        # p / xidot : both up by c
  spr   = c0 * c0,                        # p / J2 : p up by c, J2 down by c
  meas  = c0,                             # dp
  dv    = 1 / c0)                         # DV
cat("      c      p.int xi^2   p/xidot   p/J2    dp    DV\n")
for (c0 in c(0.5, 1, 2, 4)) {
  s <- scal(c0)
  cat(sprintf("   %5.2f  %10.3f %9.3f %7.3f %5.2f %5.3f\n", c0, s$jac, s$red, s$spr, s$meas, s$dv))
}
cat("   The Jacobian and the reduced Gaussian are invariant, which is the good news and is what\n")
cat("   M3's inner product was chosen to achieve. The spreading factor is not: p/J2 scales as\n")
cat("   c^2.\n")

cat("\n=== 4. and the same failure is already there with NO caustic, which is the finding ===\n")
cat("   Run the identical bookkeeping on a non-degenerate pair, where the answer is known to be\n")
cat("   a scalar. There A(p) = p / (2 pi i sqrt(J1 J2)) and G = int (dp/2pi) A e^{-i p DV}:\n")
for (c0 in c(0.5, 2, 4)) {
  A_s  <- c0 * c0        # p up by c, sqrt(J1 J2) down by c
  tot  <- c0 * A_s       # times the measure dp
  cat(sprintf("      c = %.2f:  A(p) scales as %6.3f, and with dp the integral scales as %7.3f\n",
              c0, A_s, tot))
}
cat("   A propagator between two fixed points is a scalar and cannot scale at all, so the naive\n")
cat("   assignment is wrong BEFORE any caustic enters. The error is in the weight carried by the\n")
cat("   momentum-space measure under the rescaling, not in the collective coordinate.\n")
stopifnot(abs(4 * 4 * 4 - 64) < 1e-12)

cat("\n=== 5. what that leaves, named precisely ===\n")
cat("   IMPASSE: the p-space measure's weight under the residual Brinkmann rescaling.\n")
cat("   LOAD-BEARING: that A(p) as written above is the full integrand, with no compensating\n")
cat("                 factor of p or of the transverse volume hidden in the Fourier convention.\n")
cat("   EVIDENCE:     with that assumption the NON-degenerate propagator scales as c^3 instead\n")
cat("                 of being a scalar, so the assumption is false and the missing factor is\n")
cat("                 p^{-3} or equivalent.\n")
cat("   NEAREST ROUTE: fix the convention on the non-degenerate case, where the answer is the\n")
cat("                 known flat propagator, and only then put the zero mode back in. Doing it\n")
cat("                 in the other order is what produced two withdrawn numbers.\n")
cat("   This is a better place to be than before: the structure and the power are settled, and\n")
cat("   what is left is a normalisation fixable against a case whose answer is already known.\n")

cat("\n=== 6. the plant ===\n")
cat("   The invariance test must pass something that IS invariant, or it proves nothing.\n")
for (c0 in c(0.5, 2, 4)) {
  cat(sprintf("      c = %.2f:  lambda_tot/J2 = %.6f (invariant, as it must be)\n",
              c0, (lam_tot / c0) / (J2end / c0)))
  stopifnot(abs((lam_tot/c0)/(J2end/c0) - lam_tot/J2end) < 1e-12)
}
cat("      and p.int xi^2 dlambda, the piece M3's inner product was built to protect, holds too.\n")

cat("\n=== 7. fixing the convention on the case whose answer is known ===\n")
cat("   Section 5 said to calibrate on the non-degenerate case first. Doing that settles it.\n")
cat("   Flat space in Brinkmann form, ds^2 = 2 dU dV + dx^2. The transverse Schrodinger kernel\n")
cat("   with mass p in two dimensions is p/(2 pi i DU) exp(i p Dx^2 / 2 DU), so with a measure\n")
cat("   dp/(2 pi) the p integral gives int dp p e^{i p Xi} = -1/Xi^2 and\n")
cat("        G ~ 1/(DU Xi^2),   Xi = Dx^2/2DU + DV.\n")
cat("   But DU Xi is half the squared interval, so that is (interval^2)^{-2} times a stray DU,\n")
cat("   and the flat Hadamard form is (interval^2)^{-1}. The naive measure is wrong by p^{-1}.\n\n")
cat("   With the LIGHT-FRONT invariant measure dp/(4 pi p) instead:\n")
cat("        G = int dp/(4 pi p) . p/(2 pi i DU) e^{i p Xi} = (1/8 pi^2 i DU)(i/Xi) = 1/(8 pi^2 DU Xi),\n")
cat("   and with DU Xi = interval^2 / 2 that is 1/(4 pi^2 interval^2), the flat Hadamard form\n")
cat("   exactly. So the measure is dp/(4 pi p), which is the standard light-front one.\n")
DU <- 1.7; Dx <- 0.9; DV <- -0.35
Xi <- Dx^2 / (2 * DU) + DV     # the phase carries +i p DV in this convention
int2 <- 2 * DU * DV + Dx^2
cat(sprintf("      check on numbers: DU Xi = %.8f and interval^2/2 = %.8f\n", DU * Xi, int2 / 2))
stopifnot(abs(DU * Xi - int2 / 2) < 1e-12)
cat(sprintf("      so G = 1/(8 pi^2 DU Xi) = %.8f and 1/(4 pi^2 interval^2) = %.8f\n",
            1 / (8 * pi^2 * DU * Xi), 1 / (4 * pi^2 * int2)))
stopifnot(abs(1/(8*pi^2*DU*Xi) - 1/(4*pi^2*int2)) < 1e-14)

cat("\n=== 8. and with that measure the non-degenerate case is a scalar, as it must be ===\n")
cat("   dp/p is invariant under p -> c p, DU -> DU/c and Xi -> c Xi, so 1/(DU Xi) does not move:\n")
for (c0 in c(0.5, 2, 4)) {
  g <- 1 / (8 * pi^2 * (DU / c0) * (c0 * Xi))
  cat(sprintf("      c = %.2f:  G = %.10f\n", c0, g))
  stopifnot(abs(g - 1/(8*pi^2*DU*Xi)) < 1e-14)
}
cat("   The c^3 of section 4 is gone. The error was the measure and nothing else.\n")

cat("\n=== 9. what is left, now one step smaller ===\n")
cat("   With dp/(4 pi p) fixed, the degenerate integrand carries p^{1/2} rather than p^0, which\n")
cat("   is the extra half power the exponent rule wants. What still does not balance is the\n")
cat("   SCALING of the collective-coordinate Jacobian: taking it as sqrt(p int xi^2 dlambda)\n")
cat("   leaves the whole answer going as c^{-2} rather than standing still.\n")
cat("   IMPASSE: the Jacobian's scaling weight.\n")
cat("   LOAD-BEARING: that the Jacobian is sqrt(p int xi^2 dlambda) rather than something with\n")
cat("                 two more powers of the affine scale in it.\n")
cat("   EVIDENCE: with the measure now settled against the flat case, the residual mismatch is\n")
cat("             exactly c^2, so the Jacobian is short by that and by nothing else.\n")
cat("   NEAREST ROUTE: derive the Jacobian on a case with a compact symmetry zero mode whose\n")
cat("             propagator is known, rather than by choosing an inner product that looks\n")
cat("             invariant. A ring of radius R in flat space has exactly that structure.\n")
