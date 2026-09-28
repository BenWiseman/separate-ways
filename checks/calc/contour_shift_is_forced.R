#!/usr/bin/env Rscript
# The closed time path has one free parameter. The fold fixes it.
#
# In the real-time thermal formalism the second leg of the contour is displaced downward in
# imaginary time by an amount sigma, anywhere in [0, beta]. Every choice returns the same physical
# in-in correlators, so sigma is a convention: Landsman and van Weert call it exactly that, and
# sigma = beta/2 is the symmetric choice that Niemi and Semenoff and later Herzog and Son work
# with. Section 2 hypothesises that the fold's image is the second leg. That hypothesis is only
# as narrow as the contour it lands on, so it matters whether beta/2 is chosen or forced.
#
# It is forced, and by the one property the fold has to have. Theta is antilinear with Theta^2 = 1,
# so swapping the two sheets is a symmetry of the pair and the cross-sheet kernel cannot tell the
# sheets apart. There are two cross kernels on the contour, one for each order of the legs. They
# agree for all separations at exactly one sigma.
#
# Written out, with W the thermal Wightman function of a single mode and the branch-2 field at
# t - i sigma,
#     G12(t) = W(-t - i sigma),   G21(t) = W(t - i sigma),
#     G12(t) - G21(t) = i sin(omega t) sinh(omega (beta/2 - sigma)) / (omega sinh(beta omega / 2)),
# which vanishes identically in t if and only if sigma = beta/2. Nothing was assumed about the
# temperature or the mode: the Bose factor does the work, so this is the KMS condition choosing the
# contour rather than a convenient convention.
#
# Two further things follow and both are checked here. The fold-even and fold-odd amplitudes are
# coth(beta omega / 4) and tanh(beta omega / 4), whose product is 1 at every frequency and every
# temperature, so the two halves of the fold are reciprocal about the vacuum. And the
# quantum-quantum correlator of the symmetric contour does NOT vanish: it equals
# tanh(beta omega / 4) cos(omega t) / omega, and the familiar statement that it vanishes belongs to
# sigma -> 0. Section 3 must not claim otherwise, so the limit is checked here too.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) {
  if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 }
  invisible(ok)
}

nbose <- function(bw) 1 / expm1(bw)

# W(t) for one mode, complex, at inverse temperature beta = bw/omega. Argument may be complex.
W <- function(t, om, bw) {
  n <- nbose(bw)
  ((1 + n) * exp(-1i * om * t) + n * exp(1i * om * t)) / (2 * om)
}

cat("=== 1. the Wightman function obeys KMS, as a check on the algebra above ===\n")
cat("   W(t - i beta) = W(-t) is the KMS condition. Tested on a grid of real t.\n")
cat("      omega   beta*omega      max |W(t - i beta) - W(-t)|\n")
for (om in c(0.7, 1, 2.5)) for (bw in c(0.4, 2 * log(3), 5)) {
  beta <- bw / om
  ts <- seq(-4, 4, length.out = 81)
  d <- max(Mod(W(ts - 1i * beta, om, bw) - W(-ts, om, bw)))
  cat(sprintf("   %7.2f %12.5f %32.3e\n", om, bw, d))
  note(d < 1e-11, "KMS")
}

cat("\n=== 2. the two cross kernels agree at exactly one sigma ===\n")
cat("   G12(t) = W(-t - i sigma) and G21(t) = W(t - i sigma). Their difference is measured\n")
cat("   against the closed form i sin(om t) sinh(om (beta/2 - sigma)) / (om sinh(beta om / 2)).\n")
cat("      beta*omega   sigma/beta     max|G12-G21|      closed form      agree to\n")
for (bw in c(0.4, 2 * log(3), 5)) {
  om <- 1.3; beta <- bw / om
  for (frac in c(0, 0.25, 0.5, 0.75, 1)) {
    sg <- frac * beta
    ts <- seq(-3, 3, length.out = 241)
    dif <- W(-ts - 1i * sg, om, bw) - W(ts - 1i * sg, om, bw)
    cf <- 1i * sin(om * ts) * sinh(om * (beta / 2 - sg)) / (om * sinh(bw / 2))
    cat(sprintf("   %11.5f %12.3f %16.6e %16.6e %13.3e\n",
                bw, frac, max(Mod(dif)), max(Mod(cf)), max(Mod(dif - cf))))
    note(max(Mod(dif - cf)) < 1e-13, "closed form for the cross-kernel gap")
    if (abs(frac - 0.5) < 1e-14) note(max(Mod(dif)) < TOL, "gap closes at sigma = beta/2")
    if (abs(frac - 0.5) > 0.2)    note(max(Mod(dif)) > 1e-3, "gap open away from beta/2")
  }
}
cat("   The gap is odd in t and pure imaginary, so no choice of overall phase hides it. It is\n")
cat("   proportional to sinh(om (beta/2 - sigma)): one zero, and it is where the fold sits.\n")

cat("   At sigma = 0 the gap is the commutator. i sin(om t) sinh(om beta / 2) / (om sinh(beta om\n")
cat("   / 2)) is i sin(om t) / om whatever the temperature, which is [phi(t), phi(0)] itself:\n")
for (bw in c(0.4, 2 * log(3), 5)) {
  om <- 1.3; ts <- seq(-3, 3, length.out = 241)
  g0 <- W(-ts, om, bw) - W(ts, om, bw)
  cat(sprintf("      beta*om = %7.4f   max|G12 - G21 - i sin(om t)/om| = %10.3e\n",
              bw, max(Mod(g0 - 1i * sin(om * ts) / om))))
  note(max(Mod(g0 - 1i * sin(om * ts) / om)) < 1e-13, "sigma = 0 gap is the commutator")
}
cat("   So the two ends of the family are the two extremes of the same quantity: at sigma = 0 the\n")
cat("   cross orderings differ by the whole commutator, and at beta/2 they do not differ at all.\n")

cat("\n=== 3. how sharply beta/2 is selected, mode by mode ===\n")
cat("   The gap relative to the direct correlator at coincidence, with sigma moved off beta/2 by\n")
cat("   one per cent of beta. The selection is sharpest for modes at the temperature itself and\n")
cat("   weakens at both ends, since a cold mode has nothing thermal to be asymmetric about and a\n")
cat("   hot one has its asymmetry swamped by the size of the correlator:\n")
cat("      beta*omega    relative gap at sigma = 0.49 beta\n")
rg <- c()
for (bw in c(0.5, 1, 2, 5, 10, 20)) {
  om <- 1; beta <- bw / om; sg <- 0.49 * beta
  g <- Mod(sinh(om * (beta / 2 - sg)) / (om * sinh(bw / 2)))
  rg <- c(rg, g / Mod(W(0, om, bw)))
  cat(sprintf("   %11.3f %38.6e\n", bw, g / Mod(W(0, om, bw))))
}
note(which.max(rg) > 1 && which.max(rg) < length(rg), "the selection peaks in the middle")
cat(sprintf("   Peak of those six at beta*omega = %.1f, and every entry is nonzero, so no mode is\n",
            c(0.5, 1, 2, 5, 10, 20)[which.max(rg)]))
cat("   indifferent to the displacement.\n")

cat("\n=== 4. the fold-even and fold-odd amplitudes, and their product ===\n")
cat("   At sigma = beta/2 the equal-time matrix is [[W(0), K], [K, W(0)]] with K = W(-i beta/2),\n")
cat("   so the fold-even and fold-odd combinations are its eigenvectors with eigenvalues\n")
cat("   W(0) + K and W(0) - K. Those are coth(beta om / 4) / 2 om and tanh(beta om / 4) / 2 om.\n")
cat("      beta*omega     even (num)     coth/2om      odd (num)      tanh/2om     2om^2*product\n")
for (bw in c(0.3, 1, 2 * log(3), 4, 9)) {
  om <- 1.1; beta <- bw / om
  K <- Re(W(-1i * beta / 2, om, bw)); W0 <- Re(W(0, om, bw))
  ev <- W0 + K; od <- W0 - K
  cat(sprintf("   %11.5f %14.8f %13.8f %14.8f %13.8f %15.10f\n",
              bw, ev, 1 / tanh(bw / 4) / (2 * om), od, tanh(bw / 4) / (2 * om),
              ev * od * 4 * om^2))
  note(abs(ev - 1 / tanh(bw / 4) / (2 * om)) < 1e-12, "fold-even eigenvalue")
  note(abs(od - tanh(bw / 4) / (2 * om)) < 1e-12, "fold-odd eigenvalue")
  note(abs(ev * od * 4 * om^2 - 1) < 1e-12, "the product is the vacuum value squared")
}
cat("   The last column is 1 to twelve figures at every temperature. The two halves of the fold\n")
cat("   are reciprocal, and their geometric mean is the zero-temperature amplitude 1/2 omega.\n")
cat("   The manuscript normalises Phi_c = (Phi + Theta Phi)/2 and Phi_q = Phi - Theta Phi, so its\n")
cat("   two correlators are half the first eigenvalue and twice the second, and the product is\n")
cat("   untouched. Their ratio is 4 tanh^2(beta om / 4), which is the manuscript's number; it is 1\n")
cat("   at beta om = 2 ln 3, where the classical and quantum halves carry equal weight. That is a\n")
cat("   different frequency from the crossing ceiling at beta om = ln 3, and both are checked:\n")
for (bw in c(log(3), 2 * log(3))) {
  om <- 1.0
  K <- Re(W(-1i * (bw / om) / 2, om, bw)); W0 <- Re(W(0, om, bw))
  ratio <- 4 * (W0 - K) / (W0 + K)
  cat(sprintf("      beta om = %8.5f   4 tanh^2(bw/4) = %10.7f   from the matrix = %10.7f   tanh r = %8.5f\n",
              bw, 4 * tanh(bw / 4)^2, ratio, exp(-bw / 2)))
  note(abs(ratio - 4 * tanh(bw / 4)^2) < 1e-12, "the manuscript ratio from the matrix")
}
note(abs(4 * tanh(2 * log(3) / 4)^2 - 1) < 1e-14, "the halves balance at beta om = 2 ln 3")
note(abs(exp(-log(3) / 2) - tanh(asinh(1 / sqrt(2)))) < 1e-14,
     "the ceiling squeeze sits at beta om = ln 3")

cat("\n=== 5. the quantum-quantum correlator does not vanish at the symmetric contour ===\n")
cat("   G^qq = G11 + G22 - G12 - G21 with G11 time-ordered and G22 anti-time-ordered. Compared\n")
cat("   against tanh(beta om / 4) cos(om t) / om, and then followed down to sigma = 0.\n")
qq <- function(ts, om, bw, sg) {
  n <- nbose(bw)
  Wt <- W(ts, om, bw); Wmt <- W(-ts, om, bw)
  g11 <- ifelse(ts >= 0, Wt, Wmt)          # time ordered
  g22 <- ifelse(ts >= 0, Wmt, Wt)          # anti-time ordered
  g12 <- W(-ts - 1i * sg, om, bw); g21 <- W(ts - 1i * sg, om, bw)
  g11 + g22 - g12 - g21
}
om <- 1.0
for (bw in c(0.5, 2 * log(3), 6)) {
  beta <- bw / om; ts <- seq(-3, 3, length.out = 121)
  got <- qq(ts, om, bw, beta / 2)
  want <- tanh(bw / 4) * cos(om * ts) / om
  cat(sprintf("   beta*om = %7.4f   max|G^qq - tanh(bw/4)cos(om t)/om| = %10.3e   size %8.5f\n",
              bw, max(Mod(got - want)), max(Mod(want))))
  note(max(Mod(got - want)) < 1e-12, "closed form for G^qq at the symmetric contour")
  note(max(Mod(want)) > 0.1, "G^qq is not small")
  down <- sapply(c(0.5, 0.1, 0.01, 1e-4, 0) * beta, function(s) max(Mod(qq(ts, om, bw, s))))
  cat(sprintf("     sigma/beta 0.5, 0.1, 0.01, 1e-4, 0 -> %s\n",
              paste(sprintf("%9.3e", down), collapse = " ")))
  note(down[5] < TOL, "G^qq vanishes at sigma = 0")
  note(all(diff(down) < 0), "and falls monotonically on the way there")
}
cat("   So the vanishing quantum-quantum correlator is a statement about sigma = 0, not about the\n")
cat("   contour the fold picks. At the fold's contour it is tanh(beta om / 4) / om at coincidence,\n")
cat("   which is the same number as the fold-odd amplitude of section 4, and it goes to zero in\n")
cat("   the classical direction beta om -> 0 rather than being zero throughout.\n")

cat("\n=== 6. plants: each of the three results is made to fail on purpose ===\n")
cat("   (a) Drop the Bose factor, keeping only the vacuum piece. Then nothing selects beta/2,\n")
cat("       because the selection came from (1+n) e^{-om sigma} = n e^{om sigma}.\n")
Wvac <- function(t, om) exp(-1i * om * t) / (2 * om)
bw <- 2 * log(3); om <- 1.3; beta <- bw / om
ts <- seq(-3, 3, length.out = 241)
gvac <- max(Mod(Wvac(-ts - 1i * beta / 2, om) - Wvac(ts - 1i * beta / 2, om)))
gvac0 <- max(Mod(Wvac(-ts - 1i * 0.2 * beta, om) - Wvac(ts - 1i * 0.2 * beta, om)))
cat(sprintf("       vacuum-only gap at beta/2 = %.6e, at 0.2 beta = %.6e  (ratio %.4f)\n",
            gvac, gvac0, gvac / gvac0))
note(gvac > 1e-3, "PLANT (a) fires: with no Bose factor beta/2 is not special")
cat("   (b) Use the wrong shift, sigma = beta/3, and ask for the gap to close.\n")
g3 <- max(Mod(W(-ts - 1i * beta / 3, om, bw) - W(ts - 1i * beta / 3, om, bw)))
cat(sprintf("       gap at beta/3 = %.6e against %.1e at beta/2\n", g3, TOL))
note(g3 > 1e-3, "PLANT (b) fires: beta/3 leaves a gap")
cat("   (c) Reciprocity with the half-angle written as beta om / 2 instead of beta om / 4.\n")
bad <- (1 / tanh(bw / 2)) * tanh(bw / 2)
cat(sprintf("       coth(bw/4) tanh(bw/4) = %.12f, coth(bw/2) tanh(bw/4) = %.12f\n",
            (1 / tanh(bw / 4)) * tanh(bw / 4), (1 / tanh(bw / 2)) * tanh(bw / 4)))
note(abs((1 / tanh(bw / 2)) * tanh(bw / 4) - 1) > 1e-3,
     "PLANT (c) fires: the wrong half-angle breaks the product")
note(abs(bad - 1) < 1e-15, "control: the identity itself is exact")

cat("\n=== 7. what this does to the ledger ===\n")
cat("   Before: the fold's image is hypothesised to be the second leg of a contour whose\n")
cat("   displacement is a convention, and beta/2 is adopted because it is the symmetric one.\n")
cat("   After:  Theta^2 = 1 admits one displacement out of the family, and it is beta/2. The\n")
cat("   hypothesis is now only that the legs are paired at all; where the second leg sits is a\n")
cat("   consequence. One fewer choice, and it was the only free parameter the contour had.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
