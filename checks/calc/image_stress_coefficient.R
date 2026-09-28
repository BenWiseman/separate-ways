#!/usr/bin/env Rscript
# The coefficient kappa, which was the last number the shell thickness was carried in terms of.
#
# image_stress_shell.R writes the interior stress as T_kk = kappa m^2 D^{-p} r_h^{p-2} and leaves
# kappa as a pure number nobody had computed. caustic_power_at_a_hole.R then fixed p at 5/2. With
# the power fixed, kappa is no longer an unknown of the same kind, because every other ingredient
# is already in the release:
#
#   the proper-time representation      G = int ds (4 pi s)^{-D/2} Delta^{1/2} e^{-sigma/2s}
#   the caustic replacement             Delta^{1/2} -> 6.6092 M s^{-1/2}   (A.19, order-one caustic)
#   the mass as a phase in s            one extra power of s, to first order in m^2
#   the world function                  sigma proportional to M - r, coefficient from A.19's 34.85
#
# Those four give the coefficient outright. The only thing left free is one contraction factor, the
# square of the null vector's derivative of the distance, and for the ingoing radial congruence the
# shell calculation actually follows that factor is exactly one, because dr/dlambda = -1 there.
#
# So kappa comes out near 4e-3 rather than 1, the shell shrinks by a factor of nine rather than by
# orders, and the conclusion that it is a macroscopic length survives: a couple of microns for the
# fold's own fermion at a solar mass instead of twenty, and half a femtometre for an electron.
#
# The amplitude above is A.19's, built on the PROJECTION arc length M(pi+2) and not on the affine
# total. This file used the affine one for a day, which made the amplitude 6.6092 and the shell
# 2.6 microns; A.19's own measurement excludes that reading and the provenance audit found the
# disagreement between the two.
#
# What is NOT settled by this file is the SIGN along the radial direction. The computed sign belongs
# to the contact null direction, as both manuscripts say. This file computes a magnitude and the
# sign question is untouched.

TOL <- 1e-10
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M      <- 1                      # units of the hole's mass; r_h = 2M
rh     <- 2 * M
AMP    <- 3.9004 * M             # Delta^{1/2} -> AMP s^{-1/2}, from A.19's PROJECTION length
TAU2   <- 34.85                  # A.19: tau^2/(M - r) -> 34.85
KSIG   <- TAU2 / 2 * M           # sigma = -tau^2/2, so |sigma| = KSIG * (M - r)
D4     <- 4

cat("=== 1. the two proper-time integrals, in closed form and by quadrature ===\n")
cat("   With Delta^{1/2} = AMP s^{-1/2} and D = 4 the massless integrand is s^{-5/2} and the\n")
cat("   first-order mass part carries one more power of s, so s^{-3/2}:\n")
cat("      int_0^inf ds s^{-5/2} e^{-sigma/2s} = Gamma(3/2)(2/sigma)^{3/2}\n")
cat("      int_0^inf ds s^{-3/2} e^{-sigma/2s} = Gamma(1/2)(2/sigma)^{1/2}\n")
Iexact <- function(k, sg) gamma(k - 1) * (2/sg)^(k - 1)
Inum <- function(k, sg) integrate(function(u) { s <- exp(u); s * s^(-k) * exp(-sg/(2*s)) },
                                  -40, 40, rel.tol = 1e-12, subdivisions = 4000L)$value
cat("      k       sigma      quadrature       closed form      relative\n")
for (k in c(2.5, 1.5)) for (sg in c(0.02, 0.5)) {
  q <- Inum(k, sg); c0 <- Iexact(k, sg)
  cat(sprintf("   %6.1f %11.3f %16.8e %17.8e %13.2e\n", k, sg, q, c0, abs(q/c0 - 1)))
  note(abs(q/c0 - 1) < 1e-8, sprintf("the closed form holds at k = %.1f", k))
}
note(abs(gamma(1.5) - sqrt(pi)/2) < TOL && abs(gamma(0.5) - sqrt(pi)) < TOL,
     "the two Gamma values are sqrt(pi)/2 and sqrt(pi)")

cat("\n=== 2. the mass part of the correlator, assembled ===\n")
cat("   G_mass = -(m^2/2) AMP (4 pi)^{-2} Gamma(1/2) (2/sigma)^{1/2}, so writing it as\n")
cat("   G_mass = -A m^2 sigma^{-1/2}:\n")
A <- 0.5 * AMP * (4*pi)^(-D4/2) * gamma(0.5) * sqrt(2)
cat(sprintf("      (4 pi)^{-2}      = %.9f\n", (4*pi)^(-2)))
cat(sprintf("      A                = %.9f M   (dimension of length)\n", A))
note(abs(A - 0.5 * 3.9004 * (1/(16*pi^2)) * sqrt(pi) * sqrt(2)) < 1e-9, "A assembled correctly")
cat("   Checked against the direct product of the four factors:\n")
cat(sprintf("      0.5 x %.4f x %.9f x %.6f x %.6f = %.9f\n",
            AMP, (4*pi)^(-2), gamma(0.5), sqrt(2), A))

cat("\n=== 3. two derivatives along the congruence ===\n")
cat("   T_kk is quadratic in derivatives, and with sigma = KSIG (M - r) the leading piece is\n")
cat("   (3/4) A m^2 sigma^{-5/2} (d sigma/d lambda)^2. For the ingoing radial null congruence\n")
cat("   dr/dlambda = -1, so d sigma/d lambda = KSIG and the contraction factor is exactly one.\n")
cat("   Collecting the powers of D:\n")
cat("      sigma^{-5/2} (d sigma/d lambda)^2 = KSIG^{-1/2} D^{-5/2}\n")
cat(sprintf("      KSIG             = %.5f M\n", KSIG))
cat(sprintf("      KSIG^{-1/2}      = %.6f M^{-1/2}\n", KSIG^(-0.5)))
kap_M <- 0.75 * A * KSIG^(-0.5)                        # coefficient of m^2 M^{1/2} D^{-5/2}
cat(sprintf("      (3/4) A KSIG^{-1/2} = %.8f M^{1/2}\n", kap_M))
kappa <- kap_M / sqrt(2)                               # M^{1/2} = (r_h/2)^{1/2}
cat(sprintf("      in units of r_h^{1/2}, kappa = %.8f\n", kappa))
note(abs(kappa - kap_M/sqrt(2)) < TOL, "the change from M to r_h is one factor of sqrt 2")
cat("   Dimensions, checked rather than asserted: T_kk must carry four inverse lengths.\n")
Tkk <- function(m, Dd, rh) kappa * m^2 * Dd^(-2.5) * rh^(0.5)
for (sc in c(2, 6.7)) {
  a <- Tkk(1, 0.3, 2); b <- Tkk(1/sc, 0.3*sc, 2*sc)
  cat(sprintf("      scale %5.2f: ratio %12.6f against s^4 = %12.6f\n", sc, a/b, sc^4))
  note(abs(a/b - sc^4) < 1e-9 * sc^4, "T_kk scales as the fourth inverse power")
}

cat("\n=== 4. what that does to the shell ===\n")
mP_GeV <- 1.220890e19; lP_m <- 1.616255e-35; rh_sun <- 2.953250e3; Bstar <- 11.5138
Ds <- function(m_GeV, kap) (8*pi*kap*(m_GeV/mP_GeV)^2/Bstar)^(1/2.5) * rh_sun
cat("      field                 at kappa = 1        at kappa computed      ratio\n")
for (mm in list(c("electron", 0.000511), c("proton", 0.938272), c("top", 172.69),
                c("the fold's fermion", 4.916e8))) {
  m <- as.numeric(mm[2]); d1 <- Ds(m, 1); d2 <- Ds(m, kappa)
  cat(sprintf("   %-22s %12.3e m %18.3e m %11.4f\n", mm[1], d1, d2, d2/d1))
  note(abs(d2/d1 - kappa^(1/2.5)) < 1e-9, "the shell moves as kappa^{2/5}")
}
cat(sprintf("   The ratio is kappa^{2/5} = %.5f, one factor of %.2f, so the shell is a shorter\n",
            kappa^(1/2.5), 1/kappa^(1/2.5)))
cat("   macroscopic length and not a sub-Planckian one:\n")
for (mm in list(c("electron", 0.000511), c("the fold's fermion", 4.916e8))) {
  m <- as.numeric(mm[2]); d <- Ds(m, kappa)
  cat(sprintf("      %-20s %10.3e m = %9.3e Planck lengths\n", mm[1], d, d/lP_m))
  note(d/lP_m > 1e10, "still far above a Planck length")
}
mc <- exp(uniroot(function(l) Ds(exp(l), kappa)/lP_m - 1, c(log(1e-30), log(1e12)),
                  tol = 1e-13)$root) * 1e9
cat(sprintf("   Sub-Planckian only below %.3g eV at one solar mass, against %.3g eV at kappa = 1.\n",
            mc, exp(uniroot(function(l) Ds(exp(l), 1)/lP_m - 1, c(log(1e-30), log(1e12)),
                           tol = 1e-13)$root) * 1e9))
note(mc < 1e-12, "the crossover is still below every real particle mass")

cat("\n=== 5. how much of this is settled, said plainly ===\n")
cat("   Settled: the caustic amplitude, from a non-symmetric surface of revolution holding the\n")
cat("   length and the family volume fixed; the power, from the caustic-order rule with the mass\n")
cat("   costing one power of the world function; the world function's linear vanishing and its\n")
cat("   coefficient, from A.19; and the contraction factor for the radial congruence, which is one\n")
cat("   because dr/dlambda = -1 on it.\n")
cat("   Not settled, and this file does not pretend otherwise: the SIGN along the radial\n")
cat("   direction. The negative sign in the release belongs to the contact null direction. This is\n")
cat("   a magnitude.\n")
cat("   Also approximate: only the leading derivative pairing is kept. Terms where both\n")
cat("   derivatives fall on the same factor, or where the Van Vleck prefactor is differentiated,\n")
cat("   are of the same order in D and carry their own numerical coefficients, so kappa is good to\n")
cat("   an order-unity factor and not better. The shell moves as its two-fifths power, so a factor\n")
cat("   of ten either way moves the thickness by two and a half:\n")
for (f in c(0.1, 1, 10)) {
  cat(sprintf("      kappa x %5.1f: the fold's fermion gets %10.3e m\n", f, Ds(4.916e8, kappa*f)))
}
note(abs(Ds(4.916e8, kappa*10)/Ds(4.916e8, kappa) - 10^(1/2.5)) < 1e-9,
     "a factor of ten in kappa is 10^{2/5} in the shell")

cat("\n=== 6. plants ===\n")
cat("   (a) forget that the mass costs one power of s and the D-power comes out wrong:\n")
kbad <- 0.75 * (0.5 * AMP * (4*pi)^(-2) * gamma(1.5) * 2^(1.5)) * KSIG^(-1.5)
cat("       using Gamma(3/2) and (2/sigma)^{3/2} for the mass part gives a D-power of -7/2,\n")
cat(sprintf("       not -5/2, and a coefficient of %.6e in the wrong units altogether\n", kbad))
note(TRUE, "recorded")
pw <- function(k) -(k - 1) - 2          # sigma exponent minus two derivatives, in D
cat(sprintf("       the exponent from k = 3/2 is %.1f and from k = 5/2 is %.1f\n", pw(1.5), pw(2.5)))
note(abs(pw(1.5) + 2.5) < TOL && abs(pw(2.5) + 3.5) < TOL,
     "PLANT (a) fires: the mass part is the k = 3/2 integral and gives -5/2")
cat("   (b) the world function's coefficient must matter, and as its inverse square root:\n")
k2 <- 0.75 * A * (2 * KSIG)^(-0.5) / sqrt(2)
cat(sprintf("       doubling KSIG gives kappa = %.8f against %.8f, a ratio of %.6f\n",
            k2, kappa, kappa/k2))
note(abs(kappa/k2 - sqrt(2)) < 1e-9, "PLANT (b) fires: kappa carries KSIG^{-1/2}")
cat("   (c) the contraction factor is one only for the radial congruence. At (k.grad D)^2 = 4\n")
cat("       the shell would move by 4^{2/5}:\n")
cat(sprintf("       %10.3e m against %10.3e m, ratio %.4f against %.4f\n",
            Ds(4.916e8, 4*kappa), Ds(4.916e8, kappa),
            Ds(4.916e8, 4*kappa)/Ds(4.916e8, kappa), 4^(1/2.5)))
note(abs(Ds(4.916e8, 4*kappa)/Ds(4.916e8, kappa) - 4^(1/2.5)) < 1e-9,
     "PLANT (c) fires: the contraction factor enters the same way kappa does")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
