#!/usr/bin/env Rscript
# How thick the shell is in which the fold's term beats the interior's own focusing, and a
# correction to the way interior_focusing_threshold.R posed the question.
#
# That file computed two thresholds on a CONSTANT R_kk: 8/r_h^2 to turn the expansion around at
# the edge of the contact region, and 11.5138/r_h^2 to bring it to zero before the congruence
# leaves a region of extent r_h/2. Both are right about a constant. Neither is the instrument the
# geometry asks for, and this file says why and replaces them.
#
# The contact locus is the sphere r = M = r_h/2 and the world function between a point and its
# image vanishes LINEARLY there, so the image term diverges as the contact boundary is approached
# from inside. An ingoing congruence therefore meets the largest stress at the moment it enters the
# contact region and a decreasing one thereafter, which is the opposite of the constant the earlier
# file integrated. Two consequences. The integrated kick int B dlambda is divergent at the
# boundary, so the vacuum's focusing is not the leading behaviour there and no finite threshold on
# B is the question. And what is worth a number instead is the THICKNESS of the shell inside which
# the fold's term dominates: outside it the interior focuses as general relativity says, and inside
# it the fold's term is the larger of the two. Which way it pushes the focusing is a separate
# question, and contact_and_radial_are_orthogonal.R shows it is not the one the computed sign
# answers. A thickness is a magnitude and none of this file turns on that.
#
# The power is 5/2, and section 6 gives the three links that fix it. The bracket this file used to
# carry, 3/2 to 7/2, is kept in the tables below because it shows how much the power decides and how
# little kappa does. The bracket came from a rule checked against an independent stress computation. antipodal_sphere_family.R
# measures the image Green function at a caustic as delta^{-(D-2+n)/2} in the DISTANCE delta to it,
# over four members of a family. Two derivatives make a stress, so the stress goes as
# delta^{-(D+2+n)/2}: on the Einstein static universe that is n = 2, D = 4 and delta^{-4}, and
# sign_through_the_mass_slot.R's coupling slot comes out at delta^{-4.010} once its fitted power is
# read as a power of 1 + cos tau, which near the caustic is the SQUARE of the distance. The rule and
# the stress therefore agree to within a per cent on the one geometry where both are available.
# At a hole the caustic is order one, n = 1, so the same rule gives delta^{-7/2}. The one thing not
# explained is that the mass slot on the same geometry is two powers softer, delta^{-2.037}, and the
# mass slot is the one that acts at a Ricci-flat hole. So the bracket carried below is
#
#     p = 7/2   the rule applied at n = 1, validated on the coupling slot
#     p = 3/2   the same, less the two powers the mass slot is observed to be softer by
#
# with 5/2 between them, and 5/2 is the answer. The form is fixed by dimensions either way:
#
#     T_kk = kappa m^2 D^{-p} r_h^{p-2},   B = -R_kk = 8 pi G T_kk,
#
# with D the distance to the contact sphere, the dimensions fixed by T_kk having four inverse
# lengths and m^2 two of them, and kappa the pure number A.18's sum on the contact geodesic would
# supply. Setting B against the earlier file's 11.5138/r_h^2 gives
#
#     D*/r_h = (8 pi kappa (m/m_P)^2 / 11.5138)^{1/p},
#
# so kappa enters as its p-th root and hardly matters, while p decides everything.

TOL <- 1e-10
fail <- 0
note <- function(ok, what) {
  if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 }
  invisible(ok)
}

mP_GeV <- 1.220890e19          # Planck mass
lP_m   <- 1.616255e-35         # Planck length
rh_sun <- 2.953250e3           # Schwarzschild radius of one solar mass, metres
Bstar  <- 11.5138              # in units of 1/r_h^2, from interior_focusing_threshold.R

cat("=== 1. the shell thickness, in closed form and in units of r_h ===\n")
cat("   D*/r_h = (8 pi kappa (m/m_P)^2 / 11.5138)^{1/p}. Checked against the defining equation\n")
cat("   8 pi G kappa m^2 D^{-p} r_h^{p-2} = 11.5138/r_h^2 at three powers:\n")
Dstar_over_rh <- function(m_GeV, p, kappa = 1) (8 * pi * kappa * (m_GeV / mP_GeV)^2 / Bstar)^(1 / p)
cat("      p        (m/m_P)^2        D*/r_h        residual of the defining equation\n")
for (p in c(1.5, 2.5, 3.5)) {
  m <- 4.916e8                                   # GeV, the paper's 491.6 PeV fermion
  x <- Dstar_over_rh(m, p); gm2 <- (m / mP_GeV)^2
  resid <- 8 * pi * gm2 * x^(-p) - Bstar         # in units of 1/r_h^2, with r_h = 1
  cat(sprintf("   %6.2f %16.4e %14.6e %36.2e\n", p, gm2, x, resid))
  note(abs(resid) < 1e-8 * Bstar, sprintf("the closed form solves its own equation at p = %.1f", p))
}

cat("\n=== 2. dimensions, since the whole result is a dimensional argument ===\n")
cat("   T_kk must carry four inverse lengths and m^2 carries two, so the remaining factor is\n")
cat("   D^{-p} r_h^{p-2} and no other combination of the two lengths does it. Checked by\n")
cat("   scaling every length by a factor and requiring T_kk to scale as the fourth power:\n")
Tkk <- function(m, D, rh, p) m^2 * D^(-p) * rh^(p - 2)
for (p in c(1.5, 3.5)) for (s in c(2, 7.3)) {
  a <- Tkk(1, 0.3, 1.7, p); b <- Tkk(1 / s, 0.3 * s, 1.7 * s, p)
  cat(sprintf("      p = %4.1f, scale = %5.2f:  T_kk goes from %12.6f to %12.6f, ratio %10.6f\n",
              p, s, a, b, a / b))
  note(abs(a / b - s^4) < 1e-9 * s^4, "T_kk scales as the fourth inverse power of length")
}
cat("   A wrong exponent on r_h fails the same test, which is the plant in section 5.\n")

cat("\n=== 3. the shell in Planck lengths, which is the question that decides anything ===\n")
cat("   A shell thinner than a Planck length is a statement about nothing. Listed at kappa = 1,\n")
cat("   for one solar mass, at the two ends of the old bracket, to show how much the power decides:\n")
masses <- list(c("electron", 0.000511), c("muon", 0.10566), c("proton", 0.938272),
               c("b quark", 4.18), c("W", 80.377), c("Higgs", 125.25), c("top", 172.69),
               c("the fold's fermion", 4.916e8))
cat("   The last row is the 491.6 PeV dark-matter fermion of the cosmology paper.\n")
cat("      field                    m (GeV) D*/l_P p=3/2 D* (m) p=3/2  D*/l_P p=7/2  D* (m) p=7/2\n")
for (mm in masses) {
  m <- as.numeric(mm[2])
  d1 <- Dstar_over_rh(m, 1.5) * rh_sun; d2 <- Dstar_over_rh(m, 3.5) * rh_sun
  cat(sprintf("   %-22s %10.3e %12.3e %12.3e %13.3e %13.3e\n",
              mm[1], m, d1 / lP_m, d1, d2 / lP_m, d2))
}
note(abs(Dstar_over_rh(4.916e8, 3.5) * rh_sun / 4.24e-3 - 1) < 0.02,
     "the fermion's p = 7/2 shell is four millimetres")
note(abs(Dstar_over_rh(0.000511, 3.5) * rh_sun / 6.02e-10 - 1) < 0.02,
     "an electron's p = 7/2 shell is six Angstrom")
cat("   Every field on the list is above a Planck length at both ends of the bracket, which is the\n")
cat("   result: a root of a very small number is not a very small number. At p = 7/2 the shell is\n")
cat("   macroscopic, four millimetres for the fold's own fermion and six Angstrom for an electron.\n")
cat("   Across the bracket, for one solar mass:\n")
cat("      p        electron          the fold's fermion\n")
for (pp in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
  cat(sprintf("   %6.3f %14.3e m %20.3e m\n", pp,
              Dstar_over_rh(0.000511, pp) * rh_sun, Dstar_over_rh(4.916e8, pp) * rh_sun))
}
cat("   The spread across the bracket was what had to be closed, and it is larger than anything\n")
cat("   kappa could do: a factor of ten in kappa moves D* by 4.6 at p = 3/2 and by 1.9 at p = 7/2.\n")
note(abs(Dstar_over_rh(4.916e8, 1.5, 10) / Dstar_over_rh(4.916e8, 1.5, 1) - 10^(1/1.5)) < 1e-9 &&
     abs(Dstar_over_rh(4.916e8, 3.5, 10) / Dstar_over_rh(4.916e8, 3.5, 1) - 10^(1/3.5)) < 1e-9,
     "kappa enters as kappa^{1/p}")

mth <- uniroot(function(lm) Dstar_over_rh(exp(lm), 1.5) * rh_sun / lP_m - 1,
               c(log(1e-24), log(1e12)), tol = 1e-13)$root
cat(sprintf("   The crossover to sub-Planckian at p = 3/2 is m = %.4g GeV for one solar mass, so\n",
            exp(mth)))
cat(sprintf("   %.3g eV: even a field at the neutrino-mass scale has a shell above a Planck length.\n",
            exp(mth) * 1e9))
note(abs(Dstar_over_rh(exp(mth), 1.5) * rh_sun / lP_m - 1) < 1e-6, "the crossover mass solves D* = l_P")

cat("\n=== 4. and how it moves with the hole, which is the direction nobody would guess ===\n")
cat("   D* is proportional to r_h and to nothing else dimensionful, so the shell measured in\n")
cat("   Planck lengths grows in direct proportion to the hole, and the crossover mass falls as\n")
cat("   r_h^{-p/2}:\n")
cat("      hole mass (M_sun)   crossover m (eV) at p = 3/2    D*/l_P for the fold's fermion\n")
for (Ms in c(1, 10, 1e6, 1e9)) {
  rh <- rh_sun * Ms
  f <- function(lm) Dstar_over_rh(exp(lm), 1.5) * rh / lP_m - 1
  mc <- exp(uniroot(f, c(log(1e-30), log(1e12)), tol = 1e-13)$root)
  cat(sprintf("   %18.0e %29.4e %32.3e\n", Ms, mc * 1e9, Dstar_over_rh(4.916e8, 1.5) * rh / lP_m))
}
cat("   So the effect is largest at the largest holes, which is the reverse of the usual\n")
cat("   expectation that quantum corrections matter most for small ones. The reason is that the\n")
cat("   threshold it has to beat, 11.5138/r_h^2, is itself weaker for a bigger hole.\n")

cat("\n=== 5. plants ===\n")
cat("   (a) the wrong power of r_h must fail the scaling test of section 2:\n")
bad <- function(m, D, rh, p) m^2 * D^(-p) * rh^(p - 1)
s <- 3.1; a <- bad(1, 0.3, 1.7, 1); b <- bad(1 / s, 0.3 * s, 1.7 * s, 1)
cat(sprintf("       r_h^{p-1} instead of r_h^{p-2}: ratio %0.6f against s^4 = %0.6f\n",
            a / b, s^4))
note(abs(a / b - s^4) > 0.1 * s^4, "PLANT (a) fires: the wrong exponent breaks the dimensions")
cat("   (b) dropping G, which is what writing m^2 where the equation wants G m^2 amounts to. The\n")
cat("       shell has to lie inside the contact region, D* <= r_h/2, for any field the semiclassical\n")
cat("       treatment applies to. With G in place it does; without it the shell is outside the hole\n")
cat("       by seventeen orders of magnitude, which is the kind of answer a missing G gives:\n")
withG <- Dstar_over_rh(4.916e8, 1)
noG    <- (8 * pi * (4.916e8)^2 / Bstar)^(1 / 1)
cat(sprintf("       with G: D*/r_h = %.4e, inside the region.  Without: %.4e\n", withG, noG))
note(withG < 0.5 && noG > 1e10, "PLANT (b) fires: a missing G puts the shell outside the hole")
cat("   (c) the crossover must move the right way with the hole: doubling r_h must LOWER it, and\n")
cat("       by 2^{1/p} in the mass squared, so by 2^{1/2p} in the mass.\n")
pp <- 1.5
f1 <- function(lm) Dstar_over_rh(exp(lm), pp) * rh_sun / lP_m - 1
f2 <- function(lm) Dstar_over_rh(exp(lm), pp) * 2 * rh_sun / lP_m - 1
m1 <- exp(uniroot(f1, c(log(1e-30), log(1e12)), tol = 1e-13)$root)
m2 <- exp(uniroot(f2, c(log(1e-30), log(1e12)), tol = 1e-13)$root)
cat(sprintf("       r_h: %0.5g eV, 2 r_h: %0.5g eV, ratio %0.6f against 2^{p/2} = %0.6f\n",
            m1 * 1e9, m2 * 1e9, m1 / m2, 2^(pp / 2)))
note(abs(m1 / m2 - 2^(pp / 2)) < 1e-4, "PLANT (c) fires the right way: the crossover falls as r_h^{-p/2}")

cat("\n=== 6. what this does to the earlier file, and what is still open ===\n")
cat("   interior_focusing_threshold.R's 11.5138/r_h^2 is retained as what it is, the requirement\n")
cat("   on a CONSTANT stress over a region of extent r_h/2, and it is used above as the level the\n")
cat("   real stress has to reach. What is withdrawn is any reading of it as the condition for the\n")
cat("   focusing to be stopped, because the stress is not constant: it diverges at the contact\n")
cat("   boundary, which is where an ingoing congruence enters, so the integrated effect is\n")
cat("   divergent there and the vacuum's focusing is not the leading behaviour. The question that\n")
cat("   survives is the thickness of the shell in which the fold's term is the larger one, and\n")
cat("   that is what this file computes.\n")
cat("   THE BRACKET IS CLOSED. caustic_power_at_a_hole.R settles p at 5/2 by measuring the three\n")
cat("   links the spread rested on: the leading stress power is -(D+2+n)/2, verified on this\n")
cat("   geometry; the mass part of every quantity is one power softer than its massless\n")
cat("   counterpart, because differentiating a mode sum with respect to m^2 brings an extra\n")
cat("   inverse power of the mode frequency; and the further power the mass slot appeared to lose\n")
cat("   belongs to the coupling, since it comes off only when 1 - 4 xi equals dQ/dP at the caustic,\n")
cat("   which on a conformally flat geometry happens to be the conformal coupling and at a hole is\n")
cat("   a fine-tuning to a number nobody has computed. At a hole n = 1, so the leading stress is\n")
cat("   -7/2 and the mass part is -5/2.\n")
cat("      field                          D* at p = 5/2      in Planck lengths\n")
for (mm in masses) {
  m <- as.numeric(mm[2]); d <- Dstar_over_rh(m, 2.5) * rh_sun
  cat(sprintf("   %-28s %14.3e m %19.3e\n", mm[1], d, d / lP_m))
}
mc <- exp(uniroot(function(l) Dstar_over_rh(exp(l), 2.5) * rh_sun / lP_m - 1,
                  c(log(1e-30), log(1e12)), tol = 1e-13)$root) * 1e9
cat(sprintf("   Sub-Planckian only below %.3g eV at one solar mass, so the shell is a real length\n", mc))
cat("   for every field there is. The fold's own fermion gets twenty microns.\n")
note(mc < 1e-15, "the crossover at p = 5/2 is below every real particle mass")
note(abs(Dstar_over_rh(4.916e8, 2.5) * rh_sun / 1.949e-5 - 1) < 0.02,
     "the fermion's shell at p = 5/2 is twenty microns")

cat("\n   AND kappa IS NO LONGER OPEN EITHER. image_stress_coefficient.R assembles it from the four\n")
cat("   ingredients the release already had: the caustic amplitude 3.9004 M s^{-1/2}, the\n")
cat("   proper-time integrals with the mass as one extra power of s, the world function's linear\n")
cat("   vanishing with A.19's coefficient, and a contraction factor that is exactly one on the\n")
cat("   ingoing radial congruence because dr/dlambda = -1 there. It comes out 0.0039329, so the\n")
cat("   shell is 0.1091 of its kappa = 1 value, a factor of 9.2 shorter and still macroscopic:\n")
KAP <- 0.0039329
for (mm in masses) {
  m <- as.numeric(mm[2]); d <- Dstar_over_rh(m, 2.5, KAP) * rh_sun
  cat(sprintf("      %-22s %12.3e m %16.3e Planck lengths\n", mm[1], d, d / lP_m))
}
note(abs(Dstar_over_rh(4.916e8, 2.5, KAP) * rh_sun / 2.127e-6 - 1) < 0.01,
     "the fermion's shell at the computed kappa is 2.1 microns")
cat("   What remains is an order-unity factor from the sub-leading derivative pairings and, for the\n")
cat("   radial direction, the SIGN, which contact_and_radial_are_orthogonal.R shows is a separate\n")
cat("   combination of the stress from the one the release computes, with conservation pointing the\n")
cat("   two the opposite way. Nothing above depends on it, a thickness being a magnitude. kappa\n")
cat("   enters as its 2/5 power, so a factor of ten either way moves the fermion's shell between\n")
cat("   0.85 and 5.3 microns:\n")
for (f in c(0.1, 1, 10))
  cat(sprintf("      kappa x %5.1f: %10.3e m\n", f, Dstar_over_rh(4.916e8, 2.5, KAP*f) * rh_sun))
note(abs(Dstar_over_rh(4.916e8, 2.5, 10*KAP) / Dstar_over_rh(4.916e8, 2.5, KAP) - 10^(1/2.5)) < 1e-9,
     "kappa enters as its p-th root")
cat("   So the thickness is settled to an order-unity factor, and the statement that the fold's\n")
cat("   departure from general relativity inside a hole is a macroscopic shell does not rest on\n")
cat("   anything uncomputed.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
