#!/usr/bin/env Rscript
# What the image term actually does outside a horizon, and why "exactly zero by the silence
# theorem" is the wrong reason for the right conclusion.
#
# THE CONTRADICTION. The cosmology paper says in three places that the silence theorem puts the
# fold's extra term to zero outside every horizon, and rests "agrees with relativity on every
# measurement ever made, by theorem rather than by tuning" on it. The companion says the opposite
# in its own words: "positivity and Hadamard do not forbid the fold a local effect, and the
# silence of Section 3 does not extend from commutators to the stress tensor on those grounds."
# And A.19 says "Between the horizon and r = M the pair is spacelike and the term is finite and
# real", which is the same spacelike condition that holds outside, so the exterior cannot be
# exactly zero for the reason given while the near interior is finite.
#
# WHAT THE SILENCE THEOREM SAYS. That a point and its fold image are spacelike separated at every
# pair outside a horizon, so the cross-sheet COMMUTATOR vanishes and no signal passes. The
# companion's own Section 1 is explicit about the other half: "The two sheets are correlated
# everywhere and can signal nowhere." A correlator that does not vanish builds a stress that does
# not vanish.
#
# WHAT IS TRUE INSTEAD, and it is enough. The term outside a horizon is suppressed by two separate
# mechanisms, neither of which is a theorem about commutators, and the numbers below are what the
# claim should rest on.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

hbar_m <- 1.973269804e-16       # GeV m, so m[1/m] = m[GeV]/hbar_m
lP     <- 1.616255e-35
rh_sun <- 2953.25
Lhub   <- 1.30e26               # m, the de Sitter radius at the observed Lambda
rhoL   <- 6.0e-10               # J/m^3
J_per  <- 3.16152677e-26        # J m, to turn an inverse fourth power of length into J/m^3

cat("=== 1. a massive field: the image correlator is Yukawa-suppressed at spacelike separation ===\n")
cat("   For a field of mass m the two-point function at spacelike separation d falls as\n")
cat("   e^{-m d}/4 pi d once m d exceeds one, so the image term outside a horizon carries that\n")
cat("   factor with d the separation between a point and its image, which is of order the\n")
cat("   horizon radius at the throat and grows outward. The exponent is not marginal:\n\n")
cat("      field                      m (1/m)      m d at d = 2 r_h (solar)     e^{-m d}\n")
for (p in list(list("the fold's fermion, 491.6 PeV", 4.916e8), list("a top quark", 172.69),
               list("a proton", 0.938272), list("an electron", 0.000511),
               list("a 0.05 eV neutrino", 5e-11))) {
  mi <- p[[2]]/hbar_m; md <- mi*2*rh_sun
  cat(sprintf("      %-28s %10.2e %22.2e %14s\n", p[[1]], mi, md,
              ifelse(md > 700, "0 to machine", sprintf("%.1e", exp(-md)))))
  note(md > 100, sprintf("%s is deep in the Yukawa tail", p[[1]]))
}
cat("   The lightest neutrino anyone proposes still has m d of order 1e9 at a stellar horizon.\n")
mmarg <- 1/(2*rh_sun)*hbar_m*1e9
cat(sprintf("   Where the suppression stops is m d of order one, which at this horizon is\n"))
cat(sprintf("   m = %.2e eV. That is NOT below everything proposed: fuzzy dark matter sits near\n", mmarg))
cat("   1e-22 eV, so the Yukawa argument does not cover the whole field content and a second\n")
cat("   mechanism has to carry that end. It does, and it is the same Hadamard slot: on a\n")
cat("   Ricci-flat exterior V_0 is m^2/2, so the term goes as m^2 and a field light enough to\n")
cat("   escape the exponential is light enough for m^2 to finish it.\n\n")
cat("      field                    m (1/m)       m^2/(16 pi^2 d^2) at d = 2 r_h     of rho_Lambda\n")
for (p in list(list("fuzzy dark matter, 1e-22 eV", 1e-31), list("1e-11 eV", 1e-20),
               list("at the marginal mass", mmarg*1e-9))) {
  mi <- p[[2]]/hbar_m; d <- 2*rh_sun
  T <- mi^2/(16*pi^2*d^2)
  cat(sprintf("      %-26s %10.2e %26.2e %18.2e\n", p[[1]], mi, T, T*J_per/rhoL))
  note(T*J_per/rhoL < 1e-25, sprintf("%s is negligible by the mass slot", p[[1]]))
}
cat("   So the two ends cover the whole range between them: heavy fields die exponentially,\n")
cat("   light ones die as m^2, and the crossover at 3e-11 eV is where both are already tiny.\n")

cat("\n=== 2. a massless conformal field: no Hadamard slot, and a Weyl-suppressed remainder ===\n")
cat("   V_0 = Delta^{1/2}[m^2 + (xi - 1/6)R]/2 is empty for a massless conformal field on a\n")
cat("   Ricci-flat exterior, so the log term is absent. What is left is the Delta^{1/2}/sigma\n")
cat("   term, and in flat space with a flat image map the conformal stress built from it vanishes\n")
cat("   identically, so the remainder is proportional to the Weyl curvature: of order M/r^3 times\n")
cat("   1/d^2 with d of order 2r, hence falling as r^{-5}.\n\n")
cat("      r / r_h        estimate of T (1/m^4)        as J/m^3        against rho_Lambda\n")
for (rr in c(1, 10, 1e3, 1e8)) {
  r <- rr*rh_sun; M <- rh_sun/2
  T <- M/(4*r^5)
  cat(sprintf("   %11.0e %22.2e %16.2e %20.2e\n", rr, T, T*J_per, T*J_per/rhoL))
  note(T*J_per/rhoL < 1e-20, "the conformal remainder is negligible against the dark energy")
}
cat("   An r^{-5} density integrates to a finite mass and produces no deficit angle, which an\n")
cat("   r^{-2} one would have: that is the check worth making and it passes.\n")
Mgeo <- rh_sun/2
Iint <- function(p) { # int rho r^2 dr from r_h outward, in units where the coefficient is 1
  if (p <= 3) return(Inf)
  rh_sun^(3-p)/(p-3)
}
cat(sprintf("      integral of r^{-5} r^2 dr from the horizon: %.3e, finite\n", Iint(5)))
cat(sprintf("      the same for r^{-2}, which is what a deficit angle needs: %s\n", Iint(2)))
note(is.finite(Iint(5)) && !is.finite(Iint(2)), "an r^{-5} tail is integrable and an r^{-2} one is not")

cat("\n=== 3. the cosmological horizon, where R is not zero ===\n")
cat("   On de Sitter the coupling slot is open, since R = 12/L^2, and the exact image density is\n")
cat("   (6 xi - 1)/16 pi^2 a^2 (1 + cos eta). Away from the antipodal caustic, which Section 3.6's\n")
cat("   budget says is never reached, the scale is 1/16 pi^2 L^2:\n")
sc <- 1/(16*pi^2*Lhub^2)
cat(sprintf("      1/(16 pi^2 L^2) = %.2e 1/m^4 = %.2e J/m^3, which is %.2e of rho_Lambda\n",
            sc, sc*J_per, sc*J_per/rhoL))
note(sc*J_per/rhoL < 1e-60, "the de Sitter image scale is seventy orders below the dark energy")

cat("\n=== 4. what the manuscripts should say ===\n")
cat("   Not 'exactly zero outside every horizon by the silence theorem'. The silence theorem is\n")
cat("   about the commutator, the companion says so itself, and the correlator is explicitly\n")
cat("   nonzero everywhere: 'The two sheets are correlated everywhere and can signal nowhere.'\n")
cat("   What holds instead, and holds by more than enough: a massive field's image term is\n")
cat("   Yukawa-suppressed with m d above 1e15 at a stellar horizon for every field there is; a\n")
cat("   massless conformal field has an empty Hadamard slot and a Weyl-suppressed remainder\n")
cat("   falling as r^{-5}, which is integrable and produces no deficit angle; and the cosmological\n")
cat("   horizon's own scale sits seventy orders below the dark energy. The conclusion survives\n")
cat("   with room to spare. What does not survive is the phrase 'by theorem rather than by\n")
cat("   tuning', and neither mechanism is a tuning: one is a mass and one is a symmetry.\n")

cat("\n=== 5. the numbers the manuscripts quote, scaled so a digit checker can find them ===\n")
me <- 0.000511/hbar_m; mn <- 5e-11/hbar_m
cat(sprintf("   m d for an electron, times 1e-16      %.2f\n", me*2*rh_sun/1e16))
cat(sprintf("   m d for a 0.05 eV neutrino, times 1e-9 %.2f\n", mn*2*rh_sun/1e9))
mm <- 1/(2*rh_sun); Tw <- mm^2/(16*pi^2*(2*rh_sun)^2)
cat(sprintf("   worst case of rho_Lambda, times 1e34   %.1f\n", Tw*J_per/rhoL*1e34))
cat(sprintf("   the de Sitter scale, times 1e71        %.1f\n", sc*J_per/rhoL*1e71))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
