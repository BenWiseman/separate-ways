#!/usr/bin/env Rscript
# when_can_a_hole_be_two_sided.R
#
# ############################################################################
# THE CENTRAL ARGUMENT OF THIS SCRIPT IS WITHDRAWN. Read this header before
# quoting anything below it.
#
# The argument was: a two-sided hole needs its two exteriors to be the two
# sheets; the sheets are only "in contact" until the branches decohere at
# t_dec = 1.417e-32 s (section 3.5); so the two-sided class has a mass ceiling
# equal to the horizon mass at that moment, about 5700 kg, and every such hole
# evaporates in microseconds. That would have closed the class by 7.5 orders of
# magnitude and given the companion a derived headline ceiling mirroring the
# main paper's 491.6 PeV.
#
# It rests on a false premise. The sheets' causal disjointness is not an epoch,
# it is geometry: elliptic de Sitter puts every point spacelike to its antipode
# at EVERY time, which is exactly why the cross-sheet commutator vanishes outside
# any horizon. There is no moment at which the sheets stop being in contact,
# because outside horizons they never were. The decoherence time of section 3.5
# is when the two BRANCHES of the wavefunction stop interfering, which is a
# different statement and does not bound the hole's causal structure.
#
# The arithmetic below is correct and is kept, because the horizon mass at t_dec
# and its Hawking lifetime are worth knowing. The INTERPRETATION is wrong and no
# bound on two-sidedness follows from it.
#
# CHECK 1's actual answer, which stands and came from reading the source rather
# than from this timing argument: a type II-B primordial hole's region beyond the
# throat is drawn as a closed FLRW patch (arXiv:2505.00366, Fig. 1, "a patchwork
# of closed FLRW, Schwarzschild and flat FLRW"), with no second asymptotic
# exterior described. So that channel supplies the bifurcation surface the
# algebra needs and does NOT supply the gas-fed second exterior the accretion
# argument needs.
#
# Recorded rather than deleted because the near miss is the useful part: the
# calculation was internally right, validated, robust to its inputs, and pointed
# at a headline result, and every one of those is compatible with resting on an
# identification that was never established.
# ############################################################################

banner <- function() {
  cat("\n")
  cat("  ############################################################################\n")
cat("  ## WITHDRAWN. The central argument below does not hold.\n")
  cat("  ############################################################################\n")
cat("  The argument was that the sheets are only in contact until the branches decohere, giving the\n")
cat("  two-sided class a mass ceiling near 5700 kg. The premise is false: elliptic de Sitter puts\n")
cat("  every point spacelike to its antipode at EVERY time, so there is no epoch at which the sheets\n")
cat("  stop being in contact. The arithmetic is correct and kept; the interpretation is not.\n")
cat("  \n")
cat("  Section 6 of the paper now answers this question a different way, and the opposite way: the\n")
cat("  only hole the fold can fix is the maximal Nariai one. See fold_map_classification.R.\n")
  cat("  ############################################################################\n\n")
}
banner()
G  <- 6.67430e-11; cc <- 2.99792458e8; hbar <- 1.054571817e-34
Msun <- 1.98892e30
t_dec  <- 1.417e-32          # s, section 3.5, c_G = 1
t_dec2 <- 2.249e-32          # s, c_G = 1/2, the other counting convention

# horizon mass at cosmic time t in radiation domination: M_H ~ c^3 t / G
M_H <- function(t) cc^3 * t / G
# PBH evaporation lifetime, Hawking: t_ev = 5120 pi G^2 M^3 / (hbar c^4)
t_evap <- function(M) 5120*pi*G^2*M^3/(hbar*cc^4)
# invert: mass whose lifetime equals the age of the universe
age_now <- 13.797e9*3.1557e7      # s
M_evap_now <- (age_now*hbar*cc^4/(5120*pi*G^2))^(1/3)

cat("=== 0. validation: the checks must be able to fail ===\n")
stopifnot(M_H(2*t_dec) > M_H(t_dec))                 # later horizon is heavier
stopifnot(t_evap(2e11) > t_evap(1e11))               # heavier holes live longer
planted <- tryCatch({ stopifnot(t_evap(1e11) > t_evap(2e11)); TRUE }, error=function(e) FALSE)
cat(sprintf("  monotonicity holds; reversed assertion fails as required: %s\n",
            if (!planted) "yes" else "NO - BLIND"))
# independent cross-check of M_evap_now against the textbook ~5e11 kg
cat(sprintf("  mass evaporating now: %.3e kg (textbook value ~5e11 kg, ratio %.2f)\n",
            M_evap_now, M_evap_now/5e11))
stopifnot(abs(log10(M_evap_now/5e11)) < 0.5)

cat("\n=== 1. how heavy is a hole that exists while the sheets are still in contact? ===\n")
for (nm in c("c_G = 1","c_G = 1/2")) {
  td <- if (nm=="c_G = 1") t_dec else t_dec2
  m  <- M_H(td)
  cat(sprintf("  %-10s t_dec = %.3e s   horizon mass = %.3e kg = %.2e Msun\n",
              nm, td, m, m/Msun))
  cat(sprintf("  %-10s that hole's Hawking lifetime = %.3e s\n", "", t_evap(m)))
}

cat("\n=== 2. what has to be true for such a hole to still be here? ===\n")
cat(sprintf("  mass whose lifetime equals the age of the universe : %.3e kg\n", M_evap_now))
cat(sprintf("  horizon mass at the decoherence time              : %.3e kg\n", M_H(t_dec)))
cat(sprintf("  shortfall                                         : %.2e  (%.1f orders)\n",
            M_evap_now/M_H(t_dec), log10(M_evap_now/M_H(t_dec))))
t_need <- G*M_evap_now/cc^3
cat(sprintf("\n  a hole heavy enough to survive to today must form at t > %.3e s,\n", t_need))
cat(sprintf("  which is %.1f orders of magnitude AFTER the sheets decohere.\n",
            log10(t_need/t_dec)))

cat("\n=== 3. WITHDRAWN, see header ===\n")
cat(sprintf("
  The window closes, and it closes by %.0f orders of magnitude rather than
  marginally. A hole that exists while the two sheets are still in contact
  weighs about %.0f kg, roughly a %s, and Hawking evaporation removes it in
  %.0e seconds. A hole heavy enough to survive to the present cannot have formed
  until %.0e s, which is %.0f orders of magnitude after the branches stop
  interfering, by which time the sheets are causally disjoint and there is no
  second exterior for the fold to reach.

  So CHECK 1 is answered, and not the way the growth mechanism needed. Type II-B
  formation supplies the bifurcation surface the ALGEBRA requires, and it cannot
  supply the second gas-fed exterior the ACCRETION argument requires, because
  anything old enough to have two sheets available is far too light to survive
  and anything heavy enough to survive formed long after the sheets parted.
  Section 3.4's kappa is therefore 1 for every hole that can actually be observed.

  THE YES-AND, and it is worth more than the mechanism it replaces. This is a
  DERIVED closure rather than an abundance argument. It uses only the paper's own
  decoherence time, which comes from M_1 with no free parameter, and the Hawking
  lifetime. It makes the two-sided class a statement about the first 1e-32 s
  rather than about JWST, and it explains WHY nobody has seen one without
  appealing to rarity. A referee cannot answer it with a better catalogue.
", log10(M_evap_now/M_H(t_dec)), M_H(t_dec),
   if (M_H(t_dec) > 1e3) "few tonnes" else "small mass",
   t_evap(M_H(t_dec)), t_need, log10(t_need/t_dec)))

cat("=== 4. robustness of the WITHDRAWN argument, kept for the record ===\n")
# the window opens only if the two-sided mass ceiling reaches the survival mass.
# M_ceiling = c^3 t_dec / G  and  t_dec ~ 1/M_1 (section 3.5), so M_ceiling ~ 1/M_1.
need <- M_evap_now / M_H(t_dec)
cat(sprintf("  the ceiling must rise by a factor %.2e to reach the survival mass\n", need))
cat(sprintf("  t_dec scales as 1/M_1, so M_1 would have to be smaller by the same factor:\n"))
M1 <- 491.6e6          # GeV
cat(sprintf("    M_1 = %.1f PeV would have to become %.3e GeV = %.2f GeV\n",
            M1/1e6, M1/need, M1/need))
cat(sprintf("  which is %.0f orders below the paper's own ceiling and far below the\n",
            log10(need)))
cat("  245.8 PeV neutrino line that the same machinery fixes at half of it.\n")

# and the collapse-efficiency caveat: a PBH forms with ~gamma * M_H, gamma ~ 0.2
cat(sprintf("\n  Collapse efficiency makes it worse, not better: a hole forms with about\n"))
cat(sprintf("  0.2 of the horizon mass, so the ceiling is nearer %.0f kg than %.0f kg.\n",
            0.2*M_H(t_dec), M_H(t_dec)))

# sensitivity to the evaporation convention
cat(sprintf("\n  The evaporation mass computed here, %.2e kg, sits a factor %.1f below the\n",
            M_evap_now, 5e11/M_evap_now))
cat("  textbook 5e11 kg, which is the usual spread over greybody factors and species\n")
cat("  counting. The gap being 7.5 orders, a factor of three either way changes nothing.\n")

cat(sprintf("
=== 5. WITHDRAWN, see header ===

  A black hole can have two sides only while the two sheets are still in contact,
  and that ends at %.3e s. The heaviest hole that exists by then weighs about
  %.0f kg, a few tonnes, and it evaporates in %.1e s.

  That is a derived ceiling with no fitted parameter in it, in the same sense as
  the main paper's 491.6 PeV: it follows from M_1 through the decoherence time of
  section 3.5 and from the Hawking lifetime, and nothing else.

  The companion's headline is therefore the mirror of the main paper's. Separate
  Ways says the dark matter weighs no more than 491.6 PeV. This says a black hole
  can touch the other sheet only if it weighs less than a few tonnes and lives
  less than a tenth of a millisecond. Both are ceilings, both are derived, and the
  second explains why the first universe never sees the second one.
", t_dec, M_H(t_dec), t_evap(M_H(t_dec))))


banner()
