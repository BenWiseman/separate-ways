# ==========================================================================================
# two_sided_accretion.R showed a two-sided hole doubles its Salpeter e-folds. The part not
# extracted there is that the doubling is OBSERVABLE, and with no free parameter.
#
# The Eddington limit balances radiation pressure against GRAVITY. Gravity is set by the
# TOTAL mass, which is shared: one geometry, one mass parameter. The radiation pushing back
# on gas in our exterior comes only from what WE accrete. So each exterior may accrete at the
# full Eddington rate for the total mass, and from our side the hole grows twice as fast as
# our own luminosity can power.
#
# That is a signature, not just a rate.
# ==========================================================================================
cat("  PER OBJECT.\n")
cat("   observed luminosity L powers Mdot_our = L/(eta c^2)\n")
cat("   true growth is       Mdot_true = 2 Mdot_our\n")
cat("   so an observer inferring efficiency from growth gets eta_inferred = eta/2,\n")
cat("   and an observer inferring an Eddington ratio from growth gets exactly 2.\n\n")
for (eta in c(0.057, 0.1, 0.2, 0.32)) {
  cat(sprintf("   true eta = %.3f (%-18s) -> inferred eta = %.3f\n", eta,
      switch(as.character(eta), "0.057"="Schwarzschild", "0.1"="standard", "0.2"="a=0.9",
             "0.32"="maximal Kerr"), eta/2))
}
cat("\n  The literature invokes super-Eddington accretion for these objects with a free\n")
cat("  factor, typically a few to tens. This predicts EXACTLY 2 and has nothing to tune.\n")

cat("\n  POPULATION, which is where it can be constrained now.\n")
cat("  A Soltan-type argument infers the accreted black hole mass density from integrated\n")
cat("  quasar light. Two-sided holes contribute only half their mass to that integral, so\n")
cat("  with f the true fraction of black hole mass density in two-sided holes,\n\n")
cat("      rho_Soltan / rho_true = 1 - f/2,   i.e.  f = 2 (1 - rho_Soltan/rho_true).\n\n")
cat("      agreement of Soltan with the true density      implied bound on f\n")
for (agree in c(0.99, 0.95, 0.9, 0.8, 0.7)) {
  cat(sprintf("   %36.0f%% %26.2f\n", 100*agree, 2*(1-agree)))
}
cat("\n  So a Soltan argument agreeing to ten per cent bounds the two-sided share at 0.2.\n")
cat("  That is a weak bound and it is the honest one: this does not yet cost the framework\n")
cat("  anything, and it is not yet evidence for it either.\n")

cat("\n  WHAT WOULD MAKE IT A TEST. The factor is exactly 2 per object, so the discriminating\n")
cat("  measurement is per object and not on the population: a source whose mass growth is\n")
cat("  twice what its own luminosity can power, with no other super-Eddington signature. The\n")
cat("  little red dots are the natural place to look because they are reported overmassive\n")
cat("  for their epoch, but reporting a mass and reporting a growth rate are different\n")
cat("  measurements and only the second tests this.\n")
cat("\n  UNCHANGED CAVEATS from two_sided_accretion.R, all three still open: that mass\n")
cat("  accreted through one exterior raises the mass seen from the other, that the two-sided\n")
cat("  class is populated at all, and that radiation genuinely cannot couple across. The\n")
cat("  third is the seam coefficient again.\n")

cat("\n  ==========================================================================\n")
cat("  COROLLARY: two-sided accretion is necessarily CHAOTIC, so the spin is lower.\n")
cat("  ==========================================================================\n")
cat("  The two exteriors are causally disjoint, so their accretion flows cannot share an\n")
cat("  angular momentum direction. Torques from the two sides therefore add in quadrature\n")
cat("  rather than linearly, which is the chaotic-accretion regime by construction rather\n")
cat("  than by assumption about the gas supply.\n\n")
# equal-time comparison: one-sided accretes dM coherently, two-sided accretes dM per side
cat("      regime                       |J| for the same elapsed time     mass gained\n")
cat("   one-sided, coherent                        dM * l                      dM\n")
cat("   two-sided, two uncorrelated streams   sqrt(2) * dM * l                2 dM\n\n")
cat("  Spin parameter a ~ J/M^2, so in the growth-dominated limit M >> M_seed:\n")
cat(sprintf("      a(two-sided)/a(one-sided) -> sqrt(2)/4 = %.4f\n", sqrt(2)/4))
cat("\n  Lower spin means lower radiative efficiency, and that COMPOUNDS the factor of two\n")
cat("  rather than offsetting it: an observer assuming a standard eta while the object runs\n")
cat("  at a Schwarzschild-like value underestimates the accretion rate again.\n\n")
cat("      eta assumed   eta true (low spin)   apparent Eddington ratio = 2 * assumed/true\n")
for (ea in c(0.1, 0.2)) for (et in c(0.057, 0.08)) {
  cat(sprintf("   %11.3f %21.3f %35.2f\n", ea, et, 2*ea/et))
}
cat("\n  So the prediction is 'at least two, and plausibly three to four once the spin\n")
cat("  corollary is folded in'. The factor of two is the parameter-free part and is the\n")
cat("  one to quote; the rest depends on what efficiency an observer assumes.\n")
cat("\n  This is NOT a derivation of the chaotic-accretion spin-down literature, which\n")
cat("  already exists for one-sided holes with randomly oriented gas supply. What is\n")
cat("  specific here is that a two-sided hole cannot avoid it: the uncorrelated directions\n")
cat("  are forced by causal disjointness rather than assumed about the environment.\n")

cat("\n  ==========================================================================\n")
cat("  SCOPE CORRECTION, 2026-09-21, third pre-submission pass.\n")
cat("  ==========================================================================\n")
cat("  Two overclaims above.\n")
cat("   1. The factor is NOT 'exactly two with nothing to tune'. The enhancement is\n")
cat("      1 + Mdot_other/Mdot_ours, and it equals two only if the two exteriors accrete at\n")
cat("      EQUAL rates. Nothing here establishes that; it is an extra assumption.\n")
cat("   2. Causal disjointness does NOT force statistically uncorrelated angular momenta.\n")
cat("      Two regions can be causally disjoint and still share correlations from a common\n")
cat("      past, which a two-sided geometry has. The chaotic-accretion corollary therefore\n")
cat("      needs an independent argument and does not follow from disjointness alone.\n")
cat("  And the load-bearing assumption stated at the top remains unestablished: that mass\n")
cat("  accreted through one exterior raises the mass measured in the other.\n")
