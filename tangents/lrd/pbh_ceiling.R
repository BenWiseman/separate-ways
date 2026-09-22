# The framework's real constraint on primordial black holes, which is tighter and more useful
# than the total-budget ceiling.
#
# The dark sector is CLOSED: every gram is the heavy sterile neutrino. Any primordial black hole
# component competes for the same total, so a fraction f in PBHs lowers the sterile mass as
# M_1 = M_1(0) (1-f)^(2/5). That is not a nine-orders-of-slack statement. It is a bound on f,
# and PBH searches measure exactly that quantity in exactly the seed mass window.
M0 <- 491.6; sig <- 2.0          # PeV, and the quoted width
f_for <- function(dM) uniroot(function(f) M0 - M0*(1-f)^(2/5) - dM, c(1e-8, 0.9))$root

cat("  what PBH fraction would move the sterile mass by a given amount?\n\n")
cat("     mass shift (PeV)   required f_PBH    as a percentage\n")
for (dM in c(0.5, 1.0, 2.0, 5.0, 20.0, 49.2)) {
  f <- f_for(dM)
  cat(sprintf("   %16.1f %16.5f %17.3f%%\n", dM, f, 100*f))
}
cat(sprintf("\n  So the framework is indifferent to PBHs below f = %.4f (%.2f per cent), the point at\n",
    f_for(sig), 100*f_for(sig)))
cat("  which the induced shift equals the quoted width, and it is in real trouble above about\n")
cat(sprintf("  f = %.3f, where the shift is ten times the width.\n", f_for(10*sig)))

cat("\n  Why that is the useful number and the total-budget ceiling is not. The ceiling\n")
cat("  n*M_seed <= rho_DM says PBHs cannot exceed ALL the dark matter, which is true of any\n")
cat("  model and constrains nothing. The one-per-cent figure is specific to this framework,\n")
cat("  because the mass is pinned by the abundance: a competing dark component does not merely\n")
cat("  coexist, it eats the budget the sterile neutrino needs and drags the endpoint down.\n")
cat("\n  The seed window little red dots point to is roughly 1e4 to 1e6 solar masses. PBH limits\n")
cat("  in that window come from CMB accretion, dynamical friction and Lyman-alpha. This file\n")
cat("  does NOT quote them: they must be taken from the current literature rather than from\n")
cat("  memory. What it supplies is the threshold they have to be compared against.\n")
cat("\n  WITHDRAWN 2026-09-21. The next two lines were printed as a prediction and are\n")
cat("  not one; see the correction at the foot of this file. [withdrawn] this framework\n")
cat("  requires\n")
cat(sprintf("      f_PBH < %.4f  over the seed mass range,\n", f_for(sig)))
cat("  and is refuted if a PBH component above that is established there.\n")
cat("\n  THE COMPARISON, fetched rather than recalled (Green and Kavanagh, arXiv:2007.10722).\n")
cat("  In that window the limits are accretion-based: CMB energy injection gives f_PBH < ~3e-9\n")
cat("  at 1e4 Msun, dwarf-galaxy heating ~1e-4, Chandra/VLA source counts ~1e-3 at stellar\n")
cat("  masses and tighter above. Against our 0.0101 that is six orders of margin on the CMB\n")
cat("  limit and two on the weakest.\n")
cat("\n  So primordial black holes cannot rescue this framework, threaten it, or be tuned to move\n")
cat("  the endpoint. And the little red dots' seeds cannot be primordial black holes at any\n")
cat("  fraction it would notice, which is what 3.4 concludes from Poisson statistics alone.\n")
cat("  Two independent routes to the same place.\n")

cat("\n  ==========================================================================\n")
cat("  SCOPE CORRECTION, 2026-09-21, from the pre-submission review.\n")
cat("  ==========================================================================\n")
cat("  This file computed the fraction that displaces the endpoint by one width and then\n")
cat("  printed it as a PREDICTION, f_PBH < 0.0101, without any condition establishing an\n")
cat("  exclusion. There is none. The paper has a CEILING, not a measured mass, and a larger\n")
cat("  f simply lowers the ceiling:\n\n")
M0 <- 491.6
cat("        f        ceiling (PeV)   still a valid ceiling?\n")
for (f in c(0.0101, 0.05, 0.1, 0.3)) {
  M <- M0*(1-f)^0.4
  cat(sprintf("   %8.4f %14.3f          %s\n", f, M, ifelse(M<=M0,"yes","NO")))
}
cat("\n  Nothing is contradicted at any f, so nothing is excluded. The fraction and the\n")
cat("  endpoint are DEGENERATE, in the same way Omega_m and the endpoint are degenerate in\n")
cat("  3.2, and a measured endpoint would break it where a bound cannot.\n")
cat("\n  WHAT IS TRUE, and it runs the other way. External accretion limits bound f, and that\n")
cat("  bounds how far a primordial component can move OUR number:\n\n")
cat("        external limit on f     maximum shift in the ceiling (PeV)\n")
for (f in c(1e-4, 1e-6, 3e-9)) cat(sprintf("   %18.0e %30.4f\n", f, M0 - M0*(1-f)^0.4))
cat("\n  So the 491.6 PeV figure is safe from contamination by a dark compact-object component\n")
cat("  to a hundredth of its own width. That is a robustness result about our prediction, not\n")
cat("  a constraint we place on primordial black holes, and the difference is the whole point.\n")
