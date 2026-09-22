# ==========================================================================================
# CORRECTION, 2026-09-21. The companion ALREADY does this audit and does it better: it
# computes (kappa t_evap)^-1 ~ (M_Pl/M)^2 as 1e-77 at 3 Msun to 1e-94 at 1e9 Msun, notes the
# violation is linear rather than discontinuous, and identifies exactly which claims are
# exposed -- one needs KMS, three need the two-sided structure. This file was written without
# checking that, and its numbers carry the full 1280 pi prefactor where the companion quotes
# the bare (M_Pl/M)^2, a factor of 4021. Same quantity, two conventions. The paper now defers
# to the companion so the two do not disagree in print. What is kept below is the explicit
# prefactor and the three-way separation of obstructions, which is a presentation of the
# companion's result and not a new one.
#
# Original note it was answering: hawking_needed.R named the calculation that would WIDEN the
# domain, saying it was not done: the Unruh state departs from
# Hartle-Hawking by O(1/(kappa t_evap)), so if the fold's horizon machinery degrades smoothly
# in that parameter it extends from eternal holes to astrophysical ones as a controlled
# approximation. Compute the parameter.
#
# kappa = 1/(4GM) and t_evap = 5120 pi G^2 M^3 / (hbar c^4), so
#     kappa t_evap = 1280 pi G M^2 / (hbar c^4) = 1280 pi (M/M_Pl)^2
G <- 6.67430e-11; hbar <- 1.054571817e-34; c <- 2.99792458e8
Msun <- 1.98892e30
M_Pl <- sqrt(hbar*c/G)                       # kg
cat(sprintf("  Planck mass = %.4e kg\n\n", M_Pl))
kt <- function(M) 1280*pi*(M/M_Pl)^2
cat("      object                       M (kg)      kappa*t_evap      1/(kappa t_evap)\n")
rows <- list(c(1e12,"primordial, evaporating now"), c(Msun,"solar mass"),
             c(24*Msun,"stellar remnant seed"), c(1e5*Msun,"LRD seed, light"),
             c(1e7*Msun,"LRD seed, heavy"), c(1e9*Msun,"quasar"))
for (r in rows) {
  M <- as.numeric(r[1])
  cat(sprintf("   %-28s %.2e %16.3e %18.3e\n", r[2], M, kt(M), 1/kt(M)))
}
cat("\n  FLATLY: for every astrophysical mass the quasi-static parameter runs from 1.2e-43 down to\n")
cat("  3.0e-98. Time dependence is not the obstruction to extending the horizon machinery. A\n")
cat("  hole of a hundred thousand solar masses is stationary to eighty-odd decimal places\n")
cat("  over its own evaporation time, which is itself longer than the age of the universe by\n")
cat("  a vast margin.\n")

cat("\n  BUT THAT IS NOT THE WHOLE OBSTRUCTION, and saying otherwise would be the error this\n")
cat("  calculation is meant to avoid. Two things were being conflated:\n\n")
cat("   (a) TIME DEPENDENCE. Settled above: negligible at every astrophysical mass.\n")
cat("   (b) BOUNDARY CONDITION. Hartle-Hawking has an incoming thermal flux and Unruh does\n")
cat("       not. That is a different state, not a small perturbation of one, and no power of\n")
cat("       1/(kappa t_evap) turns one into the other. The outgoing sector agrees near the\n")
cat("       horizon; the ingoing sector does not.\n")
cat("   (c) TOPOLOGY. A hole formed by collapse has no past exterior at all. That is the\n")
cat("       one-sided/two-sided split of 4.3 and it is a statement about the manifold, so\n")
cat("       adiabaticity has nothing to say about it either.\n")

cat("\n  So the honest result is a NARROWER open problem than before, not a solved one. The\n")
cat("  domain restriction on the horizon machinery is not 'these holes are evolving'. It is\n")
cat("  entirely (b) and (c). Whoever wants to extend the machinery must show which sector\n")
cat("  the modular construction actually uses: if it reads only the outgoing near-horizon\n")
cat("  correlations, (b) costs nothing and only (c) remains; if it needs the ingoing modes,\n")
cat("  (b) is fatal for astrophysical holes however slowly they evaporate.\n")
cat("\n  That is a yes-or-no question about the existing construction rather than a new\n")
cat("  physical input, which is what makes it the next thing to do.\n")
