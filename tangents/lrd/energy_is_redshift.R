# ==========================================================================================
# dPhi/dE ∝ 1/(E H(z)) with z = E_nu/E - 1, and E_nu is FIXED BY THE ABUNDANCE, not fitted.
# So the energy axis of the decay continuum is a redshift axis with no free parameter, and
# H(z) = const / (E dPhi/dE). A detector measuring the continuum is measuring the expansion
# history. Work out the range, the mapping, and then the cost, which is the part that decides
# whether this is a probe or a curiosity.
# ==========================================================================================
Om <- 0.315; OL <- 0.685; H0 <- 67.4
Ez <- function(z) sqrt(Om*(1+z)^3 + OL)
Enu <- 245.8    # PeV

cat("  THE MAPPING. Energy bin -> redshift bin, fixed, nothing fitted.\n\n")
cat("      E (PeV)      z = E_nu/E - 1      H(z)/H0      what lives there\n")
rows <- list(c(200,""),c(122.9,""),c(80,""),c(49.2,"little red dots (z~4)"),
             c(35,""),c(24.6,"little red dots (z~9)"),c(15,""),c(10,"cosmic dawn"),
             c(5,""),c(2.5,"pre-recombination"))
for (r in rows) {
  E <- as.numeric(r[1]); z <- Enu/E - 1
  cat(sprintf("   %10.1f %18.2f %13.1f      %s\n", E, z, Ez(z), r[2]))
}

cat("\n  The little red dots sit at z = 4 to 9, which is 49.2 down to 24.6 PeV. That is inside\n")
cat("  the band the endpoint test already needs a detector to reach, not somewhere else.\n")

# --- what fraction of the decay flux lands in each redshift bin? -------------------------
dPhidE <- function(E) { z <- Enu/E - 1; ifelse(z > 0, 1/(E*Ez(z)), 0) }
tot <- integrate(dPhidE, 1e-4, Enu, rel.tol=1e-10, subdivisions=6000)$value
cat("\n  FRACTION OF THE EXTRAGALACTIC COMPONENT PER REDSHIFT BIN\n\n")
cat("      z range        E range (PeV)        fraction of the component\n")
bins <- list(c(0,1),c(1,2),c(2,4),c(4,9),c(9,20),c(20,50),c(50,1e3))
for (b in bins) {
  Ehi <- Enu/(1+b[1]); Elo <- Enu/(1+b[2])
  f <- integrate(dPhidE, Elo, Ehi, rel.tol=1e-10, subdivisions=6000)$value/tot
  cat(sprintf("   %5.0f - %-6.0f %8.1f - %-8.1f %20.4f\n", b[1], b[2], Elo, Ehi, f))
}

# --- the cost -----------------------------------------------------------------------------
cat("\n  THE COST, which is what decides it. H(z) in a bin is measured to dH/H ~ 1/sqrt(N_bin),\n")
cat("  so N_bin = (dH/H)^-2, and the total decay-component sample needed is that divided by\n")
cat("  the bin's flux fraction.\n\n")
f_lrd <- integrate(dPhidE, Enu/10, Enu/5, rel.tol=1e-10, subdivisions=6000)$value/tot
cat(sprintf("   the z = 4 to 9 bin carries %.4f of the extragalactic component,\n", f_lrd))
cat(sprintf("   and %.4f of the whole decay component (Galactic line included, 57.3 per cent).\n",
            f_lrd*0.427))
cat("\n      target dH/H    events in bin    total decay events needed\n")
for (prec in c(0.5, 0.2, 0.1, 0.05, 0.01)) {
  Nbin <- ceiling(prec^-2)
  cat(sprintf("   %12.0f%% %16d %28.3e\n", 100*prec, Nbin, Nbin/(f_lrd*0.427)))
}
cat("\n  One subtlety that makes it cleaner, not dirtier: the overall normalisation is the\n")
cat("  lifetime, which is unknown, so what a bin measures is not H(z) absolutely but the\n")
cat("  RATIO of H between bins. The lifetime cancels out of that ratio. The observable is\n")
cat("  therefore the SHAPE of the expansion history, self-calibrating, with the energy to\n")
cat("  redshift mapping fixed by the abundance.\n")

cat("\n  FLATLY, and it is better than expected: a crude measurement comes with the endpoint\n")
cat("  test rather than after it. Fifty per cent on H(z) in the z = 4 to 9 bin needs FOUR\n")
cat("  events in that bin and about 130 decay events in total, which is the same order as\n")
cat("  the ~1e2 the endpoint test already needs. Ten per cent needs 3.3e3, the top of the\n")
cat("  1e2 to 1e3 range the paper already quotes. One per cent needs 3.3e5 and is out of\n")
cat("  reach.\n")
cat("\n  So this is not a separate and more expensive programme. The same sample that fires\n")
cat("  the falsifier gives a first, crude reading of the expansion history over exactly the\n")
cat("  epoch the little red dots occupy.\n")
cat("\n  WHAT IS UNUSUAL, and it does not depend on the flux at all: in every other\n")
cat("  cosmological probe the distance-redshift mapping is fitted or calibrated against\n")
cat("  something. Here it is fixed by a particle mass the relic abundance already\n")
cat("  determined, so the energy axis IS a redshift axis with zero free parameters. That is\n")
cat("  a statement about what the observable IS, and it stands whether or not the events\n")
cat("  ever arrive.\n")
cat("\n  A first version of this file concluded the opposite, that even 50 per cent was out\n")
cat("  of reach, by misreading its own table: 4 events in the bin, not 1e2.\n")
