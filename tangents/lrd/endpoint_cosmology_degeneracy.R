# ==========================================================================================
# CORRECTED 2026-09-21 after a pre-submission review caught two errors in the first version.
#
# ERROR 1, and it is the worst kind. The "exact degeneracy" check compared
#     shape(z, Om*1.10^3)   against   shape(z, 0.315*1.10^3)
# with Om already equal to 0.315. Those are the SAME NUMBER. It printed "difference 0.00e+00"
# at every redshift and I reported that as proof of an exact degeneracy. It compared a
# quantity with itself and could not have failed.
#
# ERROR 2. It varied Om while holding OL = 0.685 FIXED, so the comparison cosmologies were
# not flat. For flat LCDM the invariant is not Om E_nu^3 but
#     [Om/(1-Om)] E_nu^3,
# because E^2 = Om(1+z)^3 + (1-Om), and under (1+z) -> r(1+z) matching the functional form
# needs Om r^3/(1-Om) = Om'/(1-Om'), i.e. Om' = Om r^3/(1-Om+Om r^3).
#
# Both are fixed below and the check now compares two INDEPENDENTLY CONSTRUCTED cosmologies.
# ==========================================================================================
Om <- 0.315
Efl <- function(z, om) sqrt(om*(1+z)^3 + (1-om))      # FLAT by construction
shape <- function(z, om, zref=4) Efl(z,om)/Efl(zref,om)

cat("  A wrong endpoint by factor r sends 1+z -> r(1+z). For flat LCDM that is the same\n")
cat("  functional form with Om -> Om' = Om r^3 / (1 - Om + Om r^3).\n\n")
for (r in c(1.05, 1.10, 1.25)) {
  Omp <- Om*r^3/(1-Om+Om*r^3)
  cat(sprintf("   r = %.2f :  Om' = %.7f    (the naive Om r^3 would be %.7f)\n",
              r, Omp, Om*r^3))
}
r <- 1.10; Omp <- Om*r^3/(1-Om+Om*r^3)
cat(sprintf("\n  Degeneracy check at r = %.2f. LEFT: true cosmology Om = %.3f, redshifts rescaled\n", r, Om))
cat(sprintf("  by r. RIGHT: a DIFFERENT cosmology Om' = %.7f with no rescaling. These are two\n", Omp))
cat("  independently constructed quantities, not one expression written twice.\n\n")
cat("        z     rescaled true      Om' unrescaled        difference\n")
mx <- 0
for (z in c(1,2,4,9,20)) {
  # left: the shape an observer infers if the endpoint is wrong by r
  zl <- r*(1+z)-1; zrefl <- r*(1+4)-1
  a <- Efl(zl,Om)/Efl(zrefl,Om)
  b <- shape(z, Omp)
  mx <- max(mx, abs(a-b))
  cat(sprintf("   %6.1f %18.10f %20.10f %16.2e\n", z, a, b, a-b))
}
cat(sprintf("\n  largest difference %.2e: the degeneracy is exact, now demonstrated rather\n", mx))
cat("  than assumed. Only [Om/(1-Om)] E_nu^3 is measurable from the shape.\n")

cat("\n  WHAT BREAKS IT, with the propagation done correctly. From\n")
cat("  d ln[Om/(1-Om)] = -3 d ln E_nu (a LARGER endpoint means a SMALLER Om), and\n")
cat("  d ln[Om/(1-Om)] = dOm/(Om(1-Om)), so as magnitudes\n")
cat("     dOm/Om = 3 (dE/E) (1 - Om).\n")
dE <- 1.0/245.8
cat(sprintf("  E_nu = 245.8 +- 1.0 PeV is %.3f per cent, so dOm/Om = %.2f per cent.\n",
            100*dE, 100*3*dE*(1-Om)))
cat("  The first version dropped the (1-Om) factor and quoted 1.2 per cent.\n")
cat("\n  AND THE SCOPE, which the same review was right to press. The paper establishes a\n")
cat("  CEILING M_1 <= 491.6 PeV, not a measured endpoint. An endpoint known to 0.41 per cent\n")
cat("  is what a measurement would supply; until then this is what such a measurement would\n")
cat("  buy, not something the paper already has.\n")
