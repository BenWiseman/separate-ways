# A CALIBRATED endpoint test, at a common false-alarm probability.
#
# The verification pass was right that the old numbers were not powers at a common test size: a
# fixed PER-EVENT threshold lets the sample-wide false alarm grow with N, reaching 32 per cent at
# the 12,000 events the weakest case needed. It also said calibrating this is straightforward and
# was not done. Do it, and quote the falsifier at a stated size instead of withdrawing it.
#
# Setup, unchanged from endpoint_power.R: events drawn from E^-2 above 50 PeV, lognormal response
# with sigma = 0.3 in ln E, model endpoint E_m = 245.8 PeV. Under the model no event is produced
# above E_m, so a reconstructed excursion is response only. Under an alternative with true
# endpoint r*E_m the spectrum extends further.
E_m <- 245.8; Ecut <- 50; sg <- 0.3
# reconstructed-energy tail probability for a source population truncated at Etrue_max
p_above <- function(thr, Emax) {
  # E^-2 between Ecut and Emax, then lognormal smearing
  f <- function(E) (1/E^2) * pnorm((log(E/thr))/sg)      # P(E_rec > thr | E)
  num <- integrate(f, Ecut, Emax, subdivisions=4000)$value
  den <- integrate(function(E) 1/E^2, Ecut, Emax)$value
  num/den
}
# calibrate: choose the threshold so the SAMPLE-WIDE false alarm is alpha under the model
thr_for <- function(N, alpha=0.05) {
  target <- 1-(1-alpha)^(1/N)                            # required per-event rate
  uniroot(function(t) p_above(t, E_m) - target, c(E_m, 40*E_m), tol=1e-10)$root
}
cat("  calibrated to a 5 per cent SAMPLE-WIDE false alarm, model endpoint 245.8 PeV\n\n")
cat("      N events   threshold (PeV)   power vs true endpoint  2x    3x    5x\n")
for (N in c(30, 100, 300, 1000, 3000)) {
  th <- thr_for(N)
  pw <- sapply(c(2,3,5), function(r) 1-(1-p_above(th, r*E_m))^N)
  cat(sprintf("   %10d %17.1f %22.3f %5.3f %5.3f\n", N, th, pw[1], pw[2], pw[3]))
}
cat("\n  Read the table as the falsifier's operating characteristic, which is what it now is.\n")
cat("  At a fixed 5 per cent chance of firing on a population that obeys the model, a few hundred\n")
cat("  reconstructed events at these energies distinguish the model's endpoint from one twice as\n")
cat("  high with high probability. The threshold rises with N because the test must stay honest\n")
cat("  as the sample grows, which is exactly what the uncalibrated version failed to do.\n")
