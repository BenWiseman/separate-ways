# 5.1 gives a hard endpoint: E_nu <= 245.8 +/- 1.0 PeV, and says an event confidently
# assigned to the decay component ABOVE the endpoint refutes the fixed-history
# implementation. Both referees noted no likelihood is computed. That is fair, and the
# missing piece is not a likelihood for KM3NeT's event - which needs their exposure and
# background - but something we CAN compute and which a referee actually wants: what it
# would take for the falsifier to fire.
#
# Setup, deliberately minimal and stated: n events assigned to the decay component, each
# with fractional energy resolution sigma_rel (neutrino telescopes reconstruct shower
# energy to tens of per cent). If the true endpoint is E0 and the model's is E_m, an event
# is "above the endpoint" when its RECONSTRUCTED energy exceeds E_m by more than the
# combined width.

E_m <- 245.8; sig_m <- 1.0

cat("=== 1. the endpoint's own width is negligible beside detector resolution\n\n")
cat("        sigma_rel      detector width at E_m (PeV)     model width (PeV)    ratio\n")
for (sr in c(0.10, 0.20, 0.30, 0.50)) 
  cat(sprintf("   %12.2f %28.1f %21.1f %8.1f\n", sr, sr*E_m, sig_m, sr*E_m/sig_m))
cat("\n  At 30 per cent resolution the detector width is 73 PeV against the model's 1.0.\n")
cat("  So the +/- 1.0 PeV that 5.1 propagates, while worth having for provenance, is\n")
cat("  irrelevant to whether the test fires. The test is resolution-limited, not\n")
cat("  theory-limited, and the paper should say so.\n")

cat("\n=== 2. false-alarm rate: how often does a TRUE endpoint look violated?\n\n")
cat("   Events drawn from a spectrum truncated at E_m, smeared by the resolution. An\n")
cat("   event is counted as a violation if reconstructed above E_m by k sigma.\n\n")
set.seed(11)
draw <- function(n, sr, gamma=2) {           # E^-gamma truncated at E_m
  u <- runif(n); Emin <- 50
  Etrue <- (Emin^(1-gamma) + u*(E_m^(1-gamma) - Emin^(1-gamma)))^(1/(1-gamma))
  Etrue * exp(rnorm(n, 0, sr)) }             # lognormal smearing
cat("        sigma_rel     k=1      k=2      k=3     fraction of a TRUE sample flagged\n")
for (sr in c(0.10, 0.20, 0.30)) {
  E <- draw(2e5, sr)
  f <- sapply(c(1,2,3), function(k) mean(E > E_m*exp(k*sr)))
  cat(sprintf("   %12.2f %8.4f %8.4f %8.4f\n", sr, f[1], f[2], f[3]))
}
cat("\n  At k = 3 the false-alarm rate is well under a per cent even at 30 per cent\n")
cat("  resolution, so a single clean over-endpoint event IS informative, provided the\n")
cat("  assignment to the decay component is secure. The assignment, not the energy, is\n")
cat("  the weak link.\n")

cat("\n=== 3. power: if the TRUE endpoint is higher, how many events to notice?\n\n")
cat("   True endpoint E0 = f * E_m. Probability that at least one of n events lands\n")
cat("   above E_m by 3 sigma:\n\n")
cat("        E0/E_m      p(one event flags)     n for 90% chance     n for 99%\n")
for (f in c(1.2, 1.5, 2.0, 3.0)) {
  sr <- 0.30; Emin <- 50; gamma <- 2; E0 <- f*E_m
  u <- runif(2e5)
  Et <- (Emin^(1-gamma) + u*(E0^(1-gamma) - Emin^(1-gamma)))^(1/(1-gamma))
  p <- mean(Et*exp(rnorm(2e5,0,sr)) > E_m*exp(3*sr))
  n90 <- if (p>0) ceiling(log(0.10)/log(1-p)) else NA
  n99 <- if (p>0) ceiling(log(0.01)/log(1-p)) else NA
  cat(sprintf("   %9.1f %21.5f %20s %12s\n", f, p,
      ifelse(is.na(n90),">1e6",format(n90,big.mark=",")),
      ifelse(is.na(n99),">1e6",format(n99,big.mark=","))))
}
cat("\n  Read the numbers rather than the hope. An endpoint THREE TIMES the model's still\n")
cat("  needs about 100 securely assigned events for a 90 per cent chance of one landing\n")
cat("  above; twice the model's needs about 400; twenty per cent above needs twelve\n")
cat("  thousand. A first draft of this line called the f = 2 case 'a handful', which is\n")
cat("  off by more than two orders of magnitude.\n")
cat("\n  So the flat statement is that the endpoint test is very weak at present\n")
cat("  statistics. Neutrino telescopes have a small number of events at these energies,\n")
cat("  so the POPULATION test cannot be run at all yet. What CAN fire is the single-event\n")
cat("  test of part 2: one event clearly above the endpoint and securely assigned to the\n")
cat("  decay component refutes it, and the false-alarm rate for that is under a per cent\n")
cat("  at three sigma even with thirty per cent resolution. The falsifier is real and it\n")
cat("  is one event wide, not a spectrum shape.\n")

cat("\n=== 3b. the numbers 3.2 quotes from this file, printed rather than rounded away\n\n")
# Section 2 shows the k=3 false-alarm rate as 0.0000 because it prints four decimals. 3.2
# quotes that probability as 3.2e-5 and three sample-wide rates derived from it, and a reader
# running this file could not check any of them. Print them.
sr3  <- 0.30; Emin3 <- 50                    # the 30 per cent resolution 3.2 quotes
thr3 <- E_m*exp(3*sr3)
p3   <- integrate(function(E) (1/E^2)*pnorm(log(E/thr3)/sr3), Emin3, E_m,
                  subdivisions=4000)$value /
        integrate(function(E) 1/E^2, Emin3, E_m)$value
cat(sprintf("   at %.0f per cent resolution, 3-sigma threshold on a %.1f PeV endpoint: %.1f PeV\n",
            100*sr3, E_m, thr3))
cat(sprintf("   per-event crossing probability under the model's OWN endpoint: %.4e\n\n", p3))
cat("      events      sample-wide false alarm for a fixed per-event rule\n")
for (n in c(100, 400, 12446))
  cat(sprintf("   %9d %18.4f   = %.1f per cent\n", n, 1-(1-p3)^n, 100*(1-(1-p3)^n)))
cat("\n   That growth with n is the whole reason calibrated_falsifier.R exists: a FIXED\n")
cat("   per-event threshold is not a test at a fixed size, and at the twelve thousand\n")
cat("   events the weakest case needs it is firing on a third of conforming samples.\n")
stopifnot(abs(p3 - 3.2e-5) < 0.05e-5)

cat("\n=== 4. what this does NOT do\n\n")
cat("  It is not a likelihood for KM3NeT's event. That needs their exposure, effective\n")
cat("  area, background model and look-elsewhere treatment, none of which we have. What\n")
cat("  it supplies is the falsifier's OPERATING CHARACTERISTIC: false-alarm rate and\n")
cat("  power as functions of resolution and sample size, computed from the endpoint\n")
cat("  alone. The spectrum index and resolution are stated inputs, not measurements.\n")
