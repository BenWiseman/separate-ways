# 5.6 says plainly: "No KM3NeT likelihood preference is calculated." Fill that in, because the
# arithmetic is short and a referee will do it. Reference [79] gives median neutrino energy
# 220 PeV with a 68% interval of 110-790 PeV. The endpoint predicted here is 245.8 PeV and the
# claim is that the event "sits inside" it. Ask what fraction of the posterior actually does.
#
# The interval is strongly asymmetric IN LOG SPACE: log(110/220) = -0.693 but log(790/220) =
# +1.278. So no symmetric lognormal fits, and a split (two-piece) form is the honest minimum.

m  <- 220; lo <- 110; hi <- 790; E0 <- 245.8
s1 <- log(m/lo); s2 <- log(hi/m)
cat(sprintf("   median %.0f PeV, 68%% interval %.0f-%.0f PeV\n", m, lo, hi))
cat(sprintf("   log-space half-widths: lower %.4f, upper %.4f  (ratio %.2f, so not symmetric)\n\n",
    s1, s2, s2/s1))

z <- log(E0/m)/s2
p_split <- 1 - pnorm(z)
cat(sprintf("   split-lognormal, each side carrying half the mass:\n"))
cat(sprintf("     z = log(%.1f/%.0f)/%.4f = %.4f, so P(E > %.1f PeV) = %.3f\n\n", E0, m, s2, z, E0, p_split))

# robustness: a plain lognormal matched to the UPPER half only, and one matched on the full 68%
s_up   <- s2
s_full <- (s1+s2)/2
for (nm in list(c("lognormal matched to the upper half", s_up),
                c("lognormal matched to the mean half-width", s_full),
                c("lognormal matched to the lower half", s1))) {
  s <- as.numeric(nm[[2]])
  cat(sprintf("   %-42s sigma = %.4f  P(E > %.1f) = %.3f\n",
      nm[[1]], s, E0, 1 - pnorm(log(E0/m)/s)))
}

cat("\n   And the same question asked of the muon energy, which is what is measured directly:\n")
mu <- 120; mlo <- 60; mhi <- 230        # 120 (+110, -60) PeV
cat(sprintf("     muon 120 (+110, -60) PeV; the endpoint maps to a muon energy well below the\n"))
cat(sprintf("     neutrino energy, so this is a consistency note rather than a second test.\n"))

cat("\n=== reading, flatly\n\n")
cat(sprintf("   About %.0f per cent of the reconstructed posterior lies ABOVE the predicted\n", 100*p_split))
cat("   endpoint, and the figure is between 44 and 47 per cent across every reasonable\n")
cat("   parametrisation of the published interval. So the event does not confirm the\n")
cat("   endpoint: the posterior straddles it almost evenly.\n")
cat("   It does not refute it either. Refutation needs an event CONFIDENTLY above the\n")
cat("   endpoint, and a posterior this wide is not confident about anything.\n")
cat("   WHAT THIS CHANGES: 'the median sits inside the bound' is true and is a weaker\n")
cat("   statement than it reads. The honest version is that the test is live and nearly\n")
cat("   balanced, which is better for the paper than a coincidence would be, because it\n")
cat("   means the next well-reconstructed event decides something.\n")
