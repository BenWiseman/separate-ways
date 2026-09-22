# A hostile referee pass found the entropy density wrong. Verified independently, twice, and it
# is: s_0 = 2891.2 cm^-3 is 2.2215e-38 GeV^3, not the 2.3e-38 the reconstruction uses. Price it.
kB <- 8.617333262e-14; T <- 2.7255*kB; gs <- 3.9091
s0 <- (2*pi^2/45)*gs*T^3
hbarc <- 1.9732698e-14                     # GeV cm
cat(sprintf("   from g_*s T^3          : %.6e GeV^3\n", s0))
cat(sprintf("   from 2891.2 cm^-3      : %.6e GeV^3\n", 2891.2*hbarc^3))
cat(sprintf("   paper uses             : 2.300000e-38 GeV^3   (high by %.2f per cent)\n\n",
    100*(2.3e-38/s0-1)))
f <- (2.3e-38/s0)^(2/5)                    # M1 ~ s0^(-2/5)
M0 <- 484.8; E0 <- 242.4
cat(sprintf("   M_1 ~ s_0^(-2/5), so the correction multiplies the mass by %.6f\n", f))
cat(sprintf("   endpoint  %.1f -> %.1f PeV     line  %.1f -> %.1f PeV\n", M0, M0*f, E0, E0*f))
cat(sprintf("   shift %.1f PeV against the propagated +-1.9 PeV: %.1f times the quoted width\n\n",
    M0*f-M0, (M0*f-M0)/1.9))

cat("=== does the KM3NeT comparison survive the higher endpoint?\n\n")
m<-220; lo<-110; hi<-790; s2<-log(hi/m)
for (E in c(E0, E0*f)) {
  p <- 1-pnorm(log(E/m)/s2)
  cat(sprintf("   endpoint %.1f PeV : P(E_nu > endpoint) = %.3f\n", E, p))
}
cat("\n   The endpoint RISES, so the measured event sits further inside it and the posterior\n")
cat("   mass above the bound FALLS. The correction moves the comparison in the model's\n")
cat("   favour, which is exactly why it has to be made rather than left.\n")

cat("\n=== what else moves\n\n")
cat(sprintf("   propagated width was 0.40%% of the mass, so +-1.9 -> +-%.1f PeV\n", 0.004*M0*f))
cat(sprintf("   half-mass width      +-1.0 -> +-%.1f PeV\n", 0.004*M0*f/2))
cat("   Sigma m_nu = 58.8 meV is untouched: it comes from the oscillation inputs and the\n")
cat("   stabilising rule, not from the abundance normalisation.\n")
cat("   t_dec ~ M_1^(2/3) shifts by the same factor to the 2/3 power, %.4f.\n")
