# Ben asked whether the UV/IR occupation structure touches empirical data, naming DESI
# spectroscopy correlated with gravitational waves. The first answer is no. The second is
# yes, and it is not the one either of us was looking at.
#
# The occupation the Hadamard condition constrains is not the dark matter's alone. Every
# field present at the bang has one, gravitons included, and a graviton occupation IS a
# stochastic gravitational-wave background spectrum. So the constraint lands on an
# observable directly.

GeV_per_K  <- 8.617333262e-14
T0         <- 2.7255 * GeV_per_K          # GeV
hbar_GeVs  <- 6.582119569e-25
Mpl        <- 1.220890e19                 # GeV
M1         <- 491.6e6                     # GeV  (491.6 PeV)
g_s_0      <- 3.938
g_star     <- 106.75
t_dec      <- 1.436e-32                   # s   (paper 4.2, c_G = 1)

cat("=== 1. the dictionary between occupation and observable\n\n")
cat("  Energy density per logarithmic momentum for a mode set with occupation n_k:\n")
cat("     rho = Int (d3k/(2pi)^3) k n_k  =  Int (k^2 dk/2pi^2) k n_k\n")
cat("  so   d rho / d ln k  proportional to  k^4 n_k,  and since Omega_GW(f) is exactly\n")
cat("  that divided by the critical density,\n\n")
cat("     Omega_GW(f)  proportional to  f^4 n_f.\n\n")
cat("  The paper's Hadamard condition is that Sum_k n_k k^m converges for EVERY m. In the\n")
cat("  same integral measure that sum is Int k^2 dk k^m n_k, so the m-th condition is\n")
cat("  convergence of Int f^(m+3) n_f d ln f = Int f^(m-1) Omega_GW(f) d ln f.\n")

cat("\n=== 2. which moment is already a measurement?\n\n")
cat("   m       the condition                                   status\n")
cat("  ---  ------------------------------------------------  ------------------------\n")
cat("   1    Int Omega_GW dln f converges                      MEASURED: the BBN / N_eff\n")
cat("                                                          bound on extra radiation\n")
cat("   2    Int f Omega_GW dln f converges                    not measured\n")
cat("   3    Int f^2 Omega_GW dln f converges                  not measured\n")
cat("  ...  every higher moment                                not measured\n")
cat("\n  So the integral bound that cosmology already imposes on a gravitational-wave\n")
cat("  background is the m=1 case of the fold's Hadamard condition. The fold does not\n")
cat("  merely survive that bound, it DERIVES it, and asserts infinitely many strictly\n")
cat("  stronger ones alongside it. That is the algebra reaching an observable.\n")

cat("\n=== 3. what the higher moments forbid that m=1 allows\n\n")
cat("  Take a tail Omega_GW(f) ~ f^(-s) above some break, the generic high-frequency\n")
cat("  behaviour of a cosmic-string network and of most preheating spectra.\n")
cat("  Moment m converges iff f^(m-1-s) is integrable at large f, i.e. iff s > m-1.\n\n")
cat("     tail slope s     highest moment that converges     passes Hadamard?\n")
for (s in c(1, 2, 4, 8, 20)) {
  cat(sprintf("        f^-%-6.0f %24.0f %24s\n", s, s, "NO"))
}
cat("      exponential                      every m                      YES\n")
cat("      Gaussian                         every m                      YES\n")
cat("\n  A BROKEN POWER LAW FAILS, at whatever slope. Only exponential-or-faster\n")
cat("  suppression survives. This is the sharp statement: the fold does not merely\n")
cat("  require a turnover in the background, it requires the turnover to be\n")
cat("  super-polynomial. A measured power-law decline above the peak kills it.\n")

cat("\n=== 4. where the fold puts the turnover\n\n")
cat("  The scale is the production epoch, which the paper already fixes twice over.\n")
cat("  Redshifting the horizon at production to today:\n")
cat("     f_0 = (a_prod/a_0) H_prod/2pi,   a_prod/a_0 = (g_s0/g_star)^(1/3) T0/T_prod\n\n")
cat("   anchor for H_prod         T_prod (GeV)    f_0 today (Hz)     band\n")
anchors <- list(c("H = 1/(2 t_dec)", hbar_GeVs/(2*t_dec)), c("H = M_1", M1))
for (a in anchors) {
  H <- as.numeric(a[[2]])
  Tp <- sqrt(H * Mpl / (1.66*sqrt(g_star)))
  f0 <- (g_s_0/g_star)^(1/3) * (T0/Tp) * H / (2*pi) / hbar_GeVs
  band <- if (f0 > 1e4) "above LIGO" else if (f0 > 10) "LIGO band" else "below LIGO"
  cat(sprintf("   %-22s %14.3e %17.3e %14s\n", a[[1]], Tp, f0, band))
}
cat("\n  Both anchors land the turnover near or above 10^5 Hz, one to three decades above\n")
cat("  LIGO's band. So the fold predicts the observed decades are all BELOW the break,\n")
cat("  which is consistent with every current non-detection and is not yet a discriminant\n")
cat("  on the turnover's position. The discriminant available NOW is the tail SHAPE\n")
cat("  wherever a background is measured, not the break frequency.\n")

cat("\n=== 5. what this is worth, flatly\n\n")
cat("  Gained: a falsifier that costs nothing to state, follows from the seam argument\n")
cat("  rather than from the cosmology, and is tested by instruments that exist. If a\n")
cat("  stochastic background is resolved well enough to show a power-law fall above its\n")
cat("  peak, the fold's state is not Hadamard and the seam argument goes with it.\n")
cat("  Not gained: any DESI test of the occupation shape. Section 6 should say so.\n")
