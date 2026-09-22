# Ben: would we see UV ripples from ringdowns, now that the seam is IR-transparent and
# UV-reflecting? The answer is a ratio, and it decides whether the scale dependence
# threatens the null prediction or confirms it.
#
# The seam switches on below s* = 1/(2 pi A0), i.e. above a frequency of order A0 in
# natural units. Ringdown sits at omega ~ (M omega)/M with M omega of order 0.4 for the
# fundamental mode. So the question is the ratio of the ringdown frequency to the seam
# scale.

G <- 6.67430e-11; c_ <- 2.99792458e8; hbar <- 1.054571817e-34; Msun <- 1.98847e30
tP <- sqrt(hbar*G/c_^5)                      # Planck time
fP <- 1/tP                                   # Planck frequency, Hz
f_rd <- function(M) RD_COEF*c_^3/(2*pi*G*M*Msun)   # fundamental l=m=2 ringdown, Hz

RD_COEF <- 0.374   # fundamental l=m=2 Schwarzschild QNM, in units of c^3/2 pi G M
cat(sprintf("  ringdown coefficient (l=m=2 fundamental): %.3f c^3/2 pi G M\n", RD_COEF))
cat(sprintf("  the sterile-neutrino seam scale used below: %.1f PeV\n\n", 4.916e8/1e6))
cat(sprintf("=== 1. where ringdown sits, against the Planck frequency %.3e Hz\n\n", fP))
cat("        M (Msun)      ringdown f (Hz)     f_ringdown / f_Planck\n")
for (M in c(3, 10, 65, 1e6, 1e9)) 
  cat(sprintf("   %12.3g %18.3e %26.2e\n", M, f_rd(M), f_rd(M)/fP))

cat("\n=== 2. and against seam scales other than Planck\n\n")
cat("        seam scale        f_seam (Hz)        ringdown/seam at 65 Msun\n")
scales <- list(c("Planck, 1.2e19 GeV", 1.22e19), c("GUT, 1e16 GeV", 1e16),
               c("sterile nu, 492 PeV", 4.916e8), c("TeV", 1e3), c("eV", 1e-9))
for (sc in scales) {
  E <- as.numeric(sc[2])*1e9*1.602176634e-19        # GeV -> J
  f <- E/(2*pi*hbar)
  cat(sprintf("   %-20s %14.3e %28.2e\n", sc[1], f, f_rd(65)/f))
}

cat("\n=== 3. the answer, and it CONFIRMS the null rather than threatening it\n\n")
cat("  Ringdown is deeply infrared against every candidate seam scale: 41 orders of\n")
cat("  magnitude below Planck, 37 below GUT, 30 below the 492 PeV sterile neutrino and\n")
cat("  24 below a TeV. (A first draft said 'seventeen to forty-two', which matches no row\n")
cat("  in the table above it; the eV row is 12 and is not a serious seam candidate.)\n")
cat("  A seam that is transparent in the\n")
cat("  infrared is therefore transparent at ringdown frequencies by an enormous margin,\n")
cat("  and no UV ripple appears in the waveform.\n")
cat("  So the scale dependence does not open a ringdown signature. It CLOSES one, and it\n")
cat("  reaches A.13's null prediction from a second direction: not 'the seam is\n")
cat("  transparent' but 'whatever the seam does at short distance, ringdown cannot see\n")
cat("  it'.\n")

cat("\n=== 4. what WOULD see it, since that is the useful form of the answer\n\n")
cat("  Only a probe at the seam scale itself. For a Planck-scale A0 that is inaccessible\n")
cat("  by any measurement, which is why the paper's observational content has to live in\n")
cat("  the LINEAR response - the mass-and-spin sorting of 5 - and not here.\n")
cat("  The one place the UV behaviour is not hidden is where the vacuum itself probes\n")
cat("  short distances, which is the Hadamard argument. That is the whole content: the\n")
cat("  seam's ultraviolet behaviour is constrained by regularity precisely BECAUSE no\n")
cat("  experiment can reach it.\n")
