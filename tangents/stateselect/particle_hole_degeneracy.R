# WHY the entropy argument failed, which is more interesting than the failure.
#
# A.18 leaves a one-parameter band whose two ends unitarity fixes at n_min and
# n_max = 1 - n_min. On a fermionic mode, n -> 1-n IS particle-hole conjugation: it
# swaps occupied and empty. So the band's two ends are charge conjugates of one another.
#
# A single fermionic mode has reduced spectrum {n, 1-n}. Particle-hole conjugation
# PERMUTES that spectrum and leaves it setwise fixed. Therefore EVERY entanglement
# measure that is a function of the spectrum alone is blind to the difference, not just
# the Shannon entropy the result above names. Check the whole Renyi family.
Sr <- function(n, a) {                     # Renyi-a of the spectrum {n, 1-n}
  if (n <= 0 || n >= 1) return(0)
  if (abs(a-1) < 1e-12) return(-n*log(n)-(1-n)*log(1-n))
  log(n^a + (1-n)^a)/(1-a)
}
alphas <- c(0.5, 1, 2, 3, 10)
cat("  P      n_min     n_max   ")
for (a in alphas) cat(sprintf("  |dS_%-4g|", a))
cat("\n")
worst <- 0
for (P in c(0.05, 0.2, 0.5, 0.8, 0.95, 0.999)) {
  nmin <- (1-sqrt(1-P))/2; nmax <- 1-nmin
  cat(sprintf(" %5.3f %9.6f %9.6f", P, nmin, nmax))
  for (a in alphas) {
    d <- abs(Sr(nmin,a) - Sr(nmax,a)); worst <- max(worst, d)
    cat(sprintf("  %9.2e", d))
  }
  cat("\n")
}
cat(sprintf("\n  Largest difference over the whole table and every Renyi index: %.3e\n", worst))
cat("  Every spectral entanglement measure is EXACTLY blind to the two ends of the band.\n")

cat("\n  What is NOT blind: anything that sees the occupation itself rather than its spectrum.\n")
cat("  The energy is linear in n, so the gap is n_max - n_min = sqrt(1-P):\n\n")
cat("      x       P = e^{-x^2}     n_min        gap = sqrt(1-P)\n")
for (x in c(0.25, 0.5, 1, 2, 3)) {
  P <- exp(-x^2); nmin <- (1-sqrt(1-P))/2
  cat(sprintf("   %5.2f   %12.3e   %10.6f   %14.6f\n", x, P, nmin, sqrt(1-P)))
}
cat("\n  So the degeneracy is exact in every entropy at every scale, and the splitting that\n")
cat("  breaks it is dynamical and grows with x. The selection HAS to be the energy. That is\n")
cat("  not a gap in the argument, it is a consequence of the fold's own charge reversal:\n")
cat("  the band's two ends are the particle and the hole, and the fold maps one to the other.\n")
