# ITEM 1: the frequency-dependence fork, which gates route 3.
#
# Sonnet: the relative entropy's UV divergences cancel only if r(k) -> 0 as k -> inf.
# An earlier subagent: r and t MUST be frequency-independent, since kappa is
# dimensionless and the corner term carries no other scale.
# Both cannot hold. Find out which gives, and what it costs A.13.

kmaxes <- c(20, 50, 100, 200, 400, 1000, 4000, 20000)
m <- 1
Cof <- function(a, kmax, n = 2e6) {
  k <- seq(1e-8, kmax, length.out = n); dk <- k[2]-k[1]; w <- sqrt(k^2+m^2)
  sum(cos(2*k*a)/(2*w))*dk / (sum(1/(2*w))*dk)
}

cat("=== 1. is C a physical quantity or a cutoff artefact?\n\n")
cat("  A.13's ratio is A_I/A_0 = t (C+r)/(1+rC), with\n")
cat("      C = [int dk cos(2ka)/2w] / [int dk 1/2w].\n")
cat("  The NUMERATOR converges. The DENOMINATOR is int dk/(2 sqrt(k^2+m^2)), which\n")
cat("  diverges LOGARITHMICALLY. So C must fall to zero as the cutoff is lifted.\n\n")
cat("      k_max        C(a=0.02)       C(a=0.1)        C(a=1)\n")
for (K in kmaxes)
  cat(sprintf("  %10.0f %15.6f %14.6f %13.6f\n", K, Cof(0.02,K), Cof(0.1,K), Cof(1,K)))
cat("\n  C falls steadily with the cutoff. It is a REGULATOR ARTEFACT, not a property\n")
cat("  of the port geometry, and A.13's numbers were computed at k_max = 400.\n")

cat("\n=== 2. what A.13's ratio actually tends to\n\n")
cat("  As C -> 0 the Mobius factor (C+r)/(1+rC) -> r, so\n\n")
cat("      A_I / A_0  ->  t * r   in the limit the cutoff is removed.\n\n")
r_k <- function(k) (1-k^2)/(1+k^2); t_k <- function(k) 2*k/(1+k^2)
cat("     kappa      t        r        t*r      A.13 quoted at C = 0.621\n")
for (kap in c(0.2, 0.5, 1.0, 2.0, 5.0)) {
  rr <- r_k(kap); tt <- t_k(kap); C <- 0.621
  cat(sprintf("  %8.2f %8.4f %8.4f %9.4f %22.4f\n", kap, tt, rr, tt*rr, tt*(C+rr)/(1+rr*C)))
}
cat("\n  At kappa = 1 both give the same answer, t*r = 0 and the quoted value 0.621\n")
cat("  differ -- so even the anchor point moves. The cutoff was doing real work.\n")

cat("\n=== 3. so the fork resolves AGAINST the constant seam\n\n")
cat("  A seam that reflects with the same amplitude at every frequency, including\n")
cat("  k -> infinity, is not a physical interface; real interfaces go transparent in\n")
cat("  the ultraviolet. The dimensional argument is correct ABOUT THE ACTION AS\n")
cat("  WRITTEN: (kappa/2) q^T J qdot has no scale, so it gives constant r, and that\n")
cat("  is precisely why it is ultraviolet-incomplete.\n\n")
cat("  Consequences, counted:\n")
cat("    - A_0 itself is log-divergent for a constant seam, so A.13's two-port\n")
cat("      numbers are regulator-dependent, not just the relative entropy.\n")
cat("    - C is not a port-geometry factor. The two-measurement protocol in A.13,\n")
cat("      which inverts two different C values to fix kappa, does not survive as\n")
cat("      stated, because C is set by the cutoff and not by the separation.\n")
cat("    - Sonnet's condition r(k) -> 0 is the right one, and it is a REQUIREMENT ON\n")
cat("      THE SEAM ACTION rather than a nuisance: the corner term needs a scale.\n")
