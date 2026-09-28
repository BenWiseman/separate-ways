#!/usr/bin/env Rscript
# bform_accuracy.R -- how close to the bifurcation surface the B-form of A.10's
# ratio is good to a per cent. A.10 quotes "roughly 0.03 of a Schwarzschild
# radius at moderate angles"; an audit doubted it. This settles which angle 0.03
# belongs to.
#
# A.10 measures the first-order coefficient of the residual boost dependence in
# N^2 at rH = 0.9999:  12.1 at 60 and 120 degrees, 61.8 at 150, 517 at 170, and
# exactly zero at 90. Near a Schwarzschild horizon the lapse and the proper
# distance d above the horizon are related by
#       r - 2M = d^2/(8M)   =>   N^2 = 1 - 2M/r = (r-2M)/r ~ d^2/(16 M^2),
# which is the relation A.10 uses. One per cent accuracy is c N^2 = 0.01, so
#       d = 4M sqrt(0.01/c) = 0.4 M / sqrt(c),
# and in Schwarzschild radii r_s = 2M that is d/r_s = 0.2/sqrt(c).

M <- 1; rs <- 2*M
dist_for <- function(c, tol = 0.01) 4*M*sqrt(tol/c)      # proper distance, in M

cat("=== 1. the relation N^2 = d^2/16M^2, checked against the exact integral ===\n")
# exact proper distance above the horizon: d = int_{2M}^{r} dr'/sqrt(1-2M/r')
# r = 2M + u^2 removes the inverse-square-root endpoint exactly:
#   dr/sqrt(1-2M/r) = 2u du sqrt(r)/u = 2 sqrt(r) du
dexact <- function(r) integrate(function(u) 2*sqrt(2*M + u*u), 0, sqrt(r - 2*M),
                                rel.tol=1e-12)$value
cat("      r/2M      N^2 exact     d exact     d^2/16M^2     ratio\n")
for (x in c(1.0001, 1.001, 1.01, 1.05)) {
  r <- x*2*M; N2 <- 1 - 2*M/r; d <- dexact(r)
  cat(sprintf("   %8.4f   %.3e   %.6f   %.3e   %.5f\n", x, N2, d, d^2/(16*M^2), (d^2/(16*M^2))/N2))
}
cat("  The ratio tends to 1 as the horizon is approached, so the relation is the\n")
cat("  right leading behaviour and A.10 uses it correctly.\n")

cat("\n=== 2. how close you must be, per angle ===\n")
cat("    angle    coefficient   d for 1 per cent (M)   in Schwarzschild radii\n")
tab <- list(c(60,12.1), c(120,12.1), c(150,61.8), c(170,517))
for (t in tab) {
  d <- dist_for(t[2])
  cat(sprintf("   %5.0f deg   %9.1f      %12.4f        %14.4f\n", t[1], t[2], d, d/rs))
}
cat(sprintf("\n  So a per cent at MODERATE angles needs d/r_s = %.4f, and 0.03 is what the\n  150-degree coefficient gives (%.4f). A.10's sentence attaches the\n", dist_for(12.1)/rs, dist_for(61.8)/rs))
cat("  near-antipodal figure to the moderate-angle case, which is the wrong way round:\n")
cat("  moderate angles are MORE forgiving, not less.\n")
cat(sprintf("  moderate: %.4f r_s   150 deg: %.4f r_s   170 deg: %.4f r_s\n",
            dist_for(12.1)/rs, dist_for(61.8)/rs, dist_for(517)/rs))

cat("\n=== 3. PLANTED FAILURES ===\n")
b1 <- dist_for(12.1, tol = 0.10)
cat(sprintf("  (a) asking for ten per cent instead of one gives %.4f r_s, not %.4f\n",
            b1/rs, dist_for(12.1)/rs))
stopifnot(abs(b1/rs - dist_for(12.1)/rs) > 0.03)
b2 <- 2*M*sqrt(0.01/12.1)                       # N^2 = d^2/4M^2, wrong factor
cat(sprintf("  (b) using N^2 = d^2/4M^2 instead of d^2/16M^2 gives %.4f r_s\n", b2/rs))
stopifnot(abs(b2/rs - dist_for(12.1)/rs) > 0.01)
cat("  Both move the answer.\n")

cat(sprintf("
=== flatly ===

  One per cent needs d = 0.4 M / sqrt(c). At the moderate-angle coefficient 12.1
  that is %.4f M, or %.3f Schwarzschild radii. At 150 degrees, coefficient 61.8, it
  is %.4f M or %.3f radii; at 170 degrees, %.4f M or %.4f radii.

  A.10 says the B-form is good to a per cent 'within roughly 0.03 of a Schwarzschild
  radius at moderate angles, and only much closer than that near the antipodal one'.
  The 0.03 is the 150-degree number. At moderate angles the tolerance is twice that,
  %.3f radii, and the sentence has the two attached the wrong way round.\n",
  dist_for(12.1), dist_for(12.1)/rs, dist_for(61.8), dist_for(61.8)/rs,
  dist_for(517), dist_for(517)/rs, dist_for(12.1)/rs))
