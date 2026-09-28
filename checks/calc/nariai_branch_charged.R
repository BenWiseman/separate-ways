#!/usr/bin/env Rscript
# nariai_branch_charged.R -- fold_map_classification.R finds exactly one fold-invariant
# hole in the uncharged Schwarzschild-de Sitter family, the Nariai one, and
# nariai_stability.R shows the fold removes its instability. Charge turns that point into a
# curve. This maps the curve and checks that both results survive along all of it.
#
# Reissner-Nordstrom-de Sitter: f = 1 - 2M/r + Q^2/r^2 - Lam r^2/3. A degenerate horizon
# needs f(r0) = f'(r0) = 0; eliminating M gives
#     M = Q^2/r0 + Lam r0^3/3        and        Q^2 = r0^2 (1 - Lam r0^2).

lam <- 1                              # sets the unit; r0 and Q scale as 1/sqrt(lam)
f   <- function(r,M,Q) 1 - 2*M/r + Q^2/r^2 - lam*r^2/3
fp  <- function(r,M,Q) 2*M/r^2 - 2*Q^2/r^3 - 2*lam*r/3
fpp <- function(r,M,Q) -4*M/r^3 + 6*Q^2/r^4 - 2*lam/3
MQ  <- function(r0) { Q2 <- r0^2*(1 - lam*r0^2); c(Q2/r0 + lam*r0^3/3, sqrt(max(Q2,0))) }

cat("=== 1. the degenerate curve, and which half of it is Nariai ===\n")
cat("   Q^2 = r0^2(1 - Lam r0^2) peaks at r0 = 1/sqrt(2 Lam), the ultracold point where all\n")
cat("   three horizons meet. Above it the merged pair is r_+ = r_c, which is Nariai. Below\n")
cat("   it the merged pair is r_- = r_+, the extremal cold branch.\n\n")
cat("      r0        Q         M     f(r0)    f'(r0)    f''(r0)   branch      l2\n")
for (r0 in c(1.00, 0.98, 0.93, 0.85, 0.75, 0.70711, 0.60, 0.45)) {
  p <- MQ(r0); M <- p[1]; Q <- p[2]
  b <- if (r0 > 1/sqrt(2*lam)+1e-9) "Nariai" else if (r0 < 1/sqrt(2*lam)-1e-9) "cold" else "ULTRACOLD"
  stopifnot(abs(f(r0,M,Q)) < 1e-12, abs(fp(r0,M,Q)) < 1e-12)
  l2 <- if (fpp(r0,M,Q) < -1e-9) sqrt(2/abs(fpp(r0,M,Q))) else NA
  cat(sprintf("  %7.5f %8.5f %8.5f %9.1e %9.1e  %+9.5f  %-10s %8.5f\n",
              r0, Q, M, f(r0,M,Q), fp(r0,M,Q), fpp(r0,M,Q), b, l2)) }
qm <- optimize(function(r) -r^2*(1-lam*r^2), c(0.05, 0.999))
cat(sprintf("\n   numerical peak of Q: %.8f at r0 = %.8f; closed form 1/(2 sqrt(Lam)) = %.8f\n",
            sqrt(-qm$objective), qm$minimum, 1/(2*sqrt(lam))))
cat(sprintf("   at r0 = 1/sqrt(2 Lam) = %.8f\n", 1/sqrt(2*lam)))

cat("\n=== 2. why only the Nariai half admits the fold ===\n")
cat("   f''(r0) < 0 means f <= 0 on both sides of the double root, so r is timelike there\n")
cat("   and the near-horizon factor is dS_2. Section 5 of fold_map_classification.R used\n")
cat("   nothing but that factor, so its involution applies verbatim along the whole branch.\n")
cat("   f''(r0) > 0 is the cold branch: f >= 0 on both sides, the near-horizon factor is\n")
cat("   AdS_2, the horizon is extremal with zero surface gravity, and there is no bifurcate\n")
cat("   Killing horizon for a wedge reflection to act on. The ultracold point sits at\n")
cat("   f''(r0) = 0 and belongs to neither.\n")

cat("\n=== 3. the perturbation equation, with charge ===\n")
cat("   Reducing on the sphere with a radial field, the angular equation picks up the\n")
cat("   Maxwell term and reads 1 - (grad R)^2 - R grad^2 R = Lam R^2 + Q^2/R^2. The\n")
cat("   background fixes Lam r0^2 + Q^2/r0^2 = 1, which is the same relation as above.\n")
cat("   Linearising R = r0(1 + phi),\n")
cat("       grad^2 phi = -(2/r0^2)(2 Lam r0^2 - 1) phi.\n")
cat("   The traceless part of the 2D equation still forces grad_a grad_b phi to be pure\n")
cat("   trace, and on dS_2 of radius l2 the only solutions of grad_a grad_b phi =\n")
cat("   -(1/l2^2) ghat_ab phi are phi = c.X. Consistency therefore DEMANDS\n")
cat("       l2^2 = r0^2/(2 Lam r0^2 - 1),\n")
cat("   which is a prediction, not an input, because l2 was already fixed by f''. Checking\n")
cat("   the two against each other along the branch:\n\n")
cat("      r0         l2 from f''      l2 from the perturbation      difference\n")
worst <- 0
for (r0 in c(1.00, 0.99, 0.95, 0.90, 0.85, 0.80, 0.75, 0.72)) {
  p <- MQ(r0); M <- p[1]; Q <- p[2]
  a <- sqrt(2/abs(fpp(r0,M,Q))); b2 <- sqrt(r0^2/(2*lam*r0^2 - 1))
  worst <- max(worst, abs(a-b2)/a)
  cat(sprintf("  %7.5f  %14.8f  %27.8f  %14.1e\n", r0, a, b2, a-b2)) }
cat(sprintf("\n   worst relative difference over the branch: %.1e. They are the same number.\n", worst))
stopifnot(worst < 1e-10)
cat("   Note where the second expression fails: 2 Lam r0^2 - 1 < 0 below r0 = 1/sqrt(2 Lam),\n")
cat("   which is exactly the cold branch, where l2^2 goes negative and the factor is AdS_2.\n")
cat("   The perturbation analysis and the metric agree on where Nariai stops.\n")

cat("\n=== 4. the stability argument along the branch ===\n")
cat("   The three modes are phi = c.X for constant c in R^{1,2}, unchanged by charge, since\n")
cat("   what selects them is the traceless constraint and not the coefficient. Their parity\n")
cat("   under the fold map Nl(X0,X1,X2) = (-X0,-X1,X2) is a property of the embedding\n")
cat("   functions and not of l2 or r0, so the two growing modes stay odd and the bounded\n")
cat("   one stays even at every charge. The fold removes the instability along the whole\n")
cat("   branch, not only at its uncharged end (`checks/calc/nariai_stability.R`).\n")

cat("
=== flatly ===

  Charge turns the single fold-invariant configuration into a one-parameter family. It
  runs from uncharged Nariai at r0 = 1/sqrt(Lam) and Q = 0 to the ultracold point at
  r0 = 1/sqrt(2 Lam) and Q = 1/(2 sqrt(Lam)), which is the largest charge the family
  admits. Every member has its black-hole horizon coincident with its cosmological one,
  every member carries the same fold involution, and every member has its instability
  removed by the same parity argument.

  The extremal cold branch does not qualify and the reason is structural rather than a
  near miss: zero surface gravity means no bifurcate Killing horizon, and A.10 and A.15
  both need one.

  Still a family, still not a survey: the count is inside Reissner-Nordstrom-de Sitter and
  says nothing about rotation or about anything outside that family.\n")
