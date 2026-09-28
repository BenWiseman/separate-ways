#!/usr/bin/env Rscript
# how_many_sheets.R -- why two sheets and not four, and what four would have cost.
#
# A.15 fixes WHICH involution the transverse map is: every free involution of S^2 is conjugate to
# the antipodal map, because the quotient is a closed surface of Euler characteristic one and
# RP^2 is the only one. It does not ask why an INVOLUTION at all. Ben did: what if there are four
# sheets?
#
# The same Euler characteristic settles it in one line. For a free action of a finite group G on
# a compact manifold M, the quotient is a manifold and chi(M) = |G| chi(M/G), so |G| must divide
# chi(M). The transverse sphere in D spacetime dimensions is S^(D-2), and
#     chi(S^n) = 1 + (-1)^n ,
# which is 2 for even n and 0 for odd. So in four dimensions the transverse sphere is S^2, chi is
# 2, and |G| divides 2: the fold is two-sheeted and nothing else is available. Four sheets is not
# excluded by a dynamical argument or by taste. It does not fit on the sphere.

chi <- function(n) 1 + (-1)^n
cat("=== 1. how many sheets the transverse sphere allows ===\n\n")
cat("    D    transverse sphere    chi    |G| must divide    sheets available\n")
for (D in 4:9) {
  n <- D - 2; c <- chi(n)
  cat(sprintf("   %2d    S^%-16d %-6d %-18s %s\n", D, n, c,
              if (c == 0) "no constraint" else as.character(c),
              if (c == 0) "any number: S^odd carries free Z_n actions, the lens spaces"
              else "two, and only two"))
}
stopifnot(chi(2) == 2, chi(3) == 0, chi(4) == 2)
cat("\n   Even D gives an even-dimensional transverse sphere and forces two. Odd D does not:\n")
cat("   S^3 carries free Z_n actions for every n, which is what a lens space is. So a\n")
cat("   four-sheeted fold is a five-dimensional possibility and not a four-dimensional one.\n")

cat("\n=== 2. the check must be able to fail ===\n\n")
cat("   The argument is that |G| divides chi. If that were wrong, S^2 would carry a free Z_3\n")
cat("   and the quotient would have Euler characteristic 2/3, which is not an integer:\n")
cat(sprintf("     chi(S^2)/3 = %.4f, not a whole number, so no such quotient is a closed surface.\n", 2/3))
stopifnot(2 %% 3 != 0)

cat("\n=== 3. what four sheets would have done to the growth, had they been available ===\n\n")
H0 <- 67.36; Om <- 0.3153; OL <- 0.6847
invH0 <- 9.77792/(H0/100)
age <- function(z) (2/3)*invH0/sqrt(OL)*asinh(sqrt(OL/Om)*(1+z)^(-1.5))
eps <- 0.1; tS <- 0.450*eps/(1-eps)
N1 <- (age(7) - age(20))/tS
cat("      sheets   kappa    e-folds    mass gain over one-sided\n")
for (k in c(2, 3, 4, 6)) cat(sprintf("      %-8d %-8.1f %-10.2f %.2e\n", k, k, k*N1, exp((k-1)*N1)))
cat(sprintf("\n   The overmassive holes sit one to two orders above the local relation. Two sheets\n"))
cat(sprintf("   already overshoot that by %.0e; four would overshoot by %.0e. More sheets move\n",
            exp(N1)/100, exp(3*N1)/100))
cat("   away from the observations, not towards them, so nothing is fixed by adding them.\n")
stopifnot(exp(3*N1) > exp(N1))
