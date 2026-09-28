#!/usr/bin/env Rscript
# restricted_exponent_ratio.R -- the companion records, for whoever revisits the
# withdrawn 19 per cent delay, that the ratio of the restricted first-tick exponent
# to the unrestricted one tends to one half at large A. That sentence had no script.
# This is the script.
#
# A.7's exponent is
#     S(A) = -(N_f/4) sum_{n>=2} g_n log[1 + A^4(A^2-1)/(n^2 (n^2-1)^2)]
# with g_n = n^2 the full S^3 degeneracy at level n = k+1.
#
# P_perp acts on the l multiplet inside level n as (-1)^l, so its invariant
# dimension is the sum of (2l+1) over even l up to k = n-1. That is one way to
# restrict. The other way in circulation is to keep whole levels of one parity and
# drop the rest. The companion's second number, 0.4998, is the ODD levels; the even
# ones give 0.5002, the complement. Both are computed.
# N_f cancels in every ratio below, so it is set to 1.

summand <- function(n, A) log1p(A^4*(A^2-1)/(n^2*(n^2-1)^2))

g_full <- function(n) n^2
g_inv  <- function(n) { k <- n-1; l <- seq(0, k, by=2); sum(2*l+1) }   # P_perp-even
g_lvl  <- function(n) ifelse(n %% 2 == 1, n^2, 0)                      # whole ODD levels

total <- function(g, A, nmax) sum(vapply(2:nmax, function(n) g(n)*summand(n,A), 0))

cat("\n  ratio of restricted exponent to unrestricted, against A\n")
cat("      A     invariant-dim    whole-level    (both -> 1/2)\n")
for (A in c(5, 10, 20, 40, 80, 160)) {
  nmax <- max(2000, ceiling(60*A))
  f <- total(g_full, A, nmax)
  cat(sprintf("  %5.0f      %10.4f     %10.4f\n", A,
              total(g_inv, A, nmax)/f, total(g_lvl, A, nmax)/f))
}

nmax <- 2400
A <- 40; f <- total(g_full, A, nmax)
r_inv <- total(g_inv, A, nmax)/f; r_lvl <- total(g_lvl, A, nmax)/f
cat(sprintf("\n  AT A = 40:  invariant-dimension %.4f,  whole-level %.4f\n", r_inv, r_lvl))

# the asymptotic halves, checked directly: both degeneracies average to n^2/2
k <- 2:4000
cat(sprintf("  mean(g_inv/n^2) over n=2..4000 = %.6f ; mean(g_lvl/n^2) = %.6f\n",
            mean(vapply(k, function(n) g_inv(n)/n^2, 0)),
            mean(vapply(k, function(n) g_lvl(n)/n^2, 0))))

# PLANTED FAILURE: a restriction that is NOT half must not come out at half.
g_third <- function(n) ifelse(n %% 3 == 0, n^2, 0)
r_third <- total(g_third, 40, nmax)/f
cat(sprintf("\n  CHECK (planted): keeping every third level gives %.4f, which must not be near 0.5\n", r_third))
stopifnot(abs(r_third - 0.5) > 0.1)
stopifnot(abs(r_inv - 0.5) < 0.02, abs(r_lvl - 0.5) < 0.02)
cat("  Planted third-level restriction is far from a half; both real ones are close to it.\n")

stopifnot(abs(r_inv - 0.4999) < 5e-5, abs(r_lvl - 0.4998) < 5e-5)
cat(sprintf("\n=== flatly ===\n
  Both ways of dropping modes give %.4f and %.4f at A = 40, which are the two
  numbers the companion quotes. They had no script until now. The point they
  support is unchanged and is the only one being made: the large-A limit is one
  half for both restrictions, so that limit never distinguished which half had
  been dropped and would not have caught the error it was being used to check.\n",
  r_inv, r_lvl))
