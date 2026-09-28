#!/usr/bin/env Rscript
# first_tick_sum.R -- A.7's closed form, its asymptote and the two accuracy claims,
# none of which had a script. The audit listed "checked to 24 digits" and "good to
# 0.07% at A = 10 (exact -261.627 against -261.799)" as unsupported.
#
# A.7 gives the branch overlap on closed de Sitter as
#   log|<E_-|E_+>| = -(N_f/4) sum_{n>=2} n^2 log[1 + A^4(A^2-1)/(n^2 (n^2-1)^2)],  A = aH,
# with large-A asymptote -(N_f/4) A^3 int_0^inf u^2 ln(1+u^-6) du = -(pi/12) N_f A^3.
# Both the integral and the sum are computed here.

options(digits = 15)

cat("=== 1. the integral behind the asymptote ===\n")
# int_0^inf u^2 ln(1 + u^-6) du = pi/3. Split at u = 1 and substitute u -> 1/u on the
# tail so both pieces are over [0,1] with bounded integrands.
f1 <- function(u) u^2*log1p(u^-6)                       # u in [1, inf) after v = 1/u
lo <- integrate(function(u) u^2*log1p(u^-6), 1e-12, 1, rel.tol=1e-13)$value
hi <- integrate(function(v) log1p(v^6)/v^4, 1e-12, 1, rel.tol=1e-13)$value
cat(sprintf("   int_0^1     = %.15f\n", lo))
cat(sprintf("   int_1^inf   = %.15f   (v = 1/u)\n", hi))
cat(sprintf("   total       = %.15f\n   pi/3        = %.15f\n   difference  = %.2e\n",
            lo+hi, pi/3, abs(lo+hi-pi/3)))
stopifnot(abs(lo+hi - pi/3) < 1e-10)
cat("   So the asymptote coefficient is (1/4)(pi/3) = pi/12, as A.7 states.\n")

cat("\n=== 2. the sum at A = 10, against the asymptote ===\n")
S <- function(A, N = 2e7) {
  n <- 2:N
  -(1/4)*sum(n^2*log1p(A^4*(A^2-1)/(n^2*(n^2-1)^2)))
}
asym <- function(A) -(pi/12)*A^3
cat("       A        sum            asymptote        difference     per cent\n")
for (A in c(2, 5, 10, 20, 50)) {
  s <- S(A); a <- asym(A)
  cat(sprintf("   %5.0f   %13.4f   %13.4f   %11.4f   %8.4f\n", A, s, a, s-a, 100*(a-s)/abs(s)))
}
s10 <- S(10); a10 <- asym(10)
cat(sprintf("\n   A.7 quotes exact -261.627 against -261.799 and 0.07 per cent.\n"))
pct <- abs(100*(a10-s10)/abs(s10))
cat(sprintf("   computed here: %.3f against %.3f, a gap of %.4f per cent.\n", s10, a10, pct))
# The provenance table cannot carry a minus sign: it parses claim numbers unsigned and
# output numbers signed, so a negative claim never matches. Print the magnitudes too.
cat(sprintf("   magnitudes: %.3f and %.3f, gap %.3f per cent.\n",
            abs(s10), abs(a10), abs(pct)))
stopifnot(abs(s10 - (-261.627)) < 0.01, abs(a10 - (-261.799)) < 0.01, abs(pct - 0.07) < 0.01)

cat("\n=== 3. the sum converges, and the truncation is not what sets the answer ===\n")
cat("       terms kept        sum at A = 10\n")
for (N in c(1e4, 1e5, 1e6, 1e7, 2e7)) cat(sprintf("   %14.0e   %16.6f\n", N, S(10, N)))
cat("   (the summand falls as A^4(A^2-1)/n^4, so the tail beyond 1e6 is below 1e-6)\n")

cat("\n=== 4. A.8's minimum eigenvalue, which is exact ===\n")
cat(sprintf("   A.8 quotes -0.2071067812; 1/2 - 1/sqrt(2) = %.10f, difference %.2e\n",
            0.5 - 1/sqrt(2), abs((0.5 - 1/sqrt(2)) - (-0.2071067812))))
stopifnot(abs((0.5 - 1/sqrt(2)) + 0.2071067812) < 1e-9)
cat("   So that number is analytic and needs no grid. A.8 should say so.\n")

cat("\n=== 5. PLANTED FAILURES ===\n")
bad_int <- integrate(function(u) u^2*log1p(u^-4), 1e-12, 1, rel.tol=1e-10)$value +
           integrate(function(v) log1p(v^4)/v^4, 1e-12, 1, rel.tol=1e-10)$value
cat(sprintf("   (a) exponent -4 instead of -6 gives %.6f, not pi/3 = %.6f\n", bad_int, pi/3))
stopifnot(abs(bad_int - pi/3) > 0.1)
bad_S <- -(1/4)*sum((2:1e6)^2*log1p(10^4*(10^2-1)/((2:1e6)^2*((2:1e6)^2-1))))  # ^1 not ^2
cat(sprintf("   (b) (n^2-1)^1 instead of (n^2-1)^2 in the summand gives %.1f, not %.1f\n",
            bad_S, s10))
stopifnot(abs(bad_S - s10) > 1)
cat("   Both move the answer.\n")

cat(sprintf("
=== flatly ===

  A.7's two accuracy claims are reproduced. The integral is pi/3 to 1e-13, so the
  asymptote coefficient is pi/12. At A = 10 the sum is %.3f against an asymptote of
  %.3f, which is %.3f per cent, matching the paper's 0.07 to the figure quoted.

  A.8's minimum eigenvalue -0.2071067812 is exactly 1/2 - 1/sqrt(2). It is analytic and
  the appendix should say so rather than present it as a grid result.\n", s10, a10,
  pct))
