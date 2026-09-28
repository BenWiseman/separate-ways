#!/usr/bin/env Rscript
# maximiser_charge_threshold.R -- REFEREE_OPEN item 4 carried a threshold, Q/M = 2 sqrt(6)/5
# = 0.9798, that appears nowhere else: not in the paper, not in any script, and not produced by
# the condition it was attached to. This works out what is actually true.
#
# The proof that the connecting curve's maximiser sits at X = 0 needs
#     C(r) = r(2 kappa - f'(r)) + 2 f(r) < 0
# inside the horizon, with 2 kappa = f'(r_+). For Reissner-Nordstrom, f = 1 - 2M/r + Q^2/r^2.

M <- 1
rpm <- function(Q) { d <- sqrt(M^2-Q^2); c(M+d, M-d) }
f   <- function(r,Q) 1 - 2*M/r + Q^2/r^2
fp  <- function(r,Q) 2*M/r^2 - 2*Q^2/r^3
C   <- function(r,Q) { rp <- rpm(Q)[1]; r*(fp(rp,Q) - fp(r,Q)) + 2*f(r,Q) }

cat("=== 1. Schwarzschild: the condition holds throughout the interior ===\n\n")
rs <- seq(1e-4, 2*M*(1-1e-9), length.out = 4000)
v <- sapply(rs, function(r) C(r,0))
cat(sprintf("   max C over (0, 2M): %.6f, attained at r/r_h = %.4f\n", max(v), rs[which.max(v)]/2))
cat("   C(2M) = 0 exactly, since 2 kappa = f'(r_+) there and f(r_+) = 0.\n")
stopifnot(max(v) < 1e-9, abs(C(2*M - 1e-9, 0)) < 1e-6)

cat("\n=== 2. charge breaks it immediately, at the Cauchy horizon ===\n")
cat("   For every Q other than zero, C is positive in a shell just outside r_-, because\n")
cat("   f'(r_-) is NEGATIVE, so the -r f'(r) term is positive and dominates. Section 4\n")
cat("   below gives the closed form. This is not a near-extremal\n")
cat("   effect and it does not switch on at any threshold:\n\n")
cat("      Q/M       r_-       max C on (r_-, r_+)    where, as r/r_+\n")
for (Q in c(0.01, 0.1, 0.3, 0.6, 0.9, 0.99)) {
  p <- rpm(Q); rr <- seq(p[2]*(1+1e-7), p[1]*(1-1e-7), length.out = 8000)
  vv <- sapply(rr, function(r) C(r,Q))
  cat(sprintf("   %7.3f  %8.5f  %20.4f  %16.5f\n", Q, p[2], max(vv), rr[which.max(vv)]/p[1]))
  stopifnot(max(vv) > 0) }
cat("\n   So the statement 'it fails near extremality' is too kind: it fails at every charge.\n")

cat("\n=== 3. the threshold that does exist, at the contact radius ===\n")
cat("   What can be asked is where the failure reaches the radius the result is about,\n")
cat("   r = M = (r_+ + r_-)/2, which Section 5 identifies as the contact boundary at every\n")
cat("   charge. Writing s = sqrt(1 - Q^2/M^2) = (r_+ - M)/M and q = Q^2/M^2 = 1 - s^2,\n")
cat("      C(M) = 2s/(1+s)^2 - 4s^2,\n")
cat("   so C(M) = 0 at s != 0 requires 2s(1+s)^2 = 1, a cubic:\n\n")
cat("      2 s^3 + 4 s^2 + 2 s - 1 = 0.\n\n")
cub <- function(s) 2*s^3 + 4*s^2 + 2*s - 1
s0 <- uniroot(cub, c(0.1, 0.5), tol=1e-14)$root
Q0 <- sqrt(1 - s0^2)
cat(sprintf("      s  = %.12f\n      Q/M = sqrt(1 - s^2) = %.12f\n", s0, Q0))
cat(sprintf("      direct solve of C(M) = 0:          %.12f\n",
            uniroot(function(Q) C(M,Q), c(0.5,0.999), tol=1e-14)$root))
stopifnot(abs(Q0 - uniroot(function(Q) C(M,Q), c(0.5,0.999), tol=1e-14)$root) < 1e-9)
cat("\n   and the closed form against the direct evaluation, at four charges:\n\n")
cat("      Q/M      C(M) direct     2s/(1+s)^2 - 4 s^2\n")
for (Q in c(0.3, 0.7, 0.9, 0.99)) {
  s <- sqrt(1-Q^2)
  cat(sprintf("   %7.3f  %14.8f  %19.8f\n", Q, C(M,Q), 2*s/(1+s)^2 - 4*s^2))
  stopifnot(abs(C(M,Q) - (2*s/(1+s)^2 - 4*s^2)) < 1e-10) }

cat("\n=== 4. and the number the referee file carried ===\n")
cat(sprintf("   2 sqrt(6)/5 = %.10f, which is Q^2/M^2 = 24/25 and r_+ = 6M/5, r_- = 4M/5.\n", 2*sqrt(6)/5))
cat(sprintf("   C(M) there is %+.6f, comfortably positive, so it is not where anything changes.\n", C(M, 2*sqrt(6)/5)))
cat("   The threshold is the cubic root above. Nothing in the repository produced 2 sqrt(6)/5\n")
cat("   and the referee file is corrected.\n")

cat("\n=== 5. the check has to be able to fail ===\n")
for (bad in c("drop the 2f", "use f' at r instead of r_+")) {
  Cb <- if (bad == "drop the 2f") function(r,Q) { rp <- rpm(Q)[1]; r*(fp(rp,Q)-fp(r,Q)) }
        else function(r,Q) 2*f(r,Q)
  d <- abs(Cb(M,0.9) - C(M,0.9))
  cat(sprintf("     %-28s changes C(M) at Q=0.9 by %.4f   <- as it must\n", bad, d))
  stopifnot(d > 1e-3) }

cat("
=== flatly ===

  Item 4's threshold does not exist. C(r) < 0 holds throughout the Schwarzschild interior
  and fails for EVERY nonzero charge in a shell outside the Cauchy horizon, so there is no
  charge below which the proof is safe.

  What does exist is a threshold at the radius the result is about. At r = M, the contact
  boundary at every charge, C(M) = 2s/(1+s)^2 - 4s^2 with s = sqrt(1 - Q^2/M^2), which
  vanishes at the root of 2s^3 + 4s^2 + 2s - 1 = 0, giving Q/M = 0.954828785515. Below that
  the condition holds at the contact boundary; above it, it does not hold even there.

  2 sqrt(6)/5 = 0.9798 is Q^2 = 24/25, where r_+ = 6M/5 and r_- = 4M/5. C(M) is positive
  there and nothing changes at it. The number appears in no script and in no version of the
  paper, only in the open-items file, and is withdrawn.\n")

# ---------------------------------------------------------------------------
# 4. WHY it fails there, stated correctly (added 2026-09-24).
#
# Sections 2 above and the paper both said the shell appears "because f'(r) diverges"
# at the Cauchy horizon. That is wrong: f' is finite at r_-. What is true is sharper and
# gives a one-line proof, so it replaces the old wording.
#
#   f(r)  = 1 - 2M/r + Q^2/r^2,     f'(r) = (2/r^3)(M r - Q^2)
#   r_+ r_- = Q^2  and  r_+ + r_- = 2M, so at r = r_-:
#       M r_- - Q^2 = M r_- - r_+ r_- = r_-(M - r_+) = -r_- sqrt(M^2 - Q^2)
#   hence
#       f'(r_-) = -2 sqrt(M^2 - Q^2) / r_-^2 ,
#   finite, and NEGATIVE, because f falls through zero as r increases past r_-.
#   Since f(r_-) = 0,
#       C(r_-) = r_-(2 kappa - f'(r_-)) = 2 kappa r_- + 2 sqrt(M^2 - Q^2)/r_- > 0
#   for every Q != 0. No threshold, no near-extremal story, and nothing diverging:
#   the sign of f' at the inner horizon does all the work. It is large only because
#   r_- is small at small charge, which is why the shell is worst as Q -> 0.

cat("\n=== 4. the closed form at the Cauchy horizon ===\n\n")
M <- 1
cat("     Q/M        r_-     f'(r_-) closed   f'(r_-) numeric      C(r_-) closed    C(r_-) direct\n")
for (Q in c(0.01, 0.1, 0.3, 0.6, 0.9, 0.99)) {
  rm_ <- M - sqrt(M^2 - Q^2); rp_ <- M + sqrt(M^2 - Q^2)
  f   <- function(r) 1 - 2*M/r + Q^2/r^2
  fp  <- function(r) 2*M/r^2 - 2*Q^2/r^3
  kap <- (rp_ - rm_) / (2 * rp_^2)                 # outer surface gravity
  fp_closed <- -2 * sqrt(M^2 - Q^2) / rm_^2
  C_closed  <- 2 * kap * rm_ + 2 * sqrt(M^2 - Q^2) / rm_
  C_direct  <- rm_ * (2 * kap - fp(rm_)) + 2 * f(rm_)
  cat(sprintf("  %6.3f  %9.6f  %15.6e  %15.6e  %15.6f  %15.6f\n",
              Q, rm_, fp_closed, fp(rm_), C_closed, C_direct))
  stopifnot(abs(fp_closed - fp(rm_)) < 1e-6 * abs(fp_closed))
  stopifnot(abs(C_closed - C_direct) < 1e-6 * abs(C_closed))
  stopifnot(C_closed > 0)
}
cat("\n  f' is finite at r_- at every charge, and negative. C(r_-) > 0 follows in one line.\n")
