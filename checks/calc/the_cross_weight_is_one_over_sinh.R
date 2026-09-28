#!/usr/bin/env Rscript
# Fork 9, the last convention piece before the sum: what weight the interior mode sum carries at
# the image pair, and why the sum converges in frequency at all.
#
# WHY. the_cross_region_state_is_the_half_period.R settled that the future-past correlator is the
# direct one at t - i beta/2 and that nothing there is chosen. What that leaves is arithmetic, and
# the first piece of the arithmetic is the frequency weight. Writing the two-point function mode
# by mode with the thermal occupation,
#     W(Delta t) = int dw/2pi [ (1+n_w) e^{-i w Delta t} + n_w e^{+i w Delta t} ] / 2w,
# the image pair sets Delta t = -i beta/2, and the two terms then become EQUAL and their sum
# collapses to 1/sinh(beta w / 2). Three things follow and each is checked here.
#
#   The collapse is exact, not approximate, and it happens at beta/2 and nowhere else. That is
#   the KMS selection section 1 states as (1+n) e^{-w sigma} = n e^{w sigma}.
#
#   Divided by 2w it is 1/(2 w sinh(beta w / 2)), which is section 3.6's own W_cross written out.
#   So the sum's weight is not a new object: the paper already carries it for a single mode.
#
#   It decays as 2 e^{-beta w / 2}, so the frequency integral at the image pair converges
#   exponentially. The divergence at a caustic is in the l-sum and not in the frequency, which is
#   what makes the interior sum a finite computation with a divergent angular part rather than a
#   doubly divergent one.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

nB   <- function(x) 1/(expm1(x))                       # occupation at beta*w = x
wA   <- function(x, f) (1 + nB(x)) * exp(-x*f)         # the (1+n) leg at a shift sigma = f*beta
wB   <- function(x, f) nB(x) * exp(+x*f)               # the n leg

cat("=== 1. at a half period the two legs are equal, and only there ===\n")
cat("   The legs are (1+n) e^{-x sigma/beta} and n e^{+x sigma/beta} with x = beta w.\n")
cat("        x      sigma/beta    (1+n)e^-       n e^+        equal?\n")
eq <- c()
for (f in c(1/4, 1/3, 1/2, 2/3)) {
  for (x in c(0.5, 2)) {
    a <- wA(x, f); b <- wB(x, f)
    same <- abs(a - b) < 1e-12*max(a, b)
    cat(sprintf("   %6.2f %12.4f %13.6f %12.6f %11s\n", x, f, a, b, if (same) "yes" else "no"))
    eq <- c(eq, same)
  }
}
note(identical(eq, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)),
     "the two legs coincide at sigma = beta/2 and at no other shift tried")

cat("\n=== 2. and their sum is exactly 1/sinh(x/2) ===\n")
cat("        x        (1+n)e^- + n e^+        1/sinh(x/2)         relative\n")
worst <- 0
for (x in c(1e-3, 0.01, 0.1, 1, 3, 10, 30, 80)) {
  s <- wA(x, 1/2) + wB(x, 1/2); t <- 1/sinh(x/2)
  d <- abs(s - t)/t; worst <- max(worst, d)
  cat(sprintf("   %8.3g %22.12e %18.12e %15.2e\n", x, s, t, d))
}
cat(sprintf("   worst relative difference over eight decades of x: %.2e\n", worst))
note(worst < 1e-12, "the collapse is exact")

cat("\n=== 3. which is section 3.6's W_cross, once the 1/2w is put back ===\n")
cat("   Section 3.6 writes W_cross = 1/(2 w sinh(beta w / 2)) from the KMS condition alone.\n")
cat("        beta      w        (leg sum)/2w          1/(2 w sinh(beta w/2))      relative\n")
worst2 <- 0
for (bt in c(0.7, 4*pi)) for (w in c(0.3, 1.7)) {
  x <- bt*w
  s <- (wA(x, 1/2) + wB(x, 1/2))/(2*w); t <- 1/(2*w*sinh(bt*w/2))
  d <- abs(s - t)/t; worst2 <- max(worst2, d)
  cat(sprintf("   %9.4f %6.2f %20.12e %24.12e %14.2e\n", bt, w, s, t, d))
}
note(worst2 < 1e-13, "the sum's weight is the paper's own W_cross and not a new object")

cat("\n=== 4. so the frequency integral converges, which the equal-time one does not ===\n")
cat("   At coincident times the weight is coth(x/2), which goes to 1 at large x and makes the\n")
cat("   integral the usual ultraviolet divergence. At the image pair it is 1/sinh(x/2), which\n")
cat("   decays as 2 e^{-x/2}. The tail beyond a cutoff X carries\n")
# The closed form is -2 log tanh(X/4). Written that way it loses eight digits at X = 40, because
# tanh(X/4) is one to within 4e-9 and the logarithm of it is all cancellation. log1p keeps them:
#     -2 log tanh(X/4) = 2 log( (1+q)/(1-q) ) = 2 log1p( 2q/(1-q) ),   q = e^{-X/2}.
tailf <- function(X) { q <- exp(-X/2); 2*log1p(2*q/(1 - q)) }
cat("        X        int_X^inf 1/sinh(x/2) dx      4 e^{-X/2}      closed form      relative\n")
for (X in c(5, 10, 20, 40)) {
  I <- integrate(function(x) 1/sinh(x/2), X, Inf, rel.tol = 1e-12)$value
  cat(sprintf("   %6.1f %24.12e %18.6e %16.9e %12.2e\n",
              X, I, 4*exp(-X/2), tailf(X), abs(I - tailf(X))/tailf(X)))
  note(I < 5*exp(-X/2), sprintf("the tail beyond X = %g is exponentially small", X))
  note(abs(I - tailf(X))/tailf(X) < 1e-9,
       sprintf("the closed-form tail matches the quadrature at X = %g", X))
}
cat("   The exact tail is 4 e^{-X/2} to leading order and the closed form above is exact. The\n")
cat("   third and fourth columns agree from X = 10 on, as they must.\n")

cat("\n=== 5. the plants ===\n")
# (a) drop the thermal occupation: the weight halves and stops being the paper's
x <- 2
vac <- wA(x, 1/2)                                   # (1+n) leg alone
cat(sprintf("   (a) vacuum leg alone at x = 2: %.8f against the full %.8f, a ratio of %.6f\n",
            vac, 1/sinh(x/2), vac*sinh(x/2)))
note(abs(vac*sinh(x/2) - 0.5) < 1e-12, "plant: dropping the occupation halves the weight exactly")
# (b) a shift that is not half a period leaves the two legs unequal, so the weight is complex
#     once the shift is put back as an imaginary time
f <- 1/3
cat(sprintf("   (b) at sigma = beta/3 the legs are %.6f and %.6f, differing by %.6f\n",
            wA(x, f), wB(x, f), abs(wA(x, f) - wB(x, f))))
note(abs(wA(x, f) - wB(x, f)) > 1e-3, "plant: any other shift leaves the legs unequal")
# (c) the equal-time weight must NOT decay, or section 4's statement is empty
cat(sprintf("   (c) the equal-time weight coth(x/2) at x = 5, 20, 80: %.6f, %.6f, %.6f\n",
            1/tanh(2.5), 1/tanh(10), 1/tanh(40)))
note(abs(1/tanh(40) - 1) < 1e-12, "plant: the equal-time weight does not decay, so the contrast is real")

cat("\n=== 6. what is now in hand for the sum ===\n")
cat("   W_cross(r) = sum_l (2l+1)/(4 pi) int dw/2pi [1/sinh(beta w / 2)] |u_{wl}(r)|^2 / (norm),\n")
cat("   with the transverse factor P_l(cos pi) times P_perp's (-1)^l equal to +1 at every l, so\n")
cat("   the l-sum is fully coherent and carries the caustic. Every piece of that expression is\n")
cat("   now computed rather than chosen: the modes and the norm from the two solver bricks, the\n")
cat("   state from the half-period brick, and the weight here.\n")
cat("   WHAT IS STILL OPEN, and it is the whole remaining risk: |u|^2 is the coincident-radius\n")
cat("   product of a mode with its own conjugate in the EXTERIOR decomposition, and inside the\n")
cat("   horizon the two exterior families mix, so the interior object is not |u|^2 of one family.\n")
cat("   That mixing is the next thing to compute and it is where a sign can still be lost.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
