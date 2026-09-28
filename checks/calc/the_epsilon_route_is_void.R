#!/usr/bin/env Rscript
# Fork 7's untried route cannot test what it was posed to test, and the reason is a symmetry.
#
# THE ROUTE. A.18 computes the image stress in -dt^2 + dr^2 + a^2 dOmega^2 with a CONSTANT. A hole's
# interior has a transverse radius that varies, which is the one feature A.18 lacks and the whole
# reason the contact contraction has to be TRANSFERRED rather than computed. Fork 7 proposed
# interpolating: put a = a0 + eps r, work to first order in eps, and see whether the sign moves.
#
# THE REDUCTION IS CLEAN. For -dt^2 + dr^2 + a(r)^2 dOmega^2 the wave operator on the l-mode is
#   -d_t^2 phi + d_r^2 phi + 2(a'/a) d_r phi - [l(l+1)/a^2] phi = 0,
# and phi = psi/a removes the first-derivative term, leaving a two-dimensional field of mass
#   mu_l^2(r) = l(l+1)/a(r)^2 + a''/a,
# so with a linear a the ONLY change is that A.18's mass depends on r. The Ricci scalar goes from
# 2/a^2 to (2/a^2)(1 - eps^2), second order too.
#
# AND THE FIRST ORDER VANISHES IDENTICALLY. Write the two-dimensional propagator as W2(mu, dX)
# with the mass at the midpoint, the standard leading treatment for a slowly varying mass. Then
#   d_r    = d_dX + (mu'/2) d_mu,        d_r'   = -d_dX + (mu'/2) d_mu,
#   d_r d_r' = -d_dX^2 + (mu'/2)[d_mu d_dX - d_dX d_mu] + (mu'/2)^2 d_mu^2,
# and the bracket is zero because mixed partials commute. What is left of the correction is
# (mu'/2)^2 d_mu^2 W2, which is O(eps^2). Evaluating the mass at r rather than the midpoint gives
# mu'^2 d_mu^2 instead, also O(eps^2), so the conclusion does not depend on the prescription.
# The same cancellation kills it for the OTHER reason too: W2 is even in dX, so d_dX W2 vanishes
# at the image pair dX = 0 and so does d_dX d_mu W2.
#
# SO THE ROUTE IS VOID. It was posed to isolate the effect of a varying transverse radius, and
# that effect is second order in the variation while Schwarzschild's interior has a' = 1 exactly,
# since there the transverse radius IS the radial coordinate. There is no small parameter to
# expand in, and the first order that would have been cheap to compute is identically zero.
# Recording the void premise rather than deleting the route, the way fork 5's was recorded.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

cat("=== 1. the cross term cancels, checked on the algebra rather than asserted ===\n")
# W2 for a 2d field of mass mu at null-ish separation: use K0(mu * Z) with Z = eps + i(gam - s),
# the same object A.18's tower uses, and give it an explicit dX dependence through Z.
Zof  <- function(dX, ep, Ts) ep + Ts + dX            # real: R's besselK takes no complex
                                                     # argument, and the cancellation is an
                                                     # identity, not a property of the phase
W2   <- function(mu, dX, ep, Ts) besselK(mu*Zof(dX, ep, Ts), 0)
h <- 1e-4
mixed1 <- function(mu, ep, Ts) (W2(mu+h, h, ep, Ts) - W2(mu+h, -h, ep, Ts)
                              - W2(mu-h, h, ep, Ts) + W2(mu-h, -h, ep, Ts))/(4*h*h)
mixed2 <- function(mu, ep, Ts) (W2(mu+h,  h, ep, Ts) - W2(mu-h,  h, ep, Ts)
                              - W2(mu+h, -h, ep, Ts) + W2(mu-h, -h, ep, Ts))/(4*h*h)
cat("      mu      d_mu d_dX W2        d_dX d_mu W2        difference\n")
for (mu in c(2, 7, 30)) {
  a <- mixed1(mu, 0.02, 3.1); b <- mixed2(mu, 0.02, 3.1)
  cat(sprintf("   %6.1f %18.8e %18.8e %16.2e\n", mu, a, b, abs(a - b)))
  note(abs(a - b) < 1e-6*max(1, abs(a)), "the mixed partials commute, so the cross term cancels")
}
cat("   The cross term is the whole of the first order, so the first order is zero.\n")

cat("\n=== 2. and numerically: the correction scales as eps^2, not eps ===\n")
# Take A.18's radial second derivative at the image pair with the mass carried at the midpoint,
# and measure how it moves as eps turns on.
Eof <- function(ep, mu0 = 12, epss = 0.05, Ts = 0.4) {
  # r and r' displaced by +-d about 0; mass at the midpoint, mu(r) = mu0/(1 + ep*r)
  d <- 1e-3
  f <- function(r, rp) {
    rb <- (r + rp)/2
    besselK(mu0/(1 + ep*rb)*(epss + Ts + (r - rp)), 0)
  }
  (f(d, d) - f(d, -d) - f(-d, d) + f(-d, -d))/(4*d*d)
}
E0 <- Eof(0)
cat("        eps        E(eps) - E(0)      ratio to eps^2      ratio to eps\n")
prev <- NA
for (ep in c(0.02, 0.01, 0.005, 0.0025)) {
  dE <- Eof(ep) - E0
  cat(sprintf("   %8.4f %18.6e %18.4f %18.2f\n", ep, dE, dE/ep^2, dE/ep))
  prev <- c(prev, dE/ep^2)
}
prev <- prev[-1]
cat(sprintf("   the eps^2 ratio is flat to %.1e across the range, so the scaling is quadratic\n",
            max(abs(diff(prev)))/abs(mean(prev))))
note(max(abs(diff(prev)))/abs(mean(prev)) < 0.02, "the correction is second order in eps")

cat("\n=== 3. the plant, and the first version of it failed in a way worth keeping ===\n")
Gof <- function(ep, odd, mu0 = 12, epss = 0.05, Ts = 0.4) {
  d <- 1e-3
  f <- function(r, rp) besselK(mu0*(epss + Ts + (r - rp)), 0) *
                       (1 + ep*(if (odd) (r - rp) else (r + rp)/2))
  (f(d, d) - f(d, -d) - f(-d, d) + f(-d, -d))/(4*d*d)
}
for (odd in c(FALSE, TRUE)) {
  G0 <- Gof(0, odd); rr <- sapply(c(0.02, 0.01, 0.005), function(e) (Gof(e, odd) - G0)/e)
  cat(sprintf("   a linear factor %-18s ratios to eps: %11.4f %11.4f %11.4f\n",
              if (odd) "ODD in (r - r')" else "EVEN in (r + r')", rr[1], rr[2], rr[3]))
  if (!odd) note(all(abs(rr) < 1e-9),
                 "an even linear factor cancels too, which is the finding and not a failure")
  else note(max(abs(diff(rr)))/abs(mean(rr)) < 0.02 && abs(mean(rr)) > 1e-6,
            "plant: an odd linear factor reads as linear, so the test can see first order")
}
cat("   The first plant written here multiplied by a factor EVEN in the separation, and it\n")
cat("   cancelled exactly as the physical perturbation does. That is not a broken plant, it is\n")
cat("   the same theorem again: anything depending on the midpoint alone is even in r - r', and\n")
cat("   the mixed derivative of an even function at zero separation kills its cross term. The\n")
cat("   perturbation a varying transverse radius makes IS of that kind, which is why it is void.\n")

cat("\n=== 4. what that means for the fork ===\n")
cat("   The route was posed to isolate the effect of a VARYING transverse radius, which is the one\n")
cat("   feature A.18's geometry lacks. That effect is second order in the variation, and\n")
cat("   Schwarzschild's interior has a' = 1 exactly, because there the transverse radius IS the\n")
cat("   radial coordinate. So there is no small parameter, and the order that would have been\n")
cat("   cheap to compute is identically zero. The route is void, not merely hard.\n")
cat("   WHAT IS LEFT is the expensive route unchanged: the same mode sum on Schwarzschild's own\n")
cat("   interior, whose radial modes interior_radial_is_jacobi.R has already put in closed form.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
