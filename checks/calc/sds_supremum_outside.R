#!/usr/bin/env Rscript
# sds_supremum_outside.R -- REFEREE_OPEN item 5. The maximiser argument, that the
# connecting curve sits at X = 0, rests on g = sqrt(F)/r falling with r, where F is the
# Kruskal conformal factor of the patch in question. That scan had been run in the
# black-hole patch only. Schwarzschild-de Sitter has a second horizon, so this runs it in
# the patch adapted to the cosmological one and reports where it holds.
#
#   d log(F/r^2)/dr = (f'(r) - f'(r_h)) / f(r)  -  2/r,
# with r_h the horizon the patch is built on, since 2 kappa = |f'(r_h)| and the sign of the
# Kruskal exponent flips between an inward-facing and an outward-facing horizon. Writing it
# with f'(r_h) rather than 2 kappa handles both at once and makes the numerator vanish at
# the horizon, which is the regularity the chart needs.

M <- 1
fz   <- function(r,lam) 1 - 2*M/r - lam*r^2/3
fpz  <- function(r,lam) 2*M/r^2 - 2*lam*r/3
fppz <- function(r,lam) -4*M/r^3 - 2*lam/3
rts  <- function(lam) { z <- sort(Re(polyroot(c(2*M,-1,0,lam/3)))); z[z>0][1:2] }
dlog <- function(r, lam, rh) (fpz(r,lam) - fpz(rh,lam))/fz(r,lam) - 2/r

cat("=== 1. the scan, inside and outside ===\n")
cat("      9 Lam M^2    r_b      r_c    max dlog inside r_b    max dlog beyond r_c\n")
for (y in c(0.045, 0.18, 0.45, 0.4999, 0.5001, 0.7, 0.9, 0.99)) {
  lam <- y/9; q <- rts(lam); rb <- q[1]; rc <- q[2]
  ri <- seq(0.005*rb, 0.9995*rb, length.out=6000)
  ro <- exp(seq(log(1.0000005*rc), log(1e7*rc), length.out=12000))
  cat(sprintf("   %10.4f  %7.4f %8.4f  %20.5f  %+21.5f\n",
              y, rb, rc, max(dlog(ri,lam,rb)), max(dlog(ro,lam,rc)))) }
cat("\n   Inside the black-hole horizon it is comfortably negative at every Lambda, which\n")
cat("   is the previously checked case. Beyond the cosmological horizon it changes sign.\n")

cat("\n=== 2. where it changes sign, in closed form ===\n")
cat("   The worst point is the horizon itself. As r -> r_c from outside, both the\n")
cat("   numerator and f vanish linearly, so\n")
cat("      dlog -> f''(r_c)/f'(r_c) - 2/r_c,\n")
cat("   and the condition is r_c f''(r_c) < 2 f'(r_c). With f' = 2M/r^2 - 2 Lam r/3 and\n")
cat("   f'' = -4M/r^3 - 2 Lam/3 that reduces to Lam r_c^3 = 12M at equality, while\n")
cat("   f(r_c) = 0 gives Lam r_c^3 = 3 r_c - 6M. Equating, r_c = 6M and Lam = 1/(18 M^2),\n")
cat("   so the threshold is exactly\n\n")
cat("        9 Lam M^2 = 1/2,      r_c = 6M.\n\n")
lam0 <- 1/18; q0 <- rts(lam0)
cat(sprintf("   Check: at Lam = 1/18, f(6M) = %.1e and r_c = %.10f against 6.\n", fz(6,lam0), q0[2]))
stopifnot(abs(fz(6,lam0)) < 1e-12, abs(q0[2]-6) < 1e-9)
lim <- function(y) { lam <- y/9; rc <- rts(lam)[2]; fppz(rc,lam)/fpz(rc,lam) - 2/rc }
cat("\n      9 Lam M^2     limit of dlog at r_c+\n")
for (y in c(0.2, 0.4, 0.49, 0.5, 0.51, 0.6, 0.9)) cat(sprintf("   %10.4f   %+22.6e\n", y, lim(y)))
root <- uniroot(lim, c(0.2, 0.9), tol=1e-13)$root
cat(sprintf("\n   numerical sign change at 9 Lam M^2 = %.12f, against 0.5 exactly\n", root))
stopifnot(abs(root - 0.5) < 1e-9)

cat("\n=== 3. what the failure does and does not cost ===\n")
cat("   Above the threshold the sufficient condition fails near r_c; it does not follow\n")
cat("   that some other curve wins, only that this argument stops proving it does not.\n")
cat("   The direction is the safe one either way. If the maximiser leaves X = 0 it\n")
cat("   delivers MORE angle, so contact reaches further out and the computed radius\n")
cat("   becomes a lower bound rather than the answer. Nothing already claimed shrinks.\n")
cat("   It also lands where it matters least. The observed Lambda gives 9 Lam M^2 below\n")
cat("   10^-40 for any astrophysical hole, and at the other end, where 9 Lam M^2 -> 1,\n")
cat("   the contact region is the whole interior already, so a lower bound that cannot be\n")
cat("   improved is the answer.\n")
cat(sprintf("\n   For a ten solar mass hole at the observed Lambda, 9 Lam M^2 = %.2e.\n",
            9*1.1056e-52*(10*1.98847e30*6.6743e-11/(2.99792458e8)^2)^2))

cat("\n=== 4. where each chart does blow up, which is not at its own horizon ===\n")
cat("   At r_c the expression is continuous, because the numerator and f both have simple\n")
cat("   zeros there and the ratio tends to f''(r_c)/f'(r_c) from both sides. The blow-up in\n")
cat("   each chart is at the OTHER horizon, where the numerator does not vanish and f does.\n")
cat("   That is the chart failing to be regular where it was never regular, not a new gap.\n\n")
cat("      9 Lam M^2   cosmo patch at 1.001 r_b   BH patch at 0.999 r_c   cosmo patch at r_c+-\n")
for (y in c(0.045, 0.45, 0.9)) { lam <- y/9; q <- rts(lam); rb <- q[1]; rc <- q[2]
  cat(sprintf("   %10.4f  %24.1f  %22.4f  %15.5f / %.5f\n", y,
              dlog(1.001*rb,lam,rc), dlog(0.999*rc,lam,rb),
              dlog(0.999999*rc,lam,rc), dlog(1.000001*rc,lam,rc))) }
cat("\n   Columns one and two are the mismatched pairings and they diverge as expected.\n")
cat("   Column three shows the matched pairing agreeing across r_c to five figures, which\n")
cat("   is the continuity the closed form in section 2 relies on.\n")
cat("
=== flatly ===

  Item 5 is closed with a boundary rather than a clean yes. The supremum argument extends
  past the cosmological horizon exactly when 9 Lam M^2 <= 1/2, which is r_c >= 6M, and
  fails above that. The threshold is exact and falls halfway to Nariai.

  Above it the result survives as a lower bound, because a maximiser that leaves X = 0
  delivers more angle and enlarges the contact region rather than shrinking it. That is the
  same direction the charged near-extremal gap runs in, and the paper should say so in both
  places rather than in neither.\n")
