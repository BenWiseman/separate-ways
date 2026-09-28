#!/usr/bin/env Rscript
# ===========================================================================
# SUPERSEDED, 2026-09-23, by pinch_4d_exact.R. DO NOT QUOTE 59.21.
#
# This script proves that exp(3 sqrt3 pi/4) = 59.2075 is exact rather than a
# bound, GIVEN that the connecting curve may sweep the transverse angle at the
# photon sphere r = 3M. It may not. Along a future-directed causal curve both
# Kruskal U and V increase, so a curve joining Theta x to x is confined to
# |UV| <= UV(r_x) <= 1, while |UV| at r = 3M is 2.241. The photon sphere is
# unreachable by any curve that does the job.
#
# The four-dimensional answer is r <= M, half the horizon radius, and the
# optimal curve never leaves the closed interior. See pinch_4d_exact.R.
# The arithmetic below is correct and is kept as the record of the argument.
# ===========================================================================
# pinch_boundary_exact.R
#
# The cross-sheet commutator is nonzero only inside a horizon AND only deep
# enough in. Lead 11 computes the depth by a sufficient construction: leave one
# region, sweep pi in angle on the photon sphere at r = 3M, and fall back. That
# costs coordinate time pi/Omega = 16.32 M, and with surface gravity 1/4M it
# multiplies the Kruskal V by exp(kappa dt) = 59.2.
#
# Lead 11 then calls 59.2 "an upper bound from an unstable circular orbit, not
# the optimal curve", and says geodesic maximisation would give the real value,
# "a half-day, not a paragraph".
#
# This script tests that. If no causal curve can sweep angle faster than the
# photon sphere does at ANY radius, then 16.32 M is not an upper bound on the
# cost, it is the minimum cost, and 59.2 is exact rather than bounding.

banner <- function() {
  cat("\n")
  cat("  ############################################################################\n")
cat("  ## SUPERSEDED by pinch_4d_exact.R, 2026-09-23. DO NOT QUOTE 59.21.\n")
  cat("  ############################################################################\n")
cat("  The result below assumes the connecting curve may sweep the transverse angle at the photon\n")
cat("  sphere r = 3M. It may not: along a future-directed causal curve both Kruskal coordinates\n")
cat("  increase, so a curve joining Theta x to x is confined to |UV| <= 1 while |UV| at r = 3M is\n")
cat("  2.241. The four-dimensional answer is r <= M, half the horizon radius.\n")
  cat("  ############################################################################\n\n")
}
banner()
M <- 1  # units of M throughout

# tangential null curve at radius r: ds^2 = 0 with dr = 0 gives
#   (1 - 2M/r) dt^2 = r^2 dphi^2   =>   dphi/dt = sqrt(1 - 2M/r)/r
omega <- function(r) sqrt(1 - 2*M/r)/r

cat("=== 1. is the photon sphere the fastest angular sweep available? ===\n")
rs <- seq(2.0001, 60, length.out = 400000)
w  <- omega(rs)
i  <- which.max(w)
cat(sprintf("  numerical maximum of dphi/dt at r = %.6f M   (analytic: 3 M)\n", rs[i]))
cat(sprintf("  value there            = %.9f /M\n", w[i]))
cat(sprintf("  analytic 1/(3 sqrt3 M) = %.9f /M\n", 1/(3*sqrt(3)*M)))
stopifnot(abs(rs[i] - 3) < 1e-3)
stopifnot(abs(w[i] - 1/(3*sqrt(3))) < 1e-9)

# analytic confirmation: d/dr log(omega) = M/(r^2(1-2M/r)) - 1/r, zero at r = 3M
dlog <- function(r) M/(r^2*(1-2*M/r)) - 1/r
cat(sprintf("  d(log omega)/dr at r=3M = %.3e  (must be 0)\n", dlog(3)))
stopifnot(abs(dlog(3)) < 1e-12)

cat("\n=== 2. validation: the check must be able to fail ===\n")
# plant: a curve that is NOT tangential-null must be slower. Add radial motion.
# for a null curve with dr != 0: (1-2M/r)dt^2 = dr^2/(1-2M/r) + r^2 dphi^2
omega_with_radial <- function(r, drdt) {
  val <- (1 - 2*M/r) - (drdt^2)/(1 - 2*M/r)
  ifelse(val > 0, sqrt(val)/r, NA)
}
w0 <- omega(3); w1 <- omega_with_radial(3, 0.10); w2 <- omega_with_radial(3, 0.30)
cat(sprintf("  tangential at r=3M          : %.9f\n", w0))
cat(sprintf("  with |dr/dt| = 0.10         : %.9f\n", w1))
cat(sprintf("  with |dr/dt| = 0.30         : %.9f\n", w2))
stopifnot(w1 < w0, w2 < w1)
cat("  radial motion strictly reduces the angular rate, as it must\n")

cat("\n=== 3. the consequence ===\n")
dt_min <- pi/(1/(3*sqrt(3)*M))          # minimum coordinate time to sweep pi
kappa  <- 1/(4*M)                        # surface gravity
ratio  <- exp(kappa*dt_min)
cat(sprintf("  minimum coordinate time to sweep pi : %.6f M   ( = 3 sqrt3 pi M )\n", dt_min))
cat(sprintf("  check against 3*sqrt(3)*pi          : %.6f\n", 3*sqrt(3)*pi))
cat(sprintf("  V multiplier exp(kappa dt)          : %.4f\n", ratio))
cat(sprintf("  closed form exp(3 sqrt3 pi / 4)     : %.4f\n", exp(3*sqrt(3)*pi/4)))
stopifnot(abs(ratio - exp(3*sqrt(3)*pi/4)) < 1e-9)

cat(sprintf("
=== 4. flatly ===

  No causal curve sweeps angle faster than a tangential null curve, and among
  tangential null curves the rate sqrt(1-2M/r)/r is maximised at exactly r = 3M,
  confirmed both analytically and on a 4e5-point grid. Radial motion strictly
  reduces the rate. So the photon sphere is not one construction among many; it
  is the optimum, and the angular sweep cannot be done faster by any route.

  That upgrades the project's own note. 59.2 is NOT an upper bound awaiting a
  geodesic maximisation. It is the value, in closed form:

        V_pinch / V_horizon  =  exp(3 sqrt(3) pi / 4)  =  %.2f

  parameter-free, independent of the hole's mass, and exact within the
  tangential-null reduction the result is stated in.

  WHAT THIS GIVES THE COMPANION. The two halves of a folded universe are
  correlated everywhere and can influence each other nowhere outside a horizon.
  Inside one they can, but not immediately: the region where they touch begins a
  definite depth in, and that depth is a pure number with no free parameter in it.
  A reader can hold that: the mirror is not merely hidden behind a horizon, it is
  hidden behind a horizon and then a further factor of sixty.

  WHAT IS STILL OWED. The s-wave reduction is exact and the full 4D statement is
  not proved here; lead 11 flags that and this script does not close it. What is
  closed is the optimisation, which was the part called a half-day of work.
", ratio))


banner()
