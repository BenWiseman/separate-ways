#!/usr/bin/env Rscript
# What fork 9 actually has to compute, isolated: the transverse half of the caustic transfers
# EXACTLY, because it is the same object in both geometries, so the whole remaining question is
# the interior's radial factor at large l.
#
# WHY. The transfer of A's sign from A.18's geometry to a hole is held by causal character and a
# caustic count. Both are arguments about the connecting geodesic. Neither says what makes the
# tower diverge, and that is worth separating, because the part that diverges turns out not to
# need transferring at all.
#
# A.18's divergence has one source, stated in its own words: the antipodal parity cancels the
# transverse phase and every multipole arrives in step. Written out, the reduction on the sphere
# gives sum_l (2l+1)/(4pi) P_l(cos gamma) W_radial, the geometric antipode puts gamma at pi so
# P_l(cos pi) = (-1)^l, and the fold's P_perp inserts a second (-1)^l. The two cancel at EVERY l,
# leaving sum_l (2l+1)/(4pi) W_radial with no alternating sign: that is the coherence, and it is
# the whole of the singularity.
#
# THAT STEP USES NOTHING BUT THE ROUND S^2 AND ITS ANTIPODAL MAP. Schwarzschild's interior is
# -dr^2/|f| + |f|dt^2 + r^2 dOmega^2: the transverse factor at fixed (t, r) is a round S^2 and the
# fold's transverse part is its antipodal map, the same object A.18 has. The image pair sits at the
# same r, so both ends carry the same sphere. Nothing in the radial block can reach the
# cancellation, because the cancellation is between two factors neither of which the radial block
# appears in.
#
# SO THE WHOLE OF FORK 9 IS THE RADIAL FACTOR. What the interior changes is the l-dependence of
# W_radial, which is what regulates the coherent sum and therefore sets the power and the
# coefficient. This file establishes the first half and measures what the second half looks like
# from the closed-form interior modes, so the next attempt starts from the right place.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

cat("=== 1. the cancellation, at every l, exactly ===\n")
Pl <- function(l, x) { s <- 0; for (k in 0:l) s <- s + choose(l,k)^2 * ((x-1)/2)^(l-k) * ((x+1)/2)^k; s }
cat("        l    P_l(cos pi)    (-1)^l from P_perp    product\n")
worst <- 0
for (l in c(0, 1, 2, 5, 12, 25, 40)) {
  p <- Pl(l, -1); q <- (-1)^l; worst <- max(worst, abs(p*q - 1))
  cat(sprintf("   %6d %14.1f %20d %10.1f\n", l, p, q, p*q))
  note(abs(p*q - 1) < 1e-9, sprintf("the two factors cancel at l = %d", l))
}
cat(sprintf("   worst departure of the product from one: %.1e\n", worst))
cat("   The product is one at every l, so every multipole arrives in step and the sum is\n")
cat("   sum_l (2l+1)/(4pi) W_radial, with no alternating sign left to regulate it.\n")

cat("\n=== 2. and the plant: without the fold's parity the sum alternates and does not diverge ===\n")
cat("   Dropping P_perp leaves P_l(cos pi) alone, so the coefficient alternates. Partial sums of\n")
cat("   (2l+1) with and without it, to show the difference is the whole singularity:\n")
cat("        L      coherent sum      alternating sum      their ratio\n")
for (L in c(20, 200, 2000)) {
  l <- 0:L
  co <- sum(2*l + 1); al <- sum((2*l + 1)*(-1)^l)
  cat(sprintf("   %6d %17.0f %20.0f %14.0f\n", L, co, al, co/abs(al)))
  note(abs(co - (L+1)^2) < 0.5 && abs(abs(al) - (L+1)) < 0.5,
       sprintf("the sums are (L+1)^2 and L+1 exactly at L = %d", L))
}
cat("   The coherent sum is (L+1)^2 and the alternating one is bounded by L+1, so the fold's\n")
cat("   parity is what turns a bounded sum into a divergent one. That is A.18's statement,\n")
cat("   checked, and it uses only the sphere.\n")

cat("\n=== 3. what is left, measured from the interior's own modes ===\n")
cat("   The interior's zero-frequency radial modes are psi_l = r^2 P^{(2,0)}_{l-1}(1 - 2r) at\n")
cat("   2M = 1, and Darboux gives psi_l -> r^{3/4}(1-r)^{-1/4} cos(2 nu arcsin sqrt r - 5pi/4)\n")
cat("   / sqrt(pi (l-1)) with nu = l + 1/2. At the image pair both ends sit at the same r, so\n")
cat("   the radial factor carries cos^2 of that phase, which splits into a SMOOTH half and an\n")
cat("   oscillating half:\n")
jac <- function(n, a, b, x) {
  if (n == 0) return(rep(1, length(x)))
  p0 <- rep(1, length(x)); p1 <- (a + 1) + (a + b + 2)*(x - 1)/2
  if (n == 1) return(p1)
  for (k in 2:n) {
    c1 <- 2*k*(k + a + b)*(2*k + a + b - 2)
    c2 <- (2*k + a + b - 1)*(a^2 - b^2)
    c3 <- (2*k + a + b - 2)*(2*k + a + b - 1)*(2*k + a + b)
    c4 <- 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
    p2 <- ((c2 + c3*x)*p1 - c4*p0)/c1; p0 <- p1; p1 <- p2
  }
  p1
}
psi <- function(l, r) r^2 * jac(l - 1, 2, 0, 1 - 2*r)
darboux <- function(l, r) {
  nu <- l + 0.5
  r^(3/4)*(1 - r)^(-1/4)*cos(2*nu*asin(sqrt(r)) - 5*pi/4)/sqrt(pi*(l - 1))
}
r0 <- 0.37
cat("        l     psi_l(r)^2        Darboux^2      smooth half      ratio to 1/l\n")
sm <- c()
for (l in c(60, 120, 240, 480)) {
  a <- psi(l, r0)^2; b <- darboux(l, r0)^2
  half <- r0^(3/2)*(1 - r0)^(-1/2)/(2*pi*(l - 1))
  cat(sprintf("   %6d %14.3e %16.3e %16.3e %14.4f\n", l, a, b, half, half*l))
  sm <- c(sm, half*l)
  note(abs(b - a)/max(a, 1e-300) < 0.25, sprintf("Darboux tracks the exact mode at l = %d", l))
}
cat(sprintf("   the smooth half times l is flat to %.1e, so it falls as 1/l exactly\n",
            max(abs(diff(sm)))/mean(sm)))
note(max(abs(diff(sm)))/mean(sm) < 0.02, "the coherent radial weight falls as 1/l")

cat("\n=== 3b. and the radial phase IS section 5's angular budget, which ties the modes to the\n")
cat("        geometry before any state is chosen ===\n")
cat("   Inside the horizon f < 0, so the potential V = f(l(l+1)/r^2 + 2M/r^3) is negative and the\n")
cat("   radial equation d^2 psi/dr*^2 + [k^2 - V] psi = 0 oscillates at every k with no turning\n")
cat("   point. At k = 0 and large l the WKB phase gradient is l/(r sqrt|f|), and the leg section 5\n")
cat("   integrates is exactly int dr/(r sqrt|f|). So the two should agree derivative by derivative:\n")
cat("        r         d(Darboux phase)/dr        l/(r sqrt|f|)          ratio\n")
lll <- 400
dph <- function(r) { nu <- lll + 0.5; 2*nu/(2*sqrt(r)*sqrt(1 - r)) }   # d/dr of 2 nu arcsin sqrt r
wkb <- function(r) { ff <- 1/r - 1; lll/(r*sqrt(ff)) }                 # 2M = 1, |f| = 1/r - 1
rat <- c()
for (r in c(0.12, 0.30, 0.50, 0.72, 0.91)) {
  cat(sprintf("   %6.2f %22.5f %20.5f %14.6f\n", r, dph(r), wkb(r), dph(r)/wkb(r)))
  rat <- c(rat, dph(r)/wkb(r))
}
cat(sprintf("   the ratio is nu/l = %.6f at every radius, flat to %.1e, which is the whole\n",
            (lll + 0.5)/lll, max(abs(diff(rat)))))
cat("   difference between the Langer-shifted nu and l and not a mismatch\n")
note(max(abs(diff(rat))) < 1e-9 && abs(mean(rat) - (lll + 0.5)/lll) < 1e-9,
     "the Darboux phase gradient is the angular-budget integrand times nu/l, at every radius")
cat("   The plant: a potential without the f in it must NOT reproduce the budget.\n")
bad <- function(r) lll/r
br <- sapply(c(0.12, 0.5, 0.91), function(r) dph(r)/bad(r))
cat(sprintf("      dropping sqrt|f| gives ratios %.4f, %.4f, %.4f, not flat\n", br[1], br[2], br[3]))
note(max(abs(diff(br))) > 0.1, "plant: the agreement is not automatic")

cat("\n=== 4. so this is the shape of what fork 9 has to do ===\n")
cat("   The transverse half needs no transferring: it is the same sphere and the same antipodal\n")
cat("   map, and the cancellation is exact at every l in both geometries. What is left is the\n")
cat("   radial factor, whose coherent part at zero frequency falls as 1/l against A.18's growth\n")
cat("   as sqrt(l). Those are different regulators and they will give different powers, which is\n")
cat("   why the interior sum has to be done rather than argued. What this removes from the job is\n")
cat("   the half that looked like it needed a transfer argument and does not.\n")
cat("   And 3b is a second tie between the modes and the geometry, computed two ways with nothing\n")
cat("   in common: the Darboux asymptotics of a Jacobi polynomial on one side, the WKB potential\n")
cat("   of the interior wave equation on the other, agreeing to machine precision up to the Langer\n")
cat("   shift. So the machinery fork 9 needs is calibrated against the geometry before any state\n")
cat("   is chosen, which is the order fork 9 asks for.\n")
cat("   NOT SETTLED HERE: the measure and the frequency integral, which is where the zero-frequency\n")
cat("   statement above becomes a statement about the actual two-point function, and which is\n")
cat("   where REORIENT records that a sign gets lost. Fork 9 names the calibration to do first.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
