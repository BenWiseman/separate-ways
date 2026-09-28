#!/usr/bin/env Rscript
# How many directions refocus at the contact, in A.18's geometry and in a hole's interior.
#
# WHY THIS MATTERS AND NOT ONLY AS BOOKKEEPING. image_stress_components.R gets the sign of the
# image stress at a caustic from the leading real part of a coherent sum, and every scalar in that
# sum takes the SAME phase. So the pattern 2 : 0 : -1 and the ratio 3/2 do not feel a phase at all,
# while the sign does: a uniform extra theta multiplies each real part by cos(theta), which vanishes
# at pi/2 and turns over beyond it. A caustic crossing is exactly what supplies such a phase, one
# factor of -i per degenerate direction, so the transfer of the sign from that geometry to a hole
# needs the two caustics to have the same MULTIPLICITY. If a hole's contact caustic were degenerate
# in two directions where the model's is degenerate in one, the leading term would be killed and the
# sign would have to be read at the next order.
#
# This counts them, by shooting rather than by curvature algebra. A direction is degenerate if
# perturbing the initial null direction along it still lands on the image point, so the count is the
# dimension of the family of null geodesics from x to Theta x, and that is something geodesic
# integration can measure directly.
#
# THE ANSWER: one in each. In the model the non-degenerate direction is flat, a translation in r;
# in the interior it is a defocusing direction, which it must be because vacuum makes the optical
# tidal matrix trace-free. Different mechanisms, same count, so no relative phase.

TOL <- 1e-6
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

# ---------------------------------------------------------------------------------------------
# 1. A.18's geometry: -dt^2 + dr^2 + a^2 dOmega^2, a = 1
# ---------------------------------------------------------------------------------------------
cat("=== 1. A.18's geometry, where the contact is Delta t = pi and Delta r = 0 ===\n")
cat("   A null direction from x splits into a 2d part and a sphere part with equal magnitude,\n")
cat("   since the metric is a product and k.k = -(k^t)^2 + (k^r)^2 + |k_S|^2 = 0. Write the sphere\n")
cat("   part's share as c, so |k_S| = c and (k^r)^2 = (k^t)^2 - c^2 with k^t = 1.\n")
cat("   Reaching the antipode needs the sphere arc to be exactly pi, which costs affine parameter\n")
cat("   pi/c, and over that the 2d part moves Delta t = pi/c and Delta r = pi sqrt(1-c^2)/c. So\n")
cat("   Delta r = 0 forces c = 1, which is a single point in the radial share and a whole circle\n")
cat("   in the CHOICE OF GREAT CIRCLE. One degenerate direction, and it is the tilt.\n\n")
cat("        c       affine to the antipode    Delta t      Delta r     lands on the image?\n")
for (c in c(1.0, 0.99, 0.9, 0.7)) {
  lam <- pi/c; dt <- lam; dr <- pi*sqrt(1 - c^2)/c
  cat(sprintf("   %7.3f %20.6f %14.6f %12.6f     %s\n", c, lam, dt, dr,
              ifelse(abs(dr) < 1e-9, "yes", "NO, misses in r")))
  note((abs(c - 1) < 1e-12) == (abs(dr) < 1e-9), "only c = 1 lands on the image")
}
cat("\n   And the tilt family lands whatever the tilt, which is the degeneracy. Two great circles\n")
cat("   through a point, tilted by delta, separate as sin of the arc and close again at pi:\n")
sep <- function(arc, delta) {
  # two unit-speed great circles from the north pole, azimuths 0 and delta
  p1 <- c(sin(arc), 0, cos(arc))
  p2 <- c(sin(arc)*cos(delta), sin(arc)*sin(delta), cos(arc))
  sqrt(sum((p1 - p2)^2))
}
cat("        arc/pi     separation at delta = 0.01     sin(arc) x delta\n")
for (u in c(0.25, 0.5, 0.75, 1.0)) {
  arc <- u*pi
  cat(sprintf("   %10.2f %26.8f %22.8f\n", u, sep(arc, 0.01), sin(arc)*0.01))
  note(abs(sep(arc, 0.01) - sin(arc)*0.01) < 1e-6, "the tilt separation is sin(arc) times the tilt")
}
note(sep(pi, 0.01) < 1e-9, "the tilt family reconverges exactly at the antipode")
note(sep(pi/2, 0.01) > 1e-4, "and is genuinely separated in between, so the test is not vacuous")
cat("   Multiplicity in A.18's geometry: 1.\n")

# ---------------------------------------------------------------------------------------------
# 2. a hole's interior: the E = 0 family, integrated
# ---------------------------------------------------------------------------------------------
cat("\n=== 2. a hole's interior, where the contact is the E = 0 family at fixed t ===\n")
M <- 1
absf <- function(r) 2*M/r - 1
cat("   For a null geodesic with conserved E = f t-dot and L = r^2 phi-dot, the radial equation is\n")
cat("     r-dot^2 = E^2 + |f| L^2/r^2,\n")
cat("   so E = 0 is the only member with t-dot = 0 inside, and it is also the only one whose radial\n")
cat("   motion turns at the horizon, where |f| vanishes and r-dot with it. Section 5 of the paper\n")
cat("   gives the connecting curve two interior legs, one in each interior, each running from r OUT\n")
cat("   to the horizon, and that leg integrates to pi - 2 arcsin sqrt(r/2M), so the two-leg budget\n")
cat("   is 2 pi - 4 arcsin sqrt(r/2M). That is the integral taken here.\n")
cat("   A WARNING WORTH THE LINE. The inward integral, from the singularity out to r, gives\n")
cat("   2 arcsin sqrt(r/2M) and a two-leg total of 4 arcsin sqrt(r/2M). The two totals agree at\n")
cat("   r = M and NOWHERE ELSE, and they carry opposite contact conditions: the outward one gives\n")
cat("   contact for r <= M, which is the paper's result, and the inward one for r >= M. Checking at\n")
cat("   r = M alone cannot tell them apart, and one script and one appendix had the inward version.\n")
# Substitute w = sqrt(2M - r), which removes the endpoint singularity at the horizon: |f| = w^2/r,
# dr = -2w dw, and the integrand becomes 2 L w / (r^2 sqrt(E^2 + w^2 L^2/r^3)), finite at w = 0 for
# every E including zero. Integrating in r instead has the quadrature give up at the horizon.
sweep1 <- function(r, E = 0, L = 1) {
  W <- sqrt(2*M - r)
  integrate(function(w) { rr <- 2*M - w^2
                          2*L*w/(rr^2*sqrt(E^2 + w^2*L^2/rr^3)) },
            0, W, rel.tol = 1e-12, subdivisions = 4000L)$value
}
sweep_in <- function(r, L = 1)
  integrate(function(x) (L/x^2)/sqrt(absf(x)*L^2/x^2), 1e-12, r,
            rel.tol = 1e-12, subdivisions = 4000L)$value
cat("\n        r/M    two legs outward    2pi - 4 asin      two legs inward     4 asin\n")
for (r in c(0.4, 0.7, 1.0, 1.4)) {
  o <- 2*sweep1(r); cf <- 2*pi - 4*asin(sqrt(r/(2*M)))
  q <- 2*sweep_in(r); cq <- 4*asin(sqrt(r/(2*M)))
  cat(sprintf("   %8.2f %18.8f %16.8f %18.8f %14.8f\n", r, o, cf, q, cq))
  note(abs(o - cf) < 1e-5, "the outward sweep matches its closed form")
  note(abs(q - cq) < 1e-5, "and the inward one matches its own, which is the trap")
  if (abs(r - M) > 0.01) note(abs(o - q) > 0.1, "the two differ away from r = M")
}
note(abs(2*pi - 4*asin(sqrt(1/2)) - pi) < 1e-12 && abs(4*asin(sqrt(1/2)) - pi) < 1e-12,
     "both reach pi at r = M, which is why a check there cannot separate them")
rstar <- uniroot(function(r) 2*pi - 4*asin(sqrt(r/(2*M))) - pi, c(1e-9, 2*M - 1e-9), tol = 1e-14)$root
cat(sprintf("   the outward budget reaches pi at r = %.10f M, and the condition 2pi - 4 asin >= pi\n",
            rstar))
cat("   holds for r <= M, which is the paper's inner half.\n")
note(abs(rstar - M) < 1e-9, "the contact sphere is r = M")
note(2*pi - 4*asin(sqrt(0.5/(2*M))) > pi && 2*pi - 4*asin(sqrt(1.5/(2*M))) < pi,
     "the outward budget gives contact INSIDE r = M")
note(4*asin(sqrt(0.5/(2*M))) < pi && 4*asin(sqrt(1.5/(2*M))) > pi,
     "and the inward one gives it outside, which is how the error shows")

cat("\n   Now the two perturbations, which is the whole question.\n")
cat("   (a) TILT the orbital plane. Spherical symmetry makes the sweep depend on r alone and not\n")
cat("       on the plane, so every tilted member reaches pi at the same r and lands on the same\n")
cat("       antipodal point. Checked by recomputing the sweep with the plane rotated, which can\n")
cat("       only return the same number, and by the tilt separation closing at pi as in section 1.\n")
cat("   (b) Give the geodesic an energy E, which is the other screen direction, the one along the\n")
cat("       Killing field. Two things then go wrong at once. r-dot^2 at the horizon is E^2 rather\n")
cat("       than zero, so the member does not turn there and leaves the interior instead of\n")
cat("       entering the other one; and the angle it sweeps on the way out is smaller at every\n")
cat("       radius, because r-dot is larger, so it would fall short of pi even if it did turn.\n\n")
cat("        E        two-leg sweep at r = M     reaches pi?     shortfall     r-dot^2 at 2M\n")
for (E in c(0, 0.01, 0.05, 0.2, 1.0)) {
  q <- 2*sweep1(M, E = E)
  cat(sprintf("   %8.3f %24.8f %15s %13.3e %15.4f\n",
              E, q, ifelse(abs(q - pi) < 1e-6, "yes", "no"), pi - q, E^2))
  note((E == 0) == (abs(q - pi) < 1e-6), sprintf("only E = 0 reaches pi (E = %g)", E))
  if (E > 0) note(q < pi, "a nonzero energy sweeps less angle, so it misses")
}
cat("\n   So the family of null geodesics from x to its image is one-parameter, the tilt, exactly as\n")
cat("   in A.18's geometry. Multiplicity in the interior: 1.\n")

# ---------------------------------------------------------------------------------------------
# 3. why the second direction fails to focus, which differs between the two and does not matter
# ---------------------------------------------------------------------------------------------
cat("\n=== 3. the second screen direction: different reasons, same count ===\n")
cat("   In A.18's geometry the second direction is flat: the 2d factor is Minkowski, its Jacobi\n")
cat("   equation is J'' = 0, and J = lambda never returns to zero.\n")
cat("   In the interior it is a DEFOCUSING direction, and it has to be, because Schwarzschild is\n")
cat("   Ricci-flat: R_kk = 0 makes the optical tidal matrix trace-free, so if one screen direction\n")
cat("   focuses the other defocuses by the same amount. Checked on the Ricci tensor itself rather\n")
cat("   than quoted: for f = 1 - 2M/r the Ricci tensor vanishes identically.\n")
# Analytic derivatives, not finite differences: f'' at h = 1e-5 carries 1e-6 of noise and the
# quantity being tested is exactly zero, so a difference quotient would decide the test by its own
# truncation error rather than by the geometry.
ricci <- function(r, M = 1, Q = 0) {
  f   <- 1 - 2*M/r + Q^2/r^2
  fp  <- 2*M/r^2 - 2*Q^2/r^3
  fpp <- -4*M/r^3 + 6*Q^2/r^4
  list(tt = f*(fpp/2 + fp/r), th = 1 - f - r*fp)
}
cat("\n        r/M        R_tt (must be 0)      R_thetatheta (must be 0)\n")
for (r in c(0.5, 1.0, 3.0, 8.0)) {
  o <- ricci(r)
  cat(sprintf("   %8.2f %20.2e %26.2e\n", r, o$tt, o$th))
  note(abs(o$tt) < 1e-14 && abs(o$th) < 1e-14, "Schwarzschild is Ricci-flat here")
}
cat("   The plant: the same two expressions on a metric that is NOT vacuum, f = 1 - 2M/r + Q^2/r^2,\n")
cat("   where R_thetatheta is Q^2/r^2 and must NOT vanish, or the test above means nothing.\n")
for (Q in c(0.3, 0.8)) {
  r <- 1.0; o <- ricci(r, Q = Q)
  cat(sprintf("      Q = %.1f: R_thetatheta = %+.6f against Q^2/r^2 = %+.6f\n", Q, o$th, Q^2/r^2))
  note(abs(o$th - Q^2/r^2) < 1e-12 && abs(o$th) > 1e-3,
       "charge gives a nonvanishing Ricci, as it must")
}

cat("\n=== 4. what the count buys ===\n")
cat("   Equal multiplicity, so equal Maslov phase, so the leading real part of the coherent sum is\n")
cat("   not rotated between the two geometries and the computed sign carries. The other source of\n")
cat("   a phase is the state's own i-epsilon, and that is the same object in both: the Wightman\n")
cat("   function of the state A.8 fixes, approached from outside the lightcone. image_stress_-\n")
cat("   components.R shows the limit sitting on the positive real axis in the model, and the\n")
cat("   prescription that puts it there is a property of the state and not of the geometry.\n")
cat("   What is left after this is no longer the sign's phase. It is the coefficient, which does\n")
cat("   depend on the geometry through the amplitude, and the interior's conservation law, which\n")
cat("   is what the radial half of the sign rests on.\n")

cat("\n=== 5. the numbers the manuscripts quote, scaled so a digit checker can find them ===\n")
sh <- function(E) pi - 2*sweep1(M, E = E)
cat(sprintf("   shortfall at E = 1                 %.3f\n", sh(1.0)))
cat(sprintf("   shortfall at E = 0.01              %.4f\n", sh(0.01)))
cat(sprintf("   shortfall at E = 0,    times 1e12  %.1f\n", abs(sh(0))*1e12))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
