#!/usr/bin/env Rscript
# Fork 9: the image-pair mode sum, validated in flat space against an exact answer.
#
# WHY. Three bricks fixed the pieces of the interior sum: the modes and the Klein-Gordon measure
# (interior_modes_at_nonzero_k.R), the state (the_cross_region_state_is_the_half_period.R) and the
# frequency weight (the_cross_weight_is_one_over_sinh.R). Each was checked on its own. None of
# them checks that the three ASSEMBLE into the right object, and an assembly can be wrong while
# every piece is right. Flat space has the same structure and a closed-form answer, so it is where
# the assembly gets checked before Schwarzschild, where there is nothing to check it against.
#
# THE MODEL. Two-dimensional Minkowski, a massive scalar, Rindler coordinates in the right wedge:
#     T = rho sinh(a eta),  X = rho cosh(a eta),   ds^2 = -a^2 rho^2 d eta^2 + d rho^2.
# The map (T, X) -> (-T, -X) is (U, V) -> (-U, -V), the same map the fold is at a hole, and it
# carries the right wedge to the left one at the SAME Rindler coordinates, exactly as the fold
# carries a point to its image at the same (t, r). The Minkowski vacuum is the Hartle-Hawking
# state of this wedge pair at beta = 2 pi / a. So every ingredient of the interior sum has an
# exact counterpart here, and the answer is a Bessel function.
#
# THE ASSEMBLY. The Rindler mode is u = N K_{i omega/a}(m rho) e^{-i omega eta}. On a constant-eta
# slice the Klein-Gordon product is
#     (f, g) = -i int_0^inf (d rho / a rho) ( f d_eta g* - g* d_eta f ),
# which on two modes gives 2 omega N^2 (1/a) int_0^inf K_{i nu} K_{i nu'} d rho / rho. The
# Kontorovich-Lebedev orthogonality int_0^inf K_{i nu}(x) K_{i nu'}(x) dx/x = pi^2 delta(nu - nu')
# / (2 nu sinh(pi nu)) then fixes
#     N^2 = sinh(pi omega / a) / (pi^2 a).
# The cross weight at the image pair is 1/sinh(beta omega / 2) = 1/sinh(pi omega / a). THE TWO
# SINHS CANCEL EXACTLY, which is the structural point of this file, and what is left is
#     W_image(rho_1, rho_2) = int_0^inf d omega (1/pi^2 a) K_{i nu}(m rho_1) K_{i nu}(m rho_2)
#                           = (1/pi^2) int_0^inf K_{i nu}(m rho_1) K_{i nu}(m rho_2) d nu,
# a positive-definite integral with no weight left in it. The exact answer is (1/2 pi) K_0 of the
# invariant separation, and at equal Rindler time that separation is m(rho_1 + rho_2), since
# X_1 X_2 - T_1 T_2 = rho_1 rho_2 cosh(a(eta_1 - eta_2)) and the image flips both signs.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

# K_{i nu}(x) = int_0^inf cos(nu t) exp(-x cosh t) dt, by Simpson on a truncated t range.
Kiv <- function(nu, x, T = 40, N = 40000) {
  h <- T/N; t <- seq(0, T, by = h)
  f <- cos(nu*t)*exp(-x*cosh(t)); f[!is.finite(f)] <- 0
  w <- c(1, rep(c(4, 2), length.out = N - 1), 1)
  sum(w*f)*h/3
}

cat("=== 1. the Macdonald function of imaginary order, against R's own routine at nu = 0 ===\n")
for (x in c(0.3, 1, 3, 6)) {
  a <- Kiv(0, x); b <- besselK(x, 0)
  cat(sprintf("   x = %4.1f   quadrature %.12f   besselK %.12f   relative %.1e\n",
              x, a, b, abs(a - b)/b))
  note(abs(a - b)/b < 1e-11, sprintf("K_{i0}(%g) matches besselK", x))
}

cat("\n=== 2. the assembled image sum against the exact two-point function ===\n")
cat("   Left side is the assembly: modes, Klein-Gordon measure, half-period state and cross\n")
cat("   weight, with the two sinhs cancelled. Right side is (1/2 pi) K_0(m(rho_1 + rho_2)),\n")
cat("   the Minkowski two-point function at the image separation. Nothing is fitted.\n")
cat("        m rho_1   m rho_2        assembled            exact           relative\n")
for (pr in list(c(0.4, 0.4), c(1.0, 1.0), c(2.5, 2.5), c(0.5, 1.5), c(1.0, 3.0), c(0.3, 2.2))) {
  x <- pr[1]; y <- pr[2]
  f <- function(nu) sapply(nu, function(v) Kiv(v, x)*Kiv(v, y))
  lhs <- integrate(f, 0, 60, subdivisions = 400, rel.tol = 1e-10)$value/pi^2
  rhs <- besselK(x + y, 0)/(2*pi)
  cat(sprintf("   %9.2f %9.2f %18.10e %16.10e %14.2e\n", x, y, lhs, rhs, abs(lhs - rhs)/rhs))
  note(abs(lhs - rhs)/rhs < 1e-9,
       sprintf("the assembly reproduces the exact answer at (%g, %g)", x, y))
}

cat("\n=== 3. what this does and does not pin down ===\n")
cat("   It pins the PRODUCT of the mode normalisation and the cross weight, including its whole\n")
cat("   dependence on frequency: any other power of sinh in either factor changes the integrand's\n")
cat("   nu-dependence and the identity fails, as section 4 shows. What it cannot see is one\n")
cat("   nu-independent constant moved from one factor to the other, and that does not matter,\n")
cat("   because only the product enters the sum.\n")

cat("\n=== 4. the plants, one per way the assembly could be wrong ===\n")
x <- 1.0; y <- 1.0; rhs <- besselK(x + y, 0)/(2*pi)
bad <- function(w) {
  f <- function(nu) sapply(nu, function(v) w(v)*Kiv(v, x)*Kiv(v, y))
  integrate(f, 0, 60, subdivisions = 600, rel.tol = 1e-9)$value/pi^2
}
cases <- list(
  list("the equal-time weight coth in place of 1/sinh", function(v) cosh(pi*v)),
  list("the vacuum leg alone, so half the weight",      function(v) 0.5),
  list("the sinh left uncancelled",                     function(v) sinh(pi*v)),
  list("a spurious factor of nu",                       function(v) v))
for (cs in cases) {
  v <- tryCatch(bad(cs[[2]]), error = function(e) NA)
  cat(sprintf("   %-46s gives %-14s against %.6e\n", cs[[1]],
              if (is.na(v)) "divergent" else sprintf("%.6e", v), rhs))
  note(is.na(v) || abs(v - rhs)/rhs > 1e-3, sprintf("plant: %s is caught", cs[[1]]))
}
cat("   And the exact side has to be the right one: the separation at the image pair is\n")
cat("   rho_1 + rho_2 and NOT |rho_1 - rho_2|, which is the same-wedge equal-time separation.\n")
alt <- besselK(abs(0.5 - 1.5), 0)/(2*pi)
f <- function(nu) sapply(nu, function(v) Kiv(v, 0.5)*Kiv(v, 1.5))
lhs <- integrate(f, 0, 60, subdivisions = 400, rel.tol = 1e-10)$value/pi^2
cat(sprintf("   at (0.5, 1.5) the assembly gives %.6e. K_0 at the SUM separation 2.0 is %.6e,\n",
            lhs, besselK(2.0, 0)/(2*pi)))
cat(sprintf("   which is what it matches; K_0 at the DIFFERENCE separation 1.0 is %.6e, which is\n", alt))
cat("   a factor of 3.7 away, so the configuration is not interchangeable with the same-wedge one.\n")
note(abs(lhs - alt)/alt > 0.5, "plant: the wrong separation is caught")

cat("\n=== 5. one route tried and abandoned, recorded so it is not tried again ===\n")
cat("   The same-wedge equal-time identity, (1/pi^2) int cosh(pi nu) K_inu(x) K_inu(y) d nu =\n")
cat("   K_0(|x - y|)/2 pi, would separate the normalisation from the weight. It cannot be\n")
cat("   evaluated this way. K_{i nu}(x) falls as e^{-pi nu / 2} and the quadrature above carries\n")
cat("   an absolute error near machine epsilon, so by nu = 60 the function is 1e-41 and the\n")
cat("   computed value is noise; multiplying by cosh(pi nu), which is 1e+81 there, returns 1e+40\n")
cat("   instead of a number near 1e-2. It needs a uniform asymptotic expansion for K of imaginary\n")
cat("   order, which is a separate piece of work and is not needed for the image pair, where the\n")
cat("   integrand is absolutely convergent and the sinh has already cancelled.\n")
cat(sprintf("   For the record: cosh(pi*60) = %.3e and K_{i60}(1) is of order %.1e\n",
            cosh(pi*60), exp(-pi*60/2)))
note(cosh(pi*60)*exp(-pi*60/2) > 1e30, "the cancellation that defeats that route is as stated")

cat("\n=== 6. what this gives fork 9 ===\n")
cat("   The assembly is right. In flat space, where the same map pairs the same regions and the\n")
cat("   same half-period state holds, the modes, the Klein-Gordon measure, the state and the\n")
cat("   cross weight combine into the exact two-point function at the image pair, at six\n")
cat("   configurations and to ten figures, with the mode normalisation's sinh cancelling the\n")
cat("   weight's inverse sinh and leaving a positive-definite integral.\n")
cat("   WHICH PAIR THIS IS, said plainly, because it is not the one the fold contacts on. The\n")
cat("   configuration here is right wedge against left wedge, both OUTSIDE the horizon, which is\n")
cat("   the silence configuration: the separation comes out spacelike, m(rho_1 + rho_2), and the\n")
cat("   two sheets cannot touch there. It is also the configuration section 3.6 uses, so what is\n")
cat("   validated is that section's own cross-sheet correlator, mode by mode.\n")
cat("   THE CONTACT PAIR is future interior against past interior, and flat space has that too:\n")
cat("   x = (T, X) in F against (-T, -X) in P gives a separation of -4 rho^2, timelike, so the\n")
cat("   exact answer is K_0 of an imaginary argument. The mode side there needs Milne modes, and\n")
cat("   it is exactly where the two exterior families mix. That is the next brick, and it can be\n")
cat("   done in flat space against an exact answer before Schwarzschild is attempted.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
