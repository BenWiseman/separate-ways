#!/usr/bin/env Rscript
# The interior's zero-frequency radial modes are Jacobi polynomials, and that hands over the one
# constant fork 7 is missing: the phase the coherent sum carries at a hole.
#
# WHY IT IS WANTED. image_stress_components.R gets the sign of the image stress at a caustic from
# the leading real part of a coherent sum, and every scalar in that sum takes the same phase. The
# pattern and the ratio do not feel a phase; the sign does, and a Wightman function's reality at
# spacelike separation forces any residual to be a multiple of pi, so the sign either carries from
# A.18's geometry to a hole or turns over. Deciding needs the phase on Schwarzschild, which needs
# the interior's radial modes, which are not usually available in closed form.
#
# THEY ARE HERE. At zero frequency and in the eikonal limit, where the 2M/r^3 piece of the
# potential is dropped beside l(l+1)/r^2, the Regge-Wheeler equation inside the horizon becomes,
# with 2M = 1 and r the ordinary radius,
#
#   r(1-r) psi'' - psi' + l(l+1) psi = 0.
#
# That is hypergeometric with singular points at 0, 1 and infinity, exponents {0,2} at the
# singularity, {0,0} at the horizon, which is the double exponent giving a constant and a
# logarithm there, and {l, -l-1} at infinity. The branch regular in the strong sense at r = 0 is
#
#   psi_l(r) = r^2 P^{(2,0)}_{l-1}(1 - 2r),
#
# a polynomial of degree l+1. Jacobi polynomials have Darboux asymptotics with a known phase, so
# the constant this file is after can be read off rather than fitted:
#
#   P_n^{(a,b)}(cos th) ~ cos( (n + (a+b+1)/2) th - (2a+1) pi/4 ) / (...),
#
# and with th = 2 arcsin sqrt r, n = l-1, a = 2, b = 0 the phase is 2 nu arcsin sqrt r - 5 pi/4,
# nu = l + 1/2. The variable part is the WKB phase the contact geodesic already carries, since
# nu (pi - 2 arcsin sqrt r) is A.19's leg from r out to the horizon. THE CONSTANT IS -5 pi/4.
#
# WHAT THIS FILE DOES AND DOES NOT CLAIM. It establishes the identification and the constant. It
# does NOT assemble the two-point function, which also needs the measure and the saddle of the
# frequency integral, and it does not on its own settle whether the sign carries. That assembly is
# the remaining step and this is its hardest ingredient.

TOL <- 1e-6
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

jacobi <- function(n, al, be, x) {           # the standard stable three-term recurrence
  if (n == 0) return(rep(1, length(x)))
  p0 <- rep(1, length(x)); p1 <- (al + 1) + (al + be + 2)*(x - 1)/2
  if (n == 1) return(p1)
  for (k in 2:n) {
    a1 <- 2*k*(k + al + be)*(2*k + al + be - 2)
    a2 <- (2*k + al + be - 1)*(al^2 - be^2)
    a3 <- (2*k + al + be - 2)*(2*k + al + be - 1)*(2*k + al + be)
    a4 <- 2*(k + al - 1)*(k + be - 1)*(2*k + al + be)
    p0t <- p1; p1 <- ((a2 + a3*x)*p1 - a4*p0)/a1; p0 <- p0t
  }
  p1
}
psi <- function(l, r, al = 2, be = 0, nshift = -1)
  r^2 * jacobi(l + nshift, al, be, 1 - 2*r)

ode_res <- function(l, r, al = 2, be = 0, nshift = -1, h = 1e-5) {
  f <- function(x) psi(l, x, al, be, nshift)
  d1 <- (f(r+h) - f(r-h))/(2*h); d2 <- (f(r+h) - 2*f(r) + f(r-h))/h^2
  (r*(1-r)*d2 - d1 + l*(l+1)*f(r)) / max(1e-300, abs(l*(l+1)*f(r)))
}

cat("=== 1. the polynomial solves the equation, at every l tried ===\n")
cat("      l      residual r=0.2     r=0.5        r=0.8\n")
for (l in c(3, 8, 20, 50, 120, 300)) {
  o <- sapply(c(0.2, 0.5, 0.8), function(r) ode_res(l, r))
  cat(sprintf("%7d %17.2e %12.2e %12.2e\n", l, o[1], o[2], o[3]))
  note(max(abs(o)) < 1e-4, sprintf("r^2 P^{(2,0)}_{l-1}(1-2r) solves the equation at l = %d", l))
}
cat("   The residual is finite-difference noise: it grows with l because the polynomial's own\n")
cat("   scale does, and it falls as h^2 when h is reduced, which the plant below uses.\n")

cat("\n=== 2. the plants: three neighbouring candidates must all fail ===\n")
cat("   A polynomial family that nearly works is exactly what a check has to be able to reject,\n")
cat("   so the two adjacent Jacobi indices and the adjacent degree are tried on the same equation.\n")
cat("      candidate                    residual at r = 0.5, l = 20\n")
for (p in list(list("r^2 P^{(2,0)}_{l-1}  (the claim)", 2, 0, -1),
               list("r^2 P^{(1,0)}_{l-1}", 1, 0, -1),
               list("r^2 P^{(3,0)}_{l-1}", 3, 0, -1),
               list("r^2 P^{(2,1)}_{l-1}", 2, 1, -1),
               list("r^2 P^{(2,0)}_{l}  ", 2, 0, 0))) {
  v <- ode_res(20, 0.5, p[[2]], p[[3]], p[[4]])
  cat(sprintf("      %-30s %14.4e   %s\n", p[[1]], v,
              ifelse(abs(v) < 1e-4, "solves it", "REJECTED, as it must be")))
  if (p[[2]] == 2 && p[[3]] == 0 && p[[4]] == -1) note(abs(v) < 1e-4, "the claim passes")
  else note(abs(v) > 1e-3, sprintf("%s is rejected", p[[1]]))
}

cat("\n=== 3. the zeros sit where the Darboux phase says, and the constant is -5 pi/4 ===\n")
cat("   A zero of cos(2 nu arcsin sqrt r - C) sits at r = sin^2( ((k+1/2) pi + C) / (2 nu) ).\n")
cat("   Every zero in (0,1) is located and compared with that, for C = 5 pi/4 and for the two\n")
cat("   neighbouring quarter-turns, which must do worse or the constant is not being measured.\n\n")
zeros <- function(l, n = 400000) {
  g <- seq(1e-7, 1 - 1e-7, length.out = n)
  s <- sign(psi(l, g)); k <- which(diff(s) != 0); (g[k] + g[k+1])/2
}
cat("      l   zeros  l-1     max error at C = 5pi/4    at 3pi/4      at 7pi/4\n")
for (l in c(8, 20, 50, 120, 300)) {
  nu <- l + 0.5; z <- zeros(l)
  err <- function(C) {
    kk <- round((2*nu*asin(sqrt(z)) - C)/pi - 0.5)
    max(abs(sin(((kk + 0.5)*pi + C)/(2*nu))^2 - z))
  }
  e5 <- err(5*pi/4); e3 <- err(3*pi/4); e7 <- err(7*pi/4)
  cat(sprintf("%7d %6d %5d %22.3e %13.3e %13.3e\n", l, length(z), l-1, e5, e3, e7))
  note(length(z) == l - 1, "the zero count is the Jacobi degree")
  note(e5 < e3 && e5 < e7, "5 pi/4 beats both neighbouring quarter-turns")
}
cat("   And the error at 5 pi/4 falls as the SQUARE of 1/l, because a phase error of order 1/nu\n")
cat("   moves a zero by that over dPhi/dr, which is itself of order nu:\n")
e <- sapply(c(20, 50, 120, 300), function(l) {
  nu <- l + 0.5; z <- zeros(l); C <- 5*pi/4
  kk <- round((2*nu*asin(sqrt(z)) - C)/pi - 0.5)
  max(abs(sin(((kk + 0.5)*pi + C)/(2*nu))^2 - z))
})
sl <- lm(log(e) ~ log(c(20, 50, 120, 300)))$coefficients[[2]]
cat(sprintf("      log-slope of the error against l: %.3f, against -2\n", sl))
note(abs(sl + 2) < 0.25, "the zero-position error is second order in 1/l")

cat("\n=== 3b. and the sign, because zeros fix a phase only modulo pi ===\n")
cat("   cos(x - C) and cos(x - C - pi) have the SAME zeros, so section 3 cannot tell 5 pi/4 from\n")
cat("   pi/4, and the whole question this file is for is a phase modulo TWO pi. The full Darboux\n")
cat("   form settles it, amplitude and all:\n")
cat("        P_n^{(a,b)}(cos th) = cos(Phi) / ( sqrt(pi n) sin^{a+1/2}(th/2) cos^{b+1/2}(th/2) )\n")
cat("   so with psi = r^2 P, th/2 = arcsin sqrt r, a = 2 and b = 0,\n")
cat("        psi_l(r) sqrt(pi (l-1)) r^{-3/4} (1-r)^{1/4}  ->  cos(2 nu arcsin sqrt r - C).\n")
cat("   The left side is computed and the right compared at both candidates.\n\n")
amp <- function(l, r) psi(l, r)*sqrt(pi*(l-1))*r^(-3/4)*(1-r)^(1/4)
cat("      l      max |lhs - cos(Phi)| at C = 5pi/4      at C = pi/4\n")
for (l in c(50, 120, 300, 800)) {
  nu <- l + 0.5; rr <- seq(0.08, 0.92, length.out = 4000)
  lhs <- amp(l, rr)
  d5 <- max(abs(lhs - cos(2*nu*asin(sqrt(rr)) - 5*pi/4)))
  d1 <- max(abs(lhs - cos(2*nu*asin(sqrt(rr)) - 1*pi/4)))
  cat(sprintf("%7d %28.4f %20.4f\n", l, d5, d1))
  note(d5 < 0.2 && d1 > 1.0, "the amplitude picks 5 pi/4 and rejects pi/4")
}
cat("   The rejected candidate sits near 2, which is the largest a difference of two cosines can\n")
cat("   be, so it is not a near miss: it is the opposite sign everywhere.\n")

cat("\n=== 4. what the constant is, and what it is not yet ===\n")
cat("   The interior's eikonal zero-frequency radial mode carries phase\n")
cat("        2 nu arcsin sqrt r - 5 pi/4 = nu pi - Phi_h(r) - 5 pi/4,\n")
cat("   where Phi_h = nu (pi - 2 arcsin sqrt r) is A.19's leg from r out to the horizon. So the\n")
cat("   radial factor's constant is -5 pi/4, in closed form and not fitted.\n")
cat("   A.18's geometry has no radial mode structure at all: its 2d factor is Minkowski and the\n")
cat("   phase there is the free propagator's own, e^{-i pi/4} out of sqrt(1/z). The two constants\n")
cat("   are not yet comparable, because the Schwarzschild two-point function also carries the\n")
cat("   measure and the saddle of the frequency integral, whose stationary point for equal-time\n")
cat("   endpoints is exactly E = 0 and which contributes a phase of its own. Assembling those is\n")
cat("   the remaining step of the transfer, and this file supplies its hardest ingredient.\n")
cat("   Stated plainly so nobody reads more into it: a difference of pi between the assembled\n")
cat("   constants would turn the sign over, and this file does not yet say whether there is one.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
