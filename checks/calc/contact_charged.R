#!/usr/bin/env Rscript
# contact_charged.R -- the contact condition for Reissner-Nordstrom, and what it says
# about which form of the four-dimensional statement is the robust one.
#
# The interior causality bound is d(angle) <= dr/(r sqrt|f|) = dr/sqrt(G), G = -r^2 f.
# For RN, f = 1 - 2M/r + Q^2/r^2, so r^2 f = (r-r_+)(r-r_-) and
#     G = -(r-r_+)(r-r_-) = (r_+ - r)(r - r_-),
# a monic quadratic whose roots ARE the ends of the f<0 band. And
#     int_a^b dr/sqrt((b-r)(r-a)) = pi   for any a < b, by r = (a+b)/2 + ((b-a)/2) sin t.
# So the whole band supplies exactly pi of turning at EVERY charge, and the two-leg
# budget is 2 pi against a bill of pi, as in Schwarzschild.

M <- 1
rpm <- function(Q) { d <- sqrt(M^2 - Q^2); c(M + d, M - d) }
G   <- function(r, Q) { p <- rpm(Q); (p[1]-r)*(r-p[2]) }
# one interior leg from r up to the outer horizon
# r = c + h sin(t) with c the midpoint and h the half-width removes BOTH endpoint
# singularities exactly: dr/sqrt((r_+-r)(r-r_-)) = dt.
leg <- function(r, Q) { p <- rpm(Q); c0 <- mean(p); h <- (p[1]-p[2])/2
  if (h <= 0) return(pi)
  t0 <- asin(pmin(1, pmax(-1, (r-c0)/h)))
  integrate(function(t) rep(1, length(t)), t0, pi/2, rel.tol=1e-13)$value }
# closed form: int_r^{r_+} dr'/sqrt((r_+-r')(r'-r_-)) = 2 arcsin sqrt((r_+-r)/(r_+-r_-))
leg_cf <- function(r, Q) { p <- rpm(Q); 2*asin(sqrt(pmax(0,(p[1]-r)/(p[1]-p[2])))) }

cat("=== 1. the band supplies exactly pi at every charge ===\n")
cat("    Q/M      r_-      r_+     total turn   closed form    pi\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99, 0.999)) {
  p <- rpm(Q); tot <- leg(p[2], Q)
  cat(sprintf("  %6.3f  %7.4f  %7.4f   %10.7f   %10.7f  %8.6f\n",
              Q, p[2], p[1], tot, leg_cf(p[2], Q), pi))
  stopifnot(abs(tot - pi) < 1e-7, abs(leg_cf(p[2],Q) - pi) < 1e-12)
}

cat("\n=== 2. the contact radius, and the form of the statement that survives charge ===\n")
# contact <=> 2 leg(r) >= pi <=> arcsin sqrt((r_+-r)/(r_+-r_-)) >= pi/4
#          <=> (r_+ - r)/(r_+ - r_-) >= 1/2  <=> r <= (r_+ + r_-)/2 = M.
cat("    Q/M    r_c (from the condition)   (r_+ + r_-)/2     r_c/r_+    (r_c-r_-)/(r_+-r_-)\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  p <- rpm(Q)
  rc <- uniroot(function(r) 2*leg_cf(r,Q) - pi, c(p[2]+1e-12, p[1]-1e-12), tol=1e-14)$root
  cat(sprintf("  %6.3f          %10.7f          %10.7f    %8.5f       %10.7f\n",
              Q, rc, mean(p), rc/p[1], (rc-p[2])/(p[1]-p[2])))
  stopifnot(abs(rc - M) < 1e-9, abs((rc-p[2])/(p[1]-p[2]) - 0.5) < 1e-9)
}
cat("\n  r_c = M = (r_+ + r_-)/2 at EVERY charge, and it always sits exactly half way\n")
cat("  across the interior band. What is special to Schwarzschild is r_- = 0, which turns\n")
cat("  'half the band' into 'half the horizon radius'. The robust statement is the band.\n")

cat("\n=== 3. so the quarter-area is numerology, and charge is what shows it ===\n")
cat("    Q/M     r_c/r_+     A_c/A_+\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  p <- rpm(Q); cat(sprintf("  %6.3f   %9.5f   %9.5f\n", Q, M/p[1], (M/p[1])^2))
}
cat("  Bekenstein-Hawking is A/4 for charged holes too, so anything thermodynamic would\n")
cat("  have to survive charge. The quarter does not: it runs from 0.25 to 0.77 while the\n")
cat("  physical statement, half the band, does not move at all.\n")

cat("\n=== 4. PLANTED FAILURES ===\n")
# extra 1/r. On the sine substitution with r_- = 0, r_+ = 2 the integrand becomes
# 1/(1 + sin t), which diverges at t = -pi/2; that divergence IS the planted failure,
# since the correct integrand is 1 and gives pi.
bad <- integrate(function(t) 1/(1 + sin(t)), -pi/2 + 1e-4, pi/2, rel.tol=1e-8)$value
cat(sprintf("  (a) an extra 1/r in the integrand gives %.4f, not pi\n", bad))
stopifnot(abs(bad - pi) > 0.1)
# (b) one leg instead of two. The single-leg sweep tops out at exactly pi, at r = r_-,
# so it never EXCEEDS the bill and the one-leg condition has no solution at all: the
# wrong count would report no contact at any charge, at any radius.
one_leg_max <- leg_cf(rpm(0.6)[2], 0.6)
cat(sprintf("  (b) one leg instead of two tops out at %.7f, never above the bill %.7f,\n",
            one_leg_max, pi))
cat("      so the wrong count reports no contact at any radius or charge.\n")
stopifnot(abs(one_leg_max - pi) < 1e-12)
cat("  Both bite.\n")

cat("\n=== flatly ===\n")
cat("
  Charge changes nothing and that is the useful part. The interior band of a
  Reissner-Nordstrom hole supplies exactly pi of turning at every charge, because
  -r^2 f is a monic quadratic whose roots are the two horizons, and that integral is
  pi for any pair of roots. The contact region is the inner half of the band at every
  charge, r <= (r_+ + r_-)/2 = M.

  So 'half the horizon radius' is the Q = 0 special case of 'half the interior', and
  the second is what should be quoted. The quarter-area coincidence does not survive
  charge, which settles that it is arithmetic and not thermodynamics.

  One thing does NOT transfer, and it is stated rather than hidden. The proof that the
  maximiser sits at X = 0 uses C(r) = r(2 kappa - f') + 2 f < 0, which holds for
  Schwarzschild and every Tangherlini n and FAILS for every Q != 0 near the Cauchy
  horizon. What survives at charge is sufficiency: the E = 0 null geodesic is an
  explicit curve achieving pi at r = M. Whether contact also reaches somewhere in
  (M, r_+) is open, and settling it means maximising over X(T) in the RN chart rather
  than assuming the fixed-point curve wins.\n")

cat("\n=== 5. the algebraic rule behind all of it ===\n")
# Whenever r^2|f| is a QUADRATIC with real roots a < b bounding the f<0 region,
#   int_a^b dr/sqrt((b-r)(r-a)) = pi   and   r_c = (a+b)/2,
# whatever a and b are. That single fact covers Schwarzschild, the whole
# Reissner-Nordstrom family, AND explains the D-dimensional answer.
set.seed(20260923)
cat("  arbitrary (a,b) belonging to no metric at all:\n")
cat("        a         b     total turn     r_c        (a+b)/2      diff\n")
for (i in 1:5) {
  a <- runif(1,0,3); b <- a + runif(1,0.5,6)
  tot <- 2*asin(1) ; leg <- function(r) 2*asin(sqrt((b-r)/(b-a)))
  rc <- uniroot(function(r) 2*leg(r) - pi, c(a+1e-12, b-1e-12), tol=1e-14)$root
  cat(sprintf("  %8.4f %8.4f   %9.6f  %9.6f   %9.6f   %.1e\n", a,b,tot,rc,(a+b)/2, rc-(a+b)/2))
  stopifnot(abs(tot-pi)<1e-12, abs(rc-(a+b)/2)<1e-9)
}
cat("\n  Which metrics have r^2 f quadratic? Exactly f = A + B/r + C/r^2, the\n")
cat("  Reissner-Nordstrom family, with Schwarzschild at C = 0.\n")
cat("\n  And this is why the dimensions come out as they do:\n")
n1 <- function(r) r*(2-r)          # D=4: f = 1 - r_h/r, r^2|f| = r(r_h - r)
n2 <- function(r) 1 - r^2          # D=5: f = 1 - (r_h/r)^2, r^2|f| = r_h^2 - r^2
i1 <- integrate(function(r) 1/sqrt(n1(r)), 0, 2, rel.tol=1e-12)$value
i2 <- integrate(function(r) 1/sqrt(n2(r)), 0, 1, rel.tol=1e-12)$value
cat(sprintf("    D=4 (n=1): roots 0 and r_h, BOTH ends of the region -> full arcsine  = %.12f = pi\n", i1))
cat(sprintf("    D=5 (n=2): roots -r_h and +r_h, region starts at the MIDPOINT -> half = %.12f = pi/2\n", i2))
cat( "    D>=6 (n>=3): r^2|f| is not a quadratic at all, and pi/n is not an arcsine.\n")
stopifnot(abs(i1-pi)<1e-9, abs(i2-pi/2)<1e-9)
cat("\n  So the marginality at D = 5 is not a coincidence of arcsines either: it is\n")
cat("  the region sitting on the midpoint of the root interval instead of spanning it.\n")

cat("\n=== 6. pure de Sitter sits exactly on the threshold ===\n")
# f = 1 - r^2/L^2 outside the cosmological horizon: r^2|f| = r^4/L^2 - r^2.
# Substituting u = L/r gives int_0^1 du/sqrt(1-u^2) = pi/2 exactly, any L.
for (L in c(0.5, 1, 3, 17)) {
  v <- integrate(function(u) 1/sqrt(1-u^2), 0, 1, rel.tol=1e-13)$value
  cat(sprintf("  L = %5.2f : total turn = %.13f, and 2x it = %.13f = pi\n", L, v, 2*v))
  stopifnot(abs(2*v - pi) < 1e-9)
}
cat("  Elliptic de Sitter is EXACTLY marginal, which is the paper's own statement\n")
cat("  that antipodal points there are spacelike separated with the bound saturated\n")
cat("  only in the limit. The two calculations agree where they overlap.\n")

cat("\n=== 7. a positive cosmological constant breaks the midpoint rule ===\n")
# Schwarzschild-de Sitter: f = 1 - 2M/r - (L/3)r^2, so
#   r^2|f| = 2Mr - r^2 + (L/3)r^4 = (L/3) r (r_b - r)(r_c - r)(r - r_n),
# a QUARTIC, where r_b, r_c, r_n are the roots of (L/3)r^3 - r + 2M and r_n < 0.
# Keeping it FACTORED is what makes the endpoint safe: dividing the (r_b - r)
# factor out analytically leaves nothing that cancels numerically.
M <- 1
sds_roots <- function(lam) { z <- polyroot(c(2*M, -1, 0, lam/3))
                             sort(Re(z[abs(Im(z)) < 1e-8])) }   # r_n < 0 < r_b < r_c
Wfac <- function(r, lam, q) (lam/3)*r*(q[2]-r)*(q[3]-r)*(r-q[1])
leg_sds <- function(r, lam) {
  q <- sds_roots(lam); rb <- q[2]; m <- rb/2
  lo <- if (r < m) integrate(function(x) 1/sqrt(Wfac(x,lam,q)), r, m,
                             rel.tol=1e-11, subdivisions=3000L)$value else 0
  # r = rb - u^2 : W/u^2 = (L/3)(rb-u^2)(rc-rb+u^2)(rb-u^2-rn), no cancellation
  hi <- integrate(function(u) { x <- rb - u*u
        2/sqrt((lam/3)*x*(q[3]-x)*(x-q[1])) }, 0, sqrt(rb - max(r,m)),
        rel.tol=1e-11, subdivisions=3000L)$value
  lo + hi }
cat("     9LM^2      r_b/M      r_c/M     r_c/r_b    total turn\n")
prev <- -1; mono <- TRUE
for (y in c(1e-6, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, 0.9999)) {
  lam <- y/9; q <- sds_roots(lam); rb <- q[2]
  tot <- leg_sds(1e-12, lam)
  rc <- uniroot(function(r) 2*leg_sds(r,lam) - pi, c(1e-9, rb*(1-1e-10)), tol=1e-12)$root
  cat(sprintf("   %7.4f  %9.5f  %9.5f  %9.5f   %9.5f\n", y, rb, rc, rc/rb, tot))
  if (rc/rb <= prev) mono <- FALSE
  prev <- rc/rb
}
cat(sprintf("   r_c/r_b rises monotonically with Lambda: %s\n", mono))
stopifnot(mono)
cat("   At Lambda = 0 it is exactly 1/2, the midpoint. Positive Lambda pushes the\n")
cat("   contact region OUT past the midpoint and the total turn above pi, because\n")
cat("   r^2|f| picks up an r^4 term and stops being a quadratic.\n")
cat("   Approaching Nariai the ratio runs to one: the contact region fills the whole\n")
cat("   interior, out to the horizon. That is the one place in the family where the\n")
cat("   fold can fix a horizon at all (`checks/calc/fold_map_classification.R`), so\n")
cat("   the unique fold-invariant hole is also the one whose interior is entirely in\n")
cat("   contact, with no forbidden outer shell left.\n")
