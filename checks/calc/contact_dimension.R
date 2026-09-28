#!/usr/bin/env Rscript
# contact_dimension.R -- the contact condition in D spacetime dimensions, and the
# dimension selection that falls out of it.
#
# Schwarzschild-Tangherlini in D dimensions has f(r) = 1 - (r_h/r)^n with n = D - 3,
# a simple zero at r_h and the same Penrose diagram as the four-dimensional case.
#
# Inside the horizon, causality gives d(angle) <= dr/(r sqrt|f|). Substituting
# u = r/r_h and then w = u^n turns the interior leg into a Beta integral:
#
#   Dphi_1(r) = int_r^{r_h} dr'/(r' sqrt(|f|))
#             = (1/n) int_{(r/r_h)^n}^{1} w^{-1/2}(1-w)^{-1/2} dw
#             = (2/n) [ pi/2 - arcsin( (r/r_h)^{n/2} ) ].
#
# So the WHOLE interior, r = 0 to r_h, supplies exactly pi/n of turning. In four
# dimensions that is pi: half a revolution, whatever the mass.
#
# The fold's transverse map is antipodal, so it demands pi. The connecting curve has
# two interior legs, giving 2 Dphi_1(r), and contact needs that to reach pi.
#
# Everything below is computed. Planted failures at the end.

nn  <- function(D) D - 3
f   <- function(r, rh, n) 1 - (rh/r)^n
dphi1_closed <- function(r, rh, n) (2/n)*(pi/2 - asin((r/rh)^(n/2)))
dphi1_num <- function(r, rh, n)
  integrate(function(x) 1/(x*sqrt(pmax(0, (rh/x)^n - 1))), r, rh,
            rel.tol=1e-11, stop.on.error=FALSE)$value

cat("=== 1. the interior leg, closed form against direct integration ===\n")
cat("      D    n    r/r_h     closed      numeric        diff\n")
for (D in 4:8) { n <- nn(D); for (x in c(0.05, 0.3, 0.7)) {
  a <- dphi1_closed(x, 1, n); b <- dphi1_num(x, 1, n)
  cat(sprintf("   %4d %4d    %5.2f   %9.6f   %9.6f   %9.2e\n", D, n, x, a, b, abs(a-b)))
  stopifnot(abs(a-b) < 1e-6) }}

cat("\n=== 2. the whole interior supplies pi/n ===\n")
cat("      D    n    total turn, horizon to singularity      = pi/n ?\n")
for (D in 4:9) { n <- nn(D)
  tot <- dphi1_closed(0, 1, n)
  cat(sprintf("   %4d %4d           %10.6f   (pi/n = %8.6f)\n", D, n, tot, pi/n))
  stopifnot(abs(tot - pi/n) < 1e-12) }
cat("  In D = 4 that is exactly pi, half a revolution, and it is the ONLY dimension\n")
cat("  in which the interior supplies the whole half-turn the antipodal map asks for.\n")

cat("\n=== 3. g is decreasing in r inside the horizon, in EVERY dimension ===\n")
# The maximiser argument needs g = sqrt(F)/r to fall with r, where F is the Kruskal
# conformal factor F = |f| e^{-2 kappa r*}/kappa^2. Then
#   d log g^2/dr = (f' - 2 kappa)/f - 2/r,
# and with f' = n r_h^n/r^{n+1} and 2 kappa = n/r_h this is
#   (n/r_h)[(r_h/r)^{n+1} - 1] / f  -  2/r,
# whose first term has a positive numerator and a negative denominator for r < r_h.
dlog <- function(r, rh, n) (n/rh)*((rh/r)^(n+1) - 1)/f(r,rh,n) - 2/r
cat("      D      max over r in (0.01, 0.99) r_h of dlog(g^2)/dr   (must be < 0)\n")
for (D in 4:9) { n <- nn(D)
  v <- max(vapply(seq(0.01, 0.99, by=0.005), function(x) dlog(x,1,n), 0))
  cat(sprintf("   %4d                    %12.4f\n", D, v)); stopifnot(v < 0) }
cat("  So X = 0 is the maximiser in every D, the connecting curve never leaves the\n")
cat("  closed interior, and the two interior legs are the whole of the budget.\n")

cat("\n=== 4. the contact condition, and what it selects ===\n")
# contact <=> 2 Dphi_1(r) >= pi <=> arcsin((r/r_h)^{n/2}) <= pi(2-n)/4
rc_over_rh <- function(n) { a <- pi*(2-n)/4; if (a < 0) NA_real_ else sin(a)^(2/n) }
cat("      D    n     r_c/r_h        contact region\n")
for (D in 4:9) { n <- nn(D); x <- rc_over_rh(n)
  msg <- if (is.na(x)) "NONE, anywhere" else if (x < 1e-12) "the singularity only" else
         sprintf("r <= %.4f r_h", x)
  cat(sprintf("   %4d %4d    %8s     %s\n", D, n, ifelse(is.na(x),"--",sprintf("%.6f",x)), msg)) }
stopifnot(abs(rc_over_rh(1) - 0.5) < 1e-12, abs(rc_over_rh(2)) < 1e-15, is.na(rc_over_rh(3)))
cat("\n  D = 4 gives r_c = r_h/2 = M exactly, reproducing the four-dimensional result.\n")
cat("  D = 5 gives contact only in the limit r -> 0, the singularity itself.\n")
cat("  D >= 6 gives no contact anywhere: the interior cannot supply half a turn.\n")

cat("\n=== 5. PLANTED FAILURES ===\n")
bad1 <- function(n) { a <- pi*(2-n)/4; if (a<0) NA_real_ else sin(a)^(1/n) }   # wrong exponent
cat(sprintf("  (a) exponent 1/n instead of 2/n gives r_c/r_h = %.4f at D=4, not 0.5\n", bad1(1)))
stopifnot(abs(bad1(1) - 0.5) > 0.1)
bad2 <- tryCatch(integrate(function(x) 1/(x^2*sqrt(pmax(0,1/x - 1))), 1e-9, 1)$value,
                 error=function(e) Inf)
cat(sprintf("  (b) r^2 sqrt|f| in the denominator: total turn = %s, not pi\n",
            ifelse(is.finite(bad2), sprintf("%.4f", bad2), "divergent")))
stopifnot(!is.finite(bad2) || abs(bad2 - pi) > 0.1)
bad3 <- (2/1)*(pi/2 - asin(0.5^(1/2)))*2                      # forgetting one of the two legs
cat(sprintf("  (c) counting one interior leg instead of two puts the D=4 budget at %.4f,\n", bad3/2))
cat(sprintf("      which is below pi, so it would report no contact in any dimension.\n"))
stopifnot(bad3/2 < pi)
cat("  Both slips move the answer.\n")

cat(sprintf("
=== flatly ===

  The interior of a D-dimensional Schwarzschild hole supplies exactly pi/(D-3) of
  turning to any causal curve crossing it, horizon to singularity, whatever the mass.
  The fold's transverse map is the antipode, which asks for pi. The connecting curve
  gets two interior legs, so it has 2 pi/(D-3) to spend against a bill of pi.

  That pays in D = 4 with a factor of two to spare, which is why the contact region
  is r <= r_h/2 rather than a point. It pays exactly in D = 5, where the region
  shrinks to the singularity. It does not pay in D = 6 or above, where the two sheets
  cannot touch anywhere at all.

  So this construction has something to say about a black hole only in four or five
  spacetime dimensions. That is not a modelling choice and nothing was tuned to make
  it come out: the antipodal map fixes the bill at pi because it is an involution of a
  sphere, and the geometry fixes the budget at pi/(D-3) because f has a simple zero
  and falls as r^{-(D-3)}. The two happen to match at D = 4.\n"))

cat("\n=== 9. REFEREE CHECK: is D = 5 'the singularity only', or empty? ===\n")
# At n = 2 the two-leg budget is 2*(2/2)[pi/2 - arcsin(r/r_h)] = pi - 2 arcsin(r/r_h),
# which is STRICTLY below pi for every r > 0 and reaches pi only in the limit r -> 0.
# r = 0 is not a point of the spacetime, so the D = 5 contact set is EMPTY in the manifold.
dphi_two <- function(x, n) (4/n)*(pi/2 - asin(x^(n/2)))   # two-leg budget
for (x in c(1e-1, 1e-3, 1e-6, 1e-9, 1e-12)) {
  v <- dphi_two(x, 2)
  cat(sprintf("   D=5, r/r_h = %8.1e : sweep = %.12f   (pi = %.12f)  short by %.3e\n",
              x, v, pi, pi - v))
  stopifnot(v < pi)
}
cat("   The sweep approaches pi from BELOW and never attains it at any r > 0.\n")
cat("   So D = 5 belongs in the 'nowhere' row and the theorem selects D = 4 alone.\n")

cat("\n=== 10. REFEREE CHECK: monotonicity OUTSIDE the horizon ===\n")
# The maximiser argument excludes curves that stray off X = 0, and those go OUTSIDE the
# horizon. The printed proof only covered r < r_h, which is where the excluded curves are
# NOT. Extend it: for r > r_h the numerator turns negative and f turns positive, so the
# ratio stays negative.
cat("      D    max of dlog(g^2)/dr over r in (1.001, 50) r_h    (must stay < 0)\n")
for (D in 4:9) { n <- D-3
  rr <- exp(seq(log(1.001), log(50), length.out=3000))
  v <- max(vapply(rr, function(x) dlog(x, 1, n), 0))
  cat(sprintf("   %4d                      %12.4f\n", D, v)); stopifnot(v < 0) }
cat("   Negative in the exterior too, in every D, so exterior excursions strictly lose\n")
cat("   and the maximiser really is X = 0 over ALL causal curves, not just interior ones.\n")

cat("\n=== 11. REFEREE CHECK: the two-leg count, argued directly ===\n")
cat("   A future-directed causal curve has dU >= 0 and dV >= 0, so once it leaves P it\n")
cat("   cannot return, and it cannot visit both R and L because every point of R is\n")
cat("   spacelike to every point of L. The itinerary is therefore P, then at most one\n")
cat("   exterior, then F: exactly one leg in each interior. The count does not depend on\n")
cat("   the maximiser and should not be derived from it.\n")
