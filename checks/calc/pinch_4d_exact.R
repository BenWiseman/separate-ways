#!/usr/bin/env Rscript
# pinch_4d_exact.R -- the four-dimensional contact condition, which A.15 left open.
#
# A.15 asks how deep inside a hole the two sheets can reach each other, and answers with
# a Kruskal ratio exp(3 sqrt3 pi/4) = 59.21 obtained by letting the connecting curve sweep
# the transverse angle at the photon sphere r = 3M. It says the result is exact only within
# the s-wave reduction and that the 4D version is not proved. This is the 4D version.
#
# SETUP. Eternal Schwarzschild, Kruskal:
#     ds^2 = F(r)(-dT^2 + dX^2) + r^2 dOmega^2,   F = 32 M^3 e^{-r/2M}/r,
#     T^2 - X^2 = UV = (1 - r/2M) e^{r/2M}.
# The fold is Theta = J . P_perp with J:(T,X) -> (-T,-X) and P_perp the map
# (theta,phi) -> (pi-theta, phi+pi). P_perp sends the unit vector n to -n, so it is the
# antipodal map of S^2 and the angular separation between x and Theta x is exactly pi for
# EVERY x. No generality is lost by working in the equatorial plane.
#
# THE QUESTION. Is x in the causal future of Theta x? Causality gives
#     r^2 dOmega^2 <= F (dT^2 - dX^2),   so   d(angle) <= g(r) dtau_2,  g = sqrt(F)/r,
# with dtau_2 the proper time of the 2D (T,X) part. So the maximum angular sweep of any
# causal curve joining the two points is the length of the longest timelike curve in the
# conformally rescaled 2D metric g^2(-dT^2+dX^2), and contact holds exactly when that
# maximum is at least pi.
#
# BOOST GAUGE. A boost sends (U,V) -> (U/L, LV) and fixes the angles, and it maps the pair
# (x, Theta x) to (x', Theta x'). So the answer depends only on the boost invariant UV, i.e.
# only on r. Put x at X = 0, so x = (T0, 0, pi/2, 0) and Theta x = (-T0, 0, pi/2, pi) with
# T0 = sqrt((1 - r/2M) e^{r/2M}).

M <- 1; rh <- 2*M
UV      <- function(r) (1 - r/(2*M))*exp(r/(2*M))     # = T^2 - X^2
Fk      <- function(r) 32*M^3*exp(-r/(2*M))/r
g       <- function(r) sqrt(Fk(r))/r
r_of_UV <- function(s) vapply(s, function(z)
             uniroot(function(r) UV(r)-z, c(1e-12, rh), tol=1e-14)$root, 0)

cat("=== 1. the optimal curve sits at X = 0, and this needs no numerics ===\n")
cat(sprintf("  g is decreasing in r:      dln(g^2)/dr = -1/2M - 3/r < 0 everywhere\n"))
cat(sprintf("  UV is decreasing in r:     d(UV)/dr = -(r/4M^2) e^{r/2M} < 0 on (0, 2M]\n"))
cat("  Moving off X = 0 lowers T^2-X^2, hence RAISES r, hence LOWERS g; and it costs\n")
cat("  sqrt(1-X'^2) <= 1 in the proper time. Both factors fall, so X = 0 is the unique\n")
cat("  maximiser. It is also the fixed set of X -> -X, so it is a geodesic.\n")
rr <- seq(0.02, 1.999, length.out=400)
stopifnot(all(diff(g(rr)) < 0), all(diff(UV(rr)) < 0))
cat("  Both monotonicities confirmed numerically on r in (0, 2M).\n")

cat("\n=== 2. the sweep along X = 0, in closed form ===\n")
# On X = 0, dtau_2 = dT and T = sqrt(UV(r)). The integrand collapses:
#   g(r) |dT/dr| = 1/sqrt(r(2M-r))     (every exponential cancels)
# so   Dphi(r) = 2 * int_r^{2M} dr'/sqrt(r'(2M-r')) = 2[2 asin sqrt(r'/2M)]_r^{2M}
#              = 2 pi - 4 asin sqrt(r/2M).
dphi_closed <- function(r) 2*pi - 4*asin(sqrt(r/(2*M)))

integrand <- function(r) { Tr <- sqrt(UV(r)); dTdr <- -(r*exp(r/(2*M)))/(8*M^2*Tr); g(r)*abs(dTdr) }
dphi_num <- function(r) 2*integrate(integrand, r, rh, rel.tol=1e-12)$value

cat("      r/2M     closed form   direct integral        diff\n")
for (x in c(0.05,0.2,0.4,0.5,0.6,0.8,0.95)) {
  r <- x*2*M
  cat(sprintf("     %5.2f      %9.6f       %9.6f   %9.2e\n",
              x, dphi_closed(r), dphi_num(r), abs(dphi_closed(r)-dphi_num(r))))
}
stopifnot(max(abs(vapply(c(0.05,0.2,0.4,0.5,0.6,0.8,0.95),
        function(x) dphi_closed(x*2*M)-dphi_num(x*2*M), 0))) < 1e-8)
cat("  Closed form and direct integration agree to better than 1e-8.\n")
cat(sprintf("  Check that the integrand really does collapse: g(r)|dT/dr| vs 1/sqrt(r(2M-r))\n"))
tst <- vapply(seq(0.05,1.95,by=0.05), function(r) integrand(r) - 1/sqrt(r*(2*M-r)), 0)
cat(sprintf("     max |difference| over r in [0.05, 1.95] = %.3e\n", max(abs(tst))))
stopifnot(max(abs(tst)) < 1e-10)

cat("\n=== 3. the contact radius ===\n")
# contact <=> Dphi >= pi <=> 4 asin sqrt(r/2M) <= pi <=> r <= 2M sin^2(pi/4) = M.
r_c <- uniroot(function(r) dphi_closed(r) - pi, c(1e-9, rh), tol=1e-14)$root
cat(sprintf("  solving Dphi(r) = pi gives r_c = %.12f M   (exactly M = %.1f)\n", r_c/M, M))
stopifnot(abs(r_c - M) < 1e-10)
cat(sprintf("  2M sin^2(pi/4) = %.12f M.  The contact region is the inner HALF of the hole.\n",
            2*M*sin(pi/4)^2/M))
cat(sprintf("  In Kruskal terms contact needs UV >= UV(M) = sqrt(e)/2 = %.6f, with the\n", UV(M)))
cat(sprintf("  singularity at UV = 1. So the s-wave claim 'nonzero whenever UV > 0' is too\n"))
cat(sprintf("  generous: the true band is %.4f <= UV <= 1.\n", UV(M)))
stopifnot(abs(UV(M) - sqrt(exp(1))/2) < 1e-12)

cat("\n=== 4. the photon sphere is not reachable, which is what breaks the old argument ===\n")
# Along a future-directed causal curve U and V both increase, so between Theta x = (-U0,-V0)
# and x = (U0,V0) the curve is confined to |U| <= U0 and |V| <= V0, hence |UV| <= U0 V0.
UV_ext <- function(r) (r/(2*M) - 1)*exp(r/(2*M))       # |UV| outside the horizon
cat(sprintf("  |UV| at the photon sphere r = 3M: %.4f\n", UV_ext(3*M)))
cat("      r/2M   max |UV| on the curve   r_max reachable\n")
for (x in c(0.1,0.3,0.5,0.9,0.999)) {
  r <- x*2*M; cap <- UV(r)
  rmax <- uniroot(function(rr) UV_ext(rr)-cap, c(2*M, 10*M), tol=1e-12)$root
  cat(sprintf("     %5.3f            %9.4f          %7.4f M\n", x, cap, rmax/M))
}
cat(sprintf("  The cap is UV(r) <= UV(0) = 1, so no connecting curve ever passes r = %.4f M,\n",
            uniroot(function(rr) UV_ext(rr)-1, c(2*M,10*M), tol=1e-12)$root/M))
cat("  well short of 3M. The optimal curve in fact never leaves the closed interior: it\n")
cat("  touches r = 2M once, at the bifurcation surface, and that is its outermost point.\n")
stopifnot(uniroot(function(rr) UV_ext(rr)-1, c(2*M,10*M))$root < 3*M)

cat("\n=== 5. the optimum is attained, not approached ===\n")
# The bound is saturated by a null curve with dX = 0 and dphi = g dT. Check ds^2 = 0.
rt <- 0.7*2*M; Tt <- sqrt(UV(rt))
ds2 <- Fk(rt)*(-1) + rt^2*g(rt)^2
cat(sprintf("  ds^2 for dT=1, dX=0, dphi=g(r): %.3e  (null)\n", ds2))
stopifnot(abs(ds2) < 1e-12)
cat("  So the maximiser is a genuine causal curve: a null spiral sitting at X = 0.\n")

cat("\n=== 6. PLANTED FAILURE ===\n")
# Two deliberate slips, each of which must move the answer.
# (a) mis-power the conformal factor, g^2 = F/r^2 written as F/r^3. A plausible typo.
g_bad <- function(r) sqrt(Fk(r))/r^1.5
int_bad <- function(r) { Tr <- sqrt(UV(r)); g_bad(r)*abs(-(r*exp(r/(2*M)))/(8*M^2*Tr)) }
dphi_bad <- function(r) 2*integrate(int_bad, r, rh, rel.tol=1e-10)$value
lo <- 0.25*M; hi <- 0.999*rh          # bad integrand is log-divergent at r -> 0
stopifnot(dphi_bad(lo) > pi, dphi_bad(hi) < pi)
r_bad <- uniroot(function(r) dphi_bad(r) - pi, c(lo, hi), tol=1e-12)$root
cat(sprintf("  (a) g^2 = F/r^3 instead of F/r^2 gives r_c = %.4f M, not %.4f M\n", r_bad/M, r_c/M))
stopifnot(abs(r_bad - r_c) > 0.02*M)

# (b) the old argument's assumption, that the sweep happens at the photon sphere. If the
# curve could sit at r = 3M the required coordinate time would be pi/Omega(3M) = 3 sqrt3 pi M
# and the answer would be a V ratio, not a radius. Check the two are not the same statement.
Omega <- function(r) sqrt(1 - 2*M/r)/r
r_opt <- optimize(Omega, c(2.0001*M, 20*M), maximum=TRUE)$maximum
cat(sprintf("  (b) Omega is maximised at r = %.4f M, and exp(kappa pi/Omega) = %.2f,\n",
            r_opt/M, exp((1/(4*M))*pi/Omega(3*M))))
cat(sprintf("      but the curve cannot reach %.4f M, so that optimum is unavailable.\n", r_opt/M))
stopifnot(abs(r_opt - 3*M) < 1e-3, abs(exp((1/(4*M))*pi/Omega(3*M)) - exp(3*sqrt(3)*pi/4)) < 1e-6)
cat("  Both slips move the answer, so neither check is blind.\n")

cat(sprintf("
=== flatly ===

  The four-dimensional contact condition is r <= M: the two sheets can reach each
  other only in the inner half of the hole, and nowhere else in the universe. It is
  exact, it is boost invariant, and the mass cancels because it is a fraction of the
  horizon radius, not a length.

  It also supersedes A.15's 59.21. That number came from letting the connecting
  curve sweep the transverse angle at the photon sphere. Along a future-directed
  causal curve U and V both increase, so a curve joining Theta x to x is confined to
  |UV| <= UV(r_x) <= 1, while |UV| at r = 3M is %.3f. The photon sphere is not
  reachable by any curve that does the job, and the optimal curve never leaves the
  closed interior at all: it touches r = 2M once, at the bifurcation surface.

  What survives from A.15 is the shape of the claim and its motivation, which were
  right: contact needs the transverse angle covered as well as the fall, it does not
  open at the horizon, and how deep it opens is a pure number free of the mass. The
  number is 1/2, not 59.21.\n", UV_ext(3*M)))

cat("\n=== 7. the contact window, which r <= M shortens ===\n")
# Maximal proper time from radius r to the singularity, for a radial infaller:
#   tau_max(r) = int_0^r dr'/sqrt(2M/r' - 1) = 2M[ asin sqrt(r/2M) - sqrt(r/2M) sqrt(1-r/2M) ]
tau_max <- function(r) { u <- sqrt(r/(2*M)); 2*M*(asin(u) - u*sqrt(1-u^2)) }
tau_num <- function(r) integrate(function(rr) sqrt(rr/(2*M-rr)), 0, r, rel.tol=1e-12)$value
cat(sprintf("  closed form vs numeric at r=2M: %.10f  %.10f\n", tau_max(2*M), tau_num(2*M)))
stopifnot(abs(tau_max(2*M) - pi*M) < 1e-10, abs(tau_max(2*M) - tau_num(2*M)) < 1e-9,
          abs(tau_max(M) - tau_num(M)) < 1e-9)
cat(sprintf("  horizon to singularity : tau = pi M          = %.4f M\n", tau_max(2*M)/M))
cat(sprintf("  r = M to singularity   : tau = (pi/2 - 1) M  = %.4f M   (%.1f%% of it)\n",
            tau_max(M)/M, 100*tau_max(M)/tau_max(2*M)))
stopifnot(abs(tau_max(M) - M*(pi/2 - 1)) < 1e-12)

G <- 6.67430e-11; cc <- 2.99792458e8; Msun <- 1.98892e30
geo <- function(Msol) G*Msol*Msun/cc^3          # one mass in seconds
cat("\n            object          M/Msun      pi M        (pi/2-1) M\n")
for (nm in list(c("stellar remnant","3"), c("Sgr A*","4.297e6"), c("M87*","6.5e9"))) {
  Ms <- as.numeric(nm[2]); gs <- geo(Ms)
  fmt <- function(x) if (x < 1) sprintf("%8.1f us", x*1e6) else
                     if (x < 3600) sprintf("%8.1f s ", x) else sprintf("%8.2f h ", x/3600)
  cat(sprintf("  %16s  %9.3g  %s  %s\n", nm[1], Ms, fmt(pi*gs), fmt((pi/2-1)*gs)))
}
cat("\n  The paper quoted pi M as 50 us at 3 Msun, which is right (46.4 us), and as
  four hours at 6.5e9 Msun, which is not: it is 27.9 hours. And with contact confined
  to r <= M the window that matters is (pi/2 - 1) M, eighteen per cent of it.\n")

cat("\n=== 8. is contact mediated through B, or only the optimal curve? ===\n")
# The maximiser sits at X = 0 and so passes through the bifurcation surface. But for r
# strictly less than M there is slack, so curves that avoid B should also work. Test an
# explicit family with X(T) = a sin(pi T / (2 T0)) ... which vanishes at both ends only for
# the odd extension; use X(T) = a cos(pi T /(2 T0)) instead, which is a at T = 0 and 0 at
# the endpoints, so it misses B by a.
sweep <- function(r, a) {
  T0 <- sqrt(UV(r))
  Xf  <- function(t) a*cos(pi*t/(2*T0))
  Xd  <- function(t) -a*(pi/(2*T0))*sin(pi*t/(2*T0))
  f <- function(t) {
    s <- t^2 - Xf(t)^2
    if (any(s <= 0)) return(rep(0, length(t)))
    rr <- r_of_UV(s); xd <- Xd(t)
    ifelse(abs(xd) >= 1, 0, g(rr)*sqrt(pmax(0, 1 - xd^2)))
  }
  2*integrate(Vectorize(function(t) f(t)), 1e-9, T0, rel.tol=1e-8, subdivisions=400L)$value
}
cat("      r/M     a = 0 (through B)   a = 0.05   a = 0.15   still >= pi?\n")
for (rx in c(0.2, 0.5, 0.8, 0.99, 1.0)) {
  r <- rx*M
  v0 <- dphi_closed(r); v1 <- sweep(r, 0.05); v2 <- sweep(r, 0.15)
  cat(sprintf("    %5.2f          %8.4f     %8.4f   %8.4f   %s\n", rx, v0, v1, v2,
              ifelse(v2 >= pi, "yes", ifelse(v1 >= pi, "at a=0.05 only", "no"))))
}
cat("
  So contact is NOT mediated through B in general. The maximising curve passes through
  the bifurcation surface, and at r = M it is the only curve that reaches pi, but deeper
  in there is slack and routes that miss B by a finite amount work too. The correct
  statement is about the boundary of the region, not about the region.\n")

cat("\n=== 12. how the light cone is crossed at r = M ===\n")
# Dphi_max(r) - pi is the angular deficit. Its zero at r = M is TRANSVERSAL, and the
# slope is exactly -2/M:
#   d/dr [2 pi - 4 asin sqrt(r/2M)] = -4 * (1/sqrt(1 - r/2M)) * (1/(2 sqrt(r/2M))) * (1/2M)
# which at r = M is -4 / (2M * 2 * (1/sqrt2)(1/sqrt2)) = -2/M.
slope_cf <- -2/M
slope_nd <- { h <- 1e-6; (dphi_closed(M+h) - dphi_closed(M-h))/(2*h) }
cat(sprintf("  closed form  d(Dphi)/dr at r=M = %.10f\n", slope_cf))
cat(sprintf("  central difference            = %.10f\n", slope_nd))
stopifnot(abs(slope_cf - slope_nd) < 1e-6, abs(slope_cf + 2/M) < 1e-12)
cat("  So the deficit has a SIMPLE zero: the identified pair crosses its own light cone\n")
cat("  transversally, not tangentially. That is the input a self-censoring argument would\n")
cat("  need, since an image sum on a transversal crossing diverges as a power.\n")

# planted: a tangential crossing must give slope zero, so the check can tell them apart
tang <- { h <- 1e-4; f <- function(r) (r - M)^2; (f(M+h)-f(M-h))/(2*h) }
cat(sprintf("  CHECK (planted): a tangential crossing gives slope %.3e, not %.1f\n", tang, slope_cf))
stopifnot(abs(tang) < 1e-8, abs(slope_cf) > 1)
cat("  The two are distinguishable, so the statement is not vacuous.\n")
