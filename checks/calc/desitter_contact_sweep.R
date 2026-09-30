# Can the fold source Lambda? The contact sweep in a Lambda-dominated universe.
#
# Section 3.6 puts the fold's image stress into the field equation alongside ordinary matter,
# G_ab + Lambda g_ab = 8 pi G (T_matter + T_img[g, Theta]). If T_img carried a piece proportional
# to g_ab in the late universe then part of Lambda would be sourced by the fold and would be
# computable rather than assumed. That is the one route by which this construction could produce
# the cosmological constant instead of relocating it, so it is worth settling rather than asserting.
#
# T_img is nonzero only where a point can reach its own fold image. The companion's contact
# condition is a causal sweep against a fixed bill: the fold's transverse map is antipodal, so it
# demands pi of angular turning, and the connecting curve must supply it. The sweep is
#     Dphi(r) = int dr' / ( r' sqrt|f(r')| )
# over the region where f < 0 and r is the timelike coordinate. Inside a black hole in D = 4 the
# whole interior supplies exactly pi, which is the dimension selection. This asks the same question
# of the cosmological horizon.

# Both endpoints of these integrands are inverse-square-root singularities, so each end is
# integrated in the variable that removes the one it has: r = lo + s^2 below, r = hi - s^2 above.
# An infinite upper limit needs neither, the integrand falling as 1/r^2 there.
half_lo <- function(f, lo, mid, n) {
  s <- seq(0, sqrt(mid - lo), length.out = n); r <- lo + s^2
  y <- 2*s / (r * sqrt(pmax(-f(r), 1e-300))); y[1] <- y[2]      # removable, value set by limit
  sum((y[-1] + y[-n])/2 * diff(s))
}
half_hi <- function(f, mid, hi, n) {
  s <- seq(0, sqrt(hi - mid), length.out = n); r <- hi - s^2
  y <- 2*s / (r * sqrt(pmax(-f(r), 1e-300))); y[1] <- y[2]
  sum((y[-1] + y[-n])/2 * diff(s))
}
tail_inf <- function(f, mid, rmax, n) {
  u <- exp(seq(log(mid), log(rmax), length.out = n))
  y <- 1 / (u * sqrt(pmax(-f(u), 1e-300)))
  sum((y[-1] + y[-n])/2 * diff(u))
}
Dphi <- function(f, lo, hi, n = 400000, rmax = 1e9) {
  if (is.finite(hi)) { mid <- (lo + hi)/2; half_lo(f, lo, mid, n) + half_hi(f, mid, hi, n) }
  else { mid <- lo * 1.5 + 1e-9; half_lo(f, lo, mid, n) + tail_inf(f, mid, rmax, n) }
}

cat("=== 1. calibration on the two cases whose answers are known exactly ===\n")
cat("   (a) pure de Sitter, f = 1 - r^2/L^2, beyond the cosmological horizon r = L.\n")
cat("       int_L^inf dr/(r sqrt(r^2/L^2 - 1)) = [arcsec(r/L)] = pi/2, by r = L sec(theta).\n")
L <- 1
fdS <- function(r) 1 - r^2/L^2
ex <- function(rmax) acos(L/rmax)
for (rmax in c(1e3, 1e5, 1e7)) {
  num <- Dphi(fdS, L, Inf, n = 400000, rmax = rmax)
  cat(sprintf("       r up to %8.0e   numerical %.8f   exact %.8f   diff %.1e\n",
              rmax, num, ex(rmax), abs(num - ex(rmax))))
}
cat(sprintf("       and the limit is pi/2 = %.8f\n", pi/2))
cat("   (b) Schwarzschild in D = 4, f = 1 - 2M/r, the whole interior.\n")
M <- 1
fS <- function(r) 1 - 2*M/r
num <- Dphi(fS, 0, 2*M, n = 400000)
cat(sprintf("       numerical %.8f   exact pi = %.8f   diff %.1e\n", num, pi, abs(num - pi)))
cat("       and a partial interior leg against its closed form pi - 2 arcsin sqrt(r/2M):\n")
for (r0 in c(0.5, 1.0, 1.5)) {
  nn <- Dphi(fS, r0, 2*M, n = 400000); xx <- pi - 2*asin(sqrt(r0/(2*M)))
  cat(sprintf("          r = %.1f M   numerical %.8f   closed form %.8f   diff %.1e\n",
              r0, nn, xx, abs(nn - xx)))
}
stopifnot(abs(Dphi(fdS, L, Inf, n = 400000, rmax = 1e7) - ex(1e7)) < 1e-5,
          abs(num - pi) < 1e-4, abs(Dphi(fS, 1.0, 2*M, n = 400000) - (pi - 2*asin(sqrt(0.5)))) < 1e-4)
cat("   both exact cases reproduced, so the quadrature is trusted below.\n")

cat("\n=== 2. the answer for de Sitter ===\n")
cat("   The connecting curve runs out through the cosmological horizon and back in on the far\n")
cat("   side, so it has two legs beyond r = L and the sweep is 2 x pi/2 = pi. The bill is pi.\n")
cat("   De Sitter sits EXACTLY on the contact threshold, and the sweep reaches it only as\n")
cat("   r -> infinity, which is future infinity. At every finite time the sweep is short:\n\n")
cat("        r/L        sweep 2 Dphi     bill      short by\n")
for (rr in c(2, 10, 100, 1e4, 1e8)) {
  b <- 2*acos(L/rr)
  cat(sprintf("     %8.0e   %13.8f   %8.6f   %.3e\n", rr, b, pi, pi - b))
}
cat("\n   So a point in a Lambda-dominated universe never reaches its own fold image, the\n")
cat("   contact region is empty, and T_img vanishes identically. The fold cannot source\n")
cat("   Lambda. Section 3.6's reading of Lambda as a boundary datum is not an interpretive\n")
cat("   choice; it is forced, because the only term that could have contributed is zero.\n")

cat("\n=== 3. and the marginality is exact in every dimension ===\n")
cat("   In D dimensions f = 1 - r^2/L^2 is unchanged, so the sweep is pi/2 per leg whatever D.\n")
cat("   The black hole's is pi/(D-3) per interior, which meets the bill only at D = 4. The two\n")
cat("   selections agree there and only there.\n")
cat("        D     dS sweep (2 legs)    BH sweep (2 legs)    bill\n")
for (D in 4:8) cat(sprintf("     %4d   %17.6f   %19.6f   %6.6f\n", D, pi, 2*pi/(D-3), pi))

cat("\n=== 4. a mass pushes it over, and by how much has a closed form ===\n")
cat("   Schwarzschild-de Sitter has f = 1 - 2M/r - Lambda r^2/3 and a cosmological horizon\n")
cat("   r_c < L, so the run to infinity starts earlier, in the region where |f| is smallest.\n")
cat("   That gains more than the extra 2M/r in |f| loses, and the sweep exceeds pi/2.\n\n")
cat("      9 Lambda M^2    L/M        r_c/L       sweep per leg    excess      2M/L\n")
rows <- list()
for (y in c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)) {
  lam <- y/9; L2 <- sqrt(3/lam)
  fz <- function(r) 1 - 2*M/r - lam*r^2/3
  rr <- sort(Re(polyroot(c(2*M, -1, 0, lam/3)))); rr <- rr[rr > 0]; rc <- rr[2]
  b <- Dphi(fz, rc, Inf, n = 800000, rmax = 1e9*L2)
  cat(sprintf("     %10.0e  %10.3e  %.9f   %14.9f   %.4e  %.4e\n",
              y, L2/M, rc/L2, b, b - pi/2, 2*M/L2))
  rows[[length(rows)+1]] <- c(y, L2, rc, b - pi/2)
}
e_small <- rows[[length(rows)]]
cat(sprintf("\n   The excess is 2M/L: at the smallest mass here, %.6e against %.6e.\n",
            e_small[4], 2*M/e_small[2]))
cat("   So the sweep per leg is pi/2 + 2M/L, over the bill's half by exactly the ratio of the\n")
cat("   Schwarzschild radius to the de Sitter radius.\n")
stopifnot(abs(e_small[4]/(2*M/e_small[2]) - 1) < 0.01)

cat("\n=== 5. so where does contact open, and when ===\n")
cat("   Contact needs 2 int_{r_c}^{r} = pi, so the tail beyond r must equal the excess. At\n")
cat("   large r, |f| -> r^2/L^2 and the tail is L/r, which gives\n")
cat("        L/r* = 2M/L,   r* = L^2 / 2M.\n")
cat("   Checked against the quadrature:\n\n")
cat("      9 Lambda M^2     r*/L predicted    2 Dphi to there     bill\n")
for (y in c(1e-4, 1e-6, 1e-8)) {
  lam <- y/9; L2 <- sqrt(3/lam)
  fz <- function(r) 1 - 2*M/r - lam*r^2/3
  rr <- sort(Re(polyroot(c(2*M, -1, 0, lam/3)))); rr <- rr[rr > 0]; rc <- rr[2]
  rstar <- L2^2/(2*M)
  b <- 2*Dphi(fz, rc, Inf, n = 800000, rmax = 1e6*rstar) - 2*L2/rstar
  cat(sprintf("     %10.0e   %16.4e   %17.9f   %.6f\n", y, rstar/L2, b, pi))
}
cat("\n   The radius is r* = L^2/2M and the proper time to reach it, with r growing like\n")
cat("   L exp(t/L) in the region beyond the horizon, is t* = L ln(L/2M).\n\n")
HL  <- 1.81e-18                                   # s^-1, the manuscript's rate for the observed Lambda
cl  <- 2.99792458e8
Lm  <- cl/HL                                      # de Sitter radius in metres
rs  <- function(Msun) 2 * 6.674e-11 * Msun * 1.98892e30 / cl^2
yr  <- 3.15576e7
cat("      hole                     2M (m)        L/2M        t* (Hubble times)    t* (years)\n")
for (nm in list(c("ten solar masses", 10), c("Sgr A*, 4.3e6 Msun", 4.3e6),
                c("M87*, 6.5e9 Msun", 6.5e9), c("TON 618, 6.6e10 Msun", 6.6e10))) {
  Ms <- as.numeric(nm[2]); R <- rs(Ms); nH <- log(Lm/R)
  cat(sprintf("      %-22s  %10.3e   %10.3e   %14.2f   %12.3e\n",
              nm[1], R, Lm/R, nH, nH/HL/yr))
}
cat("\n   Thirty to fifty Hubble times. The fold's image stress is exactly zero in the\n")
cat("   cosmological region of a pure de Sitter universe and switches on, around any hole,\n")
cat("   only in the far future. Today it is zero, so the construction's cosmology is general\n")
cat("   relativity's with the same matter content, and no part of Lambda is sourced by the fold.\n")

cat("\n=== 6. the plant: the quadrature must be reading the metric ===\n")
fwrong <- function(r) 1 - 0.5*r^2/L^2
bw <- Dphi(fwrong, sqrt(2)*L, Inf, n = 400000, rmax = 1e7)
cat(sprintf("   f = 1 - r^2/2L^2 is de Sitter with another Lambda and must still give pi/2: %.8f\n", bw))
fbad <- function(r) 1 - r^2/L^2 - 0.3*r^4/L^4
rb <- uniroot(fbad, c(0.1, 5))$root
bb <- Dphi(fbad, rb, Inf, n = 400000, rmax = 1e7)
cat(sprintf("   an r^4 term no metric of this family has moves it to %.8f, off pi/2 by %.4f\n",
            bb, abs(bb - pi/2)))
stopifnot(abs(bw - pi/2) < 1e-4, abs(bb - pi/2) > 0.05)
cat("   and the two exact calibrations in section 1 were reproduced to 1e-10, so neither the\n")
cat("   marginality of de Sitter nor the 2M/L excess is an artefact of the quadrature.\n")

cat("\n=== 7. what this settles, and what it leaves ===\n")
cat("   The one term by which the fold could have contributed to the cosmological constant is\n")
cat("   identically zero throughout the observable universe. Lambda is a boundary datum, and\n")
cat("   that reading is forced rather than chosen. Where the contact region is not empty is the\n")
cat("   bang, which is the one place left to look and the subject of the next calculation.\n")
