# horizon_universality.R -- how far the A.6 horizon ratio carries, App. A.10.
# Base R only. Units M = 1 for Kerr.
#
# THE QUESTION. A.6 obtains G_alpha/G_J -> tan^2(gamma/2) at the de Sitter horizon
# for the conformally coupled scalar in Bunch-Davies. Does that value hold at any
# bifurcate Killing horizon, or is it particular to that spacetime and that state?
#
# THE SETUP. Let B be the bifurcation surface. The wedge reflection J fixes B
# pointwise; P_perp is a free involutive isometry of B; alpha = J o P_perp. Taking
# both points to B, the boost separation drops out of both kernels and
#
#     G_alpha / G_J  ->  W_B(x, P_perp y) / W_B(x, y),
#
# the state's own correlator restricted to B, at the reflected separation over the
# direct one. That much needs no metric on B and no choice of state.
#
# THE ANSWER, in two parts.
#   (1) On a ROUND B with an isotropic state the ratio is G(pi-g)/G(g). Three of its
#       features are common to every such state and three only: it vanishes like g^2,
#       it equals exactly 1 at g = pi/2, and it has a (pi-g)^-2 pole. The VALUE at
#       intermediate angles is not universal -- a 1/arclength^2 profile gives 4.00 at
#       120 degrees where the chordal profile gives 3.00.
#   (2) tan^2(g/2) is the chordal profile's answer. It is exact for the de Sitter
#       conformal scalar in Bunch-Davies (A.6) and is NOT established for the
#       Hartle-Hawking state on Schwarzschild. An earlier transverse-flat treatment
#       that obtained it there used a flat chord on a curved sphere, which is valid
#       only at small angle, and is withdrawn as a derivation. It remains consistent
#       with all three universal features.
#
# KERR. P_perp (theta -> pi-theta, phi -> phi+pi) is still an exact free involutive
# isometry of the bifurcation surface, so the structure survives. But B is squashed
# by 1 + a^2/r_+^2, so isotropy fails and with it both the unit value at pi/2 and any
# closed form. The squashing vanishes identically at a = 0.

cat("=== 1. ROUND B: the ratio is G(pi-g)/G(g) for any isotropic correlator\n\n")

# Three profiles sharing the Hadamard 1/d^2 short-distance behaviour at BOTH
# coincidences but differing at finite separation.
G_chord <- function(g) 1/(2*(1-cos(g)))                    # de Sitter conformal form
G_geod  <- function(g) 1/(g^2)                             # 1/arclength^2
G_asym  <- function(g) (1 + 0.3*g/pi)/(2*(1-cos(g)))       # same poles, asymmetric middle

gs  <- c(1,5,20,45,60,90,120,150,175,179)*pi/180
tab <- data.frame(gamma_deg = round(gs*180/pi),
                  chordal   = G_chord(pi-gs)/G_chord(gs),
                  arclength = G_geod(pi-gs)/G_geod(gs),
                  asymmetric= G_asym(pi-gs)/G_asym(gs),
                  tan2_half = tan(gs/2)^2)
print(format(tab, digits = 6), row.names = FALSE)

cat("\n=== 2. What is common to all three, and what is not\n\n")
for (nm in c("chord","geod","asym")) {
  G <- get(paste0("G_", nm))
  g1 <- 1e-3; g2 <- 2e-3
  s0 <- log((G(pi-g2)/G(g2)) / (G(pi-g1)/G(g1))) / log(g2/g1)
  sp <- log((G(g2)/G(pi-g2)) / (G(g1)/G(pi-g1))) / log(g2/g1)
  cat(sprintf("  %-6s  value at pi/2 = %.12f   small-g slope = %+.4f   near-pi slope = %+.4f\n",
              nm, G(pi/2)/G(pi/2), s0, sp))
}
cat("\n  universal: zero like g^2, exactly 1 at g = pi/2, pole like (pi-g)^-2.\n")
cat(sprintf("  not universal: at 120 deg the three give %.2f, %.2f, %.2f.\n",
            G_chord(pi/3)/G_chord(2*pi/3), G_geod(pi/3)/G_geod(2*pi/3),
            G_asym(pi/3)/G_asym(2*pi/3)))

cat("\n=== 3. Kerr: P_perp is still a free involutive isometry of B\n\n")
rp   <- function(a) 1 + sqrt(1 - a^2)
g_tt <- function(th, a) rp(a)^2 + a^2*cos(th)^2
g_pp <- function(th, a) (rp(a)^2 + a^2)^2*sin(th)^2 / (rp(a)^2 + a^2*cos(th)^2)
ths  <- c(0.1, 0.7, 1.2, pi/2, 2.0, 2.7, 3.0)
for (a in c(0.5, 0.9, 0.99))
  cat(sprintf("  a=%.2f   max|g_thth(pi-th) - g_thth(th)| = %.2e   same for g_phph = %.2e\n",
              a, max(abs(g_tt(pi-ths,a) - g_tt(ths,a))), max(abs(g_pp(pi-ths,a) - g_pp(ths,a)))))
cat("  free: the theta = pi/2 circle is moved by phi -> phi+pi. Involutive: phi -> phi+2pi.\n")

cat("\n=== 4. Kerr: but B is not round, so isotropy fails\n\n")
for (a in c(0, 0.3, 0.5, 0.7, 0.9, 0.99))
  cat(sprintf("  a/M=%4.2f   r_+ = %.5f   g_thth(pole)/g_thth(equator) = 1 + a^2/r_+^2 = %.4f\n",
              a, rp(a), 1 + a^2/rp(a)^2))
cat("\n  The departure from the round-sphere result is switched on by spin alone and\n")
cat("  vanishes identically at a = 0, the same signature the l -> l+/-2 mixing carries.\n")
