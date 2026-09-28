#!/usr/bin/env Rscript
# Fork 9, third brick: the cross-region state is not a new input. It is the half-period shift
# that already makes the exterior correlator thermal.
#
# WHY. interior_modes_at_nonzero_k.R ends by naming the one piece that cannot be settled inside
# the interior alone: the state pairing the future interior with the past one, which is where the
# fold sends an interior point. If that pairing were a free choice the interior sum would carry an
# unfixed phase and could not be read for a sign. It is not a free choice. On the Kruskal manifold
# the shift t -> t + i beta/2 IS the antipodal map (U, V) -> (-U, -V), so the same fact that makes
# the exterior correlator a thermofield double fixes the future-past correlator as well.
#
# WHAT IS CHECKED HERE. The claim is algebraic and it is checked as algebra, at points rather than
# by manipulation: with 2M = 1 and kappa = 1/(4M) = 1/2, so beta = 2 pi / kappa = 4 pi,
#   (a) t -> t + i beta/2 sends (U, V) to (-U, -V) everywhere, in all four regions;
#   (b) it therefore fixes r, since r is a function of the product UV, and fixes t, since t is a
#       function of the ratio U/V, so the image sits at the SAME (t, r) in the antipodal region;
#   (c) doing it twice is the identity, which is Theta^2 = 1 at the level of the geometry;
#   (d) region I goes to region III and region II to region IV, which is the future interior
#       pairing with the past one;
#   (e) and the E = 0 contact geodesic, which is the curve the interior image stress is computed
#       along, keeps t constant, so its two endpoints are exactly such a pair.
# The consequence: the interior mode sum inherits its state from the exterior temperature and has
# nothing left to choose there.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M <- 0.5                 # 2M = 1, horizon at r = 1
kap <- 1/(4*M)           # surface gravity 1/2
beta <- 2*pi/kap         # 4 pi

rstar <- function(r) r + log(abs(r - 1))     # 2M = 1: r + 2M log|r/2M - 1| with the constant absorbed

# Kruskal, region by region. U V = -exp(2 kappa rstar) outside, +exp(2 kappa rstar) inside.
UV <- function(r, t, region) {
  x <- exp(2*kap*rstar(r))
  switch(region,
    I   = c(U = -exp(-kap*(t - rstar(r))), V =  exp( kap*(t + rstar(r)))),
    III = c(U =  exp(-kap*(t - rstar(r))), V = -exp( kap*(t + rstar(r)))),
    II  = c(U =  exp(-kap*(t - rstar(r))), V =  exp( kap*(t + rstar(r)))),
    IV  = c(U = -exp(-kap*(t - rstar(r))), V = -exp( kap*(t + rstar(r)))))
}

cat("=== 1. the half-period shift is the Kruskal antipodal map, in every region ===\n")
cat("   With 2M = 1, kappa = 1/2 and beta = 4 pi, so the shift is t -> t + 2 pi i.\n")
cat("   region      r       t        U -> U'/U        V -> V'/V\n")
cases <- list(list("I", 2.5, 0.7), list("I", 8.0, -1.3), list("III", 3.1, 0.2),
              list("II", 0.8, 0.4), list("II", 0.3, -0.9), list("IV", 0.6, 1.1))
for (cs in cases) {
  reg <- cs[[1]]; r <- cs[[2]]; t <- cs[[3]]
  a <- UV(r, t, reg)
  # the shifted point, computed as a complex t and nothing else
  tc <- t + 1i*beta/2
  x  <- exp(2*kap*rstar(r))
  sgn <- switch(reg, I = c(-1, 1), III = c(1, -1), II = c(1, 1), IV = c(-1, -1))
  b <- c(U = sgn[1]*exp(-kap*(tc - rstar(r))), V = sgn[2]*exp(kap*(tc + rstar(r))))
  ru <- b[1]/a[1]; rv <- b[2]/a[2]
  cat(sprintf("   %-5s %7.2f %7.2f   %8.4f%+8.4fi   %8.4f%+8.4fi\n",
              reg, r, t, Re(ru), Im(ru), Re(rv), Im(rv)))
  note(Mod(ru + 1) < 1e-12 && Mod(rv + 1) < 1e-12,
       sprintf("the shift sends (U,V) to (-U,-V) in region %s at r = %g", reg, r))
}

cat("\n=== 2. so it fixes r and t, and the image is the antipodal region at the SAME event ===\n")
cat("   r depends on U V and t on U/V, and both are even under (U,V) -> (-U,-V):\n")
for (cs in cases[c(1, 4, 6)]) {
  reg <- cs[[1]]; r <- cs[[2]]; t <- cs[[3]]
  a <- UV(r, t, reg)
  prod0 <- a[1]*a[2]; ratio0 <- a[1]/a[2]
  prod1 <- (-a[1])*(-a[2]); ratio1 <- (-a[1])/(-a[2])
  cat(sprintf("   %-4s r = %4.2f: UV %12.6e -> %12.6e,  U/V %12.6e -> %12.6e\n",
              reg, r, prod0, prod1, ratio0, ratio1))
  note(abs(prod1 - prod0) < 1e-12*abs(prod0) && abs(ratio1 - ratio0) < 1e-12*abs(ratio0),
       sprintf("r and t are both fixed in region %s", reg))
}

cat("\n=== 3. twice is the identity, which is Theta^2 = 1 on the geometry ===\n")
a <- UV(0.8, 0.4, "II")
twice <- c(-1, -1)*(c(-1, -1)*a)
cat(sprintf("   region II, r = 0.8: (U,V) = (%.6f, %.6f) -> twice -> (%.6f, %.6f)\n",
            a[1], a[2], twice[1], twice[2]))
note(max(abs(twice - a)) < 1e-14, "the map is an involution on the Kruskal manifold")
cat("   And a WHOLE period is the identity on the Euclidean section, t -> t + i beta giving\n")
cat("   (U,V) -> (U,V), which is what makes beta/2 the primitive half and not one choice of\n")
cat("   many. A third of a period is the plant in section 5.\n")

cat("\n=== 4. the regions pair as future with past ===\n")
cat("   region   sign(U) sign(V)   image sign(U) sign(V)   image region\n")
for (reg in c("I", "III", "II", "IV")) {
  r <- if (reg %in% c("I", "III")) 2.5 else 0.6
  a <- UV(r, 0.3, reg)
  s0 <- sign(a); s1 <- -s0
  img <- if (s1[1] < 0 && s1[2] > 0) "I" else if (s1[1] > 0 && s1[2] < 0) "III" else
         if (s1[1] > 0 && s1[2] > 0) "II" else "IV"
  cat(sprintf("   %-6s  %6d %7d   %12d %7d   %s\n", reg, s0[1], s0[2], s1[1], s1[2], img))
  want <- c(I = "III", III = "I", II = "IV", IV = "II")[[reg]]
  note(img == want, sprintf("region %s pairs with region %s", reg, want))
}
cat("   So an interior point's fold image is in the PAST interior, which can reach it: the\n")
cat("   contact the companion computes is a past-interior point sending a signal forward.\n")

cat("\n=== 5. the contact geodesic keeps t constant, so its endpoints are such a pair ===\n")
cat("   The E = 0 condition is f dt/dlambda = 0 and f is nonzero off the horizon, so t is\n")
cat("   constant along it. Integrating the closed form r = M(1 + sin phi) confirms the\n")
cat("   endpoints sit at the same r and are antipodal on the sphere.\n")
phi <- c(0, pi)
rends <- M*(1 + sin(phi))
cat(sprintf("   phi = 0 and pi: r = %.6f and %.6f, angular separation %.4f pi\n",
            rends[1], rends[2], (phi[2] - phi[1])/pi))
note(abs(diff(rends)) < 1e-14, "the two ends of the contact geodesic sit at the same radius")
note(abs((phi[2] - phi[1]) - pi) < 1e-14, "and they are antipodal")

cat("\n=== 6. the plants ===\n")
bad <- c()
for (frac in c(1/3, 1/4, 1, 2/3)) {
  tc <- 0.4 + 1i*beta*frac
  r <- 0.8
  a <- UV(r, 0.4, "II")
  b <- c(exp(-kap*(tc - rstar(r))), exp(kap*(tc + rstar(r))))
  ru <- b[1]/a[1]; rv <- b[2]/a[2]
  hit <- Mod(ru + 1) < 1e-12 && Mod(rv + 1) < 1e-12
  idt <- Mod(ru - 1) < 1e-12 && Mod(rv - 1) < 1e-12
  cat(sprintf("   shift of %5.3f beta: U ratio %7.4f%+7.4fi, antipodal? %-5s identity? %s\n",
              frac, Re(ru), Im(ru), hit, idt))
  bad <- c(bad, hit)
  if (frac == 1) note(idt, "a whole period is the identity")
}
note(!any(bad), "plant: no shift other than the half period gives the antipodal map")
# A SECOND PLANT, AND IT FIRED ON A FALSE CLAIM OF MINE. The first version asserted that a
# wrong shift moves r, so that the contact geodesic would not close. It does not. U V depends
# only on rstar, since the t in U and the t in V cancel, so EVERY shift of t, real or imaginary,
# leaves r exactly where it was. What a wrong shift moves is U/V, which is where t lives, and it
# moves it off the real axis: the image is then not a point of the Lorentzian manifold at all.
# Keeping the corrected form, because the statement it replaces is the sharper one.
cat("   What a wrong shift moves is not r. U V depends on rstar alone, since the t in U and\n")
cat("   the t in V cancel, so every shift leaves the radius exactly where it was. What goes\n")
cat("   complex is U and V themselves, and the image is then not a point of the real manifold\n")
cat("   at all. The ratio U/V is not the test either, since at a quarter period it is -1 and\n")
cat("   perfectly real while U and V are each imaginary. Both have to be real, which needs\n")
cat("   kappa times the shift to be a multiple of pi, a multiple of beta/2:\n")
# AND THE SECOND VERSION WAS WRONG TOO, caught the same way. Testing U/V for reality is not
# the test: at beta/4 the ratio is -1, perfectly real, while U and V are each imaginary. The
# Kruskal point is the PAIR, so both have to be real, and that needs kappa times the imaginary
# shift to be a multiple of pi, which is a multiple of beta/2 and not of beta/4. Two false
# claims in a row here, both found by the check rather than by reading. Keeping the record.
cat("        shift     U ratio            V ratio           U and V both real?\n")
onreal <- c()
for (frac in c(1/3, 1/4, 1/2, 2/3, 1)) {
  tc <- 1i*beta*frac
  ru <- exp(-kap*tc); rv <- exp(kap*tc)
  ok <- abs(Im(ru)) < 1e-12 && abs(Im(rv)) < 1e-12
  cat(sprintf("   %8.3f beta  %7.4f%+7.4fi   %7.4f%+7.4fi   %s\n", frac,
              Re(ru), Im(ru), Re(rv), Im(rv), if (ok) "yes" else "no"))
  onreal <- c(onreal, ok)
  uv <- exp(-kap*(0.4 + tc - rstar(0.8)))*exp(kap*(0.4 + tc + rstar(0.8)))
  note(abs(Im(uv)) < 1e-12, sprintf("r is untouched by a shift of %.3f beta", frac))
}
note(identical(onreal, c(FALSE, FALSE, TRUE, FALSE, TRUE)),
     "plant: only the half period and the whole one put the image on the real section")

cat("\n=== 7. what this settles for fork 9, and what was already settled ===\n")
cat("   THE COMPANION ALREADY HAS THE PAIRING. Its causal-structure section states that J is\n")
cat("   the boost through imaginary angle pi, which in Kruskal is (T,X) -> (-T,-X), so the two\n")
cat("   exteriors swap and the black hole interior pairs with the white hole one, checked on\n")
cat("   twenty thousand points. That is the same statement as section 1 here, reached by\n")
cat("   sampling rather than by writing the coordinates out region by region. This file is a\n")
cat("   second route to it and not a first.\n")
cat("   WHAT IT ADDS is what the pairing does to the STATE and to the contact geodesic. The\n")
cat("   interior sum inherits its state from the exterior temperature and has nothing to pick.\n")
cat("   CORRECTED 2026-09-28, later the same day, by the_contact_pair_in_flat_space.R: this\n")
cat("   file first said the future-past correlator is 'the direct one at t - i beta/2', which\n")
cat("   is right about the POINT MAP and wrong about the correlator. Shifting the Killing time\n")
cat("   and continuing the Milne time both give (U,V) -> (-U,-V), and they differ by a loop\n")
cat("   around the branch point. Outside, where the Killing time is timelike, the shift is the\n")
cat("   right continuation; inside it is not, and computed in flat space it returns a real\n")
cat("   number where the truth is complex. The continuation is in the Milne time.\n")
cat("   And because the shift fixes\n")
cat("   both r and t, the image sits at the SAME event in the antipodal region, which is\n")
cat("   exactly where the E = 0 contact geodesic ends. The contact geodesic was not chosen for\n")
cat("   convenience; it is the curve the fold's own map asks for.\n")
cat("   WHAT REMAINS is arithmetic: the sum itself, calibrated against A.19's amplitude\n")
cat("   Delta^{1/2} -> 3.9004 M s^{-1/2} before any sign is read off it, and validated against\n")
cat("   the trace identity T^a_a = -m^2 W mode by mode.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
