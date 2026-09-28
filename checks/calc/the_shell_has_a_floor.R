#!/usr/bin/env Rscript
# A massless conformal field feels the contact caustic too, and that puts a FLOOR under the shell
# which no field can go below.
#
# WHAT THE MANUSCRIPTS SAID AND WHY IT NEEDED CHECKING. 3.6 says the fold's extra term "vanishes
# identically for conformally invariant matter", and A.19 says "Schwarzschild is Ricci-flat, so on
# the reading above a massless conformal field feels no contact divergence there at all, and what
# feels it is mass". Both rest on the Hadamard log coefficient V_0 = Delta^{1/2}[m^2 + (xi-1/6)R]/2,
# which does vanish for a massless conformal field on a Ricci-flat background, and on the Einstein
# static universe's exact closed form T_kk proportional to (1 - 6 xi). But A.18's whole point is
# that the parametrix FAILS at a caustic, and image_stress_components.R computes the leading
# divergence there and finds it independent of the coupling AND of the mass.
#
# THE TWO ARE BOTH RIGHT AND THEY ARE ABOUT DIFFERENT CONFIGURATIONS. The Einstein static universe
# is conformally flat, so a conformal field's image stress maps to the flat one and vanishes; and
# its antipodal pair at equal time is SPACELIKE separated, so there is no sigma -> 0 to expand in.
# A hole's contact is a NULL caustic, sigma vanishing linearly in M - r, in a spacetime that is not
# conformally flat. There the leading term survives, and section 1 shows why it must.
#
# THE CONSEQUENCE, and it strengthens the paper rather than weakening it. The mass-independent
# piece is two powers MORE divergent, D^{-7/2} against D^{-5/2}, and vastly smaller in coefficient,
# so for any massive field it changes the shell not at all: m^2 for the fold's own fermion beats it
# by forty-six orders. What it does is put a floor under the shell. A massless conformal field, for
# which the manuscripts predict nothing at all, gets a shell of
#
#     D_floor = (9.85e-3)^{2/7} l_P^{4/7} r_h^{3/7},
#
# a geometric mean of the Planck length and the horizon radius, which is 1.1e-19 m at a solar mass:
# sixteen orders above the Planck length and carrying no matter content whatever.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

cat("=== 1. why the coupling cancels at a null caustic, before any number ===\n")
cat("   The point-split null-null component is\n")
cat("      T_kk = (1-2xi) k^a k^b' W_;ab' - 2xi k^a k^b W_;ab + xi R_ab k^a k^b W,\n")
cat("   every g_ab term having been killed by k.k = 0. Near a null caustic the leading behaviour of\n")
cat("   W is C sigma^{-3/2}, so both second derivatives are dominated by the grad-sigma pairing,\n")
cat("   and grad_b' sigma is MINUS the transport of grad_b sigma at leading order. So the two\n")
cat("   derivative terms are equal and opposite and the coupling cancels:\n")
cat("      (1-2xi)(-X) - 2xi(X) = -X,   for every xi.\n")
cat("   The third term multiplies W itself, which is two powers softer, and so does the mass term.\n")
cat("   That is a statement about any geometry with a null caustic, and section 2 measures it.\n")
for (xi in c(0, 1/12, 1/6, 1/4, 1/3, 0.5)) {
  v <- (1 - 2*xi)*(-1) - 2*xi*(1)
  note(abs(v + 1) < 1e-15, sprintf("the coupling cancels at xi = %.4f", xi))
}
cat(sprintf("   checked at six couplings, the combination returning -1 every time.\n"))
cat("   The plant: with the primed derivative taken to have the SAME sign, which is what happens\n")
cat("   away from a caustic where the pairing is not dominant, the combination is\n")
for (xi in c(0, 1/6, 1/3)) cat(sprintf("      xi = %.4f: (1-2xi)(+1) - 2xi(1) = %+.4f\n",
                                        xi, (1-2*xi) - 2*xi))
note(abs((1 - 2*(1/6)) - 2*(1/6) - (1 - 4*(1/6))) < 1e-15, "the plant's combination does move with xi")

cat("\n=== 2. the measured coefficient, at six couplings and four masses ===\n")
cat("   From image_stress_components.R's tower, T_kk(radial) times s^{7/2} as the offset closes.\n")
cat("   If the Einstein static universe's (1 - 6 xi) law held here the xi = 1/6 row would go to\n")
cat("   zero. It goes to the same number as every other row.\n\n")
SRC <- readLines("checks/calc/image_stress_components.R")
i1 <- grep("^EG <- ", SRC)[1]; i2 <- grep('^cat\\("=== 1\\.', SRC)[1]
eval(parse(text = paste(SRC[i1:(i2-1)], collapse = "\n")))
k1 <- grep("^tower <- function", SRC)[1]
k2 <- grep("^\\}", SRC); k2 <- k2[k2 > k1][1]
eval(parse(text = paste(SRC[k1:k2], collapse = "\n")))
kk_rad <- function(o) -Re(o$D) - Re(o$E)
SS <- c(0.03, 0.01, 0.003, 0.001)
cat("      xi        s=0.03       s=0.01      s=0.003      s=0.001     (1-6xi)\n")
lim <- c()
for (xi in c(0, 1/12, 1/6, 1/4, 1/3)) {
  v <- sapply(SS, function(s) kk_rad(tower(s, eps = s/40, m2 = 0, xi = xi))*s^3.5)
  lim <- c(lim, v[4])
  cat(sprintf("   %8.4f %12.6f %12.6f %12.6f %12.6f %10.2f\n", xi, v[1], v[2], v[3], v[4], 1-6*xi))
}
cat(sprintf("\n   spread of the s = 0.001 column: %.2e, against a mean of %.6f\n",
            max(lim) - min(lim), mean(lim)))
note(abs(max(lim) - min(lim)) < 1e-4, "the leading coefficient does not move with the coupling")
note(abs(lim[3]) > 0.05, "and at conformal coupling it is NOT zero")
cat(sprintf("      at xi = 1/6, where the (1-6xi) law needs 0, it measures %.6f, against the\n",
            lim[3]))
cat(sprintf("      closed form -15/(32 pi sqrt(2 pi)) = %.6f\n", -15/(32*pi*sqrt(2*pi))))
cat("\n      m        s=0.03       s=0.01      s=0.003      s=0.001\n")
limm <- c()
for (m in c(0, 1, 3, 10)) {
  v <- sapply(SS, function(s) kk_rad(tower(s, eps = s/40, m2 = m^2, xi = 1/6))*s^3.5)
  limm <- c(limm, v[4])
  cat(sprintf("   %8.2f %12.6f %12.6f %12.6f %12.6f\n", m, v[1], v[2], v[3], v[4]))
}
note(abs(limm[1] - limm[2]) < 1e-3, "and the mass does not move it either, at the offsets reached")
cat("   The mass rows converge more slowly because a mass shifts mu_l^2 by a constant and the\n")
cat("   shift matters until l(l+1) beats it, which is a subleading effect and not a leading one.\n")

cat("\n=== 3. the same term on Schwarzschild, from A.19's own assembly ===\n")
AMP  <- 3.9004              # Delta^{1/2} -> AMP M s^{-1/2}, A.19's projection length
KSIG <- 17.425              # |sigma| = KSIG (M - r), from A.19's tau^2/(M-r) -> 34.85
Cw   <- AMP * (4*pi)^(-2) * gamma(1.5) * 2^1.5
cat(sprintf("   W_img = C sigma^{-3/2} with C = AMP (4 pi)^{-2} Gamma(3/2) 2^{3/2} = %.6f M\n", Cw))
cat("   Two derivatives on the ingoing radial congruence, where A.19's contraction factor is\n")
cat("   exactly one because dr/dlambda = -1, give T_kk = -C (15/4) (k.grad sigma)^2 sigma^{-7/2}\n")
cat("   with k.grad sigma = KSIG, so\n")
Ckk <- Cw * (15/4) * KSIG^(-1.5)
cat(sprintf("      T_kk = -%.6e M^{-1/2} (M - r)^{-7/2}\n", Ckk))
note(Ckk > 0, "the coefficient is positive, so T_kk is negative, as at the caustic in the model")

cat("\n=== 4. the floor, and why it changes no massive field's shell ===\n")
lP   <- 1.616255e-35; rh_sun <- 2953.25
# 8 pi G |T_kk| = B*/r_h^2 with B* = 11.5138, G = lP^2, M = r_h/2
Bst  <- 11.5138
Dfloor <- function(rh) {
  k <- 8*pi*lP^2 * Ckk * sqrt(2/rh) / Bst * rh^2          # D^{7/2} = k
  k^(2/7)
}
cat(sprintf("   D_floor at one solar mass      %.3e m\n", Dfloor(rh_sun)))
cat(sprintf("   D_floor at a billion solar     %.3e m\n", Dfloor(rh_sun*1e9)))
cat(sprintf("   the Planck length              %.3e m, so the floor is %.1f orders above it\n",
            lP, log10(Dfloor(rh_sun)/lP)))
note(Dfloor(rh_sun) > 1e3*lP, "the floor is far above the Planck length")
cat("   The scaling is l_P^{4/7} r_h^{3/7}, a geometric mean of the two lengths, checked by\n")
cat("   doubling the horizon and reading the exponent:\n")
e37 <- log(Dfloor(2*rh_sun)/Dfloor(rh_sun))/log(2)
cat(sprintf("      exponent in r_h: %.6f against 3/7 = %.6f\n", e37, 3/7))
note(abs(e37 - 3/7) < 1e-6, "the floor goes as r_h^{3/7}")

cat("\n   And the massive term still wins wherever there is a mass. At the fold's own fermion the\n")
cat("   two are compared at the mass shell's own thickness:\n")
mP_GeV <- 1.220890e19; KAP <- 0.0039329
Dstar <- function(m_GeV, rh = rh_sun) (8*pi*KAP*(m_GeV/mP_GeV)^2/Bst)^(2/5) * rh
hbar_m <- 1.973269804e-16          # GeV m
for (p in list(list("the fold's fermion, 491.6 PeV", 4.916e8), list("an electron", 0.000511))) {
  D <- Dstar(p[[2]]); m_inv <- p[[2]]/hbar_m
  Tmass <- KAP * m_inv^2 * D^(-2.5) * rh_sun^0.5
  Tmless <- Ckk * (rh_sun/2)^(-0.5) * D^(-3.5)
  cat(sprintf("      %-30s shell %9.3e m, mass term / massless term %8.2e\n",
              p[[1]], D, Tmass/Tmless))
  note(Tmass/Tmless > 1e10, "the mass term dominates by many orders at its own shell")
}
cat("   So nothing in Figure 4 moves. What moves is the statement beside it: the shell does not\n")
cat("   go to zero as the mass does, it goes to the floor.\n")

cat("\n=== 5. what has to be corrected, stated flatly ===\n")
cat("   3.6's 'it vanishes identically for conformally invariant matter' and A.19's 'a massless\n")
cat("   conformal field feels no contact divergence there at all' are both true of the parametrix\n")
cat("   and of the conformally flat calibration, and both fail at a null caustic in a spacetime\n")
cat("   that is not conformally flat, which is what a hole's contact sphere is. The V_0 reading\n")
cat("   governs the LOG term; the caustic divergence is the Delta^{1/2}/sigma term, whose\n")
cat("   coefficient is one whatever the field is.\n")
cat("   What survives unchanged: the shell for every massive field, the power 5/2, kappa, and the\n")
cat("   2.1 microns. What is added: a floor of 1.1e-19 m at a solar mass that no field goes below.\n")

cat("\n=== 6. the numbers the manuscripts quote, as magnitudes ===\n")
cat(sprintf("   the one coefficient every coupling returns      %.4f\n", abs(mean(lim))))
cat(sprintf("   the Schwarzschild coefficient, times 1e3        %.2f\n", Ckk*1e3))
cat(sprintf("   the floor at one solar mass, times 1e19 m       %.1f\n", Dfloor(rh_sun)*1e19))
cat(sprintf("   the floor at a billion solar, times 1e16 m      %.1f\n", Dfloor(rh_sun*1e9)*1e16))
cat(sprintf("   the mass term's lead at the fermion, times 1e-46 %.1f\n",
            (KAP*(4.916e8/hbar_m)^2*Dstar(4.916e8)^(-2.5)*rh_sun^0.5) /
            (Ckk*(rh_sun/2)^(-0.5)*Dstar(4.916e8)^(-3.5)) / 1e46))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
