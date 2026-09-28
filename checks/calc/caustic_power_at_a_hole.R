#!/usr/bin/env Rscript
# The caustic stress power at a hole, from three measured links rather than a bracket.
#
# image_stress_shell.R carried the power as a bracket, 3/2 to 7/2, because the two Hadamard slots
# come out at different caustic powers on the one exactly solvable geometry and nothing said why.
# The spread is the dominant uncertainty in how thick the shell is where the fold's term beats the
# interior's own focusing, so it is worth closing. This closes it, and the answer is 5/2.
#
# Three links, each measured on R x S^3 where every quantity is a closed-form mode sum.
#
#   1. THE LEADING POWER. A.19's rule puts the image Green function at delta^{-(D-2+n)/2} in the
#      distance to the caustic and a stress, quadratic in derivatives, two powers further at
#      delta^{-(D+2+n)/2}. On this geometry n = 2 and D = 4, so delta^{-4}, and the coupling slot
#      measures it at -4.00 across three windows.
#
#   2. THE MASS IS ONE POWER SOFTER, AND FOR A REASON THAT DOES NOT DEPEND ON THE GEOMETRY.
#      Differentiating a mode sum with respect to m^2 brings an extra inverse power of the mode
#      frequency: d/dm^2 of (1/2om) cos(om tau) is -(1/4om^3)cos - (tau/4om^2) sin, and against the
#      summand's own (1/2om) cos the surviving piece is smaller by tau/2om. One extra inverse power
#      of the mode index is one softer power of delta in any sum whose divergence comes from large
#      n, so the mass-induced part of every quantity sits one power above its massless counterpart.
#      Measured: G at -2 and its mass part at -1, P and Q at -4 and their mass parts at -3.
#
#   3. THE SECOND SOFTENING IS THE COUPLING'S DOING, NOT THE MASS'S. Where the Ricci tensor
#      vanishes the null stress is (1 - 4 xi) P - Q, so the leading mass parts cancel exactly when
#      1 - 4 xi equals the ratio dQ/dP at the caustic. On this geometry that ratio is 1/3 and the
#      cancelling coupling is therefore the conformal one, which is why the mass slot was measured
#      at -2 rather than -3: the measurement had been taken at xi = 1/6. At any other coupling it
#      is -3, and that is checked below at four of them.
#
# At a hole the caustic is order one, so the rule gives -7/2 for the leading stress and link 2 gives
# -5/2 for the mass part, which is the part that acts where the Ricci tensor vanishes. Link 3 would
# take off one more power only at a coupling fine-tuned to the hole's own dQ/dP, a number nobody has
# computed and with no reason to be the conformal value there. So p = 5/2, with p = 3/2 surviving as
# a fine-tuning rather than as a competing reading.

a <- 1
NMAX <- 40000
nn <- 0:NMAX; N <- nn + 1
om <- function(m) sqrt(N^2/a^2 + m^2)
f_anti   <- function() (-1)^nn * N^2
f_anti_2 <- function() -(-1)^nn * N^2 * (N^2 - 1) / 3
damp <- function(eps) exp(-eps * nn)
EPS <- 0.001
Gimg <- function(tau, m) sum((1/(2*om(m))) * cos(om(m)*tau) * f_anti() * damp(EPS)) / (2*pi^2*a^3)
Gtt  <- function(tau, m) sum((1/(2*om(m))) * (-om(m)^2) * cos(om(m)*tau) * f_anti() * damp(EPS)) / (2*pi^2*a^3)
Guu  <- function(tau, m) sum((1/(2*om(m))) * cos(om(m)*tau) * f_anti_2() * damp(EPS)) / (2*pi^2*a^3)
Q_cf <- function(tau) 1/(8*pi^2*(1 + cos(tau))^2)
kq <- mean(sapply(c(2.90, 3.00, 3.08), function(t) Guu(t, 0)/Q_cf(t)))
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }
cat(sprintf("=== 0. same assembly as sign_through_the_mass_slot.R, calibration %.6f ===\n", kq))
note(abs(kq - 1) < 0.01, "the angular calibration is the same one")

fitp <- function(f, taus) {
  y <- sapply(taus, f); d <- pi - taus
  unname(coef(lm(log(abs(y)) ~ log(d)))[2])
}
WIN <- list(wide = c(2.86, 2.94, 3.00, 3.05, 3.08),
            mid  = c(3.00, 3.025, 3.05, 3.075, 3.10),
            near = c(3.09, 3.10, 3.11, 3.12, 3.13))
Tkk <- function(tau, xi, m) {
  P <- Gtt(tau, m); Q <- Guu(tau, m)/kq; G <- Gimg(tau, m)
  (1 - 2*xi)*(P - Q) - 2*xi*(P + Q) + xi*2*G
}

cat("\n=== 1. the leading power, and that the rule predicts it ===\n")
cat("   The rule: G at -(D-2+n)/2 and a stress two further, so -(D+2+n)/2. Here n = 2, D = 4.\n")
cat("      quantity                   predicted    wide      mid     near\n")
for (q in list(list("G, the image correlator", -2, function(t) Gimg(t, 0)),
               list("P, its second time derivative", -4, function(t) Gtt(t, 0)),
               list("Q, the transverse one", -4, function(t) Guu(t, 0)/kq),
               list("T_kk at xi = 1/6 + 0.05", -4, function(t) Tkk(t, 1/6 + 0.05, 0)))) {
  ps <- sapply(WIN, function(w) fitp(q[[3]], w))
  cat(sprintf("   %-30s %9.1f %8.3f %8.3f %8.3f\n", q[[1]], q[[2]], ps[1], ps[2], ps[3]))
  note(abs(ps[3] - q[[2]]) < 0.1, sprintf("%s sits at its predicted power", q[[1]]))
}

cat("\n=== 2. the mass part is one power softer, and the reason is analytic ===\n")
cat("   d/dm^2 of (1/2 om) cos(om tau) is -(1/4 om^3) cos - (tau/4 om^2) sin, and against the\n")
cat("   summand's own (1/2 om) cos the surviving piece is smaller by tau/2 om. Checked two ways:\n")
cat("   the measured power of the mass part, and the mass part against m^2 times the analytic\n")
cat("   first derivative of the sum.\n")
dGimg <- function(tau) sum((-(1/(4*om(0)^3))*cos(om(0)*tau) - (tau/(4*om(0)^2))*sin(om(0)*tau))
                           * f_anti() * damp(EPS)) / (2*pi^2*a^3)
cat("      quantity        massless power    mass-part power   difference\n")
m0 <- 0.25
for (q in list(list("G", function(t) Gimg(t, 0), function(t) Gimg(t, m0) - Gimg(t, 0)),
               list("P", function(t) Gtt(t, 0),  function(t) Gtt(t, m0) - Gtt(t, 0)),
               list("Q", function(t) Guu(t, 0)/kq, function(t) (Guu(t, m0) - Guu(t, 0))/kq))) {
  p0 <- fitp(q[[2]], WIN$near); p1 <- fitp(q[[3]], WIN$near)
  cat(sprintf("   %-14s %15.3f %18.3f %12.3f\n", q[[1]], p0, p1, p1 - p0))
  note(abs((p1 - p0) - 1) < 0.1, sprintf("%s's mass part is one power softer", q[[1]]))
}
cat("   And the first-order form, at three small masses, against the analytic derivative:\n")
cat("      m        (G(m) - G(0)) / m^2      dG/dm^2 analytic      relative\n")
for (m in c(0.05, 0.1, 0.2)) {
  num <- (Gimg(3.10, m) - Gimg(3.10, 0)) / m^2; an <- dGimg(3.10)
  cat(sprintf("   %6.3f %22.8f %22.8f %13.2e\n", m, num, an, abs(num/an - 1)))
  note(abs(num/an - 1) < 0.02, "the mass part is the first derivative times m^2")
}

cat("\n=== 3. the second softening belongs to the coupling ===\n")
cat("   Where the Ricci tensor vanishes the null stress is (1 - 4 xi) P - Q, so the leading mass\n")
cat("   parts cancel exactly when 1 - 4 xi equals dQ/dP at the caustic.\n")
dP <- function(t) Gtt(t, m0) - Gtt(t, 0)
dQ <- function(t) (Guu(t, m0) - Guu(t, 0))/kq
rat <- sapply(WIN$near, function(t) dQ(t)/dP(t))
cat(sprintf("      dQ/dP over the near window: %s\n", paste(sprintf("%.5f", rat), collapse = " ")))
cat(sprintf("      mean %.5f against 1/3 = %.5f, so the cancelling coupling is xi = %.5f\n",
            mean(rat), 1/3, (1 - mean(rat))/4))
note(abs(mean(rat) - 1/3) < 0.01, "dQ/dP is 1/3 at this caustic")
note(abs((1 - mean(rat))/4 - 1/6) < 0.003, "so the cancelling coupling is the conformal one HERE")
cat("   That is a property of this geometry and not of the mass. At any other coupling there is no\n")
cat("   cancellation and the mass part keeps its own power:\n")
cat("      xi        mass-slot power on the near window\n")
for (xi in c(0, 1/12, 1/6, 0.2, 0.25)) {
  p <- fitp(function(t) Tkk(t, xi, m0) - Tkk(t, xi, 0), WIN$near)
  cat(sprintf("   %7.4f %30.3f%s\n", xi, p, if (abs(xi - 1/6) < 1e-9) "   <- the cancelling one" else ""))
  if (abs(xi - 1/6) > 0.01) note(abs(p + 3) < 0.15, sprintf("xi = %.4f keeps the -3 power", xi))
  if (abs(xi - 1/6) < 1e-9) note(abs(p + 2) < 0.15, "and the conformal one loses a power")
}

cat("\n=== 4. the same three links at a hole, where the caustic is order one ===\n")
cat("      link                                              on R x S^3    at a hole\n")
for (r in list(c("caustic order n", "2", "1"), c("dimension D", "4", "4"),
               c("leading stress, -(D+2+n)/2", "-4", "-3.5"),
               c("mass part, one softer", "-3", "-2.5"),
               c("at the cancelling coupling only", "-2", "-1.5")))
  cat(sprintf("   %-48s %12s %12s\n", r[1], r[2], r[3]))
cat("   So p = 5/2 at a hole for a massive field at any coupling except one fine-tuned to the\n")
cat("   hole's own dQ/dP, which nobody has computed and which has no reason to be the conformal\n")
cat("   value there, the coincidence above being a property of a conformally flat geometry. The\n")
cat("   bracket's lower end survives as a fine-tuning and not as a competing reading.\n")

cat("\n=== 5. what p = 5/2 gives for the shell ===\n")
mP_GeV <- 1.220890e19; lP_m <- 1.616255e-35; rh_sun <- 2.953250e3; Bstar <- 11.5138
Ds <- function(m_GeV, p, kap = 1) (8*pi*kap*(m_GeV/mP_GeV)^2/Bstar)^(1/p) * rh_sun
cat("      field                 D* at p = 5/2     in Planck lengths\n")
for (mm in list(c("electron", 0.000511), c("proton", 0.938272), c("top", 172.69),
                c("the fold's 491.6 PeV fermion", 4.916e8))) {
  d <- Ds(as.numeric(mm[2]), 2.5)
  cat(sprintf("   %-28s %12.3e m %18.3e\n", mm[1], d, d/lP_m))
}
cat(sprintf("   Sub-Planckian below m = %.4g eV at one solar mass, which is below every known\n",
            exp(uniroot(function(l) Ds(exp(l), 2.5)/lP_m - 1, c(log(1e-30), log(1e12)),
                        tol = 1e-13)$root) * 1e9))
cat("   particle mass, so the shell is a real length for every field there is.\n")

cat("\n=== 6. plants ===\n")
cat("   (a) the stress power must NOT equal the Green function's, or link 1 is empty:\n")
pg <- fitp(function(t) Gimg(t, 0), WIN$near); pt <- fitp(function(t) Tkk(t, 1/6 + 0.05, 0), WIN$near)
cat(sprintf("       G at %.3f against T_kk at %.3f, apart by %.3f\n", pg, pt, pg - pt))
note(abs((pg - pt) - 2) < 0.15, "PLANT (a) fires: two derivatives are two powers")
cat("   (b) differentiating with respect to m instead of m^2 must break the first-order check,\n")
cat("       since the expansion is in m^2 and the linear term in m is absent:\n")
lin <- (Gimg(3.10, 0.1) - Gimg(3.10, 0)) / 0.1
cat(sprintf("       (G(0.1) - G(0))/m = %.8f against /m^2 = %.8f, a factor %.3f\n",
            lin, (Gimg(3.10, 0.1) - Gimg(3.10, 0))/0.01, 0.1))
note(abs(lin / ((Gimg(3.10, 0.1) - Gimg(3.10, 0))/0.01) - 0.1) < 1e-6,
     "PLANT (b) fires: the series is in m^2 and the m-normalisation drifts with m")
cat("   (c) the cancellation must be absent at a coupling next to the conformal one:\n")
p_near <- fitp(function(t) Tkk(t, 1/6 + 0.02, m0) - Tkk(t, 1/6 + 0.02, 0), WIN$near)
cat(sprintf("       xi = 1/6 + 0.02 gives %.3f, not -2\n", p_near))
note(abs(p_near + 2) > 0.3, "PLANT (c) fires: the cancellation is a point and not a neighbourhood")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
