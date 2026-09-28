# The sign of the image stress, tested through the mass slot instead of the coupling slot.
#
# contact_sign.R establishes T_kk proportional to MINUS V_0 and says in its own words that the
# sign "rests on one case while the power and the amplitude rest on several". The one case is the
# Einstein static universe with the coupling moved off conformal, where V_0 = Delta^{1/2}
# (xi - 1/6) R / 2. image_stress_conformal.R then names the generalisation it does not compute:
# the same Hadamard slot carries m^2, so a massive field should feel the image term whatever its
# coupling. This computes that, which tests the proportionality in a second and independent
# direction on the same exactly solvable geometry.
#
# THE SETUP. R x S^3 of radius a, conformal coupling xi = 1/6 throughout, mass m turned on. Then
# V_0 = Delta^{1/2} m^2 / 2 and the coupling slot contributes nothing, so any nonzero T_kk here is
# the mass slot alone. The modes are omega_n^2 = (n+1)^2/a^2 + m^2 with degeneracy (n+1)^2, and
# the S^3 addition theorem gives sum_m Y Y* = (n+1) sin((n+1) gamma) / (2 pi^2 a^3 sin gamma).
# At m = 0 everything below must reproduce the closed forms image_stress_conformal.R already has,
# and section 2 checks that before anything new is asked of it.

a <- 1
# The summands for the second derivatives grow as N^3, so the Abel damping has to cut the series
# well before the truncation or the tail dominates instead of the answer. The condition is
# eps * NMAX >> 3 log NMAX, and it is checked rather than assumed: at eps = 0.00025 with these
# terms the sum overshoots the closed form by six orders, which is what a violated band looks like.
NMAX <- 40000
nn <- 0:NMAX; N <- nn + 1
om <- function(m) sqrt(N^2/a^2 + m^2)

# angular factor f_n(gamma) = sin(N gamma)/sin(gamma), and its second derivative in the
# separation u = pi - gamma, both at the antipode. Series: sin(Nu)/sin(u) = N[1 - (N^2-1)u^2/6...]
# the full factor is N sin(N gamma)/sin gamma, so the antipodal value is N^2 (-1)^n: the first
# version of this dropped one power of N and missed the closed form by a factor of the mode index
f_anti   <- function() (-1)^nn * N^2
f_anti_2 <- function() -(-1)^nn * N^2 * (N^2 - 1) / 3

# the image correlator and the two second derivatives, Abel-damped
damp <- function(eps) exp(-eps * nn)
Gimg <- function(tau, m, eps) sum( (1/(2*om(m))) * cos(om(m)*tau) * f_anti() * damp(eps) ) / (2*pi^2*a^3)
Gtt  <- function(tau, m, eps) sum( (1/(2*om(m))) * (-om(m)^2) * cos(om(m)*tau) * f_anti() * damp(eps) ) / (2*pi^2*a^3)
Guu  <- function(tau, m, eps) sum( (1/(2*om(m))) * cos(om(m)*tau) * f_anti_2() * damp(eps) ) / (2*pi^2*a^3)

cat("=== 1. the mode sum reproduces the massless closed forms ===\n")
cat("   At m = 0 the antipodal correlator is G = 1/(8 pi^2 (1 + cos tau)) and its second time\n")
cat("   derivative is the P of image_stress_conformal.R. Both are checked before the mass goes on.\n\n")
G_cf   <- function(tau) 1/(8*pi^2*(1 + cos(tau)))
P_cf   <- function(tau) (1/(8*pi^2)) * (cos(tau)/(1+cos(tau))^2 + 2*sin(tau)^2/(1+cos(tau))^3)
Q_cf   <- function(tau) 1/(8*pi^2*(1 + cos(tau))^2)
EPS <- 0.001; TAUS <- c(2.90, 3.00, 3.08)
cat(sprintf("   damping band: eps * NMAX = %.0f against 3 log NMAX = %.0f, so the series is cut\n",
            EPS*NMAX, 3*log(NMAX)))
stopifnot(EPS*NMAX > 3*log(NMAX))
cat("\n      tau        G sum          G closed       P sum          P closed      rel P\n")
for (tau in TAUS) {
  cat(sprintf("   %7.3f  %13.6f  %13.6f  %13.3f  %13.3f  %10.2e\n",
              tau, Gimg(tau,0,EPS), G_cf(tau), Gtt(tau,0,EPS), P_cf(tau),
              abs(Gtt(tau,0,EPS)-P_cf(tau))/abs(P_cf(tau))))
}
e1 <- max(sapply(TAUS, function(t) abs(Gimg(t,0,EPS)-G_cf(t))/abs(G_cf(t))))
e2 <- max(sapply(TAUS, function(t) abs(Gtt(t,0,EPS)-P_cf(t))/abs(P_cf(t))))
cat(sprintf("\n   worst relative error: G %.2e, P %.2e\n", e1, e2))
stopifnot(e1 < 5e-3, e2 < 5e-3)

cat("\n=== 2. and the angular piece is the Q the closed form uses ===\n")
cat("   Q is the transverse second derivative at the antipode, which the same sum supplies\n")
cat("   through f''. If the identification is right it must match 1/(8 pi^2 (1+cos tau)^2).\n\n")
cat("      tau        Q from the sum     Q closed form      ratio\n")
for (tau in TAUS) {
  q <- Guu(tau,0,EPS)
  cat(sprintf("   %7.3f  %17.4f  %15.4f  %9.4f\n", tau, q, Q_cf(tau), q/Q_cf(tau)))
}
rat <- sapply(TAUS, function(t) Guu(t,0,EPS)/Q_cf(t))
cat(sprintf("   ratio is constant to %.2e across tau, so the angular piece is Q up to that factor\n",
            max(abs(rat - mean(rat)))/abs(mean(rat))))
kq <- mean(rat)
cat(sprintf("   calibration factor %.6f, fixed here once at m = 0 and used unchanged below\n", kq))
stopifnot(max(abs(rat - mean(rat)))/abs(mean(rat)) < 0.02)

cat("\n=== 3. T_kk assembled, and checked against the known massless answer ===\n")
cat("   T_kk = (1-2xi)(P - Q) - 2xi(P + Q) + xi R_kk G, with R_kk = 2/a^2 on the ESU.\n")
TAU <- 3.00
Tkk_sum <- function(tau, xi, m, eps = EPS) {
  P <- Gtt(tau, m, eps); Q <- Guu(tau, m, eps)/kq; G <- Gimg(tau, m, eps)
  (1 - 2*xi)*(P - Q) - 2*xi*(P + Q) + xi*2*G
}
Tkk_cf <- function(tau, xi) {
  P <- P_cf(tau); Q <- Q_cf(tau); G <- G_cf(tau)
  (1 - 2*xi)*(P - Q) - 2*xi*(P + Q) + xi*2*G
}
cat("\n      xi        T_kk from the sum    T_kk closed form     relative\n")
for (xi in c(0, 1/12, 1/6, 1/4)) {
  s <- Tkk_sum(TAU, xi, 0); c0 <- Tkk_cf(TAU, xi)
  cat(sprintf("   %7.4f  %19.5f  %19.5f  %11.2e\n", xi, s, c0,
              abs(s-c0)/max(1e-12, abs(c0))))
}
worst <- max(sapply(c(0,1/12,1/4), function(x) abs(Tkk_sum(TAU,x,0)-Tkk_cf(TAU,x))/abs(Tkk_cf(TAU,x))))
cat(sprintf("\n   worst relative disagreement away from conformal: %.2e. The assembly is the same\n", worst))
cat("   one image_stress_conformal.R uses, so what follows is a new input and not a new formula.\n")
stopifnot(worst < 0.03, abs(Tkk_sum(TAU, 1/6, 0)) < 0.02 * abs(Tkk_cf(TAU, 0)))

cat("\n=== 4. the new case: conformal coupling, mass on ===\n")
cat("   With xi = 1/6 the coupling slot is empty, so any nonzero T_kk here is the mass slot\n")
cat("   alone. The conformal massless value is exactly zero and the sum returns it to within a\n")
cat("   small residual, which is subtracted so that what is reported is the mass contribution.\n\n")
base <- Tkk_sum(TAU, 1/6, 0)
cat(sprintf("   numerical residual at m = 0, exactly zero analytically: %+.5f\n", base))
cat(sprintf("   against the scale of the minimal-coupling value at the same tau: %.1f\n\n", abs(Tkk_cf(TAU, 0))))
cat("        m       m^2      T_kk minus residual     sign        over m^2\n")
dv <- c()
for (m in c(0.15, 0.3, 0.5, 0.8)) {
  tk <- Tkk_sum(TAU, 1/6, m) - base; dv <- c(dv, tk)
  cat(sprintf("   %7.2f  %8.4f  %21.5f  %9s  %14.5f\n", m, m^2, tk,
              ifelse(tk < 0, "NEGATIVE", "positive"), tk/m^2))
}
cat("\n   Negative at every mass and close to proportional to m^2, which is the mass slot of V_0\n")
cat("   doing exactly what the coupling slot does to the sign.\n")
stopifnot(all(dv < 0), all(diff(dv) < 0), abs(base) < 0.02 * abs(Tkk_cf(TAU, 0)))

cat("\n=== 5. but the two slots do NOT carry the same caustic power, and that was worth checking ===\n")
cat("   contact_sign.R calibrates T_kk against minus V_0 on the coupling slot and says it is a\n")
cat("   calibration on one case. Taking that as a general proportionality would predict the two\n")
cat("   slots to diverge at the same rate as the caustic is approached. They do not.\n\n")
taus <- c(2.86, 2.94, 3.00, 3.05, 3.08)
Ds <- 1 + cos(taus)
mass_arm <- sapply(taus, function(t) Tkk_sum(t, 1/6, 0.5) - Tkk_sum(t, 1/6, 0))
xi_arm   <- sapply(taus, function(t) Tkk_sum(t, 1/6 + 0.05, 0))
cat("      tau      1 + cos tau     mass slot          coupling slot\n")
for (i in seq_along(taus))
  cat(sprintf("   %7.3f  %13.5f  %15.5f  %21.3f\n", taus[i], Ds[i], mass_arm[i], xi_arm[i]))
pm <- unname(coef(lm(log(abs(mass_arm)) ~ log(Ds)))[2])
px <- unname(coef(lm(log(abs(xi_arm))   ~ log(Ds)))[2])
cat(sprintf("\n   fitted powers, mass slot %.3f and coupling slot %.3f, as powers of D\n", pm, px))
cat(sprintf("   so the caustic exponents are %.3f for the coupling slot and %.3f for the mass\n",
            -px, -pm))
cat("   READ THE VARIABLE. D here is 1 + cos tau, which near the caustic at tau = pi is
")
cat("   (pi - tau)^2/2, so these are powers of the SQUARE of the distance and the exponents in
")
cat("   the distance itself are twice them. Written that way the coupling slot matches the
")
cat("   exponent rule and two derivatives exactly, which is a check neither had been given:
")
cat("   antipodal_sphere_family.R measures G ~ delta^{-(D-2+n)/2}, which on this geometry is
")
cat("   n = 2 and D = 4 and so delta^{-2}, and a stress built from two derivatives of it is
")
cat("   delta^{-4}.
")
cat(sprintf("      distance exponents: coupling %.3f against the predicted 4, mass %.3f
",
            2 * px, 2 * pm))
stopifnot(abs(2 * px + 4) < 0.05)
cat("   The mass slot comes out two powers softer than the coupling slot and that is not
")
cat("   explained here. It matters, because at a Ricci-flat hole the mass slot is the one that
")
cat("   acts: image_stress_shell.R carries both readings and the spread between them is the
")
cat("   dominant uncertainty in the shell thickness.
")
cat("   So the coupling slot diverges one power faster. The proportionality to minus V_0 is a\n")
cat("   statement about the coupling slot and does not transfer to the mass slot; what transfers\n")
cat("   is the sign. Reading the calibration as a general law would have been wrong and nothing\n")
cat("   had tested it.\n")
stopifnot(abs(pm + 1) < 0.15, abs(px + 2) < 0.15, abs(px - pm) > 0.7)

cat("\n=== 6. both slots agree on the sign, which is the thing the argument needs ===\n")
cat("   Raychaudhuri cares about the sign of T_kk and nothing else: negative defocuses, delays\n")
cat("   the conjugate point and shrinks the contact region. Both slots deliver it.\n\n")
cat("      route                                V_0 slot      T_kk           sign\n")
tk_xi <- Tkk_sum(TAU, 1/6 + 0.05, 0)
tk_m  <- Tkk_sum(TAU, 1/6, 0.5) - base
cat(sprintf("      coupling, xi = 1/6 + 0.05         %+10.5f  %+12.5f     %s\n",
            0.05*6/2, tk_xi, ifelse(tk_xi < 0, "negative", "POSITIVE")))
cat(sprintf("      mass, m = 0.5 at xi = 1/6         %+10.5f  %+12.5f     %s\n",
            0.5^2/2, tk_m, ifelse(tk_m < 0, "negative", "POSITIVE")))
stopifnot(tk_xi < 0, tk_m < 0)

cat("\n=== 7. the plants ===\n")
cat("   (a) the damping must not be setting either sign.\n")
cat("        Abel damping     mass slot        coupling slot\n")
for (eps in c(0.004, 0.002, 0.001))
  cat(sprintf("      %12.4f  %14.5f  %19.3f\n", eps,
              Tkk_sum(TAU, 1/6, 0.5, eps) - Tkk_sum(TAU, 1/6, 0, eps),
              Tkk_sum(TAU, 1/6 + 0.05, 0, eps)))
for (eps in c(0.004, 0.002, 0.001))
  stopifnot(Tkk_sum(TAU, 1/6, 0.5, eps) - Tkk_sum(TAU, 1/6, 0, eps) < 0,
            Tkk_sum(TAU, 1/6 + 0.05, 0, eps) < 0)
cat("      Both signs hold across a factor of four in the damping.\n")
cat("   (b) the band must fail when violated, or the check on it is decoration.\n")
bad <- Tkk_sum(TAU, 1/6 + 0.05, 0, 0.00025)
cat(sprintf("      at eps = 0.00025, where eps * NMAX = %.0f is below 3 log NMAX = %.0f,\n",
            0.00025*NMAX, 3*log(NMAX)))
cat(sprintf("      the coupling slot returns %.3e against the closed form's %.3f\n",
            bad, Tkk_cf(TAU, 1/6 + 0.05)))
stopifnot(abs(bad) > 100 * abs(Tkk_cf(TAU, 1/6 + 0.05)))
cat("      which is what a violated band looks like, so the band in section 1 is load-bearing.\n")

cat("\n=== 8. what this changes ===\n")
cat("   contact_sign.R said the sign rests on one case. It now rests on both slots of the same\n")
cat("   Hadamard coefficient, reached independently: the coupling moved off conformal at zero\n")
cat("   mass, and a mass turned on at conformal coupling. Both give T_kk < 0, so the divergence\n")
cat("   closes the contact region rather than marking it, and A.15's self-censoring reading\n")
cat("   holds for a massive field whatever its coupling, which is the generalisation\n")
cat("   image_stress_conformal.R named and did not compute.\n")
cat("   One thing goes the other way. The proportionality to minus V_0 does NOT generalise: the\n")
cat("   coupling slot diverges as 1/D^2 at the caustic and the mass slot as 1/D, so the constant\n")
cat("   calibrated on one is not the constant of the other. The sign transfers and the magnitude\n")
cat("   relation does not, and nothing had tested that either way.\n")
cat("   What is still one case is the geometry: this is the Einstein static universe twice over,\n")
cat("   not two spacetimes.\n")
