# The sign, which is A.15's actual question: does the divergence close the contact region or
# merely mark it?
#
# The quantity that decides it is not the energy density a static observer sees. It is T_kk, the
# null-null component, because that is what Raychaudhuri uses. Positive T_kk focuses the congruence
# harder, brings the conjugate point forward and makes contact EASIER, so the region grows.
# Negative T_kk defocuses, delays the conjugate point and makes contact harder, so the region
# shrinks and the construction censors itself.

cat("=== 1. the exact case, and what it says the stress is proportional to ===\n")
cat("   On the Einstein static universe image_stress_conformal.R gives, exactly,\n")
cat("        T_kk = (1 - 6 xi)(1 - cos eta) / ( 8 pi^2 a^4 (1 + cos eta)^2 ).\n")
cat("   The Hadamard log coefficient there is V_0 = Delta^{1/2}[m^2 + (xi - 1/6) R]/2 with R =\n")
cat("   6/a^2 and m = 0, so V_0 = 3 Delta^{1/2} (xi - 1/6)/a^2. Both carry the same factor and\n")
cat("   carry it with OPPOSITE signs, since 1 - 6 xi = -6(xi - 1/6):\n\n")
a <- 1
Tkk <- function(eta, xi) (1 - 6 * xi) * (1 - cos(eta)) / (8 * pi^2 * a^4 * (1 + cos(eta))^2)
V0  <- function(xi, Dh = 1) Dh * (0 + (xi - 1/6) * 6 / a^2) / 2
cat("      xi       T_kk at eta = pi - 0.1     V_0 / Delta^{1/2}      product's sign\n")
for (xi in c(0, 1/12, 1/6, 1/4, 1/3)) {
  tk <- Tkk(pi - 0.1, xi); v <- V0(xi)
  cat(sprintf("   %7.4f  %22.4f  %20.5f      %s\n", xi, tk, v,
              ifelse(abs(tk) < 1e-9, "both zero", ifelse(tk * v < 0, "opposite", "same"))))
  if (abs(xi - 1/6) > 1e-9) stopifnot(tk * v < 0)
}
cat("\n   So the caustic-divergent stress is proportional to MINUS V_0, with a positive constant.\n")
cat("   That is a calibration on one exact case and it is stated as one.\n")

cat("\n=== 2. what that gives at a hole ===\n")
cat("   Schwarzschild is Ricci-flat, so R = 0 and the (xi - 1/6) R channel is shut whatever the\n")
cat("   coupling. The same slot carries m^2, and m^2 > 0 for every field that is not massless,\n")
cat("   so V_0 = m^2 Delta^{1/2}/2 > 0 and\n")
cat("        T_kk proportional to -V_0  <  0.\n")
cat("   A NEGATIVE null-null stress at the contact surface.\n")
for (m in c(0, 0.511e-3, 105.7e-3, 172.57)) {
  v <- m^2 / 2
  cat(sprintf("      m = %10.4f GeV:  V_0 / Delta^{1/2} = %12.4e   T_kk sign: %s\n",
              m, v, ifelse(v > 0, "negative", "zero, nothing at leading order")))
}

cat("\n=== 3. and what a negative T_kk does, from Raychaudhuri rather than from words ===\n")
cat("   dtheta/dlambda = -theta^2/2 - shear^2 - R_kk with R_kk = 8 pi G T_kk. A negative T_kk\n")
cat("   makes R_kk negative, which makes dtheta/dlambda LESS negative, so the congruence focuses\n")
cat("   more slowly and the conjugate point is pushed further along. Integrate it and see where\n")
cat("   the conjugate point goes:\n")
conj <- function(Rkk, lam_max = 60, n = 400000) {
  h <- lam_max / n; J <- 0; Jp <- 1                 # J'' = -(Rkk/2) J for an isotropic congruence
  for (i in 1:n) {
    Jpp <- -(Rkk / 2) * J
    Jp2 <- Jp + h * Jpp; J2 <- J + h * Jp
    J <- J2; Jp <- Jp2
    if (i > 10 && J < 0) return(i * h)
  }
  NA
}
cat("      R_kk      first conjugate point\n")
for (Rkk in c(0.5, 0.2, 0.05, 0, -0.05)) {
  z <- conj(Rkk)
  cat(sprintf("   %8.2f  %22s\n", Rkk, ifelse(is.na(z), "none within 60", sprintf("%.3f", z))))
}
cat("   Positive R_kk brings it forward, zero and negative push it to infinity. So a negative\n")
cat("   T_kk delays the conjugate point, and since contact REQUIRES a conjugate point at the\n")
cat("   antipodal angle, delaying it makes contact harder and the region shrinks.\n")

cat("\n=== 4. the answer to A.15's question ===\n")
cat("   The divergence CLOSES the region rather than merely marking it. The stress the contact\n")
cat("   creates is negative in the null-null direction, which defocuses the very congruence whose\n")
cat("   refocusing is what creates the contact. The construction polices its own causality\n")
cat("   violation, and it does so for the fields that exist rather than for a special choice of\n")
cat("   coupling, because at a Ricci-flat hole the channel that survives is the mass one.\n")
cat("   Status, stated once: the proportionality to -V_0 is calibrated on the Einstein static\n")
cat("   universe, which is the one exactly solvable caustic available, so the SIGN rested on one\n")
cat("   case when this was written. It does not now. sign_through_the_mass_slot.R reaches it again\n")
cat("   through the OTHER slot of the same Hadamard coefficient, turning on a mass at conformal\n")
cat("   coupling rather than moving the coupling at zero mass, and T_kk comes out negative both\n")
cat("   ways. That file also finds the proportionality below does not transfer between the slots:\n")
cat("   the coupling slot diverges as 1/D^2 at the caustic and the mass slot as 1/D, so the\n")
cat("   constant calibrated here belongs to this slot and is not a general one.\n")

cat("\n=== 5. the plants ===\n")
cat("   (a) the proportionality must fail where there is no caustic, or it is not about caustics.\n")
cat("      flat inversion: R = 0 and m = 0 so V_0 = 0, yet the stress is +1/96 pi^2 R^4, nonzero.\n")
cat(sprintf("      V_0 = %.1f there and the stress is %+.6f at R = 1, so the relation does not\n",
            0, 1 / (96 * pi^2)))
cat("      apply, exactly as it should not: that configuration has no conjugate point.\n")
cat("   (b) the Raychaudhuri integration must find the conjugate point where it is known to be.\n")
cat("      For R_kk constant the equation is J'' = -(R_kk/2) J, whose first zero is at\n")
cat("      pi/sqrt(R_kk/2):\n")
for (Rkk in c(0.5, 0.2, 0.05)) {
  cat(sprintf("        R_kk = %.2f:  integrated %.3f   closed form %.3f\n",
              Rkk, conj(Rkk), pi / sqrt(Rkk / 2)))
  stopifnot(abs(conj(Rkk) / (pi / sqrt(Rkk / 2)) - 1) < 5e-3)
}
