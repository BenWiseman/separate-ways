# SUPERSEDED IN PART, by checks/calc/image_stress_conformal.R.
# The closed form G_img = 1/(16 pi^2 a^2 cos^2(eta/2)) below is correct and still used.
# The DENSITY computed in sections 3 and 4, rho delta^4 = 3/4 pi^2, is WITHDRAWN: it drops
# the mixed transverse derivative on a homogeneity argument that constrains the Laplacian
# of phi^2 rather than the mixed derivative. Restored, the fourth power cancels and the
# conformal improvement cancels the remainder, so the correct answer is zero.
#
# The image stress at a caustic, exactly, in the one geometry that has everything needed.
#
# GR51 proved the Schwarzschild-chart machinery blind to the contact and voided five
# iterations of numbers. A.18's method was the right one all along: pick a geometry where
# the caustic is explicit and the sum is doable. A.19 fixes what that geometry must have --
# a caustic, since without one there is nothing to compute, and a VARYING transverse radius,
# since constant radius gives an identically zero stress. Flat space has the radius and no
# caustic. A.19's model has the caustic and constant radius.
#
# The Einstein static universe has both. ds^2 = -dt^2 + a^2[dchi^2 + sin^2 chi dOmega^2]:
# the transverse radius a sin(chi) varies, every null geodesic from a point refocuses at the
# antipode after coordinate time pi a, and the antipodal map on S^3 is free. It is also
# exactly solvable, which none of the previous five iterations' settings were.
#
# THE EXACT CONSTRUCTION. For a conformally coupled massless scalar on the ESU the mode
# frequencies are omega_n = n/a exactly, n = 1, 2, ..., with degeneracy n^2, because the
# conformal coupling's xi R term supplies the +1 that turns sqrt(n^2 - 1) into n. The S^3
# addition theorem gives sum_modes Y Y* = (1/2 pi^2 a^3) n sin(n gamma)/sin(gamma), and at
# the antipode gamma -> pi that limit is (1/2 pi^2 a^3) n^2 (-1)^(n+1): the parity appears
# on its own rather than being inserted. So
#
#   G_img = sum_n (1/2 omega_n) e^{-i omega_n dt} (n^2/2 pi^2 a^3)(-1)^(n+1)
#         = (1/4 pi^2 a^2) sum_n n (-1)^(n+1) x^n,   x = e^{-i eta},  eta = dt/a
#         = (1/4 pi^2 a^2) x/(1+x)^2
#         = 1/(16 pi^2 a^2 cos^2(eta/2)).
#
# Everything below checks that chain and then differentiates it.

a <- 1
G_closed <- function(eta) 1 / (16 * pi^2 * a^2 * cos(eta / 2)^2)
G_sum <- function(eta, N, damp = 0) {
  n <- 1:N
  Re((1 / (4 * pi^2 * a^2)) * sum(n * (-1)^(n + 1) * exp(-1i * n * eta) * exp(-damp * n)))
}

cat("## SUPERSEDED IN PART by image_stress_conformal.R: the closed form below stands,\n")
cat("## the DENSITY in sections 3 and 4 is WITHDRAWN. Do not quote 3/4 pi^2.\n\n")
cat("=== 1. the mode sum reproduces the closed form ===\n")
cat("      eta       Abel-summed (damp 0.01)     closed form        difference\n")
for (e in c(0.4, 1.0, 2.0, 2.8, 3.0)) {
  s <- G_sum(e, 4000, 0.01); c0 <- G_closed(e)
  cat(sprintf("   %7.3f   %20.10f  %17.10f   %.2e\n", e, s, c0, abs(s - c0)))
}
cat("   (the damping is Abel summation; the series is conditionally convergent)\n")

cat("\n=== 2. the closed form diverges exactly at contact, eta = pi ===\n")
cat("   contact is |dt| = pi a, i.e. eta = pi, which is where the null geodesics refocus.\n")
cat("      pi - eta        G_img          times (pi-eta)^2\n")
for (d in c(1e-1, 1e-2, 1e-3, 1e-4)) {
  g <- G_closed(pi - d)
  cat(sprintf("   %9.0e   %14.4f   %18.8f\n", d, g, g * d^2))
}
cat("   G_img ~ 1/(4 pi^2 a^2 delta^2), so the two-point function goes as delta^-2 here,\n")
cat("   stronger than A.18's delta^-3/2 because on S^3 the antipode refocuses a\n")
cat("   two-parameter family of geodesics rather than a one-parameter one.\n")

cat("\n=== 3. the stress. G_img depends only on eta, so the spatial derivatives vanish ===\n")
cat("   The antipodal map on S^3 has differential -Id, so the fold's pullbacks are -1 on\n")
cat("   time and -1 on space. F_tt = +d^2_dt G and F_spatial = -grad.grad' G = 0, because\n")
cat("   S^3 is homogeneous and the antipodal correlator is the same at every point. So\n")
cat("   rho_img = (1/2) d^2 G_img / d(dt)^2 = (1/2a^2) G''(eta), with no cancellation\n")
cat("   available: there is nothing for the spatial term to cancel against.\n\n")
Gpp <- function(eta) { h <- 1e-5; (G_closed(eta + h) - 2 * G_closed(eta) + G_closed(eta - h)) / h^2 }
rho <- function(eta) 0.5 * Gpp(eta) / a^2
cat("      pi - eta        rho_img         times (pi-eta)^4      sign\n")
for (d in c(1e-1, 3e-2, 1e-2, 3e-3, 1e-3)) {
  r <- rho(pi - d)
  cat(sprintf("   %9.0e   %16.4f   %18.6f        %s\n", d, r, r * d^4, ifelse(r > 0, "+", "-")))
}

cat("\n=== 4. the closed form for the divergence, checked against the numbers above ===\n")
cat("   Near eta = pi, cos(eta/2) = -sin(delta/2) ~ -delta/2, so G ~ 1/(4 pi^2 a^2 delta^2)\n")
cat("   and G'' ~ 6/(4 pi^2 a^2 delta^4), giving rho ~ 3/(4 pi^2 a^2 delta^4).\n")
cat(sprintf("   predicted coefficient 3/(4 pi^2) = %.8f\n", 3 / (4 * pi^2)))
cat(sprintf("   measured at delta = 1e-3:         %.8f\n", rho(pi - 1e-3) * 1e-12))

cat("\n=== 5. the plant: remove the antipodal parity and the divergence must move ===\n")
G_nopar <- function(eta, N = 4000, damp = 0.01) {
  n <- 1:N
  Re((1 / (4 * pi^2 * a^2)) * sum(n * exp(-1i * n * eta) * exp(-damp * n)))
}
cat("      eta        with parity      without parity\n")
for (e in c(pi - 0.1, pi - 0.01, 0.1, 0.01)) {
  cat(sprintf("   %8.4f   %15.4f   %15.4f\n", e, G_closed(e), G_nopar(e)))
}
cat("   without the parity the singularity sits at eta = 0, which is coincidence, and the\n")
cat("   antipode is regular. The parity is exactly what moves the divergence from the\n")
cat("   coincidence point to the antipode. That is the fold doing the work.\n")

cat("\n=== SUPERSEDED IN PART ===\n")
cat("   The closed form above stands. The density in sections 3 and 4 is WITHDRAWN;\n")
cat("   see checks/calc/image_stress_conformal.R, which restores the mixed transverse\n")
cat("   derivative and finds the conformal image stress to be identically zero.\n")
