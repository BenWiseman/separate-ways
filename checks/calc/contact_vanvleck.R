# The magnitude at a black hole: what the contact geodesic's focusing allows.
#
# GR52 gives an exact image stress at a caustic, rho ~ (3/4pi^2) delta^-4, in the Einstein
# static universe. GR56 shows a Schwarzschild contact curve carries the ESU's signature, so
# the SIGN transfers. The magnitude does not: a signature argument carries a sign, not a
# number. What sets the number is the Hadamard amplitude,
#
#     G_img(x, Theta x) = Delta^{1/2}(x, Theta x) / (4 pi^2 sigma) + regular,
#
# whose only curvature input is the Van Vleck determinant Delta of the connecting geodesic.
# Delta = 1 in flat space and diverges at a conjugate point, so the exponent of the image
# stress is set by how many transverse directions refocus at contact. This computes it.
#
# THE CURVE IN CLOSED FORM. Inside the horizon a null geodesic has
# rdot^2 = E^2 + L^2 (2M/r - 1)/r^2 > 0 for E =/= 0, so r is MONOTONIC unless E = 0. The fold
# preserves r (it preserves UV), so x and Theta x sit at the same radius, and only an E = 0
# geodesic can return to the radius it left. That makes the connecting geodesic unique, and
# with dr/dphi = +- sqrt(r(2M - r)) it integrates to
#
#     r(phi) = M (1 + sin phi),   phi in [0, pi],
#
# which starts and ends at r = M, grazes r = 2M at phi = pi/2, and turns through exactly pi.
# So the contact boundary r = M is not an input: it is where the geodesic that returns to its
# own radius turns through the angle the antipode needs.

M <- 1; L <- 1

r_of_phi  <- function(p) M * (1 + sin(p))
rp_of_phi <- function(p) M * cos(p)

cat("=== 1. the closed-form curve solves the geodesic system ===\n")
cat("   first order: dr/dphi must equal +- sqrt(r(2M - r))\n")
worst <- 0
for (p in seq(0.01, pi - 0.01, length.out = 400)) {
  r <- r_of_phi(p); pred <- sqrt(r * (2 * M - r)) * sign(cos(p))
  worst <- max(worst, abs(rp_of_phi(p) - pred))
}
cat(sprintf("   worst residual over 400 points: %.3e\n", worst))
stopifnot(worst < 1e-12)
cat(sprintf("   endpoints  r(0) = %.6f   r(pi) = %.6f   (both M)\n", r_of_phi(0), r_of_phi(pi)))
cat(sprintf("   apex       r(pi/2) = %.6f  (the bifurcation surface 2M)\n", r_of_phi(pi/2)))
lam_tot <- 3 * pi / 2 + 4
lam_num <- integrate(function(p) r_of_phi(p)^2 / L, 0, pi, rel.tol = 1e-12)$value
cat(sprintf("   affine length  int r^2 dphi / L = %.9f, closed form 3pi/2 + 4 = %.9f\n",
            lam_num, lam_tot))
stopifnot(abs(lam_num - lam_tot) < 1e-9)

# ---------------------------------------------------------------- numerical Riemann tensor
f_of <- function(r) 1 - 2 * M / r
gmet <- function(r, th) diag(c(-f_of(r), 1 / f_of(r), r^2, r^2 * sin(th)^2))

christ <- function(r, th, h = 1e-5) {
  g <- gmet(r, th); gi <- solve(g)
  dg <- array(0, c(4, 4, 4))                       # dg[c,a,b] = d_c g_ab
  dg[2, , ] <- (gmet(r + h, th) - gmet(r - h, th)) / (2 * h)
  dg[3, , ] <- (gmet(r, th + h) - gmet(r, th - h)) / (2 * h)
  G <- array(0, c(4, 4, 4))
  for (a in 1:4) for (b in 1:4) for (c in 1:4) {
    s <- 0
    for (d in 1:4) s <- s + gi[a, d] * (dg[b, d, c] + dg[c, d, b] - dg[d, b, c])
    G[a, b, c] <- 0.5 * s
  }
  G
}

riemann_dn <- function(r, th, h = 1e-4) {
  G  <- christ(r, th)
  dG <- array(0, c(4, 4, 4, 4))                    # dG[c,a,d,b] = d_c Gamma^a_db
  dG[2, , , ] <- (christ(r + h, th) - christ(r - h, th)) / (2 * h)
  dG[3, , , ] <- (christ(r, th + h) - christ(r, th - h)) / (2 * h)
  Rud <- array(0, c(4, 4, 4, 4))                   # R^a_bcd
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4) {
    s <- dG[cc, a, dd, b] - dG[dd, a, cc, b]
    for (e in 1:4) s <- s + G[a, cc, e] * G[e, dd, b] - G[a, dd, e] * G[e, cc, b]
    Rud[a, b, cc, dd] <- s
  }
  g <- gmet(r, th); Rdn <- array(0, c(4, 4, 4, 4))
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4) {
    Rdn[a, b, cc, dd] <- sum(g[a, ] * Rud[, b, cc, dd])
  }
  list(dn = Rdn, ud = Rud, g = g)
}

cat("\n=== 2. the numerical Riemann tensor, validated before it is used ===\n")
cat("      r      max |Ricci|     Kretschmann    48 M^2/r^6     rel. err\n")
for (r in c(1.0, 1.3, 1.7, 1.95)) {
  R <- riemann_dn(r, pi / 2); gi <- solve(R$g)
  Ric <- matrix(0, 4, 4)
  for (b in 1:4) for (dd in 1:4) Ric[b, dd] <- sum(sapply(1:4, function(a) R$ud[a, b, a, dd]))
  Kre <- 0
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4) {
    up <- 0
    for (p in 1:4) for (q in 1:4) for (s in 1:4) for (t in 1:4)
      up <- up + gi[a,p]*gi[b,q]*gi[cc,s]*gi[dd,t]*R$dn[p,q,s,t]
    Kre <- Kre + R$dn[a,b,cc,dd] * up
  }
  cat(sprintf("   %5.2f   %.3e   %13.6f  %12.6f    %.2e\n",
              r, max(abs(Ric)), Kre, 48 * M^2 / r^6, abs(Kre - 48 * M^2 / r^6) / (48 / r^6)))
}
cat("   Ricci-flat to finite-difference precision and the Kretschmann scalar is the\n")
cat("   textbook 48 M^2/r^6, so the tensor is the Schwarzschild one and not an artefact.\n")

# ---------------------------------------------------------- the transverse tidal matrix
# Both transverse directions are exactly parallel-propagated along an E = 0 null geodesic:
#   grad_k (partial_theta / r) = 0  and  grad_k (partial_t / sqrt(F)) = 0,
# because grad_k partial_theta = (k^r/r) partial_theta and grad_k partial_t = (k^r f'/2f)
# partial_t, and the normalisations cancel each rescaling exactly. So no frame ODE is needed.
tidal <- function(r) {
  R <- riemann_dn(r, pi / 2)
  Fr <- 2 * M / r - 1
  kr <- sqrt((2 * M - r) / r^3) * L; kf <- L / r^2
  k  <- c(0, kr, 0, kf)
  e_th <- c(0, 0, 1 / r, 0); e_t <- c(1 / sqrt(Fr), 0, 0, 0)
  qform <- function(e) {
    s <- 0
    for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4)
      s <- s + R$dn[a, b, cc, dd] * e[a] * k[b] * e[cc] * k[dd]
    s
  }
  c(theta = qform(e_th), time = qform(e_t))
}

cat("\n=== 3. the two transverse tidal eigenvalues, and their sum ===\n")
cat("      r      T(theta)       T(t)         sum        3 M L^2 / r^5\n")
for (r in c(1.0, 1.2, 1.5, 1.8, 1.99)) {
  tt <- tidal(r)
  cat(sprintf("   %5.2f  %+11.7f  %+11.7f  %+10.2e   %11.7f\n",
              r, tt[1], tt[2], tt[1] + tt[2], 3 * M * L^2 / r^5))
}
cat("   The sum vanishes, and it must: expanding the identity in a null frame gives\n")
cat("   R_ab k^a k^b = sum over the two transverse directions of the tidal eigenvalue,\n")
cat("   the k-n terms dropping by the antisymmetry of Riemann. Schwarzschild is vacuum,\n")
cat("   so the transverse tidal matrix is TRACELESS: diag(+Psi, -Psi).\n")
cat("   Psi = 3 M L^2 / r^5, matched above to the numerical tensor.\n")

psi_num <- sapply(c(1.0, 1.2, 1.5, 1.8, 1.99), function(r) abs(tidal(r)[1]))
psi_cf  <- 3 * M * L^2 / c(1.0, 1.2, 1.5, 1.8, 1.99)^5
cat(sprintf("   worst relative departure from 3 M L^2/r^5: %.2e\n", max(abs(psi_num / psi_cf - 1))))
stopifnot(max(abs(psi_num / psi_cf - 1)) < 1e-5)

# ------------------------------------------------------------------- the Jacobi integration
# d^2 J/dlambda^2 = -T J with T = diag(+Psi, -Psi). In phi, with dlambda/dphi = r^2/L,
#     J'' - (2 r'/r) J' +- 3 M J / r = 0,
# upper sign focusing, lower defocusing. J(0) = 0, dJ/dlambda(0) = 1 so J'(0) = r(0)^2/L.
jacobi <- function(sgn, n = 400000, psi_scale = 1) {
  h <- pi / n; y <- c(0, r_of_phi(0)^2 / L); out <- numeric(n + 1); out[1] <- 0
  deriv <- function(p, y) {
    r <- r_of_phi(p)
    c(y[2], (2 * rp_of_phi(p) / r) * y[2] - sgn * psi_scale * 3 * M * y[1] / r)
  }
  for (i in 1:n) {
    p <- (i - 1) * h
    k1 <- deriv(p, y); k2 <- deriv(p + h/2, y + h/2 * k1)
    k3 <- deriv(p + h/2, y + h/2 * k2); k4 <- deriv(p + h, y + h * k3)
    y <- y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    out[i + 1] <- y[1]
  }
  list(end = y[1], trace = out, phi = seq(0, pi, length.out = n + 1))
}

cat("\n=== 4. the integrator, checked on the case with a known answer ===\n")
flat <- jacobi(+1, psi_scale = 0)
cat(sprintf("   curvature switched off: J(pi) = %.9f, and the affine length is %.9f\n",
            flat$end, lam_tot))
cat(sprintf("   these must agree (J = lambda when nothing bends it); difference %.2e\n",
            abs(flat$end - lam_tot)))
stopifnot(abs(flat$end - lam_tot) < 1e-7)

cat("\n=== 5. both branches have an exact solution, and both come from a Killing vector ===\n")
cat("   A Killing field restricted to a geodesic is a Jacobi field, and this geodesic sits\n")
cat("   in a plane of a spherically symmetric spacetime, so two of them are transverse.\n\n")
cat("   FOCUSING (the theta direction). The rotation whose axis runs through the point at\n")
cat("   phi = 0 is xi = -sin(phi) d_theta at theta = pi/2. In the parallel frame d_theta/r\n")
cat("   its amplitude is r sin(phi), so\n")
cat("        J1(phi) = M (1 + sin phi) sin phi,\n")
cat("   which already has J1(0) = 0 and dJ1/dlambda(0) = 1, the Jacobi conditions.\n")
J1  <- function(p) r_of_phi(p) * sin(p)
J1d <- function(p) M * cos(p) * (1 + 2 * sin(p))
J1dd<- function(p) M * (2 - sin(p) - 4 * sin(p)^2)
res <- function(J, Jd, Jdd, sgn, p) {
  Jdd(p) - (2 * rp_of_phi(p) / r_of_phi(p)) * Jd(p) + sgn * 3 * M * J(p) / r_of_phi(p)
}
w1 <- max(abs(sapply(seq(0.05, pi-0.05, length.out=300),
                     function(p) res(J1, J1d, J1dd, +1, p))))
cat(sprintf("   residual in the focusing equation over 300 points: %.2e (analytic derivatives)\n", w1))
stopifnot(w1 < 1e-12)

cat("\n   DEFOCUSING (the t direction). d_t is Killing, is orthogonal to k because E = 0, and\n")
cat("   has no component along k, so it is purely transverse with amplitude sqrt(F). Its\n")
cat("   smooth continuation through r = 2M, where it changes sign, is\n")
cat("        y(phi) = tan(pi/4 - phi/2),   y^2 = F = 2M/r - 1.\n")
y1  <- function(p) tan(pi/4 - p/2)
y1d <- function(p) -0.5 / cos(pi/4 - p/2)^2
y1dd<- function(p) 0.5 * tan(pi/4 - p/2) / cos(pi/4 - p/2)^2
w2 <- max(abs(sapply(seq(0.05, pi-0.05, length.out=300),
                     function(p) res(y1, y1d, y1dd, -1, p))))
cat(sprintf("   residual in the defocusing equation over 300 points: %.2e\n", w2))
cat(sprintf("   and y^2 - F at phi = 0.7: %.2e\n", y1(0.7)^2 - (2*M/r_of_phi(0.7) - 1)))
stopifnot(w2 < 1e-12)
cat("   This one does not satisfy J(0) = 0, so it is the other solution of that branch;\n")
cat("   what it does is certify the equation, and with it the tidal 3 M L^2 / r^5 that the\n")
cat("   numerical Riemann tensor produced. Two independent exact solutions, one per branch.\n")

cat("\n=== 6. the integrator against the exact solutions ===\n")
foc <- jacobi(+1); def <- jacobi(-1)
cat(sprintf("   focusing:   J1(pi) integrated = %+.3e, exact M(1+sin pi)sin pi = 0\n", foc$end))
err <- max(abs(foc$trace - J1(foc$phi)))
cat(sprintf("   worst departure from the closed form along the whole path: %.2e\n", err))
stopifnot(err < 1e-8)
cat(sprintf("   defocusing: J2(pi) = %.6f\n", def$end))

cat("\n=== 7. the contact point IS a conjugate point, and that is not a coincidence ===\n")
cat("   J1 vanishes at phi = 0 and at phi = pi, exactly, because a rotation vanishes on its\n")
cat("   own axis and the axis through x also runs through the antipodal direction. So in ANY\n")
cat("   spherically symmetric spacetime, a geodesic that turns through exactly pi ends at a\n")
cat("   point conjugate to where it started. A.15's contact condition is that the geodesic\n")
cat("   turns through pi. CONTACT AND CONJUGACY ARE THE SAME CONDITION, not two conditions\n")
cat("   that happen to coincide at r = M.\n")
cat(sprintf("   (checked: the angle budget 2pi - 4 asin sqrt(r/2M) equals pi at r/M = %.6f)\n",
            uniroot(function(x) 2*pi - 4*asin(sqrt(x/2)) - pi, c(0.1, 1.9))$root))

cat("\n=== 8. the order of the caustic, which is what fixes the exponent ===\n")
cat("   In 4D vacuum the transverse tidal matrix is traceless, so one direction focuses and\n")
cat("   the other defocuses and det J can vanish through ONE factor only. A vacuum black\n")
cat("   hole cannot carry a point caustic. The Einstein static universe can and does: it is\n")
cat("   not vacuum, R_ab k^a k^b = (2/a^2)(k^t)^2 > 0, both transverse directions focus, and\n")
cat("   its antipode is conjugate of order two.\n\n")
slope <- J1d(pi) * L / r_of_phi(pi)^2       # dJ1/dlambda = (dJ1/dphi)(dphi/dlambda)
cat(sprintf("   order one here. Near the end J1 -> |dJ1/dlambda| (lambda_tot - lambda) with\n"))
cat(sprintf("   dJ1/dlambda at contact = %+.1f, which is -L/M exactly: J1d(pi) = -M and\n", slope))
cat("   dphi/dlambda = L/M^2 there. So\n")
amp <- lam_tot / (abs(slope) * def$end)
cat(sprintf("     det J -> %.6f (lambda_tot - lambda)\n", abs(slope) * def$end))
cat(sprintf("     Delta = lambda^2/det J -> %.6f lambda_tot/(lambda_tot - lambda)\n", amp))
cat(sprintf("     Delta^{1/2} -> %.6f sqrt(lambda_tot/(lambda_tot - lambda))\n", sqrt(amp)))
cat("   The bracket is dimensionless and invariant under rescaling the affine parameter,\n")
cat("   which L does, so it is a property of the curve and not of the units.\n")

cat("\n=== 9. what this says about the magnitude ===\n")
cat("   G_img = Delta^{1/2}/(4 pi^2 sigma). At an order-one caustic Delta^{1/2} ~ delta^-1/2\n")
cat("   and sigma ~ delta, so G ~ delta^-3/2 and the stress ~ delta^-7/2. At the ESU's\n")
cat("   order-two caustic Delta^{1/2} ~ delta^-1, so G ~ delta^-2 and the stress ~ delta^-4,\n")
cat("   which is GR52's exact answer. The two exponents are the two caustic orders, and\n")
cat("   vacuum picks the lower one.\n")
cat("   A.18 computed delta^-3/2 and the companion already says the ESU is stronger because\n")
cat("   S^3 refocuses a two-parameter family. That remark is now a derivation: the family is\n")
cat("   two-parameter exactly when R_ab k^a k^b > 0, and vacuum forbids it.\n")
cat("   So a black hole's contact surface carries A.18's exponent, GR56's sign, and an\n")
cat(sprintf("   amplitude Delta^{1/2} = %.6f sqrt(lambda_tot/(lambda_tot-lambda)) relative to\n", sqrt(amp)))
cat("   the flat Hadamard value. The delta^-3/2 is the standard fold-caustic exponent and\n")
cat("   assumes both sigma and the shortfall to the conjugate point vanish linearly in the\n")
cat("   distance off r = M, which is the generic case; pinning the COEFFICIENT needs the\n")
cat("   uniform Airy treatment the parametrix is replaced by at a caustic, not this ratio.\n")
cat("   Still open with it: the index contraction from G to T_ab.\n")

cat("\n=== 10. the plant: a check that cannot fail is not a check ===\n")
cat("   (a) break vacuum by hand and the traceless test must catch it.\n")
gmet_bad <- gmet
gmet <- function(r, th) { g <- gmet_bad(r, th); g[3, 3] <- g[3, 3] * 1.03; g }
tb <- tidal(1.5)
cat(sprintf("       fake metric: T(theta) + T(t) = %+.5f  (vacuum gives 0)\n", tb[1] + tb[2]))
stopifnot(abs(tb[1] + tb[2]) > 1e-3)
gmet <- gmet_bad
cat(sprintf("       restored:    T(theta) + T(t) = %+.2e\n", sum(tidal(1.5))))
cat("   (b) detune the tidal and the conjugate point at phi = pi must move off the end.\n")
for (s in c(0.5, 0.9, 1.0, 1.1, 2.0)) {
  fs <- jacobi(+1, psi_scale = s)
  cat(sprintf("       Psi x %.1f:  J1(pi) = %+10.6f\n", s, fs$end))
  if (abs(s - 1) < 1e-9) stopifnot(abs(fs$end) < 1e-8)
  if (abs(s - 1) > 0.05) stopifnot(abs(fs$end) > 0.05)
}
cat("   Only the true tidal strength lands the zero on the endpoint. The vanishing is a\n")
cat("   property of the Schwarzschild curvature, not of an integrator that returns zero.\n")
