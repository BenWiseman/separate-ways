# Does the image stress survive the whole stress tensor? A check on GR52's own coefficient.
#
# GR52 computed the image energy density on the Einstein static universe as
# rho = (1/2) d^2 G_img/dt^2, dropping the spatial term with the argument that "S^3 is
# homogeneous and the antipodal correlator is the same at every point". Homogeneity does
# say something, but what it says is that the TOTAL derivative along the diagonal vanishes,
#     (d/dx)^2 G(x, Theta x) = 0   i.e.   G_aa + 2 G_ab + G_bb = 0,
# which is a statement about the Laplacian of phi^2 and not about the mixed derivative
# G_ab on its own. This script computes every term instead of arguing about one of them,
# validates the operator against a case whose answer is in the textbooks, and then applies
# it to the same geometry GR52 used.

N <- 1 / (4 * pi^2)

cat("=== 1. the ESU two-point function in closed form, at GENERAL angle ===\n")
# G = sum_n (1/2 omega_n) e^{-i omega_n dt} sum_modes Y Y*, conformal coupling, a = 1,
# with the S^3 addition theorem sum_modes Y Y* = n sin(n gamma)/(2 pi^2 sin gamma). Summing
# the geometric series gives 1/(8 pi^2 (cos tau - cos gamma)). GR52 only ever checked the
# gamma = pi slice of this, which is the slice where the answer it wanted lives.
G_cf  <- function(tau, gam) 1 / (8 * pi^2 * (cos(tau) - cos(gam)))
G_sum <- function(tau, gam, Nn = 8000, damp = 0.01) {
  n <- 1:Nn
  Re((1 / (4 * pi^2)) * sum(exp(-1i * n * tau) * sin(n * gam) / sin(gam) * exp(-damp * n)))
}
cat("      tau     gamma      Abel sum        closed form       difference\n")
worst <- 0
for (p in list(c(0.5, 1.0), c(1.3, 2.0), c(2.0, 2.6), c(0.7, 3.0), c(2.5, 0.9))) {
  s <- G_sum(p[1], p[2]); c0 <- G_cf(p[1], p[2]); worst <- max(worst, abs(s - c0))
  cat(sprintf("   %7.3f  %7.3f  %14.8f  %16.8f    %.2e\n", p[1], p[2], s, c0, abs(s - c0)))
}
cat(sprintf("   worst difference over the 5 angles: %.1e (the Abel damping is 0.01)\n", worst))
cat(sprintf("   coincidence limit check: near tau=gamma=0 it must go to 1/(4 pi^2 s^2).\n"))
e <- 1e-4
cat(sprintf("   G(0, %g) * (4 pi^2 %g^2) = %.8f  (must be 1)\n", e, e, G_cf(0, e) * 4 * pi^2 * e^2))
stopifnot(abs(G_cf(0, e) * 4 * pi^2 * e^2 - 1) < 1e-7)

cat("\n=== 2. the stress operator, and the form that is actually traceless ===\n")
cat("   T_ab = grad_a phi grad_b phi - (1/2) g_ab (grad phi)^2\n")
cat("          + xi ( g_ab box - grad_a grad_b + G_ab ) phi^2\n")
cat("   Its trace on a solution of box phi = xi R phi is (6 xi - 1)[ ... ], zero at xi = 1/6,\n")
cat("   which is the sign convention fixed rather than assumed. Contracting the t index and\n")
cat("   using g_tt = -1 the second time derivatives cancel and it collapses to\n")
cat("        T_tt = (1/2)[ (d_t phi)^2 + |grad_s phi|^2 ] - xi Lap_s(phi^2) + xi G_tt phi^2,\n")
cat("   whose first bracket alone is what GR52 used.\n")

cat("\n=== 3. the validation: a single Dirichlet mirror in flat space ===\n")
cat("   Textbook: a conformally coupled scalar feels NO stress from one flat plate, while a\n")
cat("   minimally coupled one feels -1/(16 pi^2 z^4). Same operator, same image structure.\n")
mirror <- function(z, xi) {
  D0 <- 4 * z^2
  Gtxty <-  2 * N / D0^2            # d_tx d_ty of -1/(4 pi^2 D)
  Gxxu  <- -2 * N / D0^2            # transverse, mixed
  Gxxxx <-  2 * N / D0^2            # transverse, same slot
  Gzxzy <- -N * (-2 / D0^2 + 2 * (4 * z)^2 / D0^3)
  Gzxzx <- Gzxzy                    # D_zxzx = D_zxzy = 2 and D_zx = D_zy = 4z
  min_part <- 0.5 * (Gtxty + (2 * Gxxu + Gzxzy))
  lap_phi2 <- 2 * (Gxxxx + 2 * Gxxu + Gxxxx) + (Gzxzx + 2 * Gzxzy + Gzxzx)
  c(minimal = min_part, conformal = min_part - xi * lap_phi2)
}
for (z in c(0.5, 1, 2)) {
  m <- mirror(z, 1/6)
  cat(sprintf("   z = %4.1f   minimal %+.8f  (-1/16 pi^2 z^4 = %+.8f)   conformal %+.3e\n",
              z, m[1], -1 / (16 * pi^2 * z^4), m[2]))
  stopifnot(abs(m[1] + 1 / (16 * pi^2 * z^4)) < 1e-12, abs(m[2]) < 1e-14)
}
cat("   Both reproduced. The operator can return a nonzero answer and does, and the zero it\n")
cat("   returns at xi = 1/6 is the one the literature returns.\n")

cat("\n=== 4. the same operator on GR52's configuration ===\n")
# Image pair: x and Theta x antipodal on S^3, separated in time by tau. Writing the image
# correlator as a function of y, D = cos tau + cos psi with psi the angle between n_x and
# n_y, because the antipodal map sends the angle gamma to pi - psi.
esu <- function(tau, xi) {
  c0 <- cos(tau); s0 <- sin(tau); D0 <- 1 + c0
  P  <- (1 / (8 * pi^2)) * ( c0 / D0^2 + 2 * s0^2 / D0^3)   # d_tx d_ty  (= d_tx d_tx here)
  Q  <- 1 / (8 * pi^2 * D0^2)
  Gmixed_s <- -Q                                            # d_alpha^i d_beta^j / delta^ij
  Gsame_s  <- +Q                                            # d_alpha^i d_alpha^j / delta^ij
  G0 <- 1 / (8 * pi^2 * D0)
  min_part <- 0.5 * (P + 3 * Gmixed_s)
  lap_phi2 <- 3 * (Gsame_s + 2 * Gmixed_s + Gsame_s)
  Gtt_einstein <- 3                                         # R_tt = 0, R = 6, g_tt = -1
  c(minimal = min_part,
    conformal = min_part - xi * lap_phi2 + xi * Gtt_einstein * G0,
    laplacian = lap_phi2, GR52 = 0.5 * P)
}
cat("   The three columns cancel against each other, so the conformal one is reported as a\n")
cat("   FRACTION of the largest term it is built from; below 1e-2 the double-precision\n")
cat("   cancellation is what limits it, and section 5 does the algebra instead.\n")
cat("     pi - tau      GR52's rho        minimal (full)   conformal / largest term\n")
for (d in c(1e-1, 1e-2, 1e-3)) {
  v <- esu(pi - d, 1/6); scale <- 0.5 * abs((1/(8*pi^2)) * 2*sin(pi-d)^2/(1+cos(pi-d))^3)
  cat(sprintf("   %9.0e  %16.4f  %18.6f  %+18.3e\n", d, v["GR52"], v["minimal"], v["conformal"]/scale))
}
cat("\n   GR52's column diverges as delta^-4. The full minimal density is two powers weaker\n")
cat("   and NEGATIVE, and the full conformal density is zero to machine precision.\n")
cat("   The Laplacian term that homogeneity really constrains:\n")
for (d in c(1e-1, 1e-3)) cat(sprintf("      Lap_s(phi^2) at delta = %.0e : %+.3e\n", d, esu(pi-d,1/6)["laplacian"]))
cat("   It vanishes identically, which is the true content of the homogeneity argument.\n")

cat("\n=== 5. both results in closed form ===\n")
cat("   Writing D0 = 1 + cos tau and using sin^2 = (1-cos)(1+cos),\n")
cat("     minimal:    (1/2)(P - 3Q) = -1/(16 pi^2 D0) = -1/(32 pi^2 cos^2(tau/2))\n")
cat("     conformal:  that plus (1/2)/(8 pi^2 D0), which cancels it exactly.\n")
for (tau in c(1.0, 2.0, 3.0, pi - 1e-2)) {
  v <- esu(tau, 1/6)
  cat(sprintf("   tau = %8.5f   minimal %+14.6f   closed form %+14.6f   conformal %+.2e\n",
              tau, v["minimal"], -1 / (32 * pi^2 * cos(tau/2)^2), v["conformal"]))
  scale <- abs((1/(8*pi^2)) * 2*sin(tau)^2/(1+cos(tau))^3)
  stopifnot(abs(v["minimal"] + 1/(32*pi^2*cos(tau/2)^2)) < 1e-11 * max(1, scale))
  stopifnot(abs(v["conformal"]) < 1e-11 * max(1, scale))
}

cat("\n=== 6. the null-null component, which is what feeds Raychaudhuri ===\n")
cat("   T_kk = (1-2xi)(k.grad phi)^2 - 2 xi phi (k.grad)^2 phi + xi G_kk phi^2, the g_ab\n")
cat("   terms dropping because k is null. For the ESU G_kk = R_kk = 2 (k^t)^2 / a^2.\n")
Tkk <- function(tau, xi) {
  c0 <- cos(tau); s0 <- sin(tau); D0 <- 1 + c0
  P <- (1 / (8 * pi^2)) * (c0 / D0^2 + 2 * s0^2 / D0^3); Q <- 1 / (8 * pi^2 * D0^2)
  A <- P - Q; B <- P + Q; G0 <- 1 / (8 * pi^2 * D0)
  (1 - 2 * xi) * A - 2 * xi * B + xi * 2 * G0
}
cat("      pi - tau     T_kk conformal      T_kk minimal    conformal / largest term\n")
for (d in c(1e-1, 1e-2, 1e-3)) {
  sc <- abs((1/(8*pi^2)) * 2*sin(pi-d)^2/(1+cos(pi-d))^3)
  cat(sprintf("   %9.0e   %+16.3e   %+15.4f   %+18.2e\n",
              d, Tkk(pi - d, 1/6), Tkk(pi - d, 0), Tkk(pi - d, 1/6)/sc))
  stopifnot(abs(Tkk(pi - d, 1/6)) < 1e-11 * sc)
}
cat("   Conformal: zero. Minimal: nonzero, POSITIVE and delta^-4, in closed form\n")
cat("        T_kk = (1 - cos tau) / (8 pi^2 (1 + cos tau)^2),\n")
cat("   so the null-null component keeps the fourth power the energy density loses. Its\n")
cat("   sign says a non-conformal image term FOCUSES the congruence whose refocusing is\n")
cat("   what creates the contact, up to the overall parity of the image term, which the\n")
cat("   fold fixes and which is the next thing to settle: the even and odd sectors carry\n")
cat("   opposite signs and the Keldysh split puts one field in each.\n")
for (d in c(1e-1, 1e-2)) {
  cf <- (1 - cos(pi - d)) / (8 * pi^2 * (1 + cos(pi - d))^2)
  cat(sprintf("     closed form at delta = %.0e: %+14.4f  against %+14.4f\n", d, cf, Tkk(pi-d, 0)))
  stopifnot(abs(cf - Tkk(pi - d, 0)) < 1e-6 * abs(cf))
}

cat("\n=== 7. the plant: move xi off 1/6 and the zero must move with it ===\n")
for (xi in c(0, 1/12, 1/6, 1/4, 1/3)) {
  v <- esu(pi - 1e-2, xi)
  cat(sprintf("   xi = %7.4f   rho = %+16.4f   T_kk = %+16.4f\n", xi, v["conformal"], Tkk(pi-1e-2, xi)))
}
cat("   Only xi = 1/6 gives zero, in both. The cancellation is conformal invariance and\n")
cat("   not an operator that returns zero whatever it is handed.\n")

cat("\n=== 8. what this means ===\n")
cat("   GR52's delta^-4 came entirely from dropping the mixed spatial derivative. Restored,\n")
cat("   it cancels the time term at leading order, and the conformal improvement cancels\n")
cat("   what is left. The image stress of a CONFORMALLY INVARIANT field vanishes, on the\n")
cat("   ESU exactly and at a flat mirror exactly, and A.19's constant-radius zero is the\n")
cat("   third instance of the same thing rather than an artefact of a product geometry.\n")
cat("   So the fold's image stress is sourced by the BREAKING of conformal invariance. That\n")
cat("   is why A.18's model gives a divergence: its Kaluza-Klein tower has masses\n")
cat("   m_l = sqrt(l(l+1))/a, and a massive field is not conformally invariant. A.15's\n")
cat("   self-censoring stress is therefore proportional to mass and to the trace anomaly,\n")
cat("   which is a statement about matter content rather than about geometry alone.\n")
cat("   GR57's geometry is untouched by any of this: contact and conjugacy are the same\n")
cat("   condition whatever field is put on the geometry, and the caustic order is fixed by\n")
cat("   the Ricci term. What changes is which fields feel it.\n")

cat("\n=== 9. the same numbers again, from finite differences on the exact embedding ===\n")
cat("   Everything above used second derivatives worked out by hand from D = cos tau + cos psi.\n")
cat("   This redoes them with no hand algebra: points on S^3 as unit 4-vectors, the exponential\n")
cat("   map for the displacements, and central differences.\n")
pt <- function(al) { n <- sqrt(sum(al^2)); if (n < 1e-14) c(1,0,0,0) else c(cos(n), sin(n)*al/n) }
Gfun <- function(tx, al, ty, be, Tt, timeflip = TRUE) {
  tau <- if (timeflip) Tt - tx - ty else ty - tx
  cg  <- -sum(pt(al) * pt(be))               # antipodal image: cos gamma = -cos psi
  1 / (8 * pi^2 * (cos(tau) - cg))
}
d2 <- function(f, i, j, h = 1e-4) {          # central second difference in coordinates i and j
  e <- function(k) { v <- numeric(8); v[k] <- 1; v }
  if (i == j) (f(h*e(i)) - 2*f(numeric(8)) + f(-h*e(i))) / h^2
  else (f(h*e(i)+h*e(j)) - f(h*e(i)-h*e(j)) - f(-h*e(i)+h*e(j)) + f(-h*e(i)-h*e(j))) / (4*h^2)
}
esu_fd <- function(tau, xi, timeflip = TRUE) {
  Tt <- tau; tx0 <- 0
  f <- function(v) Gfun(tx0 + v[1], v[2:4], v[5], v[6:8], Tt, timeflip)
  #    index 1 = t_x, 2:4 = alpha, 5 = t_y, 6:8 = beta
  Att <- d2(f, 1, 5); Asp <- sum(sapply(1:3, function(i) d2(f, 1 + i, 5 + i)))
  lap <- sum(sapply(1:3, function(i) d2(f, 1+i, 1+i) + 2*d2(f, 1+i, 5+i) + d2(f, 5+i, 5+i)))
  G0  <- f(numeric(8))
  c(minimal = 0.5 * (Att + Asp), conformal = 0.5 * (Att + Asp) - xi * lap + xi * 3 * G0)
}
cat("      tau      minimal (hand)   minimal (fd)     conformal (hand)  conformal (fd)\n")
for (tau in c(1.0, 2.0, 2.8)) {
  h <- esu(tau, 1/6); f <- esu_fd(tau, 1/6)
  cat(sprintf("   %7.3f  %15.6f  %14.6f  %17.2e  %14.2e\n",
              tau, h["minimal"], f["minimal"], h["conformal"], f["conformal"]))
  stopifnot(abs(h["minimal"] - f["minimal"]) < 1e-5 * max(1, abs(h["minimal"])))
  stopifnot(abs(f["conformal"]) < 1e-5 * max(1, abs(h["minimal"])))
}
cat("   The hand algebra and the finite differences agree, so the cancellation is in the\n")
cat("   physics and not in a mis-differentiated cosine.\n")
cat("   And without the fold's time reflection, which is a different map:\n")
for (tau in c(1.0, 2.0)) {
  f <- esu_fd(tau, 1/6, timeflip = FALSE)
  cat(sprintf("      tau = %.1f   minimal %+12.6f   conformal %+12.6f\n", tau, f[1], f[2]))
}

cat("\n=== 10. the flat case of this appendix, run through the same operator ===\n")
cat("   Theta(t,x) = (-t,-x) is an inversion through a point, not a reflection in a plane,\n")
cat("   and it is the fold at a Rindler horizon. At t = 0 the image pair is separated by\n")
cat("   2R and D = 4R^2, with every second derivative of D equal in the two slots.\n")
flat_inv <- function(R, xi) {
  D0 <- 4 * R^2
  Gtt <- N * (2 / D0^2)                                 # D_titj = -2, D_ti = 0 at t = 0
  Gsp <- function() {                                    # sum_i over the three directions
    s <- 0
    for (i in 1:3) {
      Di <- if (i == 3) 4 * R else 0                     # R along the third axis
      s <- s + N * (-2 / D0^2 + 2 * Di^2 / D0^3)
    }
    s
  }
  A <- Gsp(); lap <- 4 * A                               # all four blocks are equal here
  c(minimal = 0.5 * (Gtt + A), conformal = 0.5 * (Gtt + A) - xi * lap)
}
for (R in c(0.5, 1, 2)) {
  v <- flat_inv(R, 1/6)
  cat(sprintf("   R = %4.1f   minimal %+12.8f  (1/32 pi^2 R^4 = %+.8f)   conformal %+12.8f  (1/96 pi^2 R^4 = %+.8f)\n",
              R, v[1], 1/(32*pi^2*R^4), v[2], 1/(96*pi^2*R^4)))
  stopifnot(abs(v[1] - 1/(32*pi^2*R^4)) < 1e-12, abs(v[2] - 1/(96*pi^2*R^4)) < 1e-12)
}
cat("   NONZERO and positive, at both couplings. So this appendix's flat result survives the\n")
cat("   full operator; it is the Einstein static universe's that does not. An inversion\n")
cat("   through a point and a reflection in a plane are not the same kind of fold, and a\n")
cat("   conformal field can tell them apart.\n")

cat("\n=== 11. where that leaves the three geometries ===\n")
cat("   flat inversion, the fold at a Rindler horizon:   +1/(96 pi^2 R^4), conformal\n")
cat("   flat single mirror:                               0, conformal (textbook)\n")
cat("   Einstein static universe, antipodal fold:         0, conformal; -(1/2) G_img, minimal\n")
cat("   constant-radius model of this appendix:           0\n")
cat("   The +1-eigenvalue sorting put the Einstein static universe with the flat case. It\n")
cat("   does not belong there, so that sorting does not survive, and with it the inference\n")
cat("   that a black hole's one +1 makes its image stress positive.\n")

cat("\n=== 12. the coupling dependence in closed form, which is the real statement ===\n")
cat("   With Lap_s(phi^2) = 0 identically, rho(xi) = -1/(16 pi^2 D0) + 3 xi/(8 pi^2 D0), so\n")
cat("        rho_img  = (6 xi - 1) / (16 pi^2 a^2 (1 + cos eta))\n")
cat("        T_kk     = (1 - 6 xi)(1 - cos eta) / (8 pi^2 a^4 (1 + cos eta)^2)\n")
cat("   Both carry the factor (1 - 6 xi) and nothing else. The image stress is not merely\n")
cat("   zero at conformal coupling, it is PROPORTIONAL to the departure from it.\n")
cat("      xi        rho (assembled)   rho (closed form)    T_kk (assembled)   T_kk (closed)\n")
for (xi in c(0, 1/12, 1/6, 1/4, 1/3, 1/2)) {
  eta <- pi - 1e-2; D0 <- 1 + cos(eta)
  ra <- esu(eta, xi)["conformal"]; rc <- (6 * xi - 1) / (16 * pi^2 * D0)
  ka <- Tkk(eta, xi);              kc <- (1 - 6 * xi) * (1 - cos(eta)) / (8 * pi^2 * D0^2)
  cat(sprintf("   %7.4f  %16.4f  %18.4f  %18.2f  %15.2f\n", xi, ra, rc, ka, kc))
  sc <- abs((1/(8*pi^2)) * 2*sin(eta)^2/D0^3)
  stopifnot(abs(ra - rc) < 1e-10 * sc, abs(ka - kc) < 1e-10 * sc)
}
cat("   Both closed forms reproduce the assembled operator at six couplings including two\n")
cat("   outside [0, 1/6], so the linearity is exact rather than a local slope.\n")
cat("   The Hadamard expansion says why: the coefficient of the log term is\n")
cat("   V_0 = Delta^{1/2}[m^2 + (xi - 1/6) R]/2, which is the only place a massless\n")
cat("   conformal field in a Ricci-flat spacetime has nothing to put. R = 6/a^2 here, so\n")
cat("   (xi - 1/6) R is the whole of the breaking and (1 - 6 xi) is its natural measure.\n")
cat("   The same slot carries m^2, so a massive field should feel the image term whatever\n")
cat("   its coupling. That is the generalisation this script does not compute, and it is\n")
cat("   the one A.18 already exhibits with its tower of masses sqrt(l(l+1))/a.\n")
