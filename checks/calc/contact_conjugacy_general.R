# Contact and conjugacy at any charge and in any dimension, and what that does to the
# quarter-area question.
#
# quarter_area_jacobson.R killed the reading of A.15's quarter as Jacobson's eta, because the
# quarter is a Schwarzschild accident: the invariant across charge and dimension is the HALF
# of the interior band, not the quarter of the area. It left the question of what the half
# actually IS, since "a half of a band" is not an object Clausius knows about.
#
# contact_vanvleck.R answers that for Schwarzschild: the contact boundary is where the
# connecting geodesic ends on a point CONJUGATE to its start, because the rotation about the
# axis through x vanishes on that axis and the axis runs through the antipodal direction too.
# That argument uses spherical symmetry and nothing else, so it should survive charge and
# dimension where the quarter does not. This checks that it does, and closes the tidal in
# closed form on the way.

M <- 1; L <- 1

# ---- the connecting geodesic in a general static spherically symmetric metric.
# E = 0 gives rdot^2 = -f L^2/r^2 and phidot = L/r^2, so dr/dphi = +- r sqrt(-f).
# The first-order form stalls at the turning point r = r_+, so integrate the second-order
# one, which is regular there: r'' = -r f - r^2 f'/2, obtained by differentiating r' = r w.
orbit <- function(r0, f, fp, phimax, n = 400000) {
  h <- phimax / n; y <- c(r0, r0 * sqrt(-f(r0))); out <- matrix(0, n + 1, 3)
  out[1, ] <- c(0, y)
  d <- function(y) c(y[2], -y[1] * f(y[1]) - y[1]^2 * fp(y[1]) / 2)
  for (i in 1:n) {
    k1 <- d(y); k2 <- d(y + h/2 * k1); k3 <- d(y + h/2 * k2); k4 <- d(y + h * k3)
    y <- y + h/6 * (k1 + 2*k2 + 2*k3 + k4); out[i + 1, ] <- c(i * h, y)
  }
  out
}
# the angle turned before returning to the radius it left
budget <- function(r0, f, fp) {
  o <- orbit(r0, f, fp, 2 * pi, 200000)
  i <- which(o[, 1] > 0.5 & o[, 2] <= r0)[1]
  approx(o[(i-1):i, 2], o[(i-1):i, 1], xout = r0)$y
}

cat("=== 1. the tidal in closed form, from the Killing field rather than from Riemann ===\n")
cat("   If J = r sin(phi) is the focusing Jacobi field, the equation it satisfies fixes the\n")
cat("   tidal. Substituting, with w = sqrt(-f) and dr/dphi = r w, every phi-dependence\n")
cat("   collapses and what is left is\n")
cat("        T_theta = (L^2/r^4) ( 1 - f + r f'/2 ).\n")
cat("   Schwarzschild: f = 1 - 2M/r gives 3 M L^2 / r^5, which is the number contact_vanvleck\n")
cat("   read off a numerically built Riemann tensor. So the closed form is testable.\n")
T_closed <- function(r, f, fp) (L^2 / r^4) * (1 - f(r) + r * fp(r) / 2)

# ---- numerical Riemann, the same machinery contact_vanvleck.R validated
gmet_of <- function(f) function(r, th) diag(c(-f(r), 1 / f(r), r^2, r^2 * sin(th)^2))
tidal_num <- function(r, f) {
  g <- gmet_of(f); h <- 1e-5
  christ <- function(r, th) {
    G0 <- g(r, th); gi <- solve(G0); dg <- array(0, c(4, 4, 4))
    dg[2, , ] <- (g(r + h, th) - g(r - h, th)) / (2 * h)
    dg[3, , ] <- (g(r, th + h) - g(r, th - h)) / (2 * h)
    G <- array(0, c(4, 4, 4))
    for (a in 1:4) for (b in 1:4) for (c in 1:4) {
      s <- 0; for (d in 1:4) s <- s + gi[a, d] * (dg[b, d, c] + dg[c, d, b] - dg[d, b, c])
      G[a, b, c] <- 0.5 * s
    }
    G
  }
  hh <- 1e-4; G <- christ(r, pi/2); dG <- array(0, c(4, 4, 4, 4))
  dG[2, , , ] <- (christ(r + hh, pi/2) - christ(r - hh, pi/2)) / (2 * hh)
  dG[3, , , ] <- (christ(r, pi/2 + hh) - christ(r, pi/2 - hh)) / (2 * hh)
  Rud <- array(0, c(4, 4, 4, 4))
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4) {
    s <- dG[cc, a, dd, b] - dG[dd, a, cc, b]
    for (e in 1:4) s <- s + G[a, cc, e] * G[e, dd, b] - G[a, dd, e] * G[e, cc, b]
    Rud[a, b, cc, dd] <- s
  }
  G0 <- g(r, pi/2); Rdn <- array(0, c(4, 4, 4, 4))
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4)
    Rdn[a, b, cc, dd] <- sum(G0[a, ] * Rud[, b, cc, dd])
  kr <- sqrt(-f(r)) * L / r; kf <- L / r^2; k <- c(0, kr, 0, kf); e <- c(0, 0, 1/r, 0)
  s <- 0
  for (a in 1:4) for (b in 1:4) for (cc in 1:4) for (dd in 1:4)
    s <- s + Rdn[a, b, cc, dd] * e[a] * k[b] * e[cc] * k[dd]
  s
}

cat("\n=== 2. the closed form against the numerical Riemann tensor, uncharged and charged ===\n")
cat("      Q/M     r      T (closed form)   T (Riemann)     rel. err\n")
for (Q in c(0, 0.5, 0.9)) {
  f  <- function(r) 1 - 2 * M / r + Q^2 / r^2
  fp <- function(r) 2 * M / r^2 - 2 * Q^2 / r^3
  rp <- M * (1 + sqrt(1 - Q^2)); rm <- M * (1 - sqrt(1 - Q^2))
  for (r in c(0.55 * rp + 0.45 * rm, 0.8 * rp + 0.2 * rm)) {
    tc <- T_closed(r, f, fp); tn <- tidal_num(r, f)
    cat(sprintf("   %6.2f  %6.4f  %15.8f  %14.8f   %.2e\n", Q, r, tc, tn, abs(tc/tn - 1)))
    stopifnot(abs(tc / tn - 1) < 2e-5)
  }
}
cat("   Charged as well as not, so the closed form is the tidal and not a Schwarzschild fit.\n")
cat("   Reissner-Nordstrom: T = L^2 (3 M r - 2 Q^2)/r^6, which CHANGES SIGN inside\n")
cat(sprintf("   r = 2Q^2/3M; at Q/M = 0.9 that is r = %.4f against r_- = %.4f.\n",
            2*0.81/3, M*(1 - sqrt(1 - 0.81))))

cat("\n=== 3. J = r sin(phi) solves it, with the tidal taken from Riemann rather than assumed ===\n")
cat("   The rotation about the axis through x is -sin(phi) d_theta on the equator whatever\n")
cat("   f is, and a Killing field restricted to a geodesic is a Jacobi field, so in the\n")
cat("   parallel frame d_theta/r the focusing field is r sin(phi). It vanishes at phi = 0\n")
cat("   and at phi = pi because a rotation vanishes on its own axis. No metric enters that.\n")
cat("   The check below is not circular: the Jacobi equation is\n")
cat("        J'' - (2 r'/r) J' + (r^4/L^2) T J = 0,\n")
cat("   and T is read off the NUMERICAL Riemann tensor at each point, not from section 1.\n")
resid <- function(r, rp_, phi, f, Tnum) {
  # r' = r w, r'' = -r f - r^2 f'/2, and J = r sin(phi)
  fpv <- (f(r + 1e-6) - f(r - 1e-6)) / 2e-6
  rpp <- -r * f(r) - r^2 * fpv / 2
  J   <- r * sin(phi)
  Jp  <- rp_ * sin(phi) + r * cos(phi)
  Jpp <- rpp * sin(phi) + 2 * rp_ * cos(phi) - r * sin(phi)
  Jpp - (2 * rp_ / r) * Jp + (r^4 / L^2) * Tnum * J
}
cat("      Q/M     phi        r        residual\n")
for (Q in c(0, 0.6)) {
  f  <- function(r) 1 - 2 * M / r + Q^2 / r^2
  fp <- function(r) 2 * M / r^2 - 2 * Q^2 / r^3
  o  <- orbit(M, f, fp, pi, 20000)
  for (ph in c(0.3, 0.8, 1.5, 2.3, 2.9)) {
    i <- which.min(abs(o[, 1] - ph))
    rr <- o[i, 2]; rpv <- o[i, 3]
    rs <- resid(rr, rpv, o[i, 1], f, tidal_num(rr, f))
    cat(sprintf("   %6.2f  %6.3f  %8.5f  %13.2e\n", Q, o[i, 1], rr, rs))
    stopifnot(abs(rs) < 1e-4)
  }
}
cat("   So the closed-form field solves the equation the real curvature writes.\n")

cat("\n=== 4. the contact boundary is the conjugate locus, at every charge ===\n")
cat("   The budget is the angle an E = 0 geodesic turns before returning to the radius it\n")
cat("   left; the bill is pi, the angle to the antipode; contact is budget >= bill.\n")
cat("      Q/M      r_+       r_-     contact r    (r_+ + r_-)/2   budget there   A_c/A_+\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  f  <- function(r) 1 - 2 * M / r + Q^2 / r^2
  fp <- function(r) 2 * M / r^2 - 2 * Q^2 / r^3
  rp <- M * (1 + sqrt(1 - Q^2)); rm <- M * (1 - sqrt(1 - Q^2))
  rc <- uniroot(function(r) budget(r, f, fp) - pi,
                c(rm + 0.02 * (rp - rm), rp - 0.02 * (rp - rm)), tol = 1e-10)$root
  cat(sprintf("   %6.2f  %8.5f  %8.5f  %10.6f  %14.6f  %12.6f  %9.6f\n",
              Q, rp, rm, rc, (rp + rm) / 2, budget(rc, f, fp), (rc / rp)^2))
  stopifnot(abs(rc - (rp + rm) / 2) < 1e-4)
}
cat("   The contact radius is (r_+ + r_-)/2 = M at every charge, so the boundary is where\n")
cat("   the turn is exactly pi, which by section 3 is exactly where r sin(phi) has its\n")
cat("   second zero. CONTACT BOUNDARY = CONJUGATE LOCUS, at every charge.\n")
cat("   The area fraction in the last column is not invariant, which is why the quarter\n")
cat("   failed where the half did not: the half is a statement about the turn.\n")
cat("\n=== 5. dimension, where the same argument gives the caustic order ===\n")
cat("   In D dimensions the transverse space of a null geodesic has D-2 directions: one is\n")
cat("   d_t and D-3 are sphere directions, and every one of the D-3 carries a rotation\n")
cat("   Killing field vanishing at phi = 0 and phi = pi. So all D-3 focus together and\n")
cat("        caustic order = D - 3,\n")
cat("   which is 1 in four dimensions. Vacuum permits it: the transverse tidal trace is\n")
cat("   R_ab k^a k^b = 0, so D-3 positive eigenvalues are paid for by the single negative\n")
cat("   one in the d_t direction, whatever D is. The Einstein static universe's order two\n")
cat("   needs D = 5 to be a vacuum answer, and in four dimensions it needs matter.\n")
cat("   Tangherlini check, f = 1 - (r_h/r)^{D-3}: the same substitution gives\n")
cat("        T_theta = (L^2/r^4) (r_h/r)^{D-3} (1 + (D-3)/2),\n")
cat("   which is 3 M L^2/r^5 at D = 4 with r_h = 2M.\n")
for (D in 4:7) {
  n <- D - 3; rh <- 2
  Td <- function(r) (L^2 / r^4) * (rh / r)^n * (1 + n / 2)
  cat(sprintf("      D = %d   order %d   T at r = r_h/2: %10.6f\n", D, n, Td(rh / 2)))
}
cat(sprintf("      and at D = 4 that is 3 M L^2/r^5 = %.6f with r = M\n", 3 * M / 1^5))

cat("\n=== 6. what this hands the Jacobson route ===\n")
cat("   quarter_area_jacobson.R asked what the invariant half IS, since a half of a band is\n")
cat("   not something Clausius knows about. It is a CONJUGATE LOCUS, which is something\n")
cat("   Clausius's machinery does know about: a conjugate locus is where the Van Vleck\n")
cat("   determinant diverges, which is where the one-loop effective action diverges, which\n")
cat("   is where an induced-gravity cutoff would have to sit. The fold therefore supplies a\n")
cat("   geometrically defined surface at which the mode counting fails, at every charge and\n")
cat("   in every dimension, without a cutoff being put in by hand.\n")
cat("   That is a route and not a result. What it does not yet supply is eta itself, which\n")
cat("   needs the entanglement entropy regulated ON that surface rather than at an arbitrary\n")
cat("   proper distance, and that calculation has not been done here.\n")

cat("\n=== 7. the plant ===\n")
cat("   (a) a metric that is not spherically symmetric must break the closed-form tidal.\n")
f_bad  <- function(r) 1 - 2 * M / r
Tb <- T_closed(1.4, f_bad, function(r) 2 * M / r^2)
Tn <- tidal_num(1.4, f_bad)
gm <- gmet_of(f_bad)
gmet_of_saved <- gmet_of
gmet_of <- function(f) function(r, th) { g <- gm(r, th); g[4, 4] <- g[4, 4] * 1.05; g }
Tn2 <- tidal_num(1.4, f_bad)
cat(sprintf("       round metric: closed %.6f vs Riemann %.6f  (agree)\n", Tb, Tn))
cat(sprintf("       squashed sphere: Riemann %.6f, closed form now wrong by %.1f%%\n",
            Tn2, 100 * abs(Tb / Tn2 - 1)))
stopifnot(abs(Tb / Tn2 - 1) > 0.01)
gmet_of <- gmet_of_saved
cat("   (b) the bill must matter: move it off pi and the contact radius must move.\n")
f <- function(r) 1 - 2 * M / r
for (bill in c(0.8 * pi, pi, 1.2 * pi)) {
  rc <- uniroot(function(r) budget(r, f, function(r) 2 * M / r^2) - bill,
                c(0.2, 1.9), tol = 1e-10)$root
  cat(sprintf("       bill %.4f  ->  contact r = %.6f\n", bill, rc))
}
cat("   Only a bill of pi puts it at M, so the coincidence of the half with the antipodal\n")
cat("   angle is doing work rather than being read into an insensitive root.\n")
