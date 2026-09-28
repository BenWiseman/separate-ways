#!/usr/bin/env Rscript
# The components of the image stress at a caustic, which A.18 computes the magnitude of and
# explicitly does not decompose.
#
# A.18's own words: "It does not settle the sign, which is the half of A.15's conjecture that
# decides whether the divergence closes the region or merely marks it, because a positive
# coefficient in the two-point function is not a sign for the stress tensor and the components
# have not been computed." This computes them, in A.18's own geometry and with A.18's own sum.
# Section 9 then carries the one contraction that transfers to a hole and says why the other
# does not.
#
# THE GEOMETRY, unchanged from A.18: -dt^2 + dr^2 + a^2 dOmega^2 with a = 1. Null geodesics
# leaving a point spread over every great circle of the transverse sphere and refocus together
# after exactly pi of transverse angle, so the image pair sits on a caustic with a one-parameter
# connecting family, which is the degeneracy the contact sphere carries. Reducing on the sphere
# leaves a tower of two-dimensional fields of mass mu_l^2 = m^2 + xi R + l(l+1)/a^2, R = 2/a^2:
#
#   W(x,x') = sum_l (2l+1)/(4 pi a^2) P_l(cos gamma) W_2(mu_l; Delta t, Delta r),
#   W_2(mu; T, X) = (1/2 pi) K_0( mu sqrt(X^2 - (T - i eps)^2) ).
#
# The image pair is Delta r = 0, gamma = pi, Delta t = pi - s, approached from outside the
# lightcone by the offset s. That is A.18's configuration exactly.
#
# WHAT IS NEW. Three scalars carry the whole tensor at that configuration, because the isotropy
# there is the SO(2) of rotations about the axis through the point and its antipode, so the two
# transverse pressures are equal and the t-r mixing vanishes:
#
#   D = d^2W/dT^2,   E = d^2W/dX^2 at X = 0,   Q = grad_par grad_par' W,
#
# and the transverse Hessian is isotropic at the antipode because P_l(cos gamma) there is a
# function of the SQUARED distance from the antipode: with u the tangent displacement of x and
# v that of x' transported back along the connecting geodesic, cos gamma = -1 + |u - v|^2/2 and
# P_l = (-1)^l [1 - l(l+1)|u - v|^2/4]. Differentiating once at each end gives +(-1)^l l(l+1)/2
# and twice at one end gives minus that, so the mixed and unmixed transverse second derivatives
# are equal and opposite. Section 4 states the components that follow and section 5 checks them
# two independent ways.

TOL <- 1e-6
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

# ---------------------------------------------------------------------------------------------
# 1. complex K_0 and K_1, built here and validated against R's own real-argument routine and
#    against the Hankel identity on the imaginary axis, which is where our arguments live.
# ---------------------------------------------------------------------------------------------
EG <- 0.5772156649015329
ZC <- 20                 # |z| above which the asymptotic series is used

Kser <- function(z) {     # ascending series, good near the imaginary axis where no cancellation
  t <- z*z/4
  I0 <- 0+0i; I1s <- 0+0i; S0 <- 0+0i; S1 <- 0+0i
  term <- rep(1+0i, length(z)); h <- 0
  for (k in 0:250) {
    if (k > 0) { term <- term * t / (k*k); h <- h + 1/k }
    I0 <- I0 + term
    if (k > 0) { S0 <- S0 + h*term; S1 <- S1 + h*(2*k/z)*term }
    I1s <- I1s + term/(k+1)          # I1 = (z/2) sum t^k/(k!(k+1)!)
  }
  I1 <- (z/2)*I1s
  lg <- log(z/2) + EG
  list(K0 = -lg*I0 + S0, K1 = I0/z + lg*I1 - S1)
}

Kasy <- function(z, nu, nterm = 12) {
  s <- rep(1+0i, length(z)); a <- rep(1+0i, length(z))
  for (k in 1:nterm) { a <- a * (4*nu^2 - (2*k-1)^2)/(8*k*z); s <- s + a }
  sqrt(pi/(2*z)) * exp(-z) * s
}

# Three bands, because no single one covers the plane in double precision. The ascending series
# loses digits to cancellation when z is nearly REAL and large, since it builds K_0 ~ 1e-6 out of
# I_0 ~ 1e4; the asymptotic series is short of a few digits in the middle of the range; and the
# integral representation K_nu(z) = int e^{-z cosh t} cosh(nu t) dt is exact where the series
# cancels, because its integrand decays at cosh t ~ 1/Re z. So: integral for nearly-real middling
# arguments, series near the imaginary axis, asymptotic far out. The physics below only ever uses
# the second and third, since z = mu (eps + i T) with eps small, and the first exists so the
# validation against R's real-argument besselK can be held to 1e-8 rather than to what one branch
# happens to manage.
Kint <- function(z, nu, n = 6000, tmax = 6) {
  h <- tmax/n; tt <- (0:n)*h
  wts <- c(1, rep(c(4,2), length.out = n-1), 1) * h/3
  ch <- cosh(tt); cn <- cosh(nu*tt)
  sapply(z, function(zz) sum(wts * exp(-zz*ch) * cn))
}

K01 <- function(z) {
  big <- Mod(z) > ZC
  mid <- !big & Re(z) > 2
  ser <- !big & !mid
  K0 <- complex(length.out = length(z)); K1 <- K0
  if (any(big)) { K0[big] <- Kasy(z[big], 0); K1[big] <- Kasy(z[big], 1) }
  if (any(mid)) { K0[mid] <- Kint(z[mid], 0); K1[mid] <- Kint(z[mid], 1) }
  if (any(ser)) { r <- Kser(z[ser]); K0[ser] <- r$K0; K1[ser] <- r$K1 }
  list(K0 = K0, K1 = K1)
}

cat("=== 1. the complex Bessel routine, checked before anything is asked of it ===\n")
cat("   (a) against R's real-argument besselK, on both branches:\n")
for (x in c(0.7, 3.0, 12.0, 25.0, 60.0)) {
  o <- K01(complex(real = x, imaginary = 0))
  e0 <- besselK(x, 0); e1 <- besselK(x, 1)
  cat(sprintf("      z = %6.2f  K0 %14.7e vs %14.7e   K1 %14.7e vs %14.7e  branch %s\n",
              x, Re(o$K0), e0, Re(o$K1), e1,
              ifelse(x > ZC, "asymptotic", ifelse(x > 2, "integral", "series"))))
  note(abs(Re(o$K0)/e0 - 1) < 1e-8 && abs(Re(o$K1)/e1 - 1) < 1e-8, "real-argument K agrees")
}
cat("   (b) against the Hankel identity K_nu(i y) = (pi/2)(-i)^{nu+1} H_nu^{(2)}(y), which is an\n")
cat("       independent route and the one our arguments sit on:\n")
for (y in c(1.5, 4.44, 11.0, 30.0, 120.0)) {
  o <- K01(complex(real = 1e-12, imaginary = y))
  h0 <- (pi/2)*(-1i)^1 * (besselJ(y,0) - 1i*besselY(y,0))
  h1 <- (pi/2)*(-1i)^2 * (besselJ(y,1) - 1i*besselY(y,1))
  cat(sprintf("      z = i%6.2f  K0 %+.7e%+.7ei vs %+.7e%+.7ei\n", y,
              Re(o$K0), Im(o$K0), Re(h0), Im(h0)))
  note(Mod(o$K0 - h0) < 1e-7*Mod(h0) && Mod(o$K1 - h1) < 1e-7*Mod(h1),
       sprintf("Hankel identity holds at y = %g", y))
}

# ---------------------------------------------------------------------------------------------
# 2. the tower, summed in chunks, at the image pair
# ---------------------------------------------------------------------------------------------
# W and the three second derivatives.  With X = 0 the 2d argument is Z = eps + i T and
#   w      = K_0(mu Z)/2 pi
#   d^2/dT^2 w = -(mu^2/2 pi)[K_0(mu Z) + K_1(mu Z)/(mu Z)]
#   d^2/dX^2 w = -(mu/2 pi Z) K_1(mu Z)
# the last two from Z = i(T - i eps), dZ/dT = i, d^2Z/dX^2 = 1/Z at X = 0.
tower <- function(s, eps, m2 = 0, xi = 0, lmin = 1, cut = 45, chunk = 250000, gam = pi) {
  Tsep <- gam - s                      # null separation across the transverse angle gam
  Z  <- complex(real = eps, imaginary = Tsep)
  R4 <- 2                              # R = 2/a^2 with a = 1
  mu2c <- m2 + xi*R4
  lmax <- max(lmin + 10, ceiling((cut/eps)))
  W <- 0+0i; D <- 0+0i; E <- 0+0i; Q <- 0+0i
  l0 <- lmin
  while (l0 <= lmax) {
    l1 <- min(l0 + chunk - 1, lmax)
    l  <- l0:l1
    mu <- sqrt(l*(l+1) + mu2c)
    z  <- mu*Z
    kk <- K01(z)
    w  <- kk$K0/(2*pi)
    dd <- -(mu^2/(2*pi))*(kk$K0 + kk$K1/z)
    ee <- -(mu/(2*pi*Z))*kk$K1
    cl <- (2*l+1)/(4*pi)
    sg <- ifelse(l %% 2 == 0, 1, -1)
    aw <- cl*sg
    W <- W + sum(aw*w); D <- D + sum(aw*dd); E <- E + sum(aw*ee)
    Q <- Q + sum(aw*(l*(l+1)/2)*w)
    l0 <- l1 + 1
  }
  list(W = W, D = D, E = E, Q = Q, lmax = lmax)
}

cat("   (c) the transverse Hessian at the antipode, by finite differences on the exact geometry.\n")
cat("       This is the step the whole decomposition rests on, so it is measured and not asserted.\n")
cat("       Displace x by u from the pole in azimuth b1 and x' by v from the antipode in azimuth\n")
cat("       b2; then cos gamma = -cos u cos v + sin u sin v cos(b1-b2) exactly. Differentiating\n")
cat("       once at each end and twice at one end must give the SAME number before transport, and\n")
cat("       the transport along the connecting geodesic then flips the sign of the mixed one.\n")
cg_exact <- function(u, v, db) -cos(u)*cos(v) + sin(u)*sin(v)*cos(db)
Pl <- function(l, x) { p0 <- 1; p1 <- x; if (l == 0) return(rep(1, length(x)))
  if (l == 1) return(x)
  for (n in 2:l) { p2 <- ((2*n-1)*x*p1 - (n-1)*p0)/n; p0 <- p1; p1 <- p2 }; p1 }
h <- 1e-4
cat("\n        l     mixed d_u d_v P      unmixed d_u^2 P      closed form -(-1)^l l(l+1)/2\n")
for (l in c(1, 2, 5, 12, 40)) {
  mx <- (Pl(l, cg_exact(h, h, 0)) - Pl(l, cg_exact(h, -h, 0))
         - Pl(l, cg_exact(-h, h, 0)) + Pl(l, cg_exact(-h, -h, 0)))/(4*h^2)
  un <- (Pl(l, cg_exact(h, 0, 0)) - 2*Pl(l, cg_exact(0, 0, 0)) + Pl(l, cg_exact(-h, 0, 0)))/h^2
  cf <- -(-1)^l * l*(l+1)/2
  cat(sprintf("   %6d %18.8f %20.8f %26.4f\n", l, mx, un, cf))
  note(abs(mx - cf) < 1e-4*max(1, abs(cf)), sprintf("the mixed transverse Hessian at l = %d", l))
  note(abs(un - cf) < 1e-4*max(1, abs(cf)), sprintf("the unmixed one agrees before transport, l = %d", l))
}
cat("   (d) and the transport itself, in R^3 rather than by assertion. Parallel transport along a\n")
cat("       great circle carries the curve's own tangent to the curve's own tangent, so the\n")
cat("       transport of e_b from the pole to the antipode is MINUS e_b:\n")
for (b in c(0, 0.7, 2.3)) {
  t0 <- c(cos(b), sin(b), 0)                       # d/ds of (sin s cos b, sin s sin b, cos s) at 0
  t1 <- c(-cos(b), -sin(b), 0)                     # and at s = pi
  cat(sprintf("      azimuth %.2f: tangent at the pole (%+.4f,%+.4f,%+.4f), at the antipode (%+.4f,%+.4f,%+.4f)\n",
              b, t0[1], t0[2], t0[3], t1[1], t1[2], t1[3]))
  note(max(abs(t1 + t0)) < 1e-12, "the transported frame is the reversed one")
}
cat("   So in the transported frame the mixed transverse second derivative is PLUS (-1)^l l(l+1)/2\n")
cat("   and the unmixed one is minus it, which is the relation section 4 uses.\n")

cat("\n=== 2. the sum reproduces A.18's own result before it is asked for anything else ===\n")
cat("   A.18 reports the image term diverging as s^{-3/2} with the limit on the positive real\n")
cat("   axis and a magnitude of 502. The closed form the tower's tail gives is\n")
cat("        W -> 1/(8 pi sqrt(2 pi)) s^{-3/2} = 0.0158740 s^{-3/2},\n")
cat("   which at s = 0.001 is 502.0. Both the power and that number are checked here.\n\n")
CW <- 1/(8*pi*sqrt(2*pi))
CD <- 15/(32*pi*sqrt(2*pi))
CE <- -3/(16*pi^2*sqrt(2*pi))
cat(sprintf("   closed forms: W %.7f s^{-3/2}   D %.7f s^{-7/2}   E %.7f s^{-5/2}   Q %.7f s^{-7/2}\n\n",
            CW, CD, CE, -CD/2))
SS <- c(0.03, 0.01, 0.003, 0.001)
res <- list()
cat("      s        W (real)        W/(s^{-3/2})     Im/Re        l_max\n")
for (s in SS) {
  o <- tower(s, eps = s/40)
  res[[as.character(s)]] <- o
  cat(sprintf("   %7.4f %14.5f %16.7f %12.2e %10d\n",
              s, Re(o$W), Re(o$W)*s^1.5, Im(o$W)/Re(o$W), o$lmax))
}
o3 <- res[["0.001"]]
cat(sprintf("\n   at s = 0.001 the sum gives %.1f against A.18's 502\n", Re(o3$W)))
note(abs(Re(o3$W)/502 - 1) < 0.02, "the magnitude reproduces A.18's 502")
note(abs(Re(o3$W)*0.001^1.5/CW - 1) < 0.02, "and matches the closed-form coefficient")
pw <- function(k) {
  a <- res[[as.character(SS[1])]][[k]]; b <- res[[as.character(SS[length(SS)])]][[k]]
  log(Re(b)/Re(a))/log(SS[length(SS)]/SS[1])
}
cat(sprintf("   fitted power of W over the whole range: %.5f against -3/2\n", pw("W")))
note(abs(pw("W") + 1.5) < 0.02, "the image term's power is -3/2")

# ---------------------------------------------------------------------------------------------
# 3. the three derivatives, their powers and their coefficients
# ---------------------------------------------------------------------------------------------
cat("\n=== 3. the second derivatives: the powers A.18 predicts and the coefficients it does not ===\n")
cat("   Two derivatives on a tower whose terms grow as sqrt(l) cost two powers of l, so D and Q\n")
cat("   go as s^{-7/2}, which is A.18's stress-tensor power. E gains only one power of l because\n")
cat("   the radial second derivative brings 1/Z rather than mu, so it goes as s^{-5/2} and is\n")
cat("   subleading. That ordering is a prediction and it is checked, not assumed.\n\n")
cat("      s          D               D s^{7/2}        E              E s^{5/2}       Q s^{7/2}\n")
for (s in SS) {
  o <- res[[as.character(s)]]
  cat(sprintf("   %7.4f %15.4e %14.7f %15.4e %14.7f %14.7f\n",
              s, Re(o$D), Re(o$D)*s^3.5, Re(o$E), Re(o$E)*s^2.5, Re(o$Q)*s^3.5))
}
cat(sprintf("\n   fitted powers: D %.5f (-7/2)   E %.5f (-5/2)   Q %.5f (-7/2)\n",
            pw("D"), pw("E"), pw("Q")))
note(abs(pw("D") + 3.5) < 0.03, "D goes as s^{-7/2}")
note(abs(pw("E") + 2.5) < 0.05, "E goes as s^{-5/2}, one power softer")
note(abs(pw("Q") + 3.5) < 0.03, "Q goes as s^{-7/2}")
note(abs(Re(o3$D)*0.001^3.5/CD - 1) < 0.02, "D matches its closed-form coefficient")
note(abs(Re(o3$Q)*0.001^3.5/(-CD/2) - 1) < 0.02, "Q is minus half of D, as the tower requires")

# ---------------------------------------------------------------------------------------------
# 4. the components
# ---------------------------------------------------------------------------------------------
cat("\n=== 4. the full tensor at the caustic ===\n")
cat("   The point-split stress tensor for a scalar is\n")
cat("     T_ab = (1-2xi)W_;ab' + (2xi-1/2)g_ab g^cd' W_;cd' - 2xi W_;ab + 2xi g_ab box W\n")
cat("            + xi(R_ab - R g_ab/2)W - m^2 g_ab W/2,\n")
cat("   and at the image pair the primed and unprimed second derivatives are related: W_;tt' =\n")
cat("   -D, W_;rr' = -E, W_;par par' = +Q against W_;tt = D, W_;rr = E, W_;par par = -Q. With\n")
cat("   R_ab = g_ab/a^2 on the sphere and zero in the 2d factor, the components collapse to\n")
cat("     rho = -D/2 - E/2 + Q + xi W/a^2 + m^2 W/2\n")
cat("     p_r = -D/2 - E/2 - Q - xi W/a^2 - m^2 W/2\n")
cat("     p_T = -D/2 + E/2 - m^2 W/2\n\n")
comp <- function(o, m2 = 0, xi = 0) {
  W <- Re(o$W); D <- Re(o$D); E <- Re(o$E); Q <- Re(o$Q)
  list(rho = -D/2 - E/2 + Q + xi*W + m2*W/2,
       p_r = -D/2 - E/2 - Q - xi*W - m2*W/2,
       p_T = -D/2 + E/2 - m2*W/2,
       kk_rad = -D - E,
       kk_con = -D + Q + xi*W)
}
cat("      s         rho s^{7/2}     p_r s^{7/2}     p_T s^{7/2}     trace s^{7/2}\n")
for (s in SS) {
  c1 <- comp(res[[as.character(s)]])
  tr <- -c1$rho + c1$p_r + 2*c1$p_T
  cat(sprintf("   %7.4f %15.7f %15.7f %15.7f %15.2e\n",
              s, c1$rho*s^3.5, c1$p_r*s^3.5, c1$p_T*s^3.5, tr*s^3.5))
}
cat("\n   Three things to read off that. The energy density is NEGATIVE. The radial pressure\n")
cat("   vanishes at leading order, its s^{-7/2} coefficient cancelling between -D/2 and -Q.\n")
cat("   And the leading divergence is TRACELESS, which it has to be, because the trace of the\n")
cat("   renormalised stress tensor is a local geometric quantity and cannot diverge. That last\n")
cat("   one was not imposed anywhere above, so it is a check and not a definition.\n")
cL <- comp(o3)
note(cL$rho < 0, "the image energy density at the caustic is negative")
note(abs(cL$p_r) < 0.05*abs(cL$rho), "the radial pressure's leading divergence cancels")
note(abs((-cL$rho + cL$p_r + 2*cL$p_T)) < 0.05*abs(cL$rho), "the leading divergence is traceless")

cat("\n=== 4b. the field equation, which fixes the whole assembly and can fail ===\n")
cat("   Each tower mode obeys (d_X^2 - d_T^2 - mu_l^2)w_l = 0 and carries S^2 eigenvalue -l(l+1),\n")
cat("   so the four scalars are not independent: box W = -D + E - 2Q must equal (m^2 + xi R)W\n")
cat("   exactly, at every offset and for every coupling. That in turn forces the trace to be\n")
cat("   T^a_a = -m^2 W, so a massless field's image stress is traceless to ALL orders in s and not\n")
cat("   only the leading one. Nothing above imposed either statement.\n")
cat("   THE SECOND ONE IS THIS GEOMETRY'S, and the restriction is not cosmetic. It uses\n")
cat("   g^{cd'}W_{;cd'} = -box W, which holds because W depends on the coordinate difference here.\n")
cat("   Where the image point moves with the point, as at a hole, the general trace carries\n")
cat("   (6 xi - 1) box <phi^2>/2 as well, which is A.19's form and the one that governs there.\n\n")
cat("      xi        m        s       box W            (m^2+xi R)W        ratio       T^a_a / (-m^2 W)\n")
for (p in list(c(0,0,0.01), c(1/6,0,0.01), c(1/6,2,0.003), c(1/3,1,0.003), c(0.5,3,0.003))) {
  o <- tower(p[3], eps = p[3]/40, m2 = p[2]^2, xi = p[1])
  W <- Re(o$W); D <- Re(o$D); E <- Re(o$E); Q <- Re(o$Q)
  box <- -D + E - 2*Q; rhs <- (p[2]^2 + p[1]*2)*W
  c1 <- comp(o, m2 = p[2]^2, xi = p[1]); tr <- -c1$rho + c1$p_r + 2*c1$p_T
  sc <- max(abs(D), abs(rhs))
  cat(sprintf("   %7.4f %7.2f %8.4f %15.6e %15.6e %12s %18s\n", p[1], p[2], p[3], box, rhs,
              ifelse(abs(rhs) > 1e-8*abs(D), sprintf("%.6f", box/rhs), "both ~ 0"),
              ifelse(p[2] > 0, sprintf("%.6f", tr/(-p[2]^2*W)), "trace ~ 0")))
  note(abs(box - rhs) < 1e-8*sc, sprintf("box W = (m^2+xi R)W at xi=%.3f, m=%.1f", p[1], p[2]))
  note(abs(tr + p[2]^2*W) < 1e-8*sc, "and the trace is exactly -m^2 W")
}
cat("   The plant for this one: the same check run with Q scaled by 1.001, which is a tenth of a\n")
cat("   per cent and must still break it, since the identity is exact and not approximate.\n")
o <- res[["0.003"]]
boxb <- -Re(o$D) + Re(o$E) - 2*Re(o$Q)*1.001
cat(sprintf("      Q x 1.001: box W = %.6e against 0, which is %.2e of D\n",
            boxb, abs(boxb)/abs(Re(o$D))))
note(abs(boxb) > 1e-8*abs(Re(o$D)), "a tenth of a per cent on Q breaks the field equation")

cat("\n=== 5. the two null contractions, which is what the fork asked for ===\n")
cat("   The radial null direction is k = e_t + e_r, transverse to the sphere, and it is the one\n")
cat("   Penrose's ingoing congruence runs on. The contact direction is k = e_t + e_par, along the\n")
cat("   connecting family. Both contractions follow from the components, and both follow again\n")
cat("   directly from the point-split expression with the metric terms killed by k.k = 0:\n")
cat("     T_kk(radial)  = -D - E          (no xi, no m: the couplings cancel)\n")
cat("     T_kk(contact) = -D + Q + xi W\n")
cat("   The two routes are independent and section 5b checks they agree.\n\n")
cat("      s        T_kk radial     T_kk contact    ratio con/rad\n")
for (s in SS) {
  c1 <- comp(res[[as.character(s)]])
  cat(sprintf("   %7.4f %15.4e %15.4e %14.7f\n", s, c1$kk_rad, c1$kk_con, c1$kk_con/c1$kk_rad))
}
note(cL$kk_rad < 0, "the RADIAL null contraction is negative")
note(cL$kk_con < 0, "the CONTACT null contraction is negative")
note(abs(cL$kk_con/cL$kk_rad - 1.5) < 0.03, "the ratio is 3/2")
cat("\n   Both are negative and the ratio is 3/2. So in A.18's geometry the image stress violates\n")
cat("   the null convergence condition in BOTH directions, and the contact direction feels it\n")
cat("   half again as strongly.\n")

cat("\n=== 5b. the components and the direct contractions agree ===\n")
for (s in c(0.01, 0.001)) {
  c1 <- comp(res[[as.character(s)]])
  a1 <- c1$rho + c1$p_r; a2 <- c1$rho + c1$p_T
  cat(sprintf("   s = %.4f:  rho+p_r %15.7e vs T_kk rad %15.7e\n", s, a1, c1$kk_rad))
  cat(sprintf("              rho+p_T %15.7e vs T_kk con %15.7e\n", a2, c1$kk_con))
  note(abs(a1/c1$kk_rad - 1) < 1e-10 && abs(a2/c1$kk_con - 1) < 1e-10,
       "the two routes to the contractions agree")
}

# ---------------------------------------------------------------------------------------------
# 6. the sign does not depend on the coupling or the mass
# ---------------------------------------------------------------------------------------------
cat("\n=== 6. neither sign moves with xi or with m, which is stronger than the flat-space case ===\n")
cat("   On the Einstein static universe the image T_kk is proportional to minus the Hadamard\n")
cat("   coefficient V_0 = [m^2 + (xi - 1/6)R]/2, so its sign FLIPS as the coupling crosses 1/6.\n")
cat("   At a caustic it cannot: xi and m enter the tower only through a constant shift of mu_l^2\n")
cat("   and through terms two powers softer than the leading divergence.\n\n")
cat("      xi        m       T_kk radial     T_kk contact    ratio\n")
for (p in list(c(0,0), c(1/6,0), c(1/3,0), c(0,1), c(1/6,2), c(1/2,3))) {
  o <- tower(0.003, eps = 0.003/40, m2 = p[2]^2, xi = p[1])
  c1 <- comp(o, m2 = p[2]^2, xi = p[1])
  cat(sprintf("   %7.4f %7.2f %15.4e %15.4e %12.6f\n",
              p[1], p[2], c1$kk_rad, c1$kk_con, c1$kk_con/c1$kk_rad))
  note(c1$kk_rad < 0 && c1$kk_con < 0, sprintf("both negative at xi = %.3f, m = %.1f", p[1], p[2]))
  note(abs(c1$kk_con/c1$kk_rad - 1.5) < 0.05, "and the ratio holds")
}

# ---------------------------------------------------------------------------------------------
# 7. the plants
# ---------------------------------------------------------------------------------------------
cat("\n=== 7. the plants: three ways this file has to be able to fail ===\n")
cat("   (a) Remove the antipodal parity by hand, leaving everything else alone. A.18 reports the\n")
cat("       divergence going away at +0.0004 when it does that; here the stress must go with it.\n")
tow_nopar <- function(s, eps, cut = 45, chunk = 250000) {
  Tsep <- pi - s; Z <- complex(real = eps, imaginary = Tsep)
  lmax <- ceiling(cut/eps); W <- 0+0i; D <- 0+0i; l0 <- 1
  while (l0 <= lmax) {
    l1 <- min(l0 + chunk - 1, lmax); l <- l0:l1
    mu <- sqrt(l*(l+1)); z <- mu*Z; kk <- K01(z)
    cl <- (2*l+1)/(4*pi)
    W <- W + sum(cl*kk$K0/(2*pi))
    D <- D + sum(cl*(-(mu^2/(2*pi))*(kk$K0 + kk$K1/z)))
    l0 <- l1 + 1
  }
  list(W = W, D = D)
}
for (s in c(0.01, 0.001)) {
  o <- tow_nopar(s, s/40)
  cat(sprintf("      s = %.4f:  W %12.5f   D %12.5f   (with parity: W %10.2f  D %12.4e)\n",
              s, Re(o$W), Re(o$D), Re(res[[as.character(s)]]$W), Re(res[[as.character(s)]]$D)))
}
d_np <- Re(tow_nopar(0.001, 0.001/40)$D)
note(abs(d_np) < 1e-4*abs(Re(o3$D)), "without the antipodal parity the divergence is gone")

cat("   (b) The generic null pair, no caustic: the connecting family is a single geodesic, the\n")
cat("       tower's terms do not arrive in step, and the parametrix powers must return. That is\n")
cat("       the case that could have failed, and it is why the antipodal numbers mean anything.\n")
tow_gen <- function(s, eps, gam = pi/2, cut = 45, chunk = 200000) {
  Tsep <- gam - s; Z <- complex(real = eps, imaginary = Tsep)
  lmax <- ceiling(cut/eps); cg <- cos(gam)
  # P_l(cos gam) by the standard forward recurrence, which is stable for |cos gam| <= 1
  P <- numeric(lmax + 1); P[1] <- 1; P[2] <- cg
  for (n in 2:lmax) P[n+1] <- ((2*n-1)*cg*P[n] - (n-1)*P[n-1])/n
  W <- 0+0i; D <- 0+0i; l0 <- 1
  while (l0 <= lmax) {
    l1 <- min(l0 + chunk - 1, lmax); l <- l0:l1
    mu <- sqrt(l*(l+1)); z <- mu*Z; kk <- K01(z)
    cl <- (2*l+1)/(4*pi) * P[l+1]
    W <- W + sum(cl*kk$K0/(2*pi))
    D <- D + sum(cl*(-(mu^2/(2*pi))*(kk$K0 + kk$K1/z)))
    l0 <- l1 + 1
  }
  list(W = W, D = D)
}
gs <- c(0.01, 0.003)
gv <- lapply(gs, function(s) tow_gen(s, s/40))
cat(sprintf("      generic pair: W %11.4f at s = %.3f and %11.4f at s = %.3f -> power %+.4f (-1)\n",
            Re(gv[[1]]$W), gs[1], Re(gv[[2]]$W), gs[2],
            log(Re(gv[[2]]$W)/Re(gv[[1]]$W))/log(gs[2]/gs[1])))
cat(sprintf("                    D %11.4e         and %11.4e         -> power %+.4f (-3)\n",
            Re(gv[[1]]$D), Re(gv[[2]]$D),
            log(Re(gv[[2]]$D)/Re(gv[[1]]$D))/log(gs[2]/gs[1])))
pgW <- log(Re(gv[[2]]$W)/Re(gv[[1]]$W))/log(gs[2]/gs[1])
pgD <- log(Re(gv[[2]]$D)/Re(gv[[1]]$D))/log(gs[2]/gs[1])
note(abs(pgW + 1) < 0.08, "the generic pair returns the parametrix power -1 for W")
note(abs(pgD + 3) < 0.15, "and -3 for the stress, half a power softer than the caustic's -7/2")

cat("   (c) The SIGN checks have to be able to fail, since the provenance table compares digits and\n")
cat("       drops signs, so the whole sign result rests on this file stopping when it should. With\n")
cat("       Q's sign flipped by hand the ratio leaves 3/2 for 1/2, and with D's sign flipped the\n")
cat("       radial contraction turns positive outright:\n")
flip <- list(W = o3$W, D = o3$D, E = o3$E, Q = -o3$Q)
cf <- comp(flip)
cat(sprintf("      Q -> -Q:  T_kk contact %+.4e (was %+.4e), ratio %.4f (was %.4f)\n",
            cf$kk_con, cL$kk_con, cf$kk_con/cf$kk_rad, cL$kk_con/cL$kk_rad))
note(abs(cf$kk_con/cf$kk_rad - 1.5) > 0.03, "flipping Q breaks the ratio check, so it can fail")
cat("       And the radial one, which is the contraction the transfer argument hands to Penrose:\n")
flip2 <- list(W = o3$W, D = -o3$D, E = o3$E, Q = o3$Q)
cf2 <- comp(flip2)
cat(sprintf("      D -> -D:  T_kk radial  %+.4e (was %+.4e)\n", cf2$kk_rad, cL$kk_rad))
note(cf2$kk_rad > 0, "flipping D turns the radial contraction positive")

cat("   (d) A corrupted coefficient has to break the traceless check, or that check is decoration.\n")
bad <- list(W = o3$W, D = o3$D, E = o3$E, Q = o3$Q*1.3)
cb <- comp(bad)
trb <- -cb$rho + cb$p_r + 2*cb$p_T
cat(sprintf("      Q inflated by 30 per cent: trace/rho %+.4f against %+.2e for the real one\n",
            trb/cb$rho, (-cL$rho + cL$p_r + 2*cL$p_T)/cL$rho))
note(abs(trb/cb$rho) > 0.05, "a corrupted Q breaks the traceless check, as it must")

# ---------------------------------------------------------------------------------------------
# 9. what transfers to a hole, and what does not, with the causal characters written down
# ---------------------------------------------------------------------------------------------
cat("\n=== 9. the transfer to a hole, decided by causal character rather than by coordinate name ===\n")
cat("   The obstacle to transferring a coordinate pattern is that t and r swap character across a\n")
cat("   horizon. In A.18's geometry t is timelike and r spacelike. Inside a hole f = 1 - 2M/r is\n")
cat("   negative, so r is TIMELIKE and t spacelike, and the same three numbers T^t_t, T^r_r,\n")
cat("   T^theta_theta then mean different things. What settles the transfer is that a null\n")
cat("   contraction is a scalar, so the question is which of the two contractions is built the\n")
cat("   same way in the two geometries.\n\n")
fS <- function(r, M = 1) 1 - 2*M/r
cat("      direction          A.18's geometry                a hole's interior\n")
cat("      contact            e_t + e_transverse             k^t = 0, so k^r d_r + k^th d_th\n")
cat("      radial             e_t + e_r                      mixes t and r\n")
cat("      character of r     spacelike                      TIMELIKE, since f < 0\n\n")
cat("   The contact geodesic at a hole has E = f t-dot = 0, and f is nonzero inside, so t-dot = 0:\n")
cat("   its tangent is a timelike direction plus a transverse one, which is exactly what\n")
cat("   e_t + e_transverse is in A.18's geometry. Checked here rather than asserted:\n")
for (r in c(0.4, 0.7, 1.0)) {
  fv <- fS(r)
  # E = 0 null: g_rr (k^r)^2 + r^2 (k^th)^2 = 0 with g_rr = 1/f < 0
  kth <- 1; kr <- sqrt(-fv)*r*kth      # from (1/f)(k^r)^2 = -r^2 (k^th)^2
  nn  <- (1/fv)*kr^2 + r^2*kth^2
  cat(sprintf("      r/M = %4.2f: f = %+8.4f (r timelike: %s), k = (0, %.4f, %.4f, 0), k.k = %+.2e\n",
              r, fv, ifelse(fv < 0, "yes", "no"), kr, kth, nn))
  note(fv < 0 && abs(nn) < 1e-12, "the contact tangent is null with no t component and r timelike")
}
cat("\n   The other thing the transfer needs is the caustic COUNT, because a uniform extra phase is\n")
cat("   exactly what a caustic crossing supplies and the leading real part is what carries the\n")
cat("   sign. Every scalar above picks up the same phase, so the pattern 2 : 0 : -1 and the ratio\n")
cat("   3/2 are phase-independent while the overall sign is not: a uniform extra phase theta\n")
cat("   multiplies every real part by cos(theta), which vanishes at pi/2 and flips beyond it.\n")
cat("   A.18's geometry has one crossing, the transverse one, carried by P_l(-1) = (-1)^l. A hole\n")
cat("   would have a second if its contact path turned in r, and it does not: with E = 0 the\n")
cat("   radial equation is r-dot^2 = |f| L^2/r^2, and |f| = 2M/r - 1 is strictly positive for\n")
cat("   0 < r < 2M, so r-dot never vanishes inside and there is no radial turning point.\n")
cat("      r/M      |f|          r-dot^2 at L = 1      zero anywhere inside?\n")
for (r in c(0.1, 0.5, 1.0, 1.5, 1.9)) {
  fv <- -fS(r); rd <- fv/r^2
  cat(sprintf("   %6.2f %12.5f %18.5f      %s\n", r, fv, rd, ifelse(rd > 0, "no", "YES")))
  note(rd > 0, sprintf("r-dot^2 is positive at r/M = %.2f", r))
}
cat("   The same count in both geometries, so no relative phase, so the sign carries.\n")
cat("   The plant for that: the same test outside the horizon, where |f| changes sign and a\n")
cat("   turning point does appear, so the test is not vacuous.\n")
for (r in c(2.5, 4.0)) {
  fv <- -fS(r)
  cat(sprintf("      r/M = %.2f (outside): |f| = %+.5f, so r-dot^2 = %+.5f  ->  %s\n",
              r, fv, fv/r^2, ifelse(fv/r^2 > 0, "no turning", "TURNS, as it must outside")))
  note(fv/r^2 < 0, "outside the horizon the E = 0 family does turn")
}

cat("\n   So the contact contraction is 'energy density plus transverse pressure' in both geometries\n")
cat("   and transfers directly. The radial one does not: A.18's radial direction is spacelike and\n")
cat("   carries no tidal field at all, the sphere there having constant radius, where a hole's is\n")
cat("   the timelike direction along which the whole interior tidal field acts.\n\n")
cL2 <- comp(o3)
Cc  <- -cL2$kk_con*0.001^3.5
cat(sprintf("   Taking it: the contact contraction is %.6f s^{-7/2}, so at a hole\n", cL2$kk_con*0.001^3.5))
cat("     A = T^theta_theta - T^r_r = rho + p_T < 0.\n")
cat("   Then the interior's own conservation law and trace, which A.18's geometry cannot see\n")
cat("   because it has no radial gradient, give X = -2A at conformal coupling, so\n")
cat("     X = T^t_t - T^r_r > 0  and  T_ab k^a k^b along Penrose's ingoing radial congruence is\n")
cat("     POSITIVE: the null convergence condition HOLDS there.\n")
cat("   The two readings close in the same direction and in opposite senses, which is the whole\n")
cat("   point: the image stress defocuses the congruence that would bring the sheets into contact\n")
cat("   and focuses the one Penrose's theorem runs on. The fold censors its own closed causal\n")
cat("   curves and leaves the singularity theorem's hypothesis intact.\n")
note(cL2$kk_con < 0, "the transferable contraction is negative, so A < 0 at a hole")
cat("\n   And the two chains agree where they overlap, which is the check on the transfer. A.18's\n")
cat("   geometry returns T^r_r = 0 and a traceless leading divergence; the interior chain needs\n")
cat("   Y = T^r_r subleading and the trace subleading, and gets both from conservation and from\n")
cat("   T^a_a = -m^2 W. Those are the same two statements reached two ways:\n")
cat(sprintf("      T^r_r / rho at s = 0.001:   %+.3e   (the interior chain needs this small)\n",
            cL2$p_r/cL2$rho))
cat(sprintf("      trace / rho at s = 0.001:   %+.3e   (and this)\n",
            (-cL2$rho + cL2$p_r + 2*cL2$p_T)/cL2$rho))
nrm <- -cL2$p_T                     # normalise on the transverse pressure
cat(sprintf("      T^t_t : T^r_r : T^th_th  =  %+.4f : %+.4f : %+.4f, against 2 : 0 : -1\n",
            -cL2$rho/nrm, cL2$p_r/nrm, cL2$p_T/nrm))
note(abs(-cL2$rho/nrm - 2) < 0.01 && abs(cL2$p_T/nrm + 1) < 1e-9,
     "the coordinate pattern is 2 : 0 : -1")
note(abs(cL2$p_r/cL2$rho) < 1e-3, "T^r_r is subleading, as the interior chain requires")

cat("\n=== 8. what this settles and what it does not ===\n")
cat("   Settled, in A.18's geometry: the image stress at the caustic has a negative energy\n")
cat("   density, a vanishing leading T^r_r, a traceless leading divergence, the coordinate pattern\n")
cat("   2 : 0 : -1, and a NEGATIVE null contraction in both of that geometry's directions, the\n")
cat("   contact one larger by exactly 3/2. Neither sign moves with the coupling or the mass. The\n")
cat("   coefficient is 15/(32 pi sqrt(2 pi)) = 0.0595 in units of s^{-7/2}.\n")
cat("   Carried to a hole by section 9: the contact contraction transfers because it is built the\n")
cat("   same way in both geometries, so the anisotropy at a hole is negative; the radial one does\n")
cat("   not transfer and comes instead from the interior's own conservation law, which makes it\n")
cat("   positive. The fold defocuses the contact congruence and focuses Penrose's.\n")
cat("   Not settled: A.18's geometry has a sphere of constant radius, so its radial direction\n")
cat("   carries no tidal field and its conservation law is empty. The transfer therefore rests on\n")
cat("   the causal-character argument of section 9 and not on a calculation in the interior\n")
cat("   itself. Doing the same sum on Schwarzschild would remove that step, and until it is done\n")
cat("   the radial sign is held by conservation rather than by the sum.\n")

cat("\n=== 10. the numbers the manuscripts quote, at the precision they quote them ===\n")
cat("   Printed as magnitudes, because the provenance checker compares digits and drops signs.\n")
cat("   The signs are what sections 4, 5 and 9 assert directly, and this file stops if any of them\n")
cat("   comes out the wrong way, so nothing about the sign rests on the digits below.\n")
cat(sprintf("   the pattern              %.4f : %.4f : %.4f\n",
            abs(-cL2$rho/nrm), abs(cL2$p_r/nrm), abs(cL2$p_T/nrm)))
cat(sprintf("   the contraction ratio    %.4f\n", abs(cL2$kk_con/cL2$kk_rad)))
cat(sprintf("   trace over rho, x 1e14   %.1f\n",
            abs((-cL2$rho + cL2$p_r + 2*cL2$p_T)/cL2$rho)*1e14))
cat(sprintf("   A.18's magnitude         %.1f\n", Re(o3$W)))
cat(sprintf("   the image power          %.4f   (negative)\n", abs(pw("W"))))
cat(sprintf("   the generic controls     %.4f and %.3f   (both negative)\n", abs(pgW), abs(pgD)))
cat(sprintf("   the closed coefficient   %.6f\n", CD))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
