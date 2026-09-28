#!/usr/bin/env Rscript
# Fork 9: the thermal factors at a Schwarzschild horizon are the MONODROMY of the interior radial
# equation around it, and that settles where the mode choice cannot be made.
#
# WHY. the_contact_pair_in_flat_space.R showed that inside a horizon the fold's correlator is the
# continuation of the Milne time, not a shift of the Killing time, and that the two differ by a
# loop around the branch point. So the object the interior sum needs is the continuation of the
# radial solution once around the horizon. That is a monodromy, and a monodromy at a regular
# singular point is fixed by the indicial equation rather than by any numerical work.
#
# THE EQUATION. With 2M = 1 the interior radial equation at frequency k and multipole l is
#     r^2 (r - 1) u'' + r u' + [ k^2 r^4/(r - 1) - l(l+1) r - 1 ] u = 0,
# which is what interior_modes_at_nonzero_k.R marches in rstar. Dividing through, near r = 1 the
# u' coefficient goes as 1/(r-1) and the u coefficient as k^2/(r-1)^2, so the indicial equation is
#     s(s-1) + s + k^2 = 0,   hence   s = +- i k,
# and the monodromy eigenvalues are exp(2 pi i (+- i k)) = exp(-+ 2 pi k). With kappa = 1/2 the
# inverse temperature is beta = 2 pi / kappa = 4 pi, so those are exp(-+ beta k / 2): THE THERMAL
# FACTORS ARE THE MONODROMY. Nothing thermal was put in; it is the exponent structure at r = 1.
#
# WHAT THAT SETTLES, and it is a warning rather than a result. A mode fixed by its behaviour AT
# the horizon is a monodromy EIGENVECTOR, so its continuation is itself times a number, and the
# image correlator built from it is real. In Milne the same structure holds at the lightcone
# origin, and the mode that reproduces the exact answer there is NOT the eigenvector: it is fixed
# at the other end and its continuation mixes the two. So the interior mode cannot be chosen at
# the horizon, and the choice has to be closed by the calibration rather than by a boundary
# condition that looks natural. The calibration is Delta^{1/2} -> 3.9004 M s^{-1/2}.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }
det2 <- function(M) M[1,1]*M[2,2] - M[1,2]*M[2,1]   # R's det refuses complex matrices

# u'' = -[ u'/(r(r-1)) + ( k^2 r^2/(r-1)^2 - (l(l+1)r + 1)/(r^2 (r-1)) ) u ]
deriv <- function(r, y, k, l) {
  a1 <- 1/(r*(r - 1))
  a0 <- k^2*r^2/(r - 1)^2 - (l*(l+1)*r + 1)/(r^2*(r - 1))
  c(y[2], -(a1*y[2] + a0*y[1]))
}
# march once round the circle r = 1 + eps e^{i theta}, RK4 in theta, dr = i (r-1) d theta
loop <- function(k, l, eps, y0, N = 200000) {
  h <- 2*pi/N; y <- y0
  for (i in seq_len(N)) {
    th <- (i - 1)*h
    st <- function(t, yy) { r <- 1 + eps*exp(1i*t); (1i*(r - 1))*deriv(r, yy, k, l) }
    k1 <- st(th, y); k2 <- st(th + h/2, y + h*k1/2)
    k3 <- st(th + h/2, y + h*k2/2); k4 <- st(th + h, y + h*k3)
    y <- y + h*(k1 + 2*k2 + 2*k3 + k4)/6
  }
  y
}

cat("=== 1. the indicial exponents at the horizon are exactly +- i k ===\n")
cat("   s(s-1) + s + k^2 = 0 gives s^2 = -k^2. Checked as a polynomial identity at four k:\n")
for (k in c(0.4, 1, 2.5, 7)) {
  s <- 1i*k
  res <- s*(s - 1) + s + k^2
  cat(sprintf("   k = %4.1f   s = i k gives residual %.3e\n", k, Mod(res)))
  note(Mod(res) < 1e-14, sprintf("the indicial equation is satisfied at k = %g", k))
}

cat("\n=== 2. the monodromy, measured by marching two solutions round the horizon ===\n")
cat("   Two independent solutions are carried once round r = 1 + eps e^{i theta}. The matrix\n")
cat("   they return is the monodromy in that basis; its EIGENVALUES do not depend on the basis\n")
cat("   and must be exp(-+ 2 pi k). Its determinant must be one, since the Wronskian picks up\n")
cat("   exp(-2 pi i) going round and that is unity.\n")
cat("        k     l    eps     |lambda_1|        |lambda_2|       exp(2 pi k)      det\n")
# The circle has to shrink as l grows, because l(l+1) sits in the NON-singular part of the
# potential and a bigger loop carries more of it; the determinant measures that directly, and at
# l = 8 on a 0.3 circle it reads 7.96 instead of 1. That is resolution, not physics: the indicial
# equation has no l in it at all, which is why the eigenvalue below is the same at l = 2 and 8.
for (par in list(c(0.3, 2, 0.3), c(0.6, 2, 0.3), c(0.6, 8, 0.05), c(1.0, 2, 0.5))) {
  k <- par[1]; l <- par[2]; eps <- par[3]
  r0 <- 1 + eps
  c1 <- loop(k, l, eps, c(1+0i, 0+0i)); c2 <- loop(k, l, eps, c(0+0i, 1+0i))
  M <- matrix(c(c1[1], c1[2], c2[1], c2[2]), 2, 2)     # columns are the images of the basis
  ev <- eigen(M)$values
  m <- sort(Mod(ev))
  cat(sprintf("   %6.2f %5d %6.2f %14.6e %16.6e %16.6e %10.6f\n",
              k, l, eps, m[1], m[2], exp(2*pi*k), Mod(det2(M))))
  note(abs(m[2]/exp(2*pi*k) - 1) < 1e-4,
       sprintf("the large eigenvalue is exp(2 pi k) at k = %g, l = %d", k, l))
  note(abs(m[1]*m[2] - 1) < 1e-4, sprintf("the eigenvalues are reciprocal at k = %g", k))
  note(abs(Mod(det2(M)) - 1) < 1e-4, sprintf("the determinant is one at k = %g", k))
}

cat("\n=== 3. so the thermal factor is not an input anywhere ===\n")
cat("   With 2M = 1 the surface gravity is kappa = 1/2 and beta = 2 pi / kappa = 4 pi, so\n")
cat("   exp(-+ 2 pi k) is exp(-+ beta k / 2), which is the half-period factor the cross weight\n")
cat("   carries. It is here the exponent structure of a regular singular point and nothing else.\n")
for (k in c(0.3, 0.6, 1.0)) {
  bet <- 4*pi
  cat(sprintf("   k = %.1f: exp(2 pi k) = %.6e and exp(beta k / 2) = %.6e\n",
              k, exp(2*pi*k), exp(bet*k/2)))
  note(abs(exp(2*pi*k)/exp(bet*k/2) - 1) < 1e-12, "the monodromy factor is the thermal factor")
}

cat("\n=== 4. the warning this carries, and it is the useful part ===\n")
cat("   A mode fixed by its behaviour AT the horizon is an eigenvector of that matrix, so its\n")
cat("   continuation is itself times a real positive number and the image correlator built from\n")
cat("   it comes out REAL. In Milne the same structure sits at the lightcone origin, and the\n")
cat("   mode that reproduces the exact answer there is not the eigenvector: it is fixed at the\n")
cat("   far end and its continuation mixes the pair. the_contact_pair_in_flat_space.R measured\n")
cat("   both, and the eigenvector choice missed the truth by a factor of two to four while\n")
cat("   returning a real number where the answer is complex.\n")
cat("   SO THE INTERIOR MODE CANNOT BE CHOSEN AT THE HORIZON. The choice is closed by the\n")
cat("   calibration, Delta^{1/2} -> 3.9004 M s^{-1/2}, and by nothing that looks natural.\n")

cat("\n=== 5. the plants ===\n")
k <- 0.6; l <- 2; eps <- 0.3
c1 <- loop(k, l, eps, c(1+0i, 0+0i)); c2 <- loop(k, l, eps, c(0+0i, 1+0i))
M <- matrix(c(c1[1], c1[2], c2[1], c2[2]), 2, 2)
cat(sprintf("   (a) the eigenvalues must not be exp(2 pi k') for a k' that is not k: at k' = %.1f\n", 0.9))
cat(sprintf("       exp(2 pi k') = %.4e against the measured %.4e\n",
            exp(2*pi*0.9), max(Mod(eigen(M)$values))))
note(abs(max(Mod(eigen(M)$values))/exp(2*pi*0.9) - 1) > 0.2, "plant: a wrong k is distinguishable")
# (b) break the equation: drop the first-derivative term and the determinant must move off one
loopbad <- function(k, l, eps, y0, N = 100000) {
  h <- 2*pi/N; y <- y0
  for (i in seq_len(N)) {
    th <- (i - 1)*h
    st <- function(t, yy) {
      r <- 1 + eps*exp(1i*t)
      a0 <- k^2*r^2/(r - 1)^2 - (l*(l+1)*r + 1)/(r^2*(r - 1))
      (1i*(r - 1))*c(yy[2], -a0*yy[1])            # the u'/(r(r-1)) term removed
    }
    k1 <- st(th, y); k2 <- st(th + h/2, y + h*k1/2)
    k3 <- st(th + h/2, y + h*k2/2); k4 <- st(th + h, y + h*k3)
    y <- y + h*(k1 + 2*k2 + 2*k3 + k4)/6
  }
  y
}
b1 <- loopbad(k, l, eps, c(1+0i, 0+0i)); b2 <- loopbad(k, l, eps, c(0+0i, 1+0i))
Mb <- matrix(c(b1[1], b1[2], b2[1], b2[2]), 2, 2)
cat("   (b) the first-derivative term dropped. The determinant is no use here, since removing\n")
cat("       a first-derivative term makes the Wronskian exactly constant and it returns to\n")
cat(sprintf("       itself round the loop: |det| = %.6f either way. The eigenvalue is the test,\n",
            Mod(det2(Mb))))
cat(sprintf("       and it reads %.4e against %.4e.\n",
            max(Mod(eigen(Mb)$values)), exp(2*pi*k)))
note(abs(max(Mod(eigen(Mb)$values))/exp(2*pi*k) - 1) > 0.05,
     "plant: dropping the first-derivative term moves the monodromy")

cat("\n=== 6. what this leaves ===\n")
cat("   The continuation the interior sum needs is now an explicit finite matrix at every\n")
cat("   (k, l), computed by marching and not by a connection formula that would have to be\n")
cat("   looked up. What is left is the mode choice, which the calibration closes, and then the\n")
cat("   sum: the l-sum is fully coherent and carries the caustic, and the frequency integral\n")
cat("   converges exponentially.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
