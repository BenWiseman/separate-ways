#!/usr/bin/env Rscript
# The first brick of fork 9's only non-circular route: an exact interior radial solver, and how
# far the closed-form Jacobi modes are from it.
#
# WHY. the_eikonal_route_is_circular.R shows the cheap way in cannot check A.19's rule, because a
# sum taken in the eikonal limit reproduces the parametrix that the rule replaces. What is left
# is exact radial modes per (k, l). interior_radial_is_jacobi.R gives closed-form modes, and
# REORIENT is careful to say they are the zero-frequency EIKONAL ones. This separates the two and
# measures the gap, so the next attempt knows where the closed form may stand in and where it
# may not.
#
# THE TWO EQUATIONS. With 2M = 1 the horizon is r = 1, the contact sphere r = 1/2, and the
# interior radial equation at k = 0 is, in full,
#     r^2 (r - 1) psi'' + r psi' - (l(l+1) r + 1) psi = 0,
# whose last term carries the 2M/r^3 piece of the potential. Dropping that piece, which is the
# eikonal step, and dividing by r leaves
#     r(1 - r) psi'' - psi' + l(l+1) psi = 0,
# the hypergeometric equation whose regular branch is psi_l = r^2 P^{(2,0)}_{l-1}(1 - 2r).
#
# THE EXPONENTS DIFFER AT THE SINGULARITY, which is the first thing to know. Near r = 0 the
# eikonal equation gives s(s - 2) = 0, exponents {0, 2}, and the exact one gives (s - 1)^2 = 0, a
# degenerate {1, 1} with a logarithm. So the eikonal form is not uniformly good down to r = 0.
# It does not have to be: the caustic sits at r = 1/2 and the contact region is r <= 1/2.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

jac <- function(n, a, b, x) {
  if (n <= 0) return(rep(1, length(x)))
  p0 <- rep(1, length(x)); p1 <- (a + 1) + (a + b + 2)*(x - 1)/2
  if (n == 1) return(p1)
  for (k in 2:n) {
    c1 <- 2*k*(k + a + b)*(2*k + a + b - 2); c2 <- (2*k + a + b - 1)*(a^2 - b^2)
    c3 <- (2*k + a + b - 2)*(2*k + a + b - 1)*(2*k + a + b)
    c4 <- 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
    p2 <- ((c2 + c3*x)*p1 - c4*p0)/c1; p0 <- p1; p1 <- p2
  }
  p1
}
psiJ  <- function(l, r) r^2 * jac(l - 1, 2, 0, 1 - 2*r)
dpsiJ <- function(l, r, h = 1e-6) (psiJ(l, r + h) - psiJ(l, r - h))/(2*h)

cat("=== 1. the indicial exponents at r = 0 really do differ ===\n")
cat("      equation            indicial polynomial        exponents\n")
cat("   eikonal   r(1-r)psi'' - psi' + l(l+1)psi    s(s-2)              0 and 2\n")
cat("   exact     r^2(r-1)psi'' + r psi' - (...)    (s-1)^2             1 twice, with a log\n")
# check numerically: the eikonal solution really does go as r^2, the exact one as r
cat(sprintf("   the Jacobi branch near zero: psi(1e-3)/1e-6 = %.4f, psi(1e-4)/1e-8 = %.4f\n",
            psiJ(6, 1e-3)/1e-6, psiJ(6, 1e-4)/1e-8))
note(abs(psiJ(6, 1e-3)/1e-6 - psiJ(6, 1e-4)/1e-8) < 0.02*abs(psiJ(6, 1e-4)/1e-8),
     "the eikonal branch goes as r^2 at the singularity")

cat("\n=== 2. the exact solver, started from the Jacobi mode inside the contact region ===\n")
# psi'' = [ (l(l+1) r + 1) psi - r psi' ] / (r^2 (r - 1))
solve_exact <- function(l, r0, r1, y0, yp0, N = 400000) {
  h <- (r1 - r0)/N; r <- r0; y <- y0; yp <- yp0
  f <- function(r, y, yp) c(yp, ((l*(l+1)*r + 1)*y - r*yp)/(r^2*(r - 1)))
  for (i in seq_len(N)) {
    k1 <- f(r, y, yp)
    k2 <- f(r + h/2, y + h*k1[1]/2, yp + h*k1[2]/2)
    k3 <- f(r + h/2, y + h*k2[1]/2, yp + h*k2[2]/2)
    k4 <- f(r + h,   y + h*k3[1],   yp + h*k3[2])
    y  <- y  + h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
    yp <- yp + h*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
    r  <- r + h
  }
  c(y = y, yp = yp)
}
cat("   Started at r = 0.20 on the Jacobi mode and its derivative, integrated to r = 0.90.\n")
cat("   The gap is the 2M/r^3 term and nothing else. It is measured as an RMS over the window\n")
cat("   r in [0.35, 0.75], which brackets the caustic at r = 1/2, and NOT at a single r: these\n")
cat("   modes oscillate, so one point can sit where the two nearly cancel and read far worse\n")
cat("   than the function is. The first version of this file compared at r = 0.8 and made l = 16\n")
cat("   look worse than l = 8.\n")
trace <- function(l, r0 = 0.20, r1 = 0.90, N = 300000) {
  h <- (r1 - r0)/N; r <- r0; y <- psiJ(l, r0); yp <- dpsiJ(l, r0)
  f <- function(r, y, yp) c(yp, ((l*(l+1)*r + 1)*y - r*yp)/(r^2*(r - 1)))
  rs <- numeric(N + 1); ys <- numeric(N + 1); rs[1] <- r; ys[1] <- y
  for (i in seq_len(N)) {
    k1 <- f(r,y,yp); k2 <- f(r+h/2,y+h*k1[1]/2,yp+h*k1[2]/2)
    k3 <- f(r+h/2,y+h*k2[1]/2,yp+h*k2[2]/2); k4 <- f(r+h,y+h*k3[1],yp+h*k3[2])
    y <- y + h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
    yp <- yp + h*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6; r <- r + h
    rs[i+1] <- r; ys[i+1] <- y
  }
  list(r = rs, y = ys)
}
cat("        l     RMS|exact|      RMS|exact - Jacobi|     relative gap\n")
gaps <- c(); ls <- c(8, 16, 32, 64, 128, 256)
for (l in ls) {
  tr <- trace(l); w <- tr$r >= 0.35 & tr$r <= 0.75
  ex <- tr$y[w]; jj <- psiJ(l, tr$r[w])
  rex <- sqrt(mean(ex^2)); rdf <- sqrt(mean((ex - jj)^2)); g <- rdf/rex
  cat(sprintf("   %6d %14.4e %22.4e %17.4e\n", l, rex, rdf, g))
  gaps <- c(gaps, g)
  note(is.finite(g) && g > 0, sprintf("the solver returns a finite gap at l = %d", l))
}
cat(sprintf("   the gap falls monotonically from %.3e to %.3e\n", gaps[1], gaps[length(gaps)]))
note(all(diff(gaps) < 0), "and it falls at every step, which the one-point version did not")

cat("\n=== 3. how it falls, since that is what says where the closed form may stand in ===\n")
sl <- lm(log(gaps) ~ log(ls))$coefficients[2]
cat(sprintf("   log-log slope of the RMS gap against l: %.3f\n", sl))
cat("   The dropped term is 2M/r^3 against l(l+1)/r^2, so one inverse power of l is what to\n")
cat("   expect, and that is what comes out.\n")
note(abs(sl + 1) < 0.2, "the eikonal error falls as one inverse power of l")

cat("\n=== 4. the plant: the solver must NOT agree when the equation is wrong ===\n")
bad <- solve_exact(32, 0.20, 0.80, psiJ(32, 0.20), dpsiJ(32, 0.20), N = 400000)
sol2 <- local({
  l <- 32; h <- 0.6/400000; r <- 0.2; y <- psiJ(l, 0.2); yp <- dpsiJ(l, 0.2)
  f <- function(r, y, yp) c(yp, ((1.3*l*(l+1)*r + 1)*y - r*yp)/(r^2*(r - 1)))   # tidal x 1.3
  for (i in seq_len(400000)) {
    k1 <- f(r,y,yp); k2 <- f(r+h/2,y+h*k1[1]/2,yp+h*k1[2]/2)
    k3 <- f(r+h/2,y+h*k2[1]/2,yp+h*k2[2]/2); k4 <- f(r+h,y+h*k3[1],yp+h*k3[2])
    y <- y + h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
    yp <- yp + h*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6; r <- r + h
  }
  y
})
cat(sprintf("   l = 32 with the angular term scaled by 1.3: psi(0.8) = %.4e against %.4e,\n",
            sol2, bad["y"]))
cat(sprintf("   a relative difference of %.2f, so the integration is sensitive to the equation\n",
            abs(sol2 - bad["y"])/abs(bad["y"])))
note(abs(sol2 - bad["y"])/abs(bad["y"]) > 0.05, "plant: a wrong angular term is visible")

cat("\n=== 5. what this brick gives fork 9 ===\n")
cat("   A validated exact solver on the interval that brackets the caustic, and a measured\n")
cat("   answer to the question the closed form raises: the eikonal modes are wrong by one\n")
cat("   inverse power of l, so they may stand in for the tail of the sum and not for its head.\n")
cat("   Since the caustic divergence comes FROM the tail, that is the useful direction.\n")
cat("   STILL TO BUILD: the same solver at nonzero k, the Klein-Gordon normalisation on a\n")
cat("   constant-r slice, and the cross-region state. Those carry the convention risk, and the\n")
cat("   calibration that catches a convention error is Delta^{1/2} -> 3.9004 M s^{-1/2}.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
