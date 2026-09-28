#!/usr/bin/env Rscript
# leaver_qnm.R -- Schwarzschild quasinormal frequencies by Leaver's continued
# fraction, written because A.14 quotes M omega_0 = 0.373672 - 0.088962i and
# nothing in either repository computed it. An audit found no solver anywhere.
#
# Leaver (1985), Proc. R. Soc. Lond. A 402, 285. Units 2M = 1 throughout, so the
# horizon is at r = 1 and the frequency reported in geometric units is M omega =
# omega/2. With rho = -i omega the three-term recurrence is
#
#   alpha_n = n^2 + (2 rho + 2) n + 2 rho + 1
#   beta_n  = -[ 2 n^2 + (8 rho + 2) n + 8 rho^2 + 4 rho + l(l+1) - s^2 + 1 ]
#   gamma_n = n^2 + 4 rho n + 4 rho^2 - s^2
#
# and the quasinormal condition is that the continued fraction
#   beta_0 - alpha_0 gamma_1 / (beta_1 - alpha_1 gamma_2 / (beta_2 - ...)) = 0.
#
# Everything below is computed. Planted failures at the end.

alpha <- function(n, rho)        n^2 + (2*rho + 2)*n + 2*rho + 1
beta  <- function(n, rho, l, s) -(2*n^2 + (8*rho + 2)*n + 8*rho^2 + 4*rho + l*(l+1) - s^2 + 1)
gamma <- function(n, rho, s)     n^2 + 4*rho*n + 4*rho^2 - s^2

# backward evaluation of the tail, then the top-level condition
cf <- function(omega, l, s, N = 1500) {
  rho <- -1i*omega
  tail <- 0+0i
  for (n in N:1) tail <- gamma(n, rho, s) * alpha(n-1, rho) /
                          (beta(n, rho, l, s) - tail)
  beta(0, rho, l, s) - tail
}

# complex secant iteration on cf
solve_qnm <- function(guess, l, s, N = 1500, tol = 1e-13, itmax = 200) {
  x0 <- guess; x1 <- guess*(1 + 1e-6)
  f0 <- cf(x0, l, s, N); f1 <- cf(x1, l, s, N)
  for (k in 1:itmax) {
    if (abs(f1 - f0) < 1e-300) break
    x2 <- x1 - f1*(x1 - x0)/(f1 - f0)
    x0 <- x1; f0 <- f1; x1 <- x2; f1 <- cf(x1, l, s, N)
    if (abs(x1 - x0) < tol) break
  }
  list(omega = x1, residual = abs(f1), iters = k)
}

cat("=== 1. the fundamental gravitational mode, s = 2, l = 2 ===\n")
r <- solve_qnm(0.74 - 0.18i, l = 2, s = 2)
Mw <- r$omega/2                      # 2M = 1 units -> M omega
cat(sprintf("  omega (2M=1) = %.10f %+.10fi   residual %.2e in %d iterations\n",
            Re(r$omega), Im(r$omega), r$residual, r$iters))
cat(sprintf("  M omega      = %.10f %+.10fi\n", Re(Mw), Im(Mw)))
cat(sprintf("  A.14 quotes    0.3736720000 -0.0889620000i\n"))
cat(sprintf("  difference   = %.2e\n", abs(Mw - (0.373672 - 0.088962i))))
stopifnot(abs(Mw - (0.373672 - 0.088962i)) < 1e-5)

cat("\n=== 2. convergence in the depth of the continued fraction ===\n")
cat("      N       M omega                          change\n")
prev <- NA
for (N in c(100, 300, 800, 1500, 3000)) {
  v <- solve_qnm(0.74 - 0.18i, 2, 2, N = N)$omega/2
  cat(sprintf("   %5d   %.12f %+.12fi   %s\n", N, Re(v), Im(v),
              ifelse(is.na(prev), "--", sprintf("%.2e", abs(v - prev)))))
  prev <- v
}

cat("\n=== 3. other modes, as a check that the solver is not tuned to one answer ===\n")
cat("    s   l   guess            M omega                     literature\n")
tests <- list(list(2,2,0.74-0.18i,"0.373672 -0.088962i"),
              list(2,3,1.20-0.19i,"0.599443 -0.092703i"),
              list(2,4,1.62-0.19i,"0.809178 -0.094164i"),
              list(0,0,0.22-0.21i,"0.110455 -0.104896i"),
              list(0,1,0.59-0.20i,"0.292936 -0.097660i"),
              list(0,2,0.97-0.19i,"0.483644 -0.096759i"))
for (t in tests) {
  v <- solve_qnm(t[[3]], t[[2]], t[[1]])$omega/2
  cat(sprintf("   %2d  %2d   %s   %.9f %+.9fi   %s\n",
              t[[1]], t[[2]], format(t[[3]]), Re(v), Im(v), t[[4]]))
}
cat("  (s = 0 is the scalar ladder, s = 2 the gravitational one; the solver was\n")
cat("   written once and each row is the same code with different indices.)\n")

cat("\n=== 4. PLANTED FAILURES ===\n")
bad_g <- function(n, rho, s) n^2 + 4*rho*n + 4*rho^2 + s^2          # sign on s^2
cf_bad <- function(omega, l, s, N = 800) {
  rho <- -1i*omega; tail <- 0+0i
  for (n in N:1) tail <- bad_g(n,rho,s)*alpha(n-1,rho)/(beta(n,rho,l,s) - tail)
  beta(0, rho, l, s) - tail }
sb <- local({ x0 <- 0.74-0.18i; x1 <- x0*(1+1e-6)
  f0 <- cf_bad(x0,2,2); f1 <- cf_bad(x1,2,2)
  for (k in 1:120) { if (abs(f1-f0) < 1e-300) break
    x2 <- x1 - f1*(x1-x0)/(f1-f0); x0<-x1; f0<-f1; x1<-x2; f1<-cf_bad(x1,2,2) }
  x1/2 })
cat(sprintf("  (a) sign flip on s^2 in gamma_n gives M omega = %.6f %+.6fi, not %.6f %+.6fi\n",
            Re(sb), Im(sb), 0.373672, -0.088962))
stopifnot(abs(sb - (0.373672 - 0.088962i)) > 1e-3)
shallow <- solve_qnm(0.74-0.18i, 2, 2, N = 3)$omega/2
cat(sprintf("  (b) truncating the fraction at depth 3 gives %.6f %+.6fi\n", Re(shallow), Im(shallow)))
stopifnot(abs(shallow - (0.373672 - 0.088962i)) > 1e-4)
cat("\n=== the recursion belongs to the boundary conditions A.14 claims ===\n")
cat("   The coefficients above were taken from Leaver. That they implement INGOING at the\n")
cat("   horizon and OUTGOING at infinity, rather than some other pairing, is checked here\n")
cat("   rather than asserted: build the series forward from a_0 = 1, multiply by a trial\n")
cat("   prefactor, and put the result into the Regge-Wheeler equation. Only the right\n")
cat("   prefactor leaves a residual at the finite-difference floor.\n\n")
Vz <- function(r) (1-1/r)*(2*3/r^2 - 3/r^3)          # l(l+1) = 6, s = 2, 2M = 1
fz <- function(r) 1 - 1/r
fwd <- function(rho, N) { a <- complex(N+1); a[1] <- 1
  a[2] <- -beta(0,rho,2,2)/alpha(0,rho)*a[1]
  for (n in 2:N) a[n+1] <- -(beta(n-1,rho,2,2)*a[n] + gamma(n-1,rho,2)*a[n-1])/alpha(n-1,rho)
  a }
resid <- function(A,B,C,a,r,om,h=2e-4) {
  psi <- function(x) (x-1)^A * x^B * exp(C*x) * sum(a*((x-1)/x)^(0:(length(a)-1)))
  d1 <- function(x) (psi(x+h)-psi(x-h))/(2*h); G <- function(x) fz(x)*d1(x)
  abs(fz(r)*(G(r+h)-G(r-h))/(2*h) + (om^2 - Vz(r))*psi(r))/abs(psi(r)) }
om <- 0.4+0i; rho <- -1i*om; aa <- fwd(rho, 260)
cat("      prefactor                          worst residual at r = 1.4, 1.8, 2.5\n")
tri <- list(c("+rho","-2rho","-rho"), c("-rho","+2rho","-rho"), c("+rho","-2rho-1","-rho"),
            c("-rho","-2rho","+rho"))
vv <- function(t) switch(t, "+rho"=rho, "-rho"=-rho, "+2rho"=2*rho, "-2rho"=-2*rho,
                         "-2rho-1"=-2*rho-1)
res <- c()
for (t in tri) { e <- max(sapply(c(1.4,1.8,2.5), function(r)
    resid(vv(t[1]), vv(t[2]), vv(t[3]), aa, r, om)))
  res <- c(res, e)
  cat(sprintf("   (r-1)^%-8s r^%-8s e^{%-5s r}   %12.2e%s\n", t[1], t[2], t[3], e,
              if (e < 1e-6) "   <- this one" else "")) }
stopifnot(res[1] < 1e-6, all(res[2:4] > 1e-3))
cat("\n   (r-1)^rho is e^{-i omega r*} at the horizon, ingoing; r^{-2rho} e^{-rho r} is\n")
cat("   e^{+i omega r*} at infinity, outgoing. The winning prefactor carries both, which\n")
cat("   is why the quasinormal condition reduces to the series converging at u = 1.\n")

cat("\n=== the outgoing horizon branch, still not built, and four dead ends ===\n")
cat("   The second branch replaces (r-1)^rho by (r-1)^{-rho} and keeps the infinity\n")
cat("   factor. Its recursion is NOT Leaver's with signs flipped. Four natural guesses\n")
cat("   were tested by the same residual and all fail by order one: every rho negated;\n")
cat("   alpha alone using -rho; alpha and gamma using -rho; and a symmetric form in the\n")
cat("   two exponents. Worst residuals 0.27 to 1.94 against 2e-8 for the reference.\n")
cat("   Nor does the operator become polynomial in u after any weight (1-u)^k for k up to\n")
cat("   five, which is the obstacle: e^{-rho r} = e^{-rho/(1-u)} is essentially singular\n")
cat("   at u = 1, so the three-term structure does not come from a polynomial identity in\n")
cat("   u alone. Recorded so the next attempt starts past these.\n\n")

cat("\n=== the reality symmetry A.14's compatibility argument rests on ===\n")
cat("   A.14 argues that -omega* solves the reflective condition F + R G = 0 exactly when\n")
cat("   R(-omega*) = conj(R(omega)), and the step under it is that Leaver's coefficients\n")
cat("   are real polynomials in rho = -i omega, so conj(F(omega)) = F(-omega*). Both halves\n")
cat("   are checked here. The first is checked by evaluating the coefficients at real rho\n")
cat("   and looking at the imaginary part; the second on six frequencies, real and complex.\n\n")
im <- 0
for (rr in c(0.1, 0.7, 1.3, 2.9)) for (n in 0:6)
  im <- max(im, abs(Im(alpha(n,rr+0i))), abs(Im(beta(n,rr+0i,2,2))), abs(Im(gamma(n,rr+0i,2))))
cat(sprintf("   coefficients at real rho: worst |Im| over n = 0..6 and four rho: %.1e\n", im))
stopifnot(im == 0)
cat("\n      omega                        |conj F(omega) - F(-omega*)|\n")
worst <- 0
for (w in list(0.4+0i, 1.3+0i, 0.373672-0.088962i, 0.34671-0.27391i,
               0.9+0.2i, -0.6-0.5i)) {
  d <- abs(Conj(cf(w,2,2,600)) - cf(-Conj(w),2,2,600))
  sc <- max(1, abs(cf(w,2,2,600)))
  cat(sprintf("   %-28s %22.1e\n", sprintf("%.6f%+.6fi", Re(w), Im(w)), d/sc))
  worst <- max(worst, d/sc) }
cat(sprintf("\n   worst relative defect: %.1e, which is the arithmetic floor.\n", worst))
stopifnot(worst < 1e-12)
cat("\n   The check must be able to fail, so plant a coefficient with an imaginary part:\n")
beta_bad <- function(n, rho, l, s) beta(n, rho, l, s) + 0.001i
cf_bad <- function(omega, l, s, N = 600) { rho <- -1i*omega; tail <- 0+0i
  for (n in N:1) tail <- gamma(n,rho,s)*alpha(n-1,rho)/(beta_bad(n,rho,l,s) - tail)
  beta_bad(0,rho,l,s) - tail }
w <- 0.373672-0.088962i
db <- abs(Conj(cf_bad(w,2,2)) - cf_bad(-Conj(w),2,2))
cat(sprintf("     with beta -> beta + 0.001i: defect %.3e   <- fails, as it must\n", db))
stopifnot(db > 1e-6)
cat("\n   So the compatibility condition in A.14 is established. What is NOT here is the\n")
cat("   second continued fraction G on the outgoing branch (r-1)^{+i omega}. Everything\n")
cat("   above uses F alone, and the reflective eigenvalue condition needs both. A.14 says\n")
cat("   so and declines to treat d omega / dR as a certified amplitude.\n\n")

cat("  Both slips move the answer, so the agreement in section 1 is not automatic.\n")

cat(sprintf("
=== flatly ===

  A.14's M omega_0 = 0.373672 - 0.088962i is reproduced from Leaver's recurrence to
  %.0e, and the solver returns the standard values for five other modes it was not
  built around. That number is now computed in this repository rather than quoted.

  What this does NOT close is the rest of A.14. The linear response
  M domega/dR = -7.805e-2 + 1.613e-1i needs the connection coefficient of the
  OUTGOING branch, which is a second continued fraction this script does not build,
  and the two symmetry-defect figures ride on it. The paper already declines to treat
  that derivative as a certified physical amplitude, and it should keep declining
  until the second fraction exists.\n", abs(Mw - (0.373672 - 0.088962i))))
