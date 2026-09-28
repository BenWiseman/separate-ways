#!/usr/bin/env Rscript
# nariai_stability.R -- A.15 finds exactly one fold-invariant hole in the
# Schwarzschild-de Sitter family, the Nariai one. A configuration that does not persist is
# no object at all, so this asks whether it persists, and then whether the fold changes the
# answer. It does.
#
# Spherical reduction. For  ds^2 = ghat_ab dx^a dx^b + R(x)^2 dOmega^2  the vacuum
# equations with Lambda, R_munu = Lambda g_munu, split into
#   (A)  Rhat_ab - (2/R) grad_a grad_b R = Lambda ghat_ab
#   (B)  1 - (grad R)^2 - R grad^2 R = Lambda R^2
# Nariai is R = l constant with l^2 = 1/Lambda, so (B) gives 1/l^2 = Lambda and (A) gives
# Rhat_ab = Lambda ghat_ab, the two-dimensional de Sitter of radius l.

lam <- 1/9; l <- 1/sqrt(lam)
cat(sprintf("=== Nariai background: Lambda = %.6f, l = 1/sqrt(Lambda) = %.6f ===\n", lam, l))

cat("\n=== 1. the perturbation equations ===\n")
cat("   Write R = l(1 + phi). Linearising (B):\n")
cat("      1 - l^2 (grad phi)^2 - l^2 (1+phi) grad^2 phi = Lambda l^2 (1+phi)^2\n")
cat("   drops to   grad^2 phi = -2 Lambda phi,   using Lambda l^2 = 1.\n")
cat("   In two dimensions Rhat_ab = (Rhat/2) ghat_ab identically, so the left side of (A)\n")
cat("   minus Lambda ghat_ab is pure trace and the traceless part of grad_a grad_b phi must\n")
cat("   vanish. With the trace fixed by the line above,\n")
cat("      grad_a grad_b phi = -Lambda ghat_ab phi.\n")
cat("   That is a strong constraint: three solutions, not a function's worth.\n")

cat("\n=== 2. checking both, in global coordinates ===\n")
# global dS_2:  ds^2 = -dt^2 + l^2 cosh^2(t/l) dtheta^2
# embedding:    X0 = l sinh(t/l), X1 = l cosh(t/l) cos(th), X2 = l cosh(t/l) sin(th)
X <- function(t, th) c(l*sinh(t/l), l*cosh(t/l)*cos(th), l*cosh(t/l)*sin(th))
ghat <- function(t) diag(c(-1, l^2*cosh(t/l)^2))
# Christoffels: G^t_thth = l cosh sinh ; G^th_t th = tanh(t/l)/l ; rest zero
box <- function(f, t, th, h=2e-4) {          # grad^2 f = -f_tt - tanh/l f_t + f_thth/(l cosh)^2
  ft  <- (f(t+h,th)-f(t-h,th))/(2*h); ftt <- (f(t+h,th)-2*f(t,th)+f(t-h,th))/h^2
  fpp <- (f(t,th+h)-2*f(t,th)+f(t,th-h))/h^2
  -ftt - tanh(t/l)/l*ft + fpp/(l*cosh(t/l))^2 }
hess <- function(f, t, th, h=2e-4) {
  ft  <- (f(t+h,th)-f(t-h,th))/(2*h); fp <- (f(t,th+h)-f(t,th-h))/(2*h)
  ftt <- (f(t+h,th)-2*f(t,th)+f(t-h,th))/h^2
  fpp <- (f(t,th+h)-2*f(t,th)+f(t,th-h))/h^2
  ftp <- (f(t+h,th+h)-f(t+h,th-h)-f(t-h,th+h)+f(t-h,th-h))/(4*h^2)
  Gt_pp <- l*cosh(t/l)*sinh(t/l); Gp_tp <- tanh(t/l)/l
  matrix(c(ftt, ftp - Gp_tp*fp, ftp - Gp_tp*fp, fpp - Gt_pp*ft), 2, 2) }

modes <- list(`phi = X0/l`=function(t,th) X(t,th)[1]/l,
              `phi = X1/l`=function(t,th) X(t,th)[2]/l,
              `phi = X2/l`=function(t,th) X(t,th)[3]/l)
cat("      mode          max |grad^2 phi + 2 Lam phi|    max |Hess + Lam ghat phi|\n")
set.seed(3); ts <- runif(40,-2,2); ths <- runif(40,0,2*pi)
for (nm in names(modes)) { f <- modes[[nm]]
  e1 <- max(mapply(function(t,th) abs(box(f,t,th) + 2*lam*f(t,th)), ts, ths))
  e2 <- max(mapply(function(t,th) max(abs(hess(f,t,th) + lam*ghat(t)*f(t,th))), ts, ths))
  cat(sprintf("   %-14s  %24.1e  %24.1e\n", nm, e1, e2)) }
cat("   Residuals are the finite-difference floor at this step, not a discrepancy: halving\n")
cat("   and doubling the step moves them the way truncation error moves. All three solve\n")
cat("   both equations, so they are the whole perturbation space. The exact statement is\n")
cat("   the standard one, that phi = c.X for constant c in R^{1,2} restricted to the\n")
cat("   hyperboloid obeys grad_a grad_b phi = -(1/l^2) ghat_ab phi.\n")
cat("\n   step dependence, for phi = X0/l:\n")
for (hh in c(1e-3, 3e-4, 1e-4)) {
  f <- modes[[1]]
  cat(sprintf("     h = %.0e  ->  |grad^2 phi + 2 Lam phi| = %.1e\n", hh,
              abs(box(f, 0.7, 1.1, h=hh) + 2*lam*f(0.7,1.1)))) }

cat("\n   The check has to be able to fail, so plant modes that should not solve it:\n")
bad <- list(`phi = (X0/l)^2`  = function(t,th) (X(t,th)[1]/l)^2,
            `phi = X0 X2/l^2` = function(t,th) X(t,th)[1]*X(t,th)[3]/l^2,
            `phi = cos(t/l)`  = function(t,th) cos(t/l))
for (nm in names(bad)) { f <- bad[[nm]]
  e1 <- abs(box(f,0.7,1.1) + 2*lam*f(0.7,1.1))
  e2 <- max(abs(hess(f,0.7,1.1) + lam*ghat(0.7)*f(0.7,1.1)))
  cat(sprintf("     %-18s  wave %8.3f   hessian %8.3f   <- fails, as it must\n", nm, e1, e2))
  stopifnot(e1 > 1e-3 || e2 > 1e-3) }

cat("\n=== 3. which of them grow ===\n")
cat("   The static patch is X1 > |X0|, and there X2^2 = l^2 + X0^2 - X1^2 < l^2, so X2/l is\n")
cat("   bounded by one while X0 and X1 are not. In static coordinates, with\n")
cat("   X0 = sqrt(l^2-rho^2) sinh(t/l) and X1 = sqrt(l^2-rho^2) cosh(t/l) and X2 = rho:\n\n")
cat("        static t/l      X0/l        X1/l        X2/l    (at rho = 0.5 l)\n")
rho <- 0.5*l; s <- sqrt(l^2-rho^2)/l
for (tt in c(0,1,2,4,8)) cat(sprintf("   %12.1f  %10.3f  %10.3f  %10.3f\n",
                                     tt, s*sinh(tt), s*cosh(tt), rho/l))
cat("\n   phi = X0/l and phi = X1/l grow without bound in static time: delta R / R diverges,\n")
cat("   which is the classical Nariai instability, the sphere swelling on one side and\n")
cat("   shrinking on the other as the solution runs away toward Schwarzschild-de Sitter.\n")
cat("   phi = X2/l is bounded by one and is the static near-Nariai deformation, no runaway.\n")

cat("\n=== 4. what the fold does to them ===\n")
cat("   The fold at Nariai is Nl(X0,X1,X2) = (-X0,-X1,X2) composed with P_perp on the\n")
cat("   sphere (`checks/calc/fold_map_classification.R`). A perturbation descends to the\n")
cat("   quotient only if it agrees at identified points, and the sphere radius is a\n")
cat("   function on the dS_2 factor alone, so the condition is phi(Nl x) = phi(x).\n\n")
cat("      mode          phi(Nl x)/phi(x)     survives the quotient\n")
for (nm in names(modes)) { f <- modes[[nm]]
  t <- 0.7; th <- 1.1; v <- f(t,th)
  Y <- X(t,th); Yr <- c(-Y[1],-Y[2],Y[3]); vr <- switch(nm, `phi = X0/l`=Yr[1]/l,
                        `phi = X1/l`=Yr[2]/l, `phi = X2/l`=Yr[3]/l)
  cat(sprintf("   %-14s  %16.3f     %s\n", nm, vr/v, if (abs(vr/v - 1) < 1e-12) "yes" else "NO")) }

cat("
=== flatly ===

  Unfolded Nariai is unstable. Two of its three perturbation modes grow without bound in
  static time, which is the known decay toward Schwarzschild-de Sitter.

  Folded Nariai is not, in this sector. Both growing modes are odd under the fold map and
  neither descends to the quotient. The one mode that survives is bounded by one and is
  the static near-Nariai deformation. So the fold removes exactly the modes that destroy
  the configuration, and the single fold-invariant hole in the family is also the one the
  fold makes stable.

  What this is NOT. It is linear stability in the spherically symmetric sector, which is
  where the classical Nariai instability lives, and nothing more. Non-spherical and
  gravitational-wave modes are not treated here, and neither is the semiclassical decay of
  the Ginsparg-Perry kind, which is a different question with a different method.\n")
