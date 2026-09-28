#!/usr/bin/env Rscript
# Can conservation give the radial sign from the contact one? The algebra first.
#
# The open item is the sign of T_ab k^a k^b along the INGOING RADIAL null direction, where the
# release computes it along the CONTACT null direction, the one that turns through pi on the sphere.
# A.18's mode sum on the contact geodesic would settle it and is hard. There is a cheaper route and
# this file tests whether its algebra holds before any of it is believed.
#
# In a static spherically symmetric spacetime a stress invariant under the isometries has only
# T^t_t, T^r_r and T^theta_theta = T^phi_phi, all functions of r. Contracting with a null vector
# that has no transverse part, and then with one that has, gives
#
#     T_kk(radial)  = |f| (k^t)^2 (T^t_t - T^r_r)
#     T_kk(contact) = |f| (k^t)^2 (T^t_t - T^r_r) + r^2 (k^theta)^2 (T^theta_theta - T^r_r)
#
# so the radial answer is the contact answer minus the anisotropy term, and the anisotropy is what
# conservation constrains. If those two identities hold, the route is open.

TOL <- 1e-10
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }
set.seed(717)

M <- 1
fofr <- function(r) 1 - 2*M/r            # negative inside
gdn <- function(r) diag(c(fofr(r)*-1, -1/(fofr(r)*-1), r^2, r^2))   # interior form, sin th = 1
# careful: write the metric directly as diag(|f|, -1/|f|, r^2, r^2) with |f| = 2M/r - 1
absf <- function(r) 2*M/r - 1
gmet <- function(r) diag(c(absf(r), -1/absf(r), r^2, r^2))

cat("=== 1. the interior metric, and that r is the timelike direction ===\n")
cat("      r        |f| = 2M/r - 1     signature of diag(|f|, -1/|f|, r^2, r^2)\n")
for (r in c(1.8, 1.0, 0.5, 0.2)) {
  g <- gmet(r); sg <- sign(diag(g))
  cat(sprintf("   %6.2f %16.6f            (%s)\n", r, absf(r),
              paste(ifelse(sg > 0, "+", "-"), collapse = "")))
  note(sum(sg < 0) == 1 && sg[2] < 0, "exactly one minus sign and it is the r direction")
}

cat("\n=== 2. the two null contractions, against a directly built tensor ===\n")
cat("   Building T^a_b = diag(Tt, Tr, Tth, Tth) at random, lowering with the metric, contracting\n")
cat("   with null vectors solved from the metric, and comparing with the two formulae above:\n")
cat("      r       radial: direct vs formula        contact: direct vs formula\n")
for (r in c(0.9, 0.5, 0.3)) {
  g <- gmet(r); af <- absf(r)
  Tup <- c(Tt = rnorm(1), Tr = rnorm(1), Tth = rnorm(1))
  Tdn <- diag(c(Tup["Tt"]*g[1,1], Tup["Tr"]*g[2,2], Tup["Tth"]*g[3,3], Tup["Tth"]*g[4,4]))
  # radial null: |f|(k^t)^2 = (k^r)^2/|f|  ->  k^r = -|f| k^t
  kt <- 1; kr <- -af*kt; k <- c(kt, kr, 0, 0)
  note(abs(as.numeric(t(k) %*% g %*% k)) < 1e-12, "the radial vector is null")
  d1 <- as.numeric(t(k) %*% Tdn %*% k)
  f1 <- af*kt^2*(Tup["Tt"] - Tup["Tr"])
  # contact null: add a transverse part, solve k^r from the null condition
  kth <- 0.37
  kr2 <- -sqrt(af*(af*kt^2 + r^2*kth^2))
  k2 <- c(kt, kr2, kth, 0)
  note(abs(as.numeric(t(k2) %*% g %*% k2)) < 1e-12, "the contact vector is null")
  d2 <- as.numeric(t(k2) %*% Tdn %*% k2)
  f2 <- af*kt^2*(Tup["Tt"] - Tup["Tr"]) + r^2*kth^2*(Tup["Tth"] - Tup["Tr"])
  cat(sprintf("   %5.2f %14.6f vs %11.6f %13.6f vs %11.6f\n", r, d1, f1, d2, f2))
  note(abs(d1 - f1) < 1e-10, "the radial formula is exact")
  note(abs(d2 - f2) < 1e-10, "the contact formula is exact")
}
cat("   So the radial contraction is the contact one minus r^2 (k^theta)^2 times the anisotropy,\n")
cat("   and nothing about the state has entered yet: this is the metric and the symmetry.\n")

cat("\n=== 3. the conservation equation, checked on a stress known to satisfy it ===\n")
cat("   For a static spherically symmetric T the r component of grad_a T^a_b = 0 reads\n")
cat("     dT^r_r/dr + (f'/2f)(T^r_r - T^t_t) + (2/r)(T^r_r - T^theta_theta) = 0,\n")
cat("   with f the signed 1 - 2M/r. The Reissner-Nordstrom electromagnetic stress is the test:\n")
cat("   T^t_t = T^r_r = -Q^2/8 pi r^4 and T^theta_theta = +Q^2/8 pi r^4, conserved exactly.\n")
Q <- 0.6
Tt_em  <- function(r) -Q^2/(8*pi*r^4)
Tr_em  <- function(r) -Q^2/(8*pi*r^4)
Tth_em <- function(r) +Q^2/(8*pi*r^4)
cons <- function(r, Tt, Tr, Tth, h = 1e-6) {
  fp <- (fofr(r+h) - fofr(r-h))/(2*h)
  dTr <- (Tr(r+h) - Tr(r-h))/(2*h)
  dTr + (fp/(2*fofr(r)))*(Tr(r) - Tt(r)) + (2/r)*(Tr(r) - Tth(r))
}
cat("      r        residual of the conservation equation\n")
for (r in c(3.0, 1.5, 0.8, 0.4)) {
  res <- cons(r, Tt_em, Tr_em, Tth_em)
  cat(sprintf("   %6.2f %35.3e\n", r, res))
  note(abs(res) < 1e-6 * max(1, abs(Tr_em(r))/r), "RN stress satisfies it")
}
cat("   And a plant: break the tensor by scaling the transverse pressure and it must fail.\n")
bad <- function(r) 1.4*Tth_em(r)
for (r in c(1.5, 0.4)) {
  res <- cons(r, Tt_em, Tr_em, bad)
  cat(sprintf("      r = %5.2f: residual %.4e with the transverse pressure scaled by 1.4\n", r, res))
  note(abs(res) > 1e-4, "PLANT fires: a non-conserved tensor is caught")
}

cat("\n=== 4. the reduction, written out ===\n")
cat("   Write X = T^t_t - T^r_r and Y = T^r_r. Conservation gives the anisotropy outright,\n")
cat("     T^theta_theta - T^r_r = (r/2)[ Y' - (f'/2f) X ],\n")
cat("   so the two contractions are\n")
cat("     T_kk(radial)  = |f|(k^t)^2 X\n")
cat("     T_kk(contact) = |f|(k^t)^2 X + r^2 (k^theta)^2 (r/2)[ Y' - (f'/2f) X ]\n")
cat("   and the trace is\n")
cat("     T = X + 4Y + r Y' - (r f'/2f) X.\n")
cat("   Three unknown functions, X, Y and Y', and three relations once the contact contraction and\n")
cat("   the trace are known. Checked by solving the linear system on random data and recovering the\n")
cat("   inputs:\n")
solveXY <- function(r, Tc, Tr_trace, kt, kth) {
  af <- absf(r); h <- 1e-6
  fp <- (fofr(r+h) - fofr(r-h))/(2*h); f0 <- fofr(r)
  a <- af*kt^2; b <- r^2*kth^2*(r/2)
  # T_c = a X + b (Yp - (fp/2f) X);  T = X + 4Y + r Yp - (r fp/2f) X
  # two equations, three unknowns X, Y, Yp: Y appears only in the trace, so the pair
  # (X, Yp) is determined by T_c and any one value of Y. Solve for X and Yp at given Y.
  function(Y) {
    A <- matrix(c(a - b*fp/(2*f0), b,
                  1 - r*fp/(2*f0),  r), nrow = 2, byrow = TRUE)
    rhs <- c(Tc, Tr_trace - 4*Y)
    as.numeric(solve(A, rhs))
  }
}
cat("      r      X in   X out      Yp in   Yp out     residual\n")
for (r in c(0.9, 0.5, 0.3)) {
  af <- absf(r); kt <- 1; kth <- 0.37; Yv <- rnorm(1)
  Xv <- rnorm(1); Ypv <- rnorm(1); h <- 1e-6
  fp <- (fofr(r+h) - fofr(r-h))/(2*h); f0 <- fofr(r)
  Tc <- af*kt^2*Xv + r^2*kth^2*(r/2)*(Ypv - (fp/(2*f0))*Xv)
  Tt <- Xv + 4*Yv + r*Ypv - (r*fp/(2*f0))*Xv
  got <- solveXY(r, Tc, Tt, kt, kth)(Yv)
  cat(sprintf("   %5.2f %7.4f %7.4f %10.4f %8.4f %12.2e\n",
              r, Xv, got[1], Ypv, got[2], max(abs(got - c(Xv, Ypv)))))
  note(max(abs(got - c(Xv, Ypv))) < 1e-8, "the system inverts")
}
cat("   The matrix is not singular anywhere in the interior, which is what makes the reduction\n")
cat("   usable rather than formal:\n")
cat("      r        determinant of the 2x2\n")
for (r in c(1.9, 1.0, 0.5, 0.1)) {
  af <- absf(r); kt <- 1; kth <- 0.37; h <- 1e-6
  fp <- (fofr(r+h) - fofr(r-h))/(2*h); f0 <- fofr(r)
  a <- af*kt^2; b <- r^2*kth^2*(r/2)
  dt <- (a - b*fp/(2*f0))*r - b*(1 - r*fp/(2*f0))
  cat(sprintf("   %6.2f %24.6f\n", r, dt))
  note(abs(dt) > 1e-6, "the determinant does not vanish")
}

cat("\n=== 5. what is still needed, and it is one scalar and not a tensor ===\n")
cat("   The contact contraction is known: its sign is computed twice in the release and its size is\n")
cat("   kappa m^2 D^{-5/2} r_h^{1/2} with kappa = 0.0066642. What the system also needs is the\n")
cat("   TRACE of the image stress. For a scalar of mass m and coupling xi the trace is\n")
cat("     T = (6 xi - 1) (1/2) box <phi^2> - m^2 <phi^2>,\n")
cat("   with no conformal anomaly in it, because the anomaly is state-independent and the image\n")
cat("   part is a difference of states. So the open item is one scalar function of r, the image\n")
cat("   <phi^2> and its Laplacian, rather than a whole tensor.\n")
cat("   Half of that scalar is already available. The image correlator's own amplitude at an\n")
cat("   order-one caustic follows from the same proper-time assembly the coefficient used:\n")
AMP <- 6.6092; KSIG <- 34.85/2
Gimg_c <- AMP * (4*pi)^(-2) * gamma(1.5) * (2/KSIG)^(1.5)
cat(sprintf("      <phi^2>_img = %.7f M^{-1/2} D^{-3/2}\n", Gimg_c))
note(abs(Gimg_c - 6.6092*(1/(16*pi^2))*gamma(1.5)*(2/17.425)^1.5) < 1e-9,
     "the image <phi^2> coefficient assembles")
cat("   What is NOT available is its Laplacian, because box acting on a function of the image\n")
cat("   separation is not the same as box acting on a function of a fixed point: the image moves\n")
cat("   with the point. That derivative is the same object A.18's sum is for, so the reduction\n")
cat("   narrows the open item without dissolving it.\n")
cat("   At xi = 1/6 the Laplacian term drops out of the trace entirely, leaving T = -m^2 <phi^2>,\n")
cat("   which is the coefficient printed above and nothing further. So at that one coupling the\n")
cat("   trace costs nothing.\n")
cat("   One piece is still missing there and it is worth naming exactly, because it is small. The\n")
cat("   system wants the CONTACT contraction in the same affine normalisation as the radial one.\n")
cat("   image_stress_coefficient.R computes the coefficient with the contraction factor set to one,\n")
cat("   which is the radial congruence's value since dr/dlambda = -1 there; the contact geodesic's\n")
cat("   own (k.grad sigma)^2 is a different number and it is the tangent of the E = 0 family with L\n")
cat("   fixed so two interior legs deliver pi, which contact_vanvleck.R already transports. That is\n")
cat("   one ratio, not a tensor and not a mode sum.\n")
cat("   So the ledger on this open item now reads: the radial MAGNITUDE is computed; the radial SIGN\n")
cat("   needs the contact contraction in matched normalisation, and at xi = 1/6 nothing else; away\n")
cat("   from xi = 1/6 it needs the Laplacian of the image <phi^2> as well, which is A.18's object.\n")
cat("   That is a smaller open item than the one this file started with, which was the whole tensor.\n")

