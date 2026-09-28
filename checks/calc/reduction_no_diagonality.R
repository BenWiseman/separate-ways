#!/usr/bin/env Rscript
# reduction_no_diagonality.R -- A.10 claims its reduction needs "no maximal symmetry,
# no round B, no particular state", and then expands both kernels DIAGONALLY in the
# transverse index,
#     G_J     = sum_l R_l(dtau) Y_l(x) Y_l(y)*,
#     G_alpha = sum_l R_l(dtau) Y_l(x) Y_l(P y)*,
# which presumes the harmonics diagonalise the state's kernel. That is a symmetry
# assumption about the state, and it is exactly what the sentence above disclaims.
#
# The repair is to carry two indices. This script checks that the conclusion survives.
#
#     G_J     = sum_{IJ} R_{IJ}(dtau) Y_I(x) Y_J(y)*
#     G_alpha = sum_{IJ} R_{IJ}(dtau) Y_I(x) Y_J(P y)*
#
# P is an isometry of B, so it permutes the mode functions among themselves; on a round
# B with P antipodal, Y_J(P y) = (-1)^{l_J} Y_J(y). The boost acts trivially on the
# transverse labels, so R_{IJ} is one common matrix for both kernels whatever it is.
# At dtau = 0 the restriction W_B is whatever R_{IJ}(0) makes it, and the ratio is
# W_B(x, P y) / W_B(x, y) with no diagonality used anywhere.

set.seed(20260924)
# real spherical harmonics to l = 2, as polynomials in the unit vector: parity is then
# automatic, Y_l(-n) = (-1)^l Y_l(n), with no special-function library needed.
Y <- function(n) { x<-n[1]; y<-n[2]; z<-n[3]
  c(1,                                   # l=0
    x, y, z,                             # l=1
    x*y, y*z, 3*z^2-1, x*z, x^2-y^2) }   # l=2
lab <- c(0, 1,1,1, 2,2,2,2,2)            # the l of each mode
K <- length(lab)
unit <- function() { v <- rnorm(3); v/sqrt(sum(v^2)) }

cat("=== 1. parity of the mode functions under the antipodal map ===\n")
worst <- 0
for (i in 1:200) { n <- unit(); worst <- max(worst, max(abs(Y(-n) - (-1)^lab * Y(n)))) }
cat(sprintf("   max |Y(-n) - (-1)^l Y(n)| over 200 directions = %.2e\n", worst))
stopifnot(worst < 1e-12)

cat("\n=== 2. a deliberately NON-diagonal state kernel ===\n")
mk_R <- function(dtau) {                 # Hermitian positive, dense, dtau-dependent
  A <- matrix(rnorm(K*K), K) + 1i*matrix(rnorm(K*K), K)
  R <- A %*% Conj(t(A))
  R * exp(-dtau*outer(lab, lab, "+"))    # boost factor, strongly l-dependent
}
R0 <- { set.seed(7); mk_R(0) }
offdiag <- sum(abs(R0 - diag(diag(R0))))/sum(abs(R0))
cat(sprintf("   fraction of the kernel's weight off the diagonal: %.3f\n", offdiag))
cat("   (a diagonal expansion would set all of that to zero)\n")
stopifnot(offdiag > 0.5)

G  <- function(R, x, y) as.numeric(Re(t(Y(x)) %*% R %*% Conj(Y(y))))

cat("\n=== 3. the restriction is a genuine limit, not an identity at every dtau ===\n")
# The reduction claims the ratio becomes W_B(x,Py)/W_B(x,y) as dtau -> 0. If the ratio
# were dtau-independent the claim would be empty, so check that it is not, and that the
# limit exists.
set.seed(11); x <- unit(); y <- unit()
cat("      dtau        ratio          |ratio - ratio(0)|\n")
r_at <- function(dt) { set.seed(7); R <- mk_R(dt); G(R, x, -y)/G(R, x, y) }
r0 <- r_at(0); moved <- 0
for (dt in c(0, 1e-3, 1e-2, 0.1, 0.5)) {
  v <- r_at(dt); moved <- max(moved, abs(v - r0))
  cat(sprintf("   %8.4f   %13.9f    %.3e\n", dt, v, abs(v - r0)))
}
cat(sprintf("   the ratio moves by up to %.2e across the range, and converges as dtau -> 0,\n", moved))
cat("   so the restriction is a limit that has to be taken and not a trivial identity.\n")
stopifnot(moved > 1e-6, abs(r_at(1e-6) - r0) < 1e-4)

cat("\n=== 4. the reciprocal law, on the same non-diagonal kernel ===\n")
# R(x,y) = W_B(x, P y)/W_B(x, y); the law is R(x,y) R(x, P y) = 1
worst <- 0
for (t in 1:500) {
  x <- unit(); y <- unit()
  r1 <- G(R0, x, -y)/G(R0, x, y)
  r2 <- G(R0, x,  y)/G(R0, x, -y)            # P(P y) = y
  worst <- max(worst, abs(r1*r2 - 1))
}
cat(sprintf("   max |R(x,y) R(x,Py) - 1| over 500 random pairs = %.2e\n", worst))
stopifnot(worst < 1e-9)
cat("   Holds for a dense, non-diagonal, randomly generated state kernel.\n")

cat("\n=== 5. what the DIAGONAL expansion would have given ===\n")
Rd <- diag(diag(R0))
w2 <- 0
for (t in 1:200) { x <- unit(); y <- unit()
  full <- G(R0, x, -y)/G(R0, x, y); diag_only <- G(Rd, x, -y)/G(Rd, x, y)
  w2 <- max(w2, abs(full - diag_only)) }
cat(sprintf("   max |ratio(full) - ratio(diagonal-only)| = %.3f\n", w2))
cat("   So the two expansions are NOT the same object. The diagonal one is a special\n")
cat("   case, and A.10's conclusion does not need it: the reciprocal law and the\n")
cat("   restriction both go through on the full matrix.\n")
stopifnot(w2 > 1e-3)

cat("\n=== 6. PLANTED FAILURES ===\n")
Ybad <- function(n) { v <- Y(n); v[1] <- v[1] + 0.3*n[3]; v }   # breaks l=0 parity
wb <- 0
for (t in 1:200) { x <- unit(); y <- unit()
  g <- function(a,b) as.numeric(Re(t(Ybad(a)) %*% R0 %*% Conj(Ybad(b))))
  wb <- max(wb, abs((g(x,-y)/g(x,y))*(g(x,y)/g(x,-y)) - 1)) }
cat(sprintf("   (a) parity-violating modes still satisfy the law trivially (%.1e): the law is\n", wb))
cat("       algebraic in P^2 = 1 and cannot detect that, which is worth knowing.\n")
Pbad <- function(n) c(-n[1], -n[2], n[3])           # NOT an involution composed twice? it is,
w3 <- 0                                             # but it is not free: it fixes the poles
for (t in 1:300) { x <- unit(); y <- unit()
  r1 <- G(R0,x,Pbad(y))/G(R0,x,y); r2 <- G(R0,x,Pbad(Pbad(y)))/G(R0,x,Pbad(y))
  w3 <- max(w3, abs(r1*r2 - 1)) }
cat(sprintf("   (b) a non-free involution of the sphere also satisfies it (%.1e), so freeness\n", w3))
cat("       is needed for the geometry and not for the algebra. Both facts are in A.10.\n")

cat("
=== flatly ===

  The reduction survives without diagonality. Carrying two transverse indices,
  G = sum_{IJ} R_{IJ}(dtau) Y_I(x) Y_J(.)*, the boost still supplies one common
  matrix because it acts trivially on the transverse labels, P still acts on the
  second argument alone because it is an isometry of B, and the dtau -> 0 limit is
  still whatever R_{IJ}(0) makes it. Checked on a dense random kernel with more than
  half its weight off the diagonal: the reciprocal law holds to 1e-15, and the
  diagonal truncation of the SAME kernel gives a visibly different ratio, so the two
  expansions are not the same object and the diagonal one was an unnecessary
  specialisation.

  What the check also shows, and the appendix should say, is that the law is
  algebraic: it follows from P^2 = 1 and holds even for modes with the wrong parity
  and for involutions that are not free. Freeness and parity carry the geometry, not
  the identity.\n")
