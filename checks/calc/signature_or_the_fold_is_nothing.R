# Does the fold have content in any signature, or only in one?
#
# Lorentzian signature is the last non-numeric entry on the assumed list, and the ledger's reason
# for it is thin: "the fold is built on a spacetime that already has one". That is true and it is
# not an argument. The question worth asking is narrower and answerable: in which signatures does
# the fold have any content at all?
#
# WHAT THE FOLD NEEDS. Everything the construction gets from Theta comes from its being a DISCRETE
# symmetry with a parity. Section 2.1's classical-quantum split is the statement that the sheet
# average is the fold-even part and the difference the fold-odd part. A parity is only invariant
# if the map is not continuously deformable to the identity: if it is, the "odd" part can be
# rotated into the "even" part along the deformation and the split means nothing.
#
# So the test is whether dTheta = -Id sits in the identity component of the isometry group. In
# four Euclidean dimensions -Id has determinant +1 and SO(4) is connected, so it does, and the
# fold is a rotation. In four Lorentzian dimensions -Id is PT, which has determinant +1 but
# reverses time orientation, and every element of O(1,3) obeys |Lambda^0_0| >= 1, so no continuous
# path of isometries can carry +1 to -1. The fold is then a genuine Z_2 with a parity to give.
# Both halves are measured below.

set.seed(31415)
eta <- diag(c(-1,1,1,1)); del <- diag(4)

cat("=== 1. Euclidean: an explicit path of isometries from the identity to -Id ===\n")
cat("   Rotate by theta in the 1-2 plane and by theta in the 3-4 plane at once. At theta = pi\n")
cat("   that is -Id, and every point on the way is an isometry of delta.\n\n")
Rpath <- function(th) {
  M <- diag(4)
  M[1,1] <- cos(th); M[1,2] <- -sin(th); M[2,1] <- sin(th); M[2,2] <- cos(th)
  M[3,3] <- cos(th); M[3,4] <- -sin(th); M[4,3] <- sin(th); M[4,4] <- cos(th)
  M
}
cat("      theta      max |R^T delta R - delta|      det R      R equals -Id?\n")
worst <- 0
for (th in c(0, pi/4, pi/2, 3*pi/4, pi)) {
  M <- Rpath(th); v <- max(abs(t(M) %*% del %*% M - del)); worst <- max(worst, v)
  cat(sprintf("   %9.4f   %24.2e   %8.4f   %s\n", th, v, det(M),
              ifelse(max(abs(M + diag(4))) < 1e-12, "yes", "no")))
}
stopifnot(worst < 1e-12, max(abs(Rpath(pi) + diag(4))) < 1e-12)
cat("   So in Euclidean signature -Id is connected to the identity through isometries. A map you\n")
cat("   can rotate away has no invariant parity, and the fold-even and fold-odd parts are not\n")
cat("   separate things.\n")

cat("\n=== 2. Lorentzian: the same path is not available, and the obstruction is a number ===\n")
cat("   Every Lambda in O(1,3) satisfies |Lambda^0_0| >= 1. Sampled over random boosts and\n")
cat("   rotations, including large rapidities:\n\n")
boost <- function(r, n) {
  n <- n/sqrt(sum(n^2)); L <- diag(4); g <- cosh(r); s <- sinh(r)
  L[1,1] <- g
  for (i in 1:3) { L[1,i+1] <- -s*n[i]; L[i+1,1] <- -s*n[i] }
  for (i in 1:3) for (j in 1:3) L[i+1,j+1] <- (i==j) + (g-1)*n[i]*n[j]
  L
}
rot <- function() { Q <- qr.Q(qr(matrix(rnorm(9),3,3))); if (det(Q)<0) Q[,1] <- -Q[,1]
  M <- diag(4); M[2:4,2:4] <- Q; M }
mn <- Inf; mx <- -Inf; iso <- 0
for (i in 1:20000) {
  L <- rot() %*% boost(rnorm(1, 0, 2), rnorm(3)) %*% rot()
  # relative tolerance: a boost of rapidity 9 has entries of order 5e3, and the round-off in
  # L^T eta L scales with their square, so an absolute bar would fail on size and not on physics
  if (max(abs(t(L) %*% eta %*% L - eta)) / max(1, max(abs(L))^2) < 1e-12) iso <- iso + 1
  mn <- min(mn, abs(L[1,1])); mx <- max(mx, abs(L[1,1]))
}
cat(sprintf("   %d of 20000 sampled maps are exact isometries of eta\n", iso))
cat(sprintf("   smallest |Lambda^0_0| found: %.10f     largest: %.2f\n", mn, mx))
stopifnot(iso == 20000, mn >= 1 - 1e-9)
cat("   The entry never enters (-1, 1). It is continuous along any path of isometries, so it\n")
cat("   cannot travel from +1 at the identity to -1 at -Id: the two lie in different components.\n")
cat("   In Lorentzian signature the fold is PT and is a genuine discrete symmetry.\n")

cat("\n=== 3. the same statistic in Euclidean signature, where it is free ===\n")
mnE <- Inf
for (i in 1:20000) {
  Q <- qr.Q(qr(matrix(rnorm(16),4,4))); if (det(Q) < 0) Q[,1] <- -Q[,1]
  mnE <- min(mnE, abs(Q[1,1]))
}
cat(sprintf("   smallest |Q_11| over 20000 random elements of SO(4): %.6f\n", mnE))
cat("   It passes freely through zero, which is the same fact as section 1: nothing separates\n")
cat("   the identity from -Id there.\n")
stopifnot(mnE < 0.05)

cat("\n=== 4. and it is the signature doing it, not the dimension ===\n")
cat("   -Id has determinant (-1)^D, so in odd dimensions it is orientation-reversing and lies\n")
cat("   outside SO(D) for a reason that has nothing to do with time. The construction is at\n")
cat("   D = 4 by the causal budget, so the relevant comparison is even-dimensional:\n\n")
cat("      D    det(-Id)   Euclidean: -Id in the identity component?   Lorentzian: ?\n")
for (D in c(2,4,6)) {
  cat(sprintf("   %4d   %8d   %-42s %s\n", D, (-1)^D, "yes, SO(D) is connected and det = +1",
              "no, it reverses time orientation"))
}
cat("   At every even D the Euclidean answer is yes and the Lorentzian answer is no, so the\n")
cat("   separation is a fact about the signature.\n")

cat("\n=== 5. the plant: a map that IS connected to the identity must fail the same test ===\n")
cat("   Take a Lorentzian element that is a pure spatial rotation by pi in the 2-3 plane. It is\n")
cat("   an involution, it is an isometry of eta, and it has Lambda^0_0 = +1, so it should sit in\n")
cat("   the identity component and carry no invariant parity.\n\n")
Rsp <- diag(c(1,1,-1,-1))
cat(sprintf("   isometry of eta: %.1e   involution: %.1e   Lambda^0_0 = %+.1f\n",
            max(abs(t(Rsp) %*% eta %*% Rsp - eta)), max(abs(Rsp %*% Rsp - diag(4))), Rsp[1,1]))
path_ok <- TRUE
for (th in seq(0, pi, length.out = 40)) {
  M <- diag(4); M[3,3] <- cos(th); M[3,4] <- -sin(th); M[4,3] <- sin(th); M[4,4] <- cos(th)
  if (max(abs(t(M) %*% eta %*% M - eta)) > 1e-12) path_ok <- FALSE
}
cat(sprintf("   and a continuous path of eta-isometries reaches it: %s\n", path_ok))
stopifnot(path_ok)
cat("   So the test separates the two cases rather than rejecting everything: a spatial rotation\n")
cat("   by pi is reachable and -Id is not, in the same signature and by the same criterion.\n")

cat("\n=== 6. what this is and what it is not ===\n")
cat("   It is not a derivation of the signature from nothing. Nothing here says a Lorentzian\n")
cat("   manifold must exist. What it says is that the fold has content in one signature only:\n")
cat("   in Euclidean signature -Id is a rotation, continuously removable, with no invariant\n")
cat("   parity, so there is no classical-quantum split, no fold-even and fold-odd sectors, and\n")
cat("   nothing for the rest of the construction to be built on. The assumption earns a better\n")
cat("   statement than the ledger's: the signature is not a convenience the fold was written in,\n")
cat("   it is the condition under which the fold is a symmetry at all.\n")
