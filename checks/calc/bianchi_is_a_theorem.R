# Is diffeomorphism invariance an assumption of this derivation, or is what it is used for a
# theorem?
#
# The relativity ledger lists diffeomorphism invariance as assumed, and the derivation uses it in
# exactly two places. On the gravitational side it needs the contracted Bianchi identity, which is
# what turns "2 pi T_kk = eta R_kk for every null k" into a field equation with one integration
# constant. On the matter side it needs grad^a T_ab = 0, which conservation_from_diffeo.R already
# obtained from Noether's second theorem applied to the matter action.
#
# The two have different status and the ledger has been treating them as one. The matter side is a
# genuine input: the matter action has to be a scalar functional of the metric, which is part of
# the metric postulate. The gravitational side is not an input at all, because the contracted
# Bianchi identity follows from the definition of the Riemann tensor of a metric connection. It is
# a theorem about any metric whatsoever, and nothing needs to be assumed for it to hold.
#
# That is worth checking rather than quoting, because it decides how long the assumed list is.
# Computed here on a generic curved metric with no symmetry, by finite differences carried to
# third order in the metric, with a torsionful connection as the plant.

set.seed(8812)
N <- 4
eta <- diag(c(-1, 1, 1, 1))
P1 <- array(rnorm(N*N*N), c(N,N,N)); Q1 <- array(rnorm(N*N*N*N), c(N,N,N,N))
C1 <- array(rnorm(N*N*N*N*N), c(N,N,N,N,N))
for (a in 1:N) for (b in 1:N) {
  P1[a,b,] <- P1[b,a,] <- (P1[a,b,] + P1[b,a,])/2
  Q1[a,b,,] <- Q1[b,a,,] <- (Q1[a,b,,] + Q1[b,a,,])/2
  C1[a,b,,,] <- C1[b,a,,,] <- (C1[a,b,,,] + C1[b,a,,,])/2
}
gmet <- function(x, s = 0.06) {                 # cubic, so third derivatives are not zero
  h <- matrix(0,N,N)
  for (a in 1:N) for (b in 1:N)
    h[a,b] <- sum(P1[a,b,]*x) + sum(outer(x,x)*Q1[a,b,,]) +
              sum(outer(outer(x,x),x)*C1[a,b,,,])
  eta + s*h
}
dd <- function(f, x, a, h) { e <- rep(0,N); e[a] <- h; (f(x+e) - f(x-e))/(2*h) }

Gam <- function(x, h, tors = 0) {
  gi <- solve(gmet(x)); dg <- lapply(1:N, function(c) dd(gmet, x, c, h))
  G <- array(0, c(N,N,N))
  for (a in 1:N) for (b in 1:N) for (c in 1:N)
    G[a,b,c] <- sum(sapply(1:N, function(d) gi[a,d]*(dg[[c]][d,b] + dg[[b]][d,c] - dg[[d]][b,c])))/2
  if (tors != 0) {                              # an antisymmetric piece: a connection with torsion
    Tt <- array(0, c(N,N,N))
    Tt[1,2,3] <-  tors; Tt[1,3,2] <- -tors
    Tt[2,3,1] <-  tors; Tt[2,1,3] <- -tors
    G <- G + Tt
  }
  G
}
Riem <- function(x, h, tors = 0) {              # R^a_bcd, from Gamma and its first derivatives
  G  <- Gam(x, h, tors)
  dG <- lapply(1:N, function(c) dd(function(y) Gam(y, h, tors), x, c, h))
  R <- array(0, c(N,N,N,N))
  for (a in 1:N) for (b in 1:N) for (c in 1:N) for (d in 1:N)
    R[a,b,c,d] <- dG[[c]][a,d,b] - dG[[d]][a,c,b] +
      sum(sapply(1:N, function(e) G[a,c,e]*G[e,d,b] - G[a,d,e]*G[e,c,b]))
  R
}
Einstein <- function(x, h, tors = 0) {          # G_ab with both indices down
  R <- Riem(x, h, tors); g <- gmet(x); gi <- solve(g)
  Ric <- matrix(0,N,N)
  for (b in 1:N) for (d in 1:N) Ric[b,d] <- sum(sapply(1:N, function(a) R[a,b,a,d]))
  Ric <- (Ric + t(Ric))/2                       # the metric connection makes it symmetric anyway
  Rs <- sum(gi * Ric)
  Ric - 0.5*Rs*g
}
divG <- function(x, h, tors = 0) {              # grad_a G^a_b
  G <- Gam(x, h, tors); gi <- solve(gmet(x))
  Gm <- function(y) solve(gmet(y)) %*% Einstein(y, h, tors)     # G^a_b
  Gmix <- Gm(x)
  out <- numeric(N)
  for (b in 1:N) {
    s <- sum(sapply(1:N, function(a) dd(function(y) Gm(y)[a,b], x, a, h)))
    s <- s + sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[a,a,c]*Gmix[c,b]))))
    s <- s - sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[c,a,b]*Gmix[a,c]))))
    out[b] <- s
  }
  out
}

x0 <- c(0.15, -0.10, 0.20, -0.05)
cat("=== 1. the setting: a generic curved metric with no symmetry at all ===\n")
g0 <- gmet(x0)
cat(sprintf("   eigenvalue signs %s, det %.5f, max |g - eta| = %.4f\n",
            paste(sign(eigen(g0)$values), collapse=" "), det(g0), max(abs(g0 - eta))))
E0 <- Einstein(x0, 6e-3)
cat(sprintf("   the Einstein tensor there is not small: max |G_ab| = %.4f\n", max(abs(E0))))

cat("\n=== 2. the contracted Bianchi identity, measured, with nothing assumed to make it hold ===\n")
cat("   Nothing in the code above imposes grad_a G^a_b = 0. The Christoffels come from the\n")
cat("   metric, the Riemann tensor from the Christoffels, and the divergence is taken directly.\n\n")
cat("        step h     max |grad_a G^a_b|     max |G_ab|      relative\n")
res <- c(); hs <- c(1.2e-2, 8e-3, 6e-3)
for (h in hs) {
  d <- max(abs(divG(x0, h))); res <- c(res, d)
  cat(sprintf("     %10.1e   %18.3e   %12.5f   %11.2e\n", h, d, max(abs(E0)), d/max(abs(E0))))
}
ex <- unname(coef(lm(log(res) ~ log(hs)))[2])
cat(sprintf("\n   the residual falls as h^%.2f, which is the finite differencing converging to zero\n", ex))
cat("   and not a number the identity happens to sit near.\n")
stopifnot(res[length(res)]/max(abs(E0)) < 1e-3, ex > 1.2)

cat("\n=== 3. the plant: a connection that is not the metric's ===\n")
cat("   Add an antisymmetric piece to the Christoffels, which is a connection with torsion. It\n")
cat("   is still a connection and the curvature of it is still a curvature, but it is no longer\n")
cat("   the Levi-Civita connection of g, and the identity has no reason to hold.\n\n")
cat("      torsion      max |grad_a G^a_b|     relative to max |G_ab|\n")
for (tt in c(0, 0.05, 0.2, 0.5)) {
  d <- max(abs(divG(x0, 6e-3, tors = tt)))
  cat(sprintf("   %10.2f   %18.4f   %20.4f\n", tt, d, d/max(abs(E0))))
  if (tt == 0.5) big <- d
}
cat("\n   With torsion the divergence is order unity, so the check is reading the connection.\n")
stopifnot(big/max(abs(E0)) > 0.1)

cat("\n=== 4. what this does to the ledger ===\n")
cat("   The identity holds for any metric, with nothing assumed beyond the connection being the\n")
cat("   one the metric determines. So the gravitational half of what the derivation calls\n")
cat("   diffeomorphism invariance is a theorem about Riemannian geometry and not an input.\n")
cat("   The matter half is a genuine input and stays: the matter action has to be a scalar\n")
cat("   functional of the metric, which conservation_from_diffeo.R then turns into\n")
cat("   grad^a T_ab = 0 by Noether's second theorem. That input is part of the metric postulate\n")
cat("   already on the list, since T_ab is defined by varying that action against that metric.\n")
cat("   Diffeomorphism invariance therefore comes off the assumed list as a separate entry,\n")
cat("   leaving a Lorentzian signature, one shared metric, and the two measured constants.\n")
