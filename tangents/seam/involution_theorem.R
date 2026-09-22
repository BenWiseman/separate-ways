# The reciprocal identity needs only P_perp^2 = 1. If so it survives Kerr squashing,
# unequal radii, and nonzero time separation, and the ninety-degree equality is its
# fixed point rather than a separate fact.
#
#   R(x,y) := F(d(x, P y)) / F(d(x, y))
#   R(x, P y) = F(d(x, PPy)) / F(d(x, P y)) = F(d(x,y)) / F(d(x,Py)) = 1 / R(x,y)
#
# One line, no geometry, no state, no correlator. Test it where the round-sphere
# argument has no purchase.

set.seed(11)
Fs <- list("conformal"=function(d) 1/d^2, "heavy"=function(d) d^(-3.4),
           "screened"=function(d) exp(-0.7*d)/d^2, "log"=function(d) log(9/d),
           "bounded"=function(d) atan(3/d), "gaussian"=function(d) exp(-d^2/4)+0.1)

cat("=== 1. general de Sitter: unequal radii AND nonzero time separation\n")
cat("   max |R(x,y) R(x,Py) - 1| over 4000 random configurations\n\n")
Z <- function(r1,r2,dt,cg) -sqrt((1-r1^2)*(1-r2^2))*cosh(dt) + r1*r2*cg
for (nm in names(Fs)) {
  G <- function(z) Fs[[nm]](sqrt(pmax(2*(1-z),1e-12)))     # F as a function of chordal sep
  w <- replicate(4000, { r1<-runif(1,0,.97); r2<-runif(1,0,.97); dt<-runif(1,-2,2); cg<-runif(1,-1,1)
    Rxy <- G(Z(r1,r2,dt,-cg))/G(Z(r1,r2,dt,cg))            # P flips cos gamma
    Rxp <- G(Z(r1,r2,dt, cg))/G(Z(r1,r2,dt,-cg))
    Rxy*Rxp - 1 })
  cat(sprintf("  %-10s %.3e\n", nm, max(abs(w))))
}

cat("\n=== 2. KERR: the SQUASHED bifurcation surface, where 'pi - gamma' means nothing\n")
cat("   Embed B with the paper's induced metric, P = the free involution\n")
cat("   (theta -> pi-theta, phi -> phi+pi), measure true geodesic distance on B.\n\n")
# Induced metric on the Kerr bifurcation surface: ds^2 = rho^2 dth^2 + ((r+^2+a^2)^2 sin^2 th / rho^2) dphi^2
# with rho^2 = r+^2 + a^2 cos^2 th. Distances by discretised geodesic (relaxation on a grid).
kerr_dist <- function(a, rp, p1, p2, n=400) {
  # straight-line-in-(theta,phi) path then relax; adequate for an identity check
  rho2 <- function(th) rp^2 + a^2*cos(th)^2
  gth  <- function(th) rho2(th)
  gph  <- function(th) (rp^2+a^2)^2*sin(th)^2/rho2(th)
  dphi <- ((p2[2]-p1[2]) + pi) %% (2*pi) - pi
  s <- seq(0,1,length.out=n)
  th <- p1[1] + s*(p2[1]-p1[1]); ph <- p1[2] + s*dphi
  L <- function(th,ph) sum(sqrt(gth((th[-1]+th[-n])/2)*diff(th)^2 + gph((th[-1]+th[-n])/2)*diff(ph)^2))
  best <- L(th,ph)
  for (it in 1:600) {                                  # relax interior nodes
    i <- 2:(n-1)
    th_new <- th; th_new[i] <- (th[i-1]+th[i+1])/2
    ph_new <- ph; ph_new[i] <- (ph[i-1]+ph[i+1])/2
    tr <- 0.5*(th+th_new); pr <- 0.5*(ph+ph_new)
    v <- L(tr,pr); if (v < best) { best <- v; th <- tr; ph <- pr } else break
  }
  best
}
Pmap <- function(p) c(pi - p[1], (p[2] + pi) %% (2*pi))
for (a in c(0.3, 0.7, 0.95)) {
  rp <- 1 + sqrt(max(1-a^2,1e-9))
  cat(sprintf("  a = %.2f  (squashing (r+^2+a^2)/r+^2 = %.4f)\n", a, (rp^2+a^2)/rp^2))
  for (nm in c("conformal","screened","bounded")) {
    G <- Fs[[nm]]
    err <- replicate(12, {
      x <- c(runif(1,0.15,pi-0.15), runif(1,0,2*pi)); y <- c(runif(1,0.15,pi-0.15), runif(1,0,2*pi))
      Py <- Pmap(y)
      Rxy <- G(kerr_dist(a,rp,x,Py))/G(kerr_dist(a,rp,x,y))
      Rxp <- G(kerr_dist(a,rp,x,Pmap(Py)))/G(kerr_dist(a,rp,x,Py))
      Rxy*Rxp - 1 })
    cat(sprintf("     %-10s max |R R' - 1| = %.3e\n", nm, max(abs(err))))
  }
}
cat("\n  Holds on the squashed surface too. The identity never used roundness: it used\n")
cat("  P_perp P_perp = 1, which A.10 already verifies on Kerr to 1.8e-15.\n")

cat("\n=== 3. so where is the ratio exactly 1? The EQUIDISTANT LOCUS of P_perp.\n")
cat("   On Kerr, solve d(x,Py) = d(x,y) for x at fixed y. Check R = 1 there.\n\n")
a <- 0.7; rp <- 1 + sqrt(1-a^2); y <- c(1.1, 0.4); Py <- Pmap(y)
f <- function(th) kerr_dist(a,rp,c(th,2.0),Py) - kerr_dist(a,rp,c(th,2.0),y)
th0 <- uniroot(f, c(0.2, pi-0.2), tol=1e-10)$root
cat(sprintf("  equidistant point found at theta = %.8f (phi = 2.0)\n", th0))
for (nm in names(Fs)) {
  G <- Fs[[nm]]
  cat(sprintf("     %-10s R = %.12f\n", nm, G(kerr_dist(a,rp,c(th0,2),Py))/G(kerr_dist(a,rp,c(th0,2),y))))
}
cat("\n  R = 1 for every profile at the same point, because the two distances are equal\n")
cat("  there and the profile cancels. On a round sphere that locus is the great circle\n")
cat("  at ninety degrees; on Kerr it is a deformed curve. Same theorem either way.\n")
