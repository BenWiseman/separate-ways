# A.10 currently says we do not establish that the state's restriction to the Kerr
# bifurcation surface is P_perp-invariant. That may be too pessimistic: P_perp factors
# into two symmetries Kerr has.
#
#   P_perp : (theta, phi) -> (pi - theta, phi + pi)
#          = R_phi(pi)  o  S_theta,
# with R_phi(pi) a rotation by pi about the axis and S_theta the equatorial reflection.
# Kerr is axisymmetric and equatorially symmetric, so a stationary axisymmetric state
# should be invariant under both, hence under their composition.

a <- 0.7; rp <- 1 + sqrt(1-a^2)
rho2 <- function(th) rp^2 + a^2*cos(th)^2
gth  <- function(th) rho2(th)
gph  <- function(th) (rp^2+a^2)^2*sin(th)^2/rho2(th)

cat("=== 1. the induced metric is invariant under EACH factor separately\n\n")
cat("        theta      g_thth(th)   g_thth(pi-th)   g_phph(th)   g_phph(pi-th)\n")
for (th in c(0.3, 0.9, 1.4, 2.1)) 
  cat(sprintf("   %9.3f %13.6f %15.6f %12.6f %14.6f\n", th, gth(th), gth(pi-th), gph(th), gph(pi-th)))
cat(sprintf("\n  max |g(th) - g(pi-th)| over 2000 samples: theta-theta %.2e, phi-phi %.2e\n",
  max(abs(sapply(seq(0.01,pi-0.01,length.out=2000), function(t) gth(t)-gth(pi-t)))),
  max(abs(sapply(seq(0.01,pi-0.01,length.out=2000), function(t) gph(t)-gph(pi-t))))))
cat("  and the metric has no phi dependence at all, so R_phi(pi) is trivially an\n")
cat("  isometry. P_perp is the composition of two isometries and is therefore one.\n")

cat("\n=== 2. so a state invariant under BOTH is P_perp-invariant. Build one and check.\n\n")
cat("  Take any kernel of the invariants an axisymmetric equatorially symmetric state\n")
cat("  can depend on: K(x,y) = F(cos th_x cos th_y, sin th_x sin th_y cos(dphi),\n")
cat("  cos^2 th_x + cos^2 th_y). Test K(P x, P y) = K(x,y).\n\n")
set.seed(5)
Ks <- list(
  "F1" = function(tx,px,ty,py) exp(-(cos(tx)*cos(ty))^2) + 0.4*sin(tx)*sin(ty)*cos(px-py),
  "F2" = function(tx,px,ty,py) 1/(2 + cos(tx)*cos(ty) + sin(tx)*sin(ty)*cos(px-py)),
  "F3" = function(tx,px,ty,py) (cos(tx)^2+cos(ty)^2)*log(3 + sin(tx)*sin(ty)*cos(px-py))
)
P <- function(t,p) c(pi-t, (p+pi) %% (2*pi))
for (nm in names(Ks)) {
  K <- Ks[[nm]]
  e <- max(replicate(20000, {
    tx<-runif(1,0,pi); px<-runif(1,0,2*pi); ty<-runif(1,0,pi); py<-runif(1,0,2*pi)
    Px<-P(tx,px); Py<-P(ty,py)
    abs(K(Px[1],Px[2],Py[1],Py[2]) - K(tx,px,ty,py)) }))
  cat(sprintf("   %-4s max |K(Px,Py) - K(x,y)| over 20000 random pairs = %.3e\n", nm, e))
}
cat("\n  Invariant. So P_perp-invariance of the restricted state follows from Kerr's\n")
cat("  axisymmetry and equatorial symmetry, and does not need to be assumed.\n")

cat("\n=== 3. BUT: invariant does not mean distance-only. That is the part that fails.\n\n")
cat("  On a squashed sphere there is no two-point homogeneity, so an invariant kernel\n")
cat("  is not a function of geodesic distance. Demonstrate: find two pairs at the SAME\n")
cat("  geodesic distance with DIFFERENT kernel values.\n\n")
dist <- function(p1,p2,n=300) {                      # same relaxation as involution_theorem.R
  dphi <- ((p2[2]-p1[2]) + pi) %% (2*pi) - pi
  s <- seq(0,1,length.out=n); th <- p1[1]+s*(p2[1]-p1[1]); ph <- p1[2]+s*dphi
  L <- function(th,ph) sum(sqrt(gth((th[-1]+th[-n])/2)*diff(th)^2 + gph((th[-1]+th[-n])/2)*diff(ph)^2))
  best <- L(th,ph)
  for (it in 1:400) { i<-2:(n-1); tn<-th; tn[i]<-(th[i-1]+th[i+1])/2; pn<-ph; pn[i]<-(ph[i-1]+ph[i+1])/2
    tr<-0.5*(th+tn); pr<-0.5*(ph+pn); v<-L(tr,pr); if (v<best){best<-v; th<-tr; ph<-pr} else break }
  best }
K <- Ks[["F2"]]
# pair A: both near the equator, separated in phi. pair B: separated in theta near the pole.
A1 <- c(pi/2, 0); A2 <- c(pi/2, 0.9)
dA <- dist(A1,A2)
f <- function(dth) dist(c(pi/2 - dth, 0.2), c(pi/2 + dth, 0.2)) - dA
r <- uniroot(f, c(0.05, 1.4), tol=1e-9)$root
B1 <- c(pi/2 - r, 0.2); B2 <- c(pi/2 + r, 0.2)
cat(sprintf("   pair A: equatorial, phi-separated   geodesic length %.9f  K = %.9f\n", dA, K(A1[1],A1[2],A2[1],A2[2])))
cat(sprintf("   pair B: theta-separated             geodesic length %.9f  K = %.9f\n", dist(B1,B2), K(B1[1],B1[2],B2[1],B2[2])))
cat(sprintf("\n   same distance to %.1e, kernel values differ by %.4f\n",
    abs(dist(B1,B2)-dA), abs(K(A1[1],A1[2],A2[1],A2[2]) - K(B1[1],B1[2],B2[1],B2[2]))))

cat("\n=== 4. so the correct statement, sharper than 'not established'\n\n")
cat("  The R = 1 locus EXISTS and is well defined on Kerr: P_perp-invariance follows\n")
cat("  from axisymmetry and equatorial symmetry, so W_B(x, P y) = W_B(P x, y) and the\n")
cat("  locus is the set where W_B(P x, y) = W_B(x, y). What fails is only its\n")
cat("  identification with the GEODESIC-equidistant curve, because an invariant kernel\n")
cat("  on a squashed sphere is not a function of distance. The six-profile check used\n")
cat("  distance-only profiles and therefore tested the geometry, as stated; the locus\n")
cat("  itself does not need them.\n")
