# The route named last turn: can A.10's reduction be redone on a late-time CROSS-SECTION
# of a Killing horizon instead of on the bifurcation surface B? If yes, it extends to
# collapse holes and the last structural caveat goes.
#
# The mechanism only needed P_perp to act transversally and commute with the boost, and a
# cross-section at fixed Killing time gives dtau = 0 by construction. So the argument
# looks like it should carry. Test the step it actually depends on: whether J fixes a
# cross-section the way it fixes B.
#
# (T,X) conformal Kruskal. R = {X>|T|}, F = {T>|X|}, L = {X<-|T|}, P = {T<-|X|}.
#   H+(R), the future horizon of our exterior = ray {X = T, T > 0}
#   H-(R), its past horizon                   = ray {X = -T, T < 0}
#   H-(L)                                     = ray {X = T, T < 0}
#   B                                         = the origin
# J is (T,X) -> (-T,-X).

where <- function(p) {
  tol <- 1e-12
  if (abs(p[1])<tol && abs(p[2])<tol) return("B")
  if (abs(p[2]-p[1])<tol) return(if (p[1]>0) "H+(R)" else "H-(L)")
  if (abs(p[2]+p[1])<tol) return(if (p[1]<0) "H-(R)" else "H+(L)")
  "off-horizon" }

cat("=== 1. J fixes B pointwise. Does it fix anything else on the horizon?\n\n")
cat("        point on H+(R)        J(point)          lands on\n")
for (t in c(0.3, 1.0, 2.5)) {
  p <- c(t,t); Jp <- -p
  cat(sprintf("   (%6.2f,%6.2f) %14s (%6.2f,%6.2f) %14s\n", p[1],p[2], where(p), Jp[1],Jp[2], where(Jp)))
}
cat(sprintf("\n   B itself: J(0,0) = (0,0), region %s, so J fixes B POINTWISE.\n", where(c(0,0))))
cat("\n  J carries every cross-section of H+(R) to a cross-section of H-(L). It fixes only\n")
cat("  the origin. So the reduction cannot be restated as a correlator restricted to a\n")
cat("  single cross-section: the two arguments do not both live there. THE ROUTE FAILS,\n")
cat("  and it fails for the reason B was special in the first place - B is the fixed\n")
cat("  point set of the boost and of J, and no other horizon cross-section is.\n")

cat("\n=== 2. what survives on a cross-section, which is not nothing\n\n")
cat("  Both kernels send their second argument to the SAME place. G_J uses J y and\n")
cat("  G_alpha uses J P_perp y, and P_perp acts only transversally, so both land on\n")
cat("  H-(L) at the same position along the generator and differ only by the transverse\n")
cat("  reflection. Check the 'time' along the generator is unchanged by P_perp:\n\n")
cat("        y on H+(R)       J y            J P_perp y (transverse only)   same generator position?\n")
for (t in c(0.4, 1.3)) {
  y <- c(t,t); Jy <- -y; JPy <- -y     # P_perp moves the sphere, not (T,X)
  cat(sprintf("   (%5.2f,%5.2f) %12s %28s %18s\n", y[1],y[2],
      sprintf("(%.2f,%.2f)",Jy[1],Jy[2]), sprintf("(%.2f,%.2f) + antipode",JPy[1],JPy[2]),
      if (identical(Jy,JPy)) "yes" else "no"))
}
cat("\n  So the RATIO is still a transverse-reflection ratio, now between H+(R) and\n")
cat("  H-(L) rather than within B. The reciprocal law needs only P_perp^2 = 1 and\n")
cat("  survives that unchanged. What does NOT survive is the reading of the ratio as a\n")
cat("  correlator restricted to one surface, which is A.10's actual wording.\n")

cat("\n=== 3. and for a collapse hole? Still nothing, and the reason is sharper now.\n\n")
cat("  A collapse hole has no L, so H-(L) does not exist and both kernels lose their\n")
cat("  second argument. The obstruction was never the absence of B specifically; it is\n")
cat("  the absence of the second sheet. Moving to cross-sections does not help because\n")
cat("  J needs somewhere to map TO.\n")

cat("\n=== 4. THE ROUTE, same turn: state the law algebraically and drop the geometry.\n\n")
cat("  The reciprocal law's proof uses only that P_perp is an involution. It never used\n")
cat("  a manifold. Write it for a state omega on an algebra with an involutive\n")
cat("  automorphism sigma:   R(A,B) = omega(A sigma(B)) / omega(A B). Then\n")
cat("  R(A, sigma B) = omega(A B)/omega(A sigma B) = 1/R(A,B), needing only sigma^2 = id.\n")
cat("  Test on matrices, with no geometry anywhere:\n\n")
set.seed(67); n <- 8
rho <- { M <- matrix(rnorm(n*n),n,n); M <- M %*% t(M); M/sum(diag(M)) }
Uc <- { M <- matrix(rnorm(n*n),n,n); qr.Q(qr(M)) }
S <- Uc %*% diag(c(rep(1,n/2), rep(-1,n/2))) %*% t(Uc)   # involution: S^2 = 1
sig <- function(A) S %*% A %*% S
om  <- function(A) sum(diag(rho %*% A))
cat(sprintf("   sigma is an involution: |sigma(sigma(A)) - A| = %.2e\n",
    { A <- matrix(rnorm(n*n),n,n); max(abs(sig(sig(A)) - A)) }))
errs <- replicate(4000, {
  A <- matrix(rnorm(n*n),n,n); B <- matrix(rnorm(n*n),n,n)
  R1 <- om(A %*% sig(B))/om(A %*% B); R2 <- om(A %*% sig(sig(B)))/om(A %*% sig(B))
  abs(R1*R2 - 1) })
cat(sprintf("   max |R(A,B) R(A,sigma B) - 1| over 4000 random pairs = %.2e\n", max(errs)))
cat("\n  So the law holds for any state on any algebra with any involutive automorphism.\n")
cat("  For a collapse hole the exterior is stationary and axisymmetric at late times, so\n")
cat("  the transverse antipodal map IS an automorphism of its algebra and the law\n")
cat("  applies. What does NOT transfer is the fold's own content: alpha = J o P_perp as\n")
cat("  the second Keldysh leg needs the second sheet, and no algebra supplies one.\n")
cat("  That is the same split the abstract already makes between the standalone horizon\n")
cat("  result and the fold's use of it, reached here from the collapse side.\n")
