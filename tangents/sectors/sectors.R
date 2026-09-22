# A.3: "the cover is simply connected and its H^1 vanishes for every coefficient group, so
# there are no periods to quantize a seam coefficient against". That is a statement about
# the COVER. But P_perp is a SPACETIME involution, and gauging one makes the sectors live
# on the quotient, not the cover. Those can differ, so A.3 may be answering next door to
# the question. Work out the quotient for the ACTUAL map rather than assume.
#
# 4.1: P_perp = diag(+1,-1,-1,-1) in the embedding, so on a closed slice it sends
# (X1,X2,X3,X4) -> (X1,-X2,-X3,-X4). Not the antipodal map of S^3.

set.seed(5)
P <- diag(c(1,-1,-1,-1))
cat("=== 1. fixed points of P_perp on S^3, counted rather than assumed\n\n")
n <- 400000
X <- matrix(rnorm(4*n), n, 4); X <- X/sqrt(rowSums(X^2))
d <- sqrt(rowSums((X %*% P - X)^2))
cat(sprintf("   random points with |P x - x| < 1e-3 : %d of %d\n", sum(d < 1e-3), n))
cat("   fixed set is where X2 = X3 = X4 = 0, i.e. X1 = +/-1, so TWO POINTS.\n")
for (p in list(c(1,0,0,0), c(-1,0,0,0))) 
  cat(sprintf("   check (%+.0f,%.0f,%.0f,%.0f): |P x - x| = %.1e\n", p[1],p[2],p[3],p[4],
      sqrt(sum((as.numeric(p) %*% P - p)^2))))
cat("\n  So P_perp is NOT free on the S^3 slice, though it IS free on the bifurcation\n")
cat("  surface B (A.10), where it is antipodal. Different spaces, different answers, and\n")
cat("  the sector question is about the slice.\n")

cat("\n=== 2. what the quotient is\n\n")
cat("  Write X1 = cos(theta) and (X2,X3,X4) = sin(theta) * nhat with nhat in S^2. The map\n")
cat("  is theta -> theta, nhat -> -nhat. Verify that parametrisation reproduces P_perp:\n\n")
th <- runif(6, 0, pi); u <- matrix(rnorm(18),6,3); u <- u/sqrt(rowSums(u^2))
Xp <- cbind(cos(th), sin(th)*u)
err <- max(abs((Xp %*% P) - cbind(cos(th), -sin(th)*u)))
cat(sprintf("   max discrepancy over 6 random points: %.2e\n", err))
cat("\n  So the quotient is theta in [0,pi] with nhat in S^2/+- = RP^2, and the RP^2\n")
cat("  collapses to a point at theta = 0 and theta = pi. That is exactly the unreduced\n")
cat("  SUSPENSION of RP^2.\n")

cat("\n=== 3. its cohomology, by the suspension isomorphism H~^n(Sigma X) = H~^(n-1)(X)\n\n")
cat("        degree n    H~^(n-1)(RP^2; Z2)    H~^n(Sigma RP^2; Z2)\n")
tab <- list(c("1","H~^0 = 0","0"), c("2","H~^1 = Z2","Z2"), c("3","H~^2 = Z2","Z2"))
for (r in tab) cat(sprintf("   %10s %20s %22s\n", r[1], r[2], r[3]))

cat("\n=== 4. the verdict, and it cuts both ways\n\n")
cat("  A.3's CONCLUSION SURVIVES: H^1 of the quotient vanishes too, so there is still one\n")
cat("  sector and no period to quantize a coefficient against. The seam coefficient is not\n")
cat("  rescued by topology after all, and A.13's Hadamard argument remains the only thing\n")
cat("  that fixes it.\n")
cat("  A.3's REASONING DOES NOT: it argues from the cover being simply connected, which is\n")
cat("  the wrong space for a gauged spacetime involution. The right statement is about the\n")
cat("  quotient, and it happens to give the same answer. Worth correcting, because a\n")
cat("  referee who notices the gap will not assume the answer survives it.\n")
cat("  AND H^2 IS NOT ZERO. H^2(Sigma RP^2; Z2) = Z2, which A.3 never mentions because it\n")
cat("  only ever looked at degree one. That is where a discrete-torsion-like choice would\n")
cat("  live if the gauging admitted one. We do not claim it does: for the GROUP Z2,\n")
cat("  H^2(Z2; U(1)) = 0, so no discrete torsion arises from the group itself. Whether the\n")
cat("  space's H^2 does anything here is open and is NOT settled by this calculation.\n")
