# Two referees independently named the same decisive calculation and I ran the
# entropy route past it. Doing it now, including the part that goes against my own
# edit.
#
# PART A. How does the sheet exchange act on the corner term (kappa/2) q^T J qdot?
#   q = (q1, q2), the boundary values on the two sheets. J = [[0,1],[-1,0]].
#   Sheet exchange E = [[0,1],[1,0]]. Time reversal sends qdot -> -qdot.

J <- matrix(c(0,1,-1,0),2,2,byrow=TRUE); E <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
cat("=== A. the exchange versus the complex structure\n\n")
cat("  E J E^-1 =\n"); print(E %*% J %*% solve(E))
cat(sprintf("\n  E J E^-1 + J  =  %s   -> the exchange ANTICOMMUTES with J\n",
    paste(round(as.vector(E%*%J%*%solve(E) + J),12), collapse=" ")))
cat("\n  So under exchange AND time reversal together:\n")
cat("    (k/2) q^T J qdot  ->  (k/2) (Eq)^T J (-E qdot) = -(k/2) q^T (E J E) qdot\n")
cat("                      =  +(k/2) q^T J qdot\n")
set.seed(3); ok <- TRUE
for (i in 1:2000) {
  q <- rnorm(2); qd <- rnorm(2); k <- rnorm(1)
  a <- (k/2) * t(q) %*% J %*% qd
  b <- (k/2) * t(E%*%q) %*% J %*% (-(E%*%qd))
  if (abs(a-b) > 1e-12) ok <- FALSE
}
cat(sprintf("\n  checked on 2000 random (q, qdot, kappa): invariant = %s\n", ok))
cat("  INVARIANT FOR EVERY KAPPA. The fold's own sheet exchange therefore fixes\n")
cat("  nothing, and the rapidity route is dead. This is the sigma = diag(1,-1) case\n")
cat("  both referees flagged, realised not by assumption but by E anticommuting with J.\n")
cat("  My earlier note called this route 'a reduction to a binary choice'. It is not a\n")
cat("  reduction; the binary is settled and it settles the wrong way.\n")

cat("\n\n=== B. is there a NON-circular criterion left? Hadamard is absolute.\n\n")
cat("  With a constant reflectivity r the seam's two-point function is\n")
cat("     W(x,y) = W0(x-y) + r W0(x - ybar),   ybar = the mirror image of y.\n")
cat("  Both sides of the seam are physical here (two sheets), so x = ybar is reachable\n")
cat("  by two DISTINCT points. Check whether W is singular there, and whether those\n")
cat("  two points are null-related. A spacelike singularity is a Hadamard violation\n")
cat("  and needs no reference state, so it is not circular.\n\n")
W0 <- function(s2) 1/s2                       # massless, W ~ 1/(interval)
# 1+1D: x=(t,z) on sheet 1, y=(t',z') on sheet 2 mapped to image (t', -z')
iv  <- function(A,B) -(A[1]-B[1])^2 + (A[2]-B[2])^2
cat("     r      t sep   z(x)  z(y)   interval(x,y)  interval(x,ybar)   W direct   W image\n")
for (r in c(0, 0.3)) for (dz in c(0.5, 0.1, 0.01)) {
  x <- c(0, dz); y <- c(0, -dz); ybar <- c(0, dz)     # image of y across the seam at z=0
  cat(sprintf("  %5.2f %8.1f %6.2f %6.2f %14.4f %18.2e %10.3f %9s\n", r, 0, x[2], y[2],
      iv(x,y), iv(x,ybar), W0(iv(x,y)), if (r==0) "0 (none)" else sprintf("%.2e", r*W0(iv(x,ybar)))))
}
cat("\n  x and y sit at mirror positions, spacelike separated (interval > 0, here 4 dz^2\n")
cat("  at zero time separation), yet the IMAGE interval is exactly zero, so the image\n")
cat("  term diverges between two spacelike-separated distinct points. No Hadamard\n")
cat("  state does that. At r = 0 the term is absent and the divergence is not there.\n")

cat("\n=== C. and a frequency-dependent r removes it, which is the whole distinction\n\n")
cat("  Smear the image term over frequency with r(w) and ask whether the coincidence\n")
cat("  limit stays finite: I(eps) = integral dw r(w) cos(w eps) for small eps.\n\n")
for (nm in list(c("r constant = 0.3", function(w) 0.3+0*w),
                c("r = 0.3 exp(-w/W)", function(w) 0.3*exp(-w/20)),
                c("r = 0.3/(1+(w/W)^2)", function(w) 0.3/(1+(w/20)^2)))) {
  f <- nm[[2]]
  v <- sapply(c(1, 0.1, 0.01, 0.001), function(e)
        integrate(function(w) f(w)*cos(w*e), 0, 2000, subdivisions=20000, rel.tol=1e-9)$value)
  cat(sprintf("  %-24s eps=1:%9.3f  0.1:%9.3f  0.01:%10.3f  0.001:%11.3f\n", nm[[1]], v[1],v[2],v[3],v[4]))
}
cat("\n  The constant reflectivity's integral grows without bound as the separation\n")
cat("  closes; both decaying ones converge to a finite number. So Hadamard is exactly\n")
cat("  the condition r(w) -> 0, stated about ONE state rather than about a comparison,\n")
cat("  and the circularity both referees found in the relative-entropy version is not\n")
cat("  present in this version.\n")
cat("\n  Within the corner term r is frequency-independent, so Hadamard holds only at\n")
cat("  r = 0, kappa = 1. Same conclusion as before, reached without a reference state.\n")
