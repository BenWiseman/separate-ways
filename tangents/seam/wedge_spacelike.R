# The load-bearing claim, and both referees asked to see it. The specific challenge:
# "mirror points across a bifurcation surface are generically null- or timelike-related,
# not spacelike." If that is right the Hadamard argument fails. Test it.
#
# The fold's two sheets are the two Rindler wedges of a bifurcate Killing horizon:
#   right wedge R = { z >  |t| },   left wedge L = { z < -|t| },
# and J maps one to the other. Claim to test: EVERY pair with x in R and y in L is
# spacelike separated - not generically, always.

set.seed(17)
inR <- function(p) p[2] >  abs(p[1])
inL <- function(p) p[2] < -abs(p[1])
iv  <- function(a,b) -(a[1]-b[1])^2 + (a[2]-b[2])^2     # signature -+ : >0 spacelike

cat("=== 1. sample opposite wedges hard and look for ANY non-spacelike pair\n\n")
n <- 400000; worst <- Inf; bad <- 0
for (i in 1:n) {
  t1 <- rnorm(1,0,3); z1 <- abs(t1) + rexp(1,0.5)          # in R
  t2 <- rnorm(1,0,3); z2 <- -abs(t2) - rexp(1,0.5)         # in L
  s <- iv(c(t1,z1), c(t2,z2))
  if (s <= 0) bad <- bad + 1
  if (s < worst) worst <- s
}
cat(sprintf("  pairs drawn: %d    non-spacelike found: %d    smallest interval: %.4e\n", n, bad, worst))
cat("\n  The proof is two lines and the sampling only confirms it: z1 > |t1| and\n")
cat("  z2 < -|t2| give z1 - z2 > |t1| + |t2| >= |t1 - t2|, so (dz)^2 > (dt)^2 and the\n")
cat("  interval is positive. ALWAYS spacelike, with no genericity qualifier. So the\n")
cat("  challenge does not land: opposite wedges are spacelike by construction, which is\n")
cat("  the defining property of the wedge pair, not an accident of a chosen slice.\n")

cat("\n=== 2. and the image locus, where the divergence sits, lies in that set\n\n")
cat("  W(x,y) = W0(x-y) + r W0(x - Jy). The image term blows up at x = Jy, i.e. when x\n")
cat("  and y are mirror images. Check those pairs specifically:\n\n")
cat("        t        z        interval(x,y)   both wedges?   image interval\n")
for (t in c(0, 0.3, 1.1, 2.5)) for (z in c(0.6, 2.0)) {
  if (z <= abs(t)) next
  x <- c(t, z); y <- c(-t, -z)                              # J: (t,z) -> (-t,-z)
  cat(sprintf("   %7.2f %8.2f %14.4f %14s %16.2e\n", t, z, iv(x,y),
      if (inR(x) && inL(y)) "yes" else "NO", iv(x, c(-y[1],-y[2]))))
}
cat("\n  Every mirror pair is spacelike separated from its partner, and the image\n")
cat("  interval is exactly zero there, so the image term diverges at every one of them.\n")
cat("  This is not one contrived configuration: it is a codimension-zero family.\n")

cat("\n=== 3. so how bad is it? Compare the two terms approaching the image locus\n\n")
cat("  Move x toward Jy along a spacelike direction, parameter eps:\n\n")
cat("       eps      direct W0(x-y)     image r W0(x-Jy), r = 0.3      ratio\n")
y <- c(-0.4, -1.5); Jy <- c(0.4, 1.5)
for (e in c(0.5, 0.1, 0.01, 1e-3, 1e-4)) {
  x <- Jy + c(0, e)
  d <- 1/abs(iv(x,y)); im <- 0.3/abs(iv(x,Jy))
  cat(sprintf("   %8.0e %17.4f %26.4e %14.2e\n", e, d, im, im/d))
}
cat("\n  The image term overwhelms the direct one by an arbitrarily large factor as the\n")
cat("  locus is approached, so it cannot be absorbed into a redefinition or treated as\n")
cat("  a finite correction. A two-point function with a divergence at spacelike\n")
cat("  separation violates the microlocal spectrum condition, so the state is not\n")
cat("  Hadamard. That is the step the argument needs, and it is stronger than stated\n")
cat("  before: the mirror pairs are not merely sometimes spacelike, they always are.\n")

cat("\n=== 4. what still has to be granted\n\n")
cat("  That Hadamard regularity is required of the seam state. It is the standard\n")
cat("  admissibility condition, needed for the stress tensor to renormalise, but it is\n")
cat("  an assumption here and a reader may decline it. If declined, the argument gives\n")
cat("  nothing and kappa is free again. Everything else above is geometry.\n")
