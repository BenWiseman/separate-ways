# If J maps F to P, what does it do to the singularities? In Kruskal the curvature
# singularity r = 0 is the hyperbola T^2 - X^2 = 1: the future branch T > 0 bounds F,
# the past branch T < 0 bounds P. J is (T,X) -> (-T,-X).

onSing <- function(p) abs(p[1]^2 - p[2]^2 - 1) < 1e-12
br <- function(p) if (p[1] > 0) "future (bounds F)" else "past (bounds P)"

cat("=== 1. does J exchange the two singular branches?\n\n")
cat("        X        T = +sqrt(1+X^2)      branch            J(point)        branch\n")
for (X in c(0, 0.7, 1.8, 4.0)) {
  p <- c(sqrt(1+X^2), X); Jp <- -p
  cat(sprintf("   %7.2f %18.4f %20s %8.2f,%7.2f %18s\n", X, p[1], br(p), Jp[1], Jp[2], br(Jp)))
}
ok <- all(replicate(20000, { X <- rnorm(1,0,4); s <- sample(c(-1,1),1)
     p <- c(s*sqrt(1+X^2), X); onSing(p) && onSing(-p) && br(p) != br(-p) }))
cat(sprintf("\n  20000 random points on the singularity: all map to the other branch: %s\n", ok))

cat("\n=== 2. so what the fold identifies\n\n")
cat("  The future singularity of the black hole and the past singularity of the white\n")
cat("  hole are the SAME locus under the fold's map, as the two exteriors are the same\n")
cat("  locus and F and P are. Under that identification the singularity is not where\n")
cat("  the two sheets each end. It is where they are joined, and it is the only place\n")
cat("  other than the bifurcation surface where they touch.\n")

cat("\n=== 3. two things this does NOT establish, stated before they are assumed\n\n")
cat("  It does not resolve the singularity. The curvature still diverges on that locus\n")
cat("  in both sheets; identifying two divergences is not removing them, and nothing\n")
cat("  here supplies the finite-curvature interior a bounce would need.\n")
cat("  And it is again the maximally extended solution. A collapse hole has no past\n")
cat("  singularity to identify with, so for an astrophysical black hole this says\n")
cat("  nothing, exactly as the F <-> P statement says nothing there.\n")

cat("\n=== 4. what it does buy\n\n")
cat("  A consistency condition rather than an endpoint condition. A single-sheet theory\n")
cat("  must say what happens AT the singularity as a terminal boundary condition. The\n")
cat("  fold instead requires the data on that locus to be J-compatible with itself,\n")
cat("  which is a self-consistency requirement of the same kind as the half-advanced\n")
cat("  half-retarded condition of section 3.3, and which is the interior counterpart of\n")
cat("  alpha^2 = 1 being the smoothness condition at the bifurcation surface.\n")
cat("  Whether that requirement has solutions is not addressed here and is the obvious\n")
cat("  next calculation: the fold's two touching loci, B and the singularity, carry the\n")
cat("  same Z2, and only one of them has been worked out.\n")
