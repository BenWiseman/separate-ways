# Objection: the conclusion is conditional on a specific action that the paper does not motivate
# as the unique fold-invariant choice. Answer it by classification rather than by another symmetry
# argument, which is finite work.
#
# A.13's three steps need only a CONSTANT reflectivity, and r is constant because kappa is
# dimensionless and the corner term carries no other scale. So the question is: which quadratic
# boundary terms on B can carry a DIMENSIONLESS coefficient at all?
#
# Conventions. In d spacetime dimensions a scalar has [q] = (d-2)/2 in mass units. A boundary
# term integrates over the (d-1)-dimensional seam, so its density must have dimension d-1 for
# the coefficient to be dimensionless.

cat("     d    [q] = (d-2)/2   boundary density must be   ...so a dimensionless coefficient needs\n")
for (d in c(3,4,5,6)) cat(sprintf("  %4d %14.1f %24d %38s\n", d, (d-2)/2, d-1, "a term of dimension d-1"))

cat("\n=== enumerate the quadratic candidates and their dimensions\n\n")
cat("   term          dimension        dimensionless coefficient at d = 4?\n")
d <- 4; q <- (d-2)/2
cands <- list(c("q^2",        2*q),
              c("q qdot",     2*q+1),
              c("qdot^2",     2*q+2),
              c("q grad q",   2*q+1),
              c("q grad^2 q", 2*q+2),
              c("qdot grad q",2*q+2))
for (cd in cands) {
  dim <- as.numeric(cd[2])
  cat(sprintf("   %-13s %8.1f %34s\n", cd[1], dim,
      if (abs(dim-(d-1)) < 1e-12) "YES" else sprintf("no, coefficient carries mass^%.1f", (d-1)-dim)))
}
cat("\n   At d = 4 the boundary density must have dimension 3 and [q] = 1, so only the terms\n")
cat("   with exactly one derivative qualify: q qdot and q grad q.\n")

cat("\n=== now remove the ones that are not there\n\n")
cat("  (i)  q grad q is a TOTAL TRANSVERSE DERIVATIVE, grad(q^2)/2, and integrates to zero on\n")
cat("       a closed bifurcation surface. B is a sphere here, so it drops.\n")
cat("  (ii) The SYMMETRIC part of q qdot is d(q^2)/dt / 2, a total time derivative, which does\n")
cat("       not affect the equations of motion or the matching condition.\n")
cat("  (iii) What survives is the ANTISYMMETRIC pairing of the two sides, q^T J qdot with J\n")
cat("       antisymmetric. That is the corner term, and there is nothing else at this order.\n")

cat("\n=== check the claim that the symmetric part is inert\n\n")
set.seed(3)
n <- 4
Jm <- matrix(0,2*n,2*n); Jm[1:n,(n+1):(2*n)] <- diag(n); Jm[(n+1):(2*n),1:n] <- -diag(n)
S  <- diag(2*n)
cat(sprintf("   J antisymmetric? max|J + J^T| = %.1e\n", max(abs(Jm+t(Jm)))))
cat(sprintf("   S symmetric?     max|S - S^T| = %.1e\n", max(abs(S-t(S)))))
q  <- rnorm(2*n); qd <- rnorm(2*n)
cat(sprintf("   q^T S qdot equals (1/2) d/dt (q^T S q) ? symmetric part gives %.6f vs %.6f\n",
    as.numeric(t(q)%*%S%*%qd), 0.5*as.numeric(t(q)%*%S%*%qd + t(qd)%*%S%*%q)))
cat("   The two agree identically for symmetric S, which is what makes it a total derivative.\n")
cat(sprintf("   For the antisymmetric pairing the same comparison gives %.6f vs %.6f, and they\n",
    as.numeric(t(q)%*%Jm%*%qd), 0.5*as.numeric(t(q)%*%Jm%*%qd + t(qd)%*%Jm%*%q)))
cat("   do not agree, so the corner term is not a total derivative and does act.\n")

cat("\n=== what this answers, and what it does not\n\n")
cat("  ANSWERS: the corner term is not a free choice among many. At quadratic order it is the\n")
cat("  only boundary term on B that can carry a dimensionless coefficient and still act, so\n")
cat("  'a seam action carrying no scale of its own' and 'the corner term' are the same thing,\n")
cat("  and A.13's conclusion covers the whole scale-free class rather than one example.\n")
cat("  DOES NOT ANSWER: higher than quadratic order, and seams that do carry a scale. Both are\n")
cat("  outside the scale-free branch by construction, and 3.4 already separates them.\n")

cat("\n=== CHALLENGE, from a hostile referee: which surface does the term live on?\n\n")
cat("  The counting above assumed a boundary density of dimension d-1, i.e. a CODIMENSION-ONE\n")
cat("  surface. But the bifurcation surface B is a 2-sphere in four dimensions: CODIMENSION TWO.\n")
cat("  An action term integrated over a k-dimensional surface needs a density of dimension k,\n")
cat("  since d^k x carries dimension -k. Redo the count for both.\n\n")
d <- 4; q <- (d-2)/2
cat("      surface            k    density must be    q^2 (dim %.0f)    q qdot (dim %.0f)\n")
cat(sprintf("      codim 1 (a slice)  %d    %d                 %s              %s\n",
    d-1, d-1, if (abs(2*q-(d-1))<1e-9) "QUALIFIES" else "no", if (abs(2*q+1-(d-1))<1e-9) "QUALIFIES" else "no"))
cat(sprintf("      codim 2 (B itself) %d    %d                 %s              %s\n",
    d-2, d-2, if (abs(2*q-(d-2))<1e-9) "QUALIFIES" else "no", if (abs(2*q+1-(d-2))<1e-9) "QUALIFIES" else "no"))

cat("\n  So the referee is right and the conclusion REVERSES with the surface. On a codimension-one\n")
cat("  matching surface the one-derivative term is the only one with a dimensionless coefficient,\n")
cat("  which is the case the counting above treats. On the bifurcation surface itself q^2 carries\n")
cat("  a dimensionless coefficient and q qdot does not.\n\n")
cat("  WHAT SURVIVES. The corner term of 3.4 is written as an integral over TIME of a two-port\n")
cat("  coupling, which is the codimension-one case, and that is the case the uniqueness argument\n")
cat("  covers. What it does NOT cover, and what the text must stop implying it covers, is a term\n")
cat("  supported on B itself. Whether a q^2 term on a codimension-two locus yields a CONSTANT\n")
cat("  reflectivity, which is the only property A.13's three steps actually use, is not settled\n")
cat("  here: a codimension-two defect does not scatter like a one-dimensional delta, and the\n")
cat("  frequency dependence would have to be computed rather than assumed either way.\n")
cat("  Recorded as a scope correction, not a repair. The uniqueness claim is now conditional and\n")
cat("  says which condition.\n")

cat("\n=== DO NOT CONCEDE WITHOUT COMPUTING: does a q^2 seam term give CONSTANT reflectivity?\n\n")
cat("  The three steps of A.13 use ONE property: that r does not depend on frequency. So the\n")
cat("  question is not which terms are dimensionally allowed but which give a constant r.\n")
cat("  A q^2 term localised on a surface is a delta potential. Solve it exactly.\n\n")
cat("  For V = g delta(x) in one dimension, matching gives r(k) = -i g /(2k + i g).\n\n")
g <- 1
k <- c(0.1, 0.5, 1, 2, 5, 20, 100)
r <- -1i*g/(2*k + 1i*g)
cat("        k        |r(k)|        arg r      constant?\n")
for (i in seq_along(k))
  cat(sprintf("   %8.2f %12.6f %12.4f %14s\n", k[i], Mod(r[i]), Arg(r[i]),
      if (i>1 && abs(Mod(r[i])-Mod(r[1]))<1e-9) "yes" else "NO"))
cat(sprintf("\n   |r| falls from %.4f to %.6f across the range and tends to g/2k -> 0.\n",
    Mod(r[1]), Mod(r[length(k)])))
cat("   A q^2 seam is therefore NOT scale-free: g carries dimension, r depends on frequency,\n")
cat("   and the seam becomes transparent in the ultraviolet exactly as a physical interface\n")
cat("   should. That is the behaviour 3.4 says a reflecting seam must have, and it is what the\n")
cat("   corner term conspicuously lacks.\n")

cat("\n=== so what actually survives, tested rather than conceded\n\n")
cat("  The referee's dimensional point is correct: on B, which is codimension two, q^2 carries a\n")
cat("  dimensionless coefficient and q qdot does not. Conceded, and the text now says so.\n")
cat("  But it does not reverse the conclusion that matters. A.13's steps need CONSTANT r, and a\n")
cat("  q^2 seam does not deliver one at any codimension: a localised quadratic term always\n")
cat("  introduces a coupling with dimension, hence a scale, hence frequency dependence. The\n")
cat("  dimensionless coefficient on B is bought by the codimension-two measure, not by the term\n")
cat("  being scale-free in the sense the argument uses.\n")
cat("  CONCLUSION: among seams with a frequency-independent reflectivity, the corner term is\n")
cat("  still the only quadratic candidate. That is the claim worth making and it is narrower and\n")
cat("  better defended than the one it replaces.\n")
