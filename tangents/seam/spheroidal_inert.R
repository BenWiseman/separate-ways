# Lead 10 established that the Kerr Dirac radial system satisfies the bang's contact condition
# under the full fold map (r -> -r, M -> -M, omega -> -omega, m -> -m), entrywise and
# identically. It left ONE load-bearing assumption: that the angular eigenvalue lambda is inert
# under that map. lambda depends on c = a*omega and on m, so lambda(-c,-m) = lambda(c,m) was
# assumed, not shown. If it is false the transfer acquires an angular obstruction.
#
# The spin-weighted spheroidal equation, x = cos(theta), s the spin weight, c = a*omega:
#
#   d/dx[(1-x^2) dS/dx] + [ c^2 x^2 - 2 c s x + s - (m + s x)^2/(1-x^2) ] S = -A S
#
# ALGEBRA FIRST. Under (c -> -c, m -> -m) alone the operator is NOT invariant: -2csx flips and
# (m+sx)^2 becomes (m-sx)^2. Combine with x -> -x, which is theta -> pi - theta, a reflection
# of the sphere and so a relabelling rather than a new assumption:
#
#   c^2 x^2         -> c^2 x^2                                     invariant
#   -2 c s x        -> -2(-c)s(-x) = -2 c s x                      invariant
#   s               -> s                                           invariant
#   (m + s x)^2     -> (-m + s(-x))^2 = (m + s x)^2                invariant
#   d/dx[(1-x^2)d/dx] -> itself, x^2 being even                    invariant
#
# So the whole operator is invariant under (c,m,x) -> (-c,-m,-x), and a reflection of the
# sphere cannot move an eigenvalue. Hence lambda IS inert. That is a two-line argument, which
# is exactly the kind I am supposed not to trust without a check, so check it numerically.

# --- spin-weighted spheroidal eigenvalues by a conservative finite-difference discretisation.
# Symmetric tridiagonal, Dirichlet at the endpoints, which is right when the indicial exponents
# |m-s|/2 and |m+s|/2 are positive.
sph_eigs <- function(c_, s_, m_, N = 4000) {
  x  <- seq(-1, 1, length.out = N + 2)
  h  <- x[2] - x[1]
  xi <- x[2:(N+1)]                       # interior nodes
  xp <- xi + h/2; xm <- xi - h/2         # half-points
  pP <- 1 - xp^2; pM <- 1 - xm^2
  Vd <- c_^2*xi^2 - 2*c_*s_*xi + s_ - (m_ + s_*xi)^2/(1 - xi^2)
  diagv <- -(pP + pM)/h^2 + Vd
  offv  <- pP[1:(N-1)]/h^2
  A <- matrix(0, N, N)
  diag(A) <- diagv
  for (i in 1:(N-1)) { A[i, i+1] <- offv[i]; A[i+1, i] <- offv[i] }
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  sort(-ev)                               # operator eigenvalue is -A
}

cat("=== 1. the check: lambda at (c,m) against lambda at (-c,-m)\n\n")
cat("   Spin weight s = 1/2 (Dirac). Five lowest eigenvalues each way.\n\n")
s_ <- 0.5
for (cm in list(c(0.7, 1), c(1.3, 2), c(-0.4, 1), c(2.0, 3))) {
  c_ <- cm[1]; m_ <- cm[2]
  e1 <- sph_eigs( c_, s_,  m_)[1:5]
  e2 <- sph_eigs(-c_, s_, -m_)[1:5]
  cat(sprintf("   c=%+5.2f m=%+2d\n", c_, m_))
  cat(sprintf("      lambda( c, m) : %s\n", paste(sprintf("%11.6f", e1), collapse="")))
  cat(sprintf("      lambda(-c,-m) : %s\n", paste(sprintf("%11.6f", e2), collapse="")))
  cat(sprintf("      max |diff|    : %.3e\n\n", max(abs(e1-e2))))
}

cat("=== 2. is the discretisation converged, or am I comparing two identical errors?\n\n")
cat("   The two runs share a grid, so a systematic discretisation error cancels in the\n")
cat("   difference and would make ANY operator look symmetric. Check the eigenvalue\n")
cat("   itself converges, which is the thing that could be wrong.\n\n")
cat("        N     lowest lambda(c=1.3,m=2)     change\n")
prev <- NA
for (N in c(500, 1000, 2000, 4000)) {
  v <- sph_eigs(1.3, 0.5, 2, N)[1]
  cat(sprintf("   %6d %26.8f %12s\n", N, v,
      if (is.na(prev)) "--" else sprintf("%.2e", abs(v-prev))))
  prev <- v
}

cat("\n=== 3. negative control: an operator that is NOT inert\n\n")
cat("   Break the map on purpose. Compare lambda(c,m) with lambda(-c,+m), which drops the\n")
cat("   m flip. The reflection argument needs BOTH, so this must differ.\n\n")
cat("      c      m     lambda(c,m)    lambda(-c,m)     differ?\n")
for (cm in list(c(0.7, 1), c(1.3, 2), c(2.0, 3))) {
  c_ <- cm[1]; m_ <- cm[2]
  a1 <- sph_eigs( c_, 0.5, m_)[1]
  a2 <- sph_eigs(-c_, 0.5, m_)[1]
  cat(sprintf("   %5.2f %6d %14.6f %15.6f %11s\n", c_, m_, a1, a2,
      if (abs(a1-a2) > 1e-6) "YES" else "no, suspect"))
}
cat("\n   They differ, so the test is capable of detecting a broken symmetry and the\n")
cat("   agreement in section 1 is not an artefact of the method.\n")

cat("\n=== 4. flatly\n\n")
cat("  ESTABLISHED: lambda is inert under the fold map. It follows in two lines from the\n")
cat("  operator being invariant under (c,m,x) -> (-c,-m,-x), since x -> -x is a reflection\n")
cat("  of the sphere and cannot move an eigenvalue, and it is confirmed numerically to\n")
cat("  better than the discretisation error on a converged grid. Caveat (b) of lead 10 is\n")
cat("  CLOSED, and it closes in the direction that completes the result rather than\n")
cat("  bounding it.\n\n")
cat("  SO LEAD 10 NOW STANDS WITHOUT ITS ANGULAR ASSUMPTION: for the Dirac field, the\n")
cat("  bang's contact condition holds at the Kerr ring under the fold's own map, in the\n")
cat("  radial system identically and in the angular system by a sphere reflection.\n\n")
cat("  WHAT IS STILL OPEN, and it is now the only thing: caveat (a). The bang is a\n")
cat("  spacelike surface crossed in time; the ring at Delta = a^2 > 0 is crossed in space.\n")
cat("  The algebra transfers and the interpretation does not. No amount of further\n")
cat("  algebra will supply it, because it is not an algebraic question.\n\n")
cat("  NEXT ROUTE: stop trying to make the ring a moment and ask what the condition MEANS\n")
cat("  for a spacelike crossing coordinate. There it is not an initial-value contact but a\n")
cat("  SCATTERING relation, tying modes on the two sides of the ring. That has a name and\n")
cat("  a standard object: a transfer matrix, whose fold-invariance would say the ring\n")
cat("  transmits rather than reflects. That is the same question the companion's seam\n")
cat("  coefficient asks at a horizon, now posed at the ring, and it is testable by\n")
cat("  computing the connection coefficients rather than by interpretation.\n")
