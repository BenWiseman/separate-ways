# Is matter conservation a separate assumption, or does diffeomorphism invariance already have it?
#
# The relativity ledger lists grad^a T_ab = 0 as an input to the Clausius step, next to
# diffeomorphism invariance. If Noether's second theorem already delivers the first from the
# second then the ledger double counts and the assumed list is one line shorter than it says.
#
# THE ARGUMENT. For a matter action S_m[g, psi] invariant under diffeomorphisms, varying by a
# vector field xi gives dg_ab = 2 grad_(a xi_b) and dpsi = L_xi psi, so
#     0 = dS_m = int sqrt(-g) T^ab grad_a xi_b  +  int (dS_m/dpsi) L_xi psi.
# On shell the second term drops and the first, integrated by parts for compactly supported xi,
# leaves grad_a T^ab = 0. The step that could fail is whether T^ab defined as the metric
# variation really is the object whose divergence the matter field equation kills.
#
# WHAT IS CHECKED. That step, as an identity that holds for ARBITRARY phi, not only solutions:
#     grad^a T_ab = (box phi - m^2 phi) d_b phi.
# Being off shell, it needs no field equation solved, so it can be measured directly on a
# concrete curved metric with an arbitrary field. If it holds, then on shell the right side is
# zero by the matter equation the Clausius step already assumes, and conservation is not an
# independent input.

set.seed(4111)
N <- 4
eta <- diag(c(-1, 1, 1, 1))
Pc  <- array(rnorm(N*N*N, sd = 1), c(N, N, N))        # linear coefficients, symmetrised below
Qc  <- array(rnorm(N*N*N*N, sd = 1), c(N, N, N, N))   # quadratic coefficients
for (a in 1:N) for (b in 1:N) { Pc[a,b,] <- Pc[b,a,] <- (Pc[a,b,] + Pc[b,a,])/2
                                Qc[a,b,,] <- Qc[b,a,,] <- (Qc[a,b,,] + Qc[b,a,,])/2 }

gmet <- function(x) {                                 # a generic curved Lorentzian metric
  h <- matrix(0, N, N)
  for (a in 1:N) for (b in 1:N)
    h[a,b] <- sum(Pc[a,b,]*x) + sum(outer(x,x) * Qc[a,b,,])
  eta + 0.10*h
}
phi <- function(x) 0.7*x[1] - 0.4*x[2]*x[3] + 0.3*sin(1.3*x[4]) + 0.25*x[1]^2*x[2] - 0.15*x[3]^3
mm  <- 0.83

d1 <- function(f, x, a, h) { e <- rep(0,N); e[a] <- h; (f(x+e) - f(x-e))/(2*h) }
d2 <- function(f, x, a, b, h) {
  if (a == b) { e <- rep(0,N); e[a] <- h; (f(x+e) - 2*f(x) + f(x-e))/h^2 }
  else { ea <- rep(0,N); ea[a] <- h; eb <- rep(0,N); eb[b] <- h
         (f(x+ea+eb) - f(x+ea-eb) - f(x-ea+eb) + f(x-ea-eb))/(4*h^2) }
}

Gamma <- function(x, h) {                             # G^a_bc from first derivatives of g
  gi <- solve(gmet(x)); dg <- lapply(1:N, function(c) d1(gmet, x, c, h))
  G <- array(0, c(N,N,N))
  for (a in 1:N) for (b in 1:N) for (c in 1:N)
    G[a,b,c] <- sum(sapply(1:N, function(d) gi[a,d]*(dg[[c]][d,b] + dg[[b]][d,c] - dg[[d]][b,c])))/2
  G
}

Tdn <- function(x, h) {                               # T_ab for a massive scalar
  g <- gmet(x); gi <- solve(g)
  dp <- sapply(1:N, function(a) d1(phi, x, a, h))
  kin <- as.numeric(t(dp) %*% gi %*% dp)
  outer(dp, dp) - g*(kin + mm^2*phi(x)^2)/2
}
Tmix <- function(x, h) solve(gmet(x)) %*% Tdn(x, h)   # T^a_b

divT <- function(x, h) {                              # grad_a T^a_b
  G <- Gamma(x, h); Tm <- Tmix(x, h)
  out <- numeric(N)
  for (b in 1:N) {
    s <- sum(sapply(1:N, function(a) d1(function(y) Tmix(y, h)[a,b], x, a, h)))
    s <- s + sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[a,a,c]*Tm[c,b]))))
    s <- s - sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[c,a,b]*Tm[a,c]))))
    out[b] <- s
  }
  out
}

boxphi <- function(x, h) {                            # (1/sqrt(-g)) d_a ( sqrt(-g) g^ab d_b phi )
  sq <- function(y) sqrt(-det(gmet(y)))
  f <- function(a) function(y) sq(y) * sum(solve(gmet(y))[a,] * sapply(1:N, function(b) d1(phi, y, b, h)))
  sum(sapply(1:N, function(a) d1(f(a), x, a, h))) / sq(x)
}

x0 <- c(0.20, -0.30, 0.15, 0.25)
cat("=== 1. the setting ===\n")
g0 <- gmet(x0)
cat(sprintf("   a generic curved metric at x0; eigenvalue signs %s, det %.6f\n",
            paste(sign(eigen(g0)$values), collapse = " "), det(g0)))
cat(sprintf("   departure from flat: max |g - eta| = %.4f, so the curvature is not a perturbation\n",
            max(abs(g0 - eta))))
cat(sprintf("   arbitrary phi, NOT a solution of anything; phi(x0) = %.6f\n", phi(x0)))

cat("\n=== 1b. the connection is validated before it is used ===\n")
cat("   Metric compatibility grad_c g_ab = 0 is not assumed by the code above, so it is a\n")
cat("   genuine check on the Christoffels the whole result rests on.\n")
compat <- function(x, h) {
  G <- Gamma(x, h); dg <- lapply(1:N, function(c) d1(gmet, x, c, h)); m <- 0
  for (c in 1:N) for (a in 1:N) for (b in 1:N) {
    v <- dg[[c]][a,b] - sum(sapply(1:N, function(d) G[d,c,a]*gmet(x)[d,b] + G[d,c,b]*gmet(x)[a,d]))
    m <- max(m, abs(v))
  }
  m
}
for (h in c(4e-3, 1e-3)) cat(sprintf("      step %8.0e   max |grad_c g_ab| = %.2e\n", h, compat(x0, h)))
cat("   and a deliberately wrong connection must break it:\n")
Gamma_bad <- Gamma; Gamma <- function(x, h) { G <- Gamma_bad(x, h); G[1,2,3] <- G[1,2,3] + 0.3; G }
cat(sprintf("      one Christoffel shifted by 0.3  ->  max |grad_c g_ab| = %.4f\n", compat(x0, 1e-3)))
stopifnot(compat(x0, 1e-3) > 0.05)
Gamma <- Gamma_bad
stopifnot(compat(x0, 1e-3) < 1e-5)
cat("      the check can fail, and does not.\n")

cat("\n=== 2. the off-shell identity, measured ===\n")
cat("   claim: grad^a T_ab = (box phi - m^2 phi) d_b phi, for any phi whatsoever.\n\n")
cat("        step h      max|LHS - RHS|    max|LHS|     relative\n")
for (h in c(4e-3, 2e-3, 1e-3)) {
  lhs <- divT(x0, h)
  rhs <- (boxphi(x0, h) - mm^2*phi(x0)) * sapply(1:N, function(b) d1(phi, x0, b, h))
  cat(sprintf("     %9.0e   %14.2e   %11.4f   %11.2e\n", h, max(abs(lhs-rhs)), max(abs(lhs)),
              max(abs(lhs-rhs))/max(abs(lhs))))
  if (h == 1e-3) { LHS <- lhs; RHS <- rhs }
}
cat("\n   component by component at h = 1e-3:\n        b      grad^a T_ab        (box - m^2)d_b phi\n")
for (b in 1:N) cat(sprintf("        %d   %+15.8f   %+18.8f\n", b, LHS[b], RHS[b]))
rel <- max(abs(LHS-RHS))/max(abs(LHS))
cat(sprintf("\n   agreement to %.1e relative, limited by the finite differencing and improving\n", rel))
cat("   as the step shrinks. The identity holds off shell.\n")
stopifnot(rel < 2e-4)

cat("\n=== 3. so on shell ===\n")
cat("   The right-hand side is the Klein-Gordon operator times d_b phi. Once the matter is on\n")
cat("   shell it is zero and the divergence is zero with it. Nothing about gravity was used.\n")

cat("\n=== 4. the plant: a stress tensor no covariant action could produce ===\n")
cat("   Add a piece built on a direction fixed in the chart rather than by the geometry, which\n")
cat("   is exactly what diffeomorphism invariance forbids. Conservation must then fail.\n")
Tbad <- function(x, h) { M <- matrix(0,N,N); M[1,1] <- phi(x)^2; solve(gmet(x)) %*% M }
divbad <- function(x, h) {
  G <- Gamma(x, h); Tm <- Tbad(x, h); out <- numeric(N)
  for (b in 1:N) {
    s <- sum(sapply(1:N, function(a) d1(function(y) Tbad(y, h)[a,b], x, a, h)))
    s <- s + sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[a,a,c]*Tm[c,b]))))
    s <- s - sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[c,a,b]*Tm[a,c]))))
    out[b] <- s
  }
  out
}
db <- divbad(x0, 1e-3)
cat(sprintf("   divergence of the non-covariant piece: %s\n", paste(sprintf("%+.5f", db), collapse = "  ")))
cat(sprintf("   largest component %.5f, nowhere near zero. The plant fires.\n", max(abs(db))))
stopifnot(max(abs(db)) > 1e-2)

cat("\n=== 5. the verdict ===\n")
cat("   The divergence of the stress tensor IS the matter field equation contracted with the\n")
cat("   field, identically, in a generic curved metric with an arbitrary field. Conservation is\n")
cat("   therefore not an extra input: it is what a diffeomorphism-invariant matter action gives\n")
cat("   once the matter is on shell, and the Clausius step already needs the matter on shell.\n")
cat("   The ledger's assumed list loses a line.\n")
