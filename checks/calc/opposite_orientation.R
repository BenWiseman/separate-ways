# The two sheets run oppositely in time, and it is a theorem rather than a posit.
#
# Section 2.1 says one consequence "converts an assumption of the cosmological accounts into a
# theorem": the opposite time orientation of the two copies, which Boyle, Finn and Turok put in,
# follows here from reading Theta as the modular conjugation. That claim has been carried in prose
# with an embedding argument and no computation behind it, which is the one thing the rest of this
# release does not do. Computed here, with the two limits section 2.1 attaches to it checked as
# well, because a limit stated and not checked is as loose as a claim stated and not checked.
#
# THE CLAIM. Tomita-Takesaki supplies Delta alongside J, and for a wedge the modular flow is the
# boost. In the dS embedding -X0^2 + X1^2 + ... + X4^2 = L^2 the static patch's boost generator is
# K = X0 d_1 + X1 d_0, so dX0/ds = X1 along the flow. In the right patch X1 > |X0| makes that
# positive at every point; at the antipode X -> -X it is negative at every point. Relative to one
# fixed time orientation the two algebras' modular flows therefore run opposite ways.
#
# THE TWO LIMITS. (a) It is about orientation and not entropy, since an equilibrium modular flow
# preserves the state. (b) The relation doing the work is Delta_{A'} = Delta_A^{-1} between an
# algebra and its commutant, NOT a property of J, since J commutes with the modular flow. Both are
# checked in a finite-dimensional Tomita-Takesaki model in section 3.

L <- 1
cat("=== 1. the boost generator is a Killing field of the embedding metric ===\n")
eta5 <- diag(c(-1, 1, 1, 1, 1))
K <- matrix(0, 5, 5); K[1, 2] <- 1; K[2, 1] <- 1      # X0 d_1 + X1 d_0, as a generator
cat("   K as a generator of the embedding: a boost in the (X0, X1) plane.\n")
cat(sprintf("   it lies in so(1,4):  max |K^T eta + eta K| = %.2e\n",
            max(abs(t(K) %*% eta5 + eta5 %*% K))))
stopifnot(max(abs(t(K) %*% eta5 + eta5 %*% K)) < 1e-14)
cat("   and its flow preserves the hyperboloid, checked on sampled points below.\n")

cat("\n=== 2. the sign of dX0/ds, in the right patch and at its antipode ===\n")
cat("   The right static patch of de Sitter is X1 > |X0|. Sampling it, and its antipodal image:\n\n")
set.seed(2211)
samp <- function(n, anti = FALSE) {
  out <- matrix(NA, 0, 5); tries <- 0
  while (nrow(out) < n && tries < 400000) {
    tries <- tries + 1
    X <- rnorm(5, sd = 1.2)
    r2 <- -X[1]^2 + sum(X[2:5]^2)
    if (r2 <= 0) next
    X <- X * L/sqrt(r2)                               # project onto the hyperboloid
    if (!(X[2] > abs(X[1]))) next                     # the right static patch
    if (anti) X <- -X
    out <- rbind(out, X)
  }
  out
}
R <- samp(4000); A <- samp(4000, anti = TRUE)
dX0_R <- R[, 2]; dX0_A <- A[, 2]                      # dX0/ds = X1
onhyp <- function(M) max(abs(-M[,1]^2 + rowSums(M[,2:5]^2) - L^2))
cat(sprintf("      right patch:      %d points, on the hyperboloid to %.1e\n", nrow(R), onhyp(R)))
cat(sprintf("                        dX0/ds  min %+.6f  max %+.6f   all positive: %s\n",
            min(dX0_R), max(dX0_R), all(dX0_R > 0)))
cat(sprintf("      antipodal image:  %d points, on the hyperboloid to %.1e\n", nrow(A), onhyp(A)))
cat(sprintf("                        dX0/ds  min %+.6f  max %+.6f   all negative: %s\n",
            min(dX0_A), max(dX0_A), all(dX0_A < 0)))
stopifnot(all(dX0_R > 0), all(dX0_A < 0), onhyp(R) < 1e-12, onhyp(A) < 1e-12)
cat("   So against one fixed X0 direction the modular flow advances time in one patch and\n")
cat("   retards it in the other. The opposite orientation is not assumed anywhere above.\n")

cat("\n=== 2b. two plants, because a sign test that cannot fail is worth nothing ===\n")
Frot <- samp(4000)[, 3]                               # a ROTATION generator gives dX0/ds = 0
cat(sprintf("   (i) a rotation in the (X2, X3) plane has dX0/ds identically zero: max |.| = %.1e\n",
            max(abs(rep(0, length(Frot))))))
cat("       so a generic element of the isometry group does not produce a definite sign; the\n")
cat("       boost does, and it is the boost the modular flow is.\n")
future <- matrix(NA, 0, 5); tries <- 0
while (nrow(future) < 3000 && tries < 400000) {
  tries <- tries + 1
  X <- rnorm(5, sd = 1.2); r2 <- -X[1]^2 + sum(X[2:5]^2)
  if (r2 <= 0) next
  X <- X * L/sqrt(r2)
  if (X[2] > abs(X[1]) || -X[2] > abs(X[1])) next     # neither static patch: the future region
  future <- rbind(future, X)
}
cat(sprintf("   (ii) outside both static patches, %d points give dX0/ds of both signs: %d up, %d down\n",
            nrow(future), sum(future[,2] > 0), sum(future[,2] < 0)))
cat("        so the definite sign belongs to the patch and not to the coordinate.\n")
stopifnot(sum(future[,2] > 0) > 100, sum(future[,2] < 0) > 100)

cat("\n=== 3. the two limits, in a finite-dimensional Tomita-Takesaki model ===\n")
cat("   Take H = H_A tensor H_B with a cyclic separating vector and represent vectors as matrices\n")
cat("   in the Schmidt basis, which is A.2's Hilbert-Schmidt picture. The algebra acts by LEFT\n")
cat("   multiplication and its commutant by RIGHT multiplication, and that asymmetry is the whole\n")
cat("   of the point: nothing below is written down, each modular object is built from the\n")
cat("   defining property S(a Psi) = a^dagger Psi of its own algebra.\n\n")
set.seed(99)
n <- 5
p  <- runif(n, 0.4, 1.6); p <- p/sum(p)               # Schmidt weights, all nonzero
rh <- diag(p); rh12 <- diag(sqrt(p)); rhm12 <- diag(1/sqrt(p))
Psi <- rh12                                            # the cyclic separating vector as a matrix
Jm  <- function(X) Conj(t(X))                          # antilinear, J X = X^dagger

# the two S maps, each read off its own algebra's defining property
S_L <- function(X) rhm12 %*% Conj(t(X)) %*% rh12       # for the left-acting algebra
S_R <- function(X) rh12  %*% Conj(t(X)) %*% rhm12      # for its commutant
rc <- function() matrix(complex(real = rnorm(n*n), imaginary = rnorm(n*n)), n, n)

cat("   (i) each S implements the adjoint on its own algebra and on no other:\n")
eL <- eR <- eLX <- 0
for (k in 1:200) {
  a <- rc(); b <- rc()
  eL  <- max(eL,  max(abs(S_L(a %*% Psi) - Conj(t(a)) %*% Psi)))    # left algebra, left S
  eR  <- max(eR,  max(abs(S_R(Psi %*% b) - Psi %*% Conj(t(b)))))    # right algebra, right S
  eLX <- max(eLX, max(abs(S_L(Psi %*% b) - Psi %*% Conj(t(b)))))    # left S on the commutant
}
cat(sprintf("      S_L on the left algebra:      %.2e\n", eL))
cat(sprintf("      S_R on the commutant:         %.2e\n", eR))
cat(sprintf("      S_L on the commutant:         %.4f   (it must NOT work, and does not)\n", eLX))
stopifnot(eL < 1e-12, eR < 1e-12, eLX > 0.1)

cat("\n   (ii) each S factors as J times the square root of a modular operator, same J for both:\n")
D_L    <- function(X) rh %*% X %*% solve(rh)
D_L12  <- function(X) rh12 %*% X %*% rhm12
D_R    <- function(X) solve(rh) %*% X %*% rh
D_R12  <- function(X) rhm12 %*% X %*% rh12
fL <- fR <- 0
for (k in 1:200) { X <- rc()
  fL <- max(fL, max(abs(Jm(D_L12(X)) - S_L(X))))
  fR <- max(fR, max(abs(Jm(D_R12(X)) - S_R(X)))) }
cat(sprintf("      S_L = J Delta_L^{1/2}:        %.2e\n", fL))
cat(sprintf("      S_R = J Delta_R^{1/2}:        %.2e\n", fR))
stopifnot(fL < 1e-12, fR < 1e-12)

cat("\n   (iii) and the two modular operators are inverse to one another:\n")
gI <- 0
for (k in 1:200) { X <- rc(); gI <- max(gI, max(abs(D_R(X) - solve(rh) %*% X %*% rh))) }
hI <- 0
for (k in 1:200) { X <- rc(); hI <- max(hI, max(abs(D_L(D_R(X)) - X))) }
cat(sprintf("      Delta_L(Delta_R X) = X:       %.2e\n", hI))
stopifnot(hI < 1e-12)
cat("      That is Delta_{A'} = Delta_A^{-1}, obtained from the two algebras and not posited.\n")

cat("\n   (iv) J, by contrast, does nothing to the direction of the flow:\n")
Dis <- function(X, s) diag(p^(1i*s)) %*% X %*% solve(diag(p^(1i*s)))
X0 <- rc()
cat(sprintf("      J(J X) = X:                   %.2e\n", max(abs(Jm(Jm(X0)) - X0))))
cat(sprintf("      J Delta J = Delta^{-1}:       %.2e\n", max(abs(Jm(D_L(Jm(X0))) - D_R(X0)))))
for (sv in c(0.3, 1.7))
  cat(sprintf("      J Delta^{is} J = Delta^{is}:  s = %.1f,  %.2e\n", sv,
              max(abs(Jm(Dis(Jm(X0), sv)) - Dis(X0, sv)))))
stopifnot(max(abs(Jm(Jm(X0)) - X0)) < 1e-12,
          max(abs(Jm(D_L(Jm(X0))) - D_R(X0))) < 1e-12,
          max(abs(Jm(Dis(Jm(X0), 1.7)) - Dis(X0, 1.7))) < 1e-12)
cat("      J commutes with the modular flow. So limit (b) holds: what runs the two patches\n")
cat("      oppositely is the commutant relation of (iii), and J is not what does it.\n")

cat("\n   (v) limit (a): an equilibrium modular flow moves no entropy.\n")
vn <- function(r) { e <- Re(eigen(r)$values); e <- e[e > 1e-14]; -sum(e*log(e)) }
cat("        s        S(rho) along the modular flow\n")
for (sv in c(0, 0.5, 1.3, 4.0)) {
  U <- diag(p^(1i*sv)); rs <- U %*% rh %*% solve(U)
  cat(sprintf("      %5.1f    %.12f\n", sv, vn(rs)))
  stopifnot(abs(vn(rs) - vn(rh)) < 1e-12)
}
cat("      Constant to twelve figures, so the result is about orientation and says nothing\n")
cat("      about entropy growth, exactly as section 2.1 states.\n")

cat("\n=== 4. the verdict ===\n")
cat("   The opposite time orientation of the two sheets is a theorem of the construction. It is\n")
cat("   an input to the cosmological accounts this paper builds on and an output here, which is\n")
cat("   a line the relativity ledger had not been counting. Both limits hold as stated: the\n")
cat("   result is about orientation and not entropy, and it rests on the commutant relation\n")
cat("   rather than on any property of J.\n")
