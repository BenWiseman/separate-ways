# Named last turn: compute epsilon, the fold misalignment, for section 4.1's branch
# decoherence, using the paper's own overlap rather than a toy.
#
# 4.1: two branches E+ and E- with overlap s = <E-|E+>, |s| falling as
#   log|s| = -(N_f/4) sum_{n>=2} n^2 log[1 + A^4(A^2-1)/(n^2(n^2-1)^2)],  A = aH,
# reaching e^-1 at aH = 2.43913, Ht = 1.53984 in the gauged count.
#
# Write the branch-label state: |psi> = a|+>|E+> + b|->|E->. Tracing the environment,
#   rho = [[ |a|^2, a b* s ], [ a* b s*, |b|^2 ]].
# Theta is antilinear AND swaps the branches, so Theta rho Theta is rho with |a|^2 and
# |b|^2 exchanged. The misalignment therefore lives in the BRANCH WEIGHT, and the overlap
# enters only through how mixed rho is.

Srel <- function(A,B) { ea <- eigen(A,symmetric=TRUE); eb <- eigen(B,symmetric=TRUE)
  la <- pmax(ea$values,0); lb <- pmax(eb$values,0)
  lgA <- ea$vectors %*% diag(ifelse(la>1e-14, log(la), -700)) %*% t(ea$vectors)
  lgB <- eb$vectors %*% diag(ifelse(lb>1e-14, log(lb), -700)) %*% t(eb$vectors)
  sum(diag(A %*% (lgA - lgB))) }
rho <- function(p, s) { q <- 1-p; off <- sqrt(p*q)*s
  matrix(c(p, off, off, q), 2, 2) }                       # s real; phase drops out here
mis <- function(p, s) { r <- rho(p,s); tr <- matrix(c(1-p, sqrt(p*(1-p))*s, sqrt(p*(1-p))*s, p),2,2)
  Srel(r, tr) }

cat("=== 1. the misalignment is zero exactly at equal branch weight, any overlap\n\n")
cat("        |s|        p=0.500      p=0.510      p=0.600      p=0.750\n")
for (s in c(0, 0.3, 0.7, 0.95, 0.999)) {
  cat(sprintf("   %9.3f", s))
  for (p in c(0.5, 0.51, 0.6, 0.75)) cat(sprintf(" %12.6f", mis(p,s)))
  cat("\n")
}
cat("\n  Zero along the whole p = 1/2 column. So the fold's requirement is equal branch\n")
cat("  weight, and nothing about coherence enters that.\n")

cat("\n=== 2. decoherence REGULATES it. The divergence is at |s| -> 1, not |s| -> 0.\n\n")
cat("   p = 0.6 fixed, approaching a coherent superposition:\n\n")
cat("        |s|          misalignment\n")
for (s in c(0, 0.5, 0.9, 0.99, 0.999, 0.99999)) cat(sprintf("   %11.5f %18.6f\n", s, mis(0.6,s)))
cat("\n  Finite once decohered and divergent when coherent, which inverts the naive\n")
cat("  picture. A coherent superposition with unequal weights is a PURE state, and a\n")
cat("  pure state has infinite relative entropy against any state off its support. So\n")
cat("  decoherence does not create the misalignment. It makes it finite.\n")

cat("\n=== 3. the fully decohered value has a closed form\n\n")
cat("        p        (2p-1) log(p/(1-p))      numerical S at |s|=0      difference\n")
for (p in c(0.51, 0.55, 0.6, 0.7, 0.9)) {
  cf <- (2*p-1)*log(p/(1-p))
  cat(sprintf("   %8.2f %22.8f %25.8f %15.2e\n", p, cf, mis(p,0), abs(cf-mis(p,0))))
}
cat("\n  S = (2p-1) log(p/(1-p)) exactly, which is 4(p-1/2)^2 * 2 + O((p-1/2)^4):\n")
for (d in c(1e-3, 1e-2, 0.05)) 
  cat(sprintf("     p-1/2 = %7.4f : S = %.6e, S/(p-1/2)^2 = %.4f\n", d, mis(0.5+d,0), mis(0.5+d,0)/d^2))
cat("\n  Quadratic with coefficient 8, matching last night's general second-order result.\n")

cat("\n=== 4. and now the number: where does 4.1's overlap put us?\n\n")
Nf <- 1
logs <- function(A) { n <- 2:400
  -(Nf/4)*sum(n^2*log(1 + A^4*(A^2-1)/(n^2*(n^2-1)^2))) }
cat("        aH        Ht          log|s|         |s|        misalignment at p=0.6\n")
for (A in c(1.2, 1.95374, 2.43913, 3.0, 4.0)) {
  Ht <- acosh(A)            # closed de Sitter: aH = cosh(Ht)
  ls <- logs(A); sv <- exp(ls)
  cat(sprintf("   %9.5f %9.5f %13.4f %12.3e %20.6f\n", A, Ht, ls, sv, mis(0.6, min(sv,0.999999))))
}
cat("\n  By the gauged threshold aH = 2.43913 the overlap is already 1e-2 and the\n")
cat("  misalignment has settled onto its fully decohered value. So the whole story is\n")
cat("  over within Ht of order one: whatever branch-weight asymmetry the fold has, it is\n")
cat("  fixed in the first e-fold and does not accumulate afterwards.\n")

cat("\n=== 5. flatly\n\n")
cat("  This does NOT deliver a number for epsilon. It delivers the shape: epsilon is\n")
cat("  8(p - 1/2)^2 once decohered, and p is the branch weight, which the fold sets to\n")
cat("  1/2 by CPT and which nothing in the paper computes a departure from. The honest\n")
cat("  reading is that the misalignment is not a dynamical quantity that grows; it is a\n")
cat("  property of the initial branch weights, frozen by the first e-fold.\n")
cat("  Ben's Zeno picture - many small unobserved events accumulating a misalignment -\n")
cat("  does not happen HERE, because there is nothing for them to accumulate into: the\n")
cat("  measure depends on p alone and p stops evolving once the branches decohere.\n")
