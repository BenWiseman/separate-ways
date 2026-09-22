# Section 6 sets the paper's own test: "If the modular structure does not fix the seam
# coefficient, the algebra has not paid for a geometric quantity that is its own business."
# A.13 fixes the coefficient by HADAMARD REGULARITY, which is generic QFT admissibility and
# needs no modular theory. So the charge stands. The question nobody has run: can the
# modular FLOW pay it, where regularity did?
#
# The flow acts geometrically as the boost (Bisognano-Wichmann, Sewell). A reflecting seam
# adds an image term W_0(x - Jy) to the two-point function. The flow can only constrain the
# seam if that image term BREAKS boost covariance. So the whole question reduces to one
# commutator: does the wedge reflection commute with the boost?

B <- function(eta) matrix(c(cosh(eta), sinh(eta), sinh(eta), cosh(eta)), 2, 2, byrow=TRUE)
Jm <- diag(c(-1,-1))                       # wedge reflection (T,X) -> (-T,-X)

cat("=== 1. the commutator that decides it\n\n")
cat("        eta      max |J B(eta) - B(eta) J|\n")
for (eta in c(0.2, 0.7, 1.5, 3.0)) 
  cat(sprintf("   %9.2f %26.2e\n", eta, max(abs(Jm %*% B(eta) - B(eta) %*% Jm))))
cat("\n  Zero. The wedge reflection COMMUTES with the boost, so the image point of a boosted\n")
cat("  point is the boost of the image point.\n")

cat("\n=== 2. therefore the image term is boost-covariant, and the flow is blind to it\n\n")
set.seed(3)
iv <- function(a,b) -(a[1]-b[1])^2 + (a[2]-b[2])^2
worst <- 0
for (i in 1:20000) {
  x <- rnorm(2); y <- rnorm(2); eta <- rnorm(1)
  # direct term and image term, before and after a boost of BOTH points
  d0 <- iv(x, y);        i0 <- iv(x, Jm %*% y)
  d1 <- iv(B(eta)%*%x, B(eta)%*%y); i1 <- iv(B(eta)%*%x, Jm %*% (B(eta)%*%y))
  worst <- max(worst, abs(d0-d1), abs(i0-i1))
}
cat(sprintf("   max change in the direct OR image interval under a boost, 20000 draws: %.2e\n", worst))
cat("\n  Both terms are boost invariant. A reflecting seam therefore preserves the boost\n")
cat("  invariance of the two-point function, so the state's modular flow still acts as\n")
cat("  the boost whatever the reflectivity is.\n")

cat("\n=== 3. FLATLY: the modular flow CANNOT pay the charge\n\n")
cat("  It is not that we failed to find the constraint. The flow is structurally incapable\n")
cat("  of supplying one, because the reflection it would have to detect commutes with it.\n")
cat("  Section 6's charge cannot be settled by the modular flow, and saying so closes a\n")
cat("  stated open question rather than leaving it hanging.\n")

cat("\n=== 4. so what DOES the algebra pay for? Exactly one thing, and it is a re-derivation.\n\n")
cat("  Separate the modular structure into its two pieces and account for each:\n\n")
cat("    the modular CONJUGATION J : supplies the map. Without it there is no alpha =\n")
cat("        J o P_perp at all, so the pairing does not exist. This is real and it is the\n")
cat("        whole reason the two fields can be composed.\n")
cat("    the modular FLOW sigma_t  : supplies NOTHING about the seam, by part 1.\n")
cat("    the KMS property          : alpha^2 = 1 reads W(t - i beta) = W(t) and returns the\n")
cat("        Gibbons-Hawking temperature. An algebraic statement fixing a geometric\n")
cat("        quantity, so the algebra DOES pay here - but A.10 already records that this\n")
cat("        is the smoothness condition in disguise, which is a re-derivation of a\n")
cat("        standard result rather than new ground.\n")
cat("\n  Honest accounting: the algebra supplies the map, and pays once, for a quantity the\n")
cat("  conical argument already delivers. That is a thinner fusion than 'the algebra does\n")
cat("  geometric work' and it is what the paper can defend.\n")

cat("\n=== 5. the route, same turn\n\n")
cat("  The flow is blind because the reflection commutes with the boost. The way to make\n")
cat("  the algebra pay is therefore to find a seam property that does NOT commute with the\n")
cat("  boost. A frequency-independent r commutes trivially. A seam carrying the horizon's\n")
cat("  own scale gives r(omega/T_H), and omega is conjugate to the boost parameter, so\n")
cat("  such a seam is NOT boost-invariant. That is exactly A.13's second branch, reached\n")
cat("  from the algebra rather than from Hadamard, and it is the one place the modular\n")
cat("  structure could still do work the geometry cannot.\n")
