# Ben: "can a wave function collapse on a singularity?" A.15 found the singularity's
# consistency condition keeps the +1 eigenspace of JU, through the projector (1 + JU)/2.
# A projector is either a von Neumann measurement or it is not, and the difference is
# whether it is SELF-ADJOINT. Idempotent alone gives an oblique projection, which is not a
# measurement. So the question is decidable rather than interpretive.

set.seed(23)
n <- 6; N <- 2*n; Im_ <- diag(n); Z <- matrix(0,n,n)
J  <- rbind(cbind(Im_,Z), cbind(Z,-Im_))
Om <- rbind(cbind(Z,Im_), cbind(-Im_,Z))
sym <- function(){ M <- matrix(rnorm(n*n),n,n); (M+t(M))/2 }
A  <- rbind(cbind(Z,sym()), cbind(sym(),Z))
EV <- eigen(A); U <- Re(EV$vectors %*% diag(exp(EV$values)) %*% solve(EV$vectors))
L  <- J %*% U
P  <- (diag(N) + L)/2

cat("=== 1. the three properties that decide it\n\n")
cat(sprintf("   L is an involution      |L^2 - 1|     = %.2e\n", max(abs(L%*%L - diag(N)))))
cat(sprintf("   P is idempotent         |P^2 - P|     = %.2e\n", max(abs(P%*%P - P))))
cat(sprintf("   P is SELF-ADJOINT?      |P - P^T|     = %.4f\n", max(abs(P - t(P)))))
cat(sprintf("   L is symmetric?         |L - L^T|     = %.4f\n", max(abs(L - t(L)))))
cat(sprintf("   U is symmetric?         |U - U^T|     = %.4f\n", max(abs(U - t(U)))))

cat("\n=== 2. so the answer is NO, and the reason is specific\n\n")
cat("  P is idempotent but NOT self-adjoint, because U is symplectic rather than symmetric\n")
cat("  and JU inherits that. An idempotent that is not self-adjoint is an OBLIQUE\n")
cat("  projection: it has a well-defined range and kernel, but they are not orthogonal\n")
cat("  complements, so it does not correspond to any von Neumann measurement and there is\n")
cat("  no Born rule attached to it.\n")

cat("\n=== 3. how far from a measurement is it? Measure the obliquity.\n\n")
cat("  The angle between the range of P and the orthogonal complement of its kernel. Zero\n")
cat("  degrees would mean an orthogonal projector and a genuine measurement.\n\n")
prang <- function(P) {
  R <- svd(P)$u[, 1:qr(P)$rank, drop=FALSE]              # range
  K <- svd(diag(nrow(P)) - P)$u[, 1:qr(diag(nrow(P))-P)$rank, drop=FALSE]  # kernel
  s <- svd(t(R) %*% K)$d
  90 - acos(pmin(pmax(s,0),1))*180/pi }                   # principal angles from orthogonality
cat("        evolution time      max obliquity (deg)   0 = orthogonal\n")
for (t in c(0.01, 0.1, 0.5, 1.0, 2.0)) {
  Ut <- Re(EV$vectors %*% diag(exp(EV$values*t)) %*% solve(EV$vectors))
  Pt <- (diag(N) + J %*% Ut)/2
  cat(sprintf("   %16.2f %22.3f\n", t, max(abs(prang(Pt)))))
}
cat("\n  It is orthogonal only in the limit of NO evolution between the branches, and the\n")
cat("  departure grows with the evolution. So the collapse reading is exact precisely\n")
cat("  where A.15 already found the two projectors share a sector - at zero evolution,\n")
cat("  which is the Wheeler-DeWitt case.\n")

cat("\n=== 4. flatly, and the part that survives\n\n")
cat("  NO: the singularity does not perform a measurement in general. Its projector is\n")
cat("  oblique, has no Born rule, and cannot be read as collapse.\n")
cat("  YES, in one place: in the constraint formulation, where there is no evolution\n")
cat("  between the branches, the projector IS orthogonal and the collapse reading is\n")
cat("  exact. That is the same timeless setting section 2.5 already commits to for three\n")
cat("  independent reasons, and this is a fourth arriving from a different direction.\n")
cat("  So Ben's question has an answer and it is conditional on the paper's own\n")
cat("  commitment rather than on taste.\n")
