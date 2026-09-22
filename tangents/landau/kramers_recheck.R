# kramers_check.R reports Theta^2 = +1 on the whole Fock space. A verification pass says its
# build() writes the transformed images into swapped basis columns and that the true answer is
# (-1)^F. Rebuilt here from the operator action, with the implementing unitary written out
# explicitly rather than exponentiated, so no convention can be inherited from the old file.
#
# Basis order |00>, |01>, |10>, |11>, with c1+|00> = |10> and c2+|00> = |01>.
e <- function(i){v <- rep(0,4); v[i] <- 1; v}
ok <- function(M,tol=1e-12) max(Mod(M)) < tol
Fpar <- diag(c(1,-1,-1,1))

# MINUS exchange, the one A.18 states: c1 -> c2, c2 -> -c1.
#   |10> = c1+|00>  ->  c2+|00> = |01>
#   |01> = c2+|00>  -> -c1+|00> = -|10>
#   |11> = c1+c2+|00> -> c2+(-c1+)|00> = +c1+c2+|00> = |11>
Um <- cbind(e(1), -e(3), e(2), e(4))
# PLUS exchange: c1 -> c2, c2 -> +c1.  Then c1+c2+ -> c2+c1+ = -c1+c2+, so |11> -> -|11>.
Up <- cbind(e(1), e(3), e(2), -e(4))

for (nm in c("minus","plus")) {
  U <- if (nm=="minus") Um else Up
  cat(sprintf("\n=== %s exchange ===\n", nm))
  cat("  U unitary:", ok(U %*% t(Conj(U)) - diag(4)), "\n")
  Th2 <- U %*% Conj(U)              # Theta = U K, so Theta^2 = U conj(U)
  cat("  Theta^2 diagonal:", paste(round(Re(diag(Th2)),12), collapse=", "), "\n")
  cat(sprintf("  ||Theta^2 - (-1)^F|| = %.3e     ||Theta^2 - (+1)|| = %.3e\n",
              max(Mod(Th2-Fpar)), max(Mod(Th2-diag(4)))))
  # which paired states does Theta leave invariant?
  cat("  invariance of sqrt(1-n)|00> + i sqrt(n)|11>:\n")
  for (nn in c(0.05,0.25,0.5,0.75,0.95)) {
    v <- c(sqrt(1-nn),0,0,1i*sqrt(nn))
    cat(sprintf("     n=%.2f  ||Theta psi - psi|| = %.3e\n", nn, max(Mod(U %*% Conj(v) - v))))
  }
  cat("  invariance of sqrt(1-n)|00> + sqrt(n)|11>  (real coefficients):\n")
  for (nn in c(0.25,0.75)) {
    v <- c(sqrt(1-nn),0,0,sqrt(nn))
    cat(sprintf("     n=%.2f  ||Theta psi - psi|| = %.3e\n", nn, max(Mod(U %*% Conj(v) - v))))
  }
}
cat("\n  READING. The minus exchange gives Theta^2 = (-1)^F, not +1. Kramers still does not kill\n")
cat("  the construction, because the paired vacuum sector is EVEN and Theta^2 = +1 there; but\n")
cat("  'Theta^2 = +1 on every state' is false, and the odd sector is Kramers-degenerate.\n")
cat("  And the exchange sign is NOT forced by the existence of an invariant paired state:\n")
cat("  each sign admits one, differing only by a phase convention on the |11> amplitude.\n")
