# A survival review proposes replacing the minimisation with an OPERATOR INEQUALITY. If it holds
# it is strictly stronger than choosing a state, and it is not what the prior contact-condition
# literature does: they select a vacuum by energy, this bounds an observable over ALL states.
# Verify the algebra before believing a word of it.
#
# Pair Fock space |n_a n_b>, basis |00>,|10>,|01>,|11>. Jordan-Wigner.
I2<-diag(2); sz<-diag(c(1,-1)); sm<-matrix(c(0,0,1,0),2,2)
a <- sm %x% I2 ; b <- sz %x% sm
ad<-t(Conj(a)); bd<-t(Conj(b))
stopifnot(max(abs(a%*%ad+ad%*%a-diag(4)))<1e-12, max(abs(a%*%b+b%*%a))<1e-12)

check <- function(P) {
  c_<-sqrt(1-P); s_<-sqrt(P)
  am <- c_*a + s_*bd            # a_- = c a + s b^dagger
  bm <- c_*b - s_*ad            # b_- = c b - s a^dagger
  amd<-t(Conj(am)); bmd<-t(Conj(bm))
  Np <- (ad%*%a + bd%*%b)/2
  Nm <- (amd%*%am + bmd%*%bm)/2
  Q  <- (Np+Nm)/2
  ev <- sort(Re(eigen(Q)$values))
  nstar <- (1-sqrt(1-P))/2
  list(ev=ev, nstar=nstar, bogo=max(abs(am%*%amd+amd%*%am-diag(4))))
}
cat("      P        eigenvalues of Q (sorted)                n_* = (1-sqrt(1-P))/2   min(ev)-n_*   Bogo ok\n")
for (P in c(0.05,0.2,0.5,0.8,0.95)) {
  r<-check(P)
  cat(sprintf("  %6.2f   %s %14.6f %14.2e %10.1e\n", P,
      paste(sprintf("%8.5f", r$ev), collapse=" "), r$nstar, min(r$ev)-r$nstar, r$bogo))
}
cat("\n  The two middle eigenvalues are the odd-parity states at exactly 1/2; the outer pair are\n")
cat("  n_* and 1-n_*. So the least eigenvalue of Q is n_* and Q >= n_* holds as an operator\n")
cat("  inequality on the whole four-dimensional block, for every density operator on it.\n")

cat("\n=== does it reproduce the paper's production integral?\n\n")
nst <- function(x) (1-sqrt(1-exp(-x^2)))/2
I <- integrate(function(x) x^2*nst(x), 0, Inf)$value/pi^2
cat(sprintf("   I from the operator bound's n_* with P = exp(-x^2): %.15f\n", I))
cat(sprintf("   the paper's quoted I:                               0.012759667363447\n"))
cat(sprintf("   difference: %.2e\n", abs(I-0.012759667363447)))

cat("\n=== what this changes, stated carefully\n\n")
cat("  IF the identification of Q with an observable the CPT condition constrains is right, the\n")
cat("  mass ceiling stops being 'the state we chose' and becomes 'no admissible state does\n")
cat("  better'. That is a different KIND of claim and it survives mixed states, correlations\n")
cat("  between momenta, and any failure of the global state to factorise, none of which a\n")
cat("  minimisation over pure Gaussian states covers.\n")
cat("  WHAT IT STILL ASSUMES, and the paper must say so: free Bogoliubov evolution, the stated\n")
cat("  asymptotic particle definitions, and equal in/out number expectations on a\n")
cat("  symmetry-complete block. It is not an interacting-theory result.\n")
cat("  DEGENERACY: at P = 1 the even block is degenerate, so the zero-momentum limit must be\n")
cat("  excluded from any uniqueness statement.\n")
