# INDEPENDENT CHECK: does a separable out-pair marginal cost more production than n_*?
# If so, the dark-matter mass is an entanglement witness, and that is a strong claim, so this
# file rebuilds the operator from scratch rather than quoting a form.
#
# Out operators a,b; in operators a_- = c a + s b^dag, b_- = c b - s a^dag, c^2 = 1-P, s^2 = P.
# With N_+- the number PER PAIR MEMBER in each region, working out (N_+ + N_-)/2 gives
#     Q = [ (1-P)(n_a+n_b) + P + sqrt(P(1-P)) (a^dag b^dag + b a) ] / 2 .
Qmat <- function(P) {                       # basis |00>,|01>,|10>,|11>
  cs <- sqrt(P*(1-P)); M <- matrix(0,4,4); nab <- c(0,1,1,2)
  for (i in 1:4) M[i,i] <- ((1-P)*nab[i] + P)/2
  M[1,4] <- cs/2; M[4,1] <- cs/2; M
}
cat("=== 1. the reconstructed operator reproduces the paper's spectrum ===\n")
cat("      P     eigenvalues                              n_*         max error\n")
for (P in c(0.05,0.2,0.5,0.8,0.95)) {
  ev <- sort(eigen(Qmat(P), symmetric=TRUE)$values); ns <- (1-sqrt(1-P))/2
  cat(sprintf("  %5.2f  %s %10.6f %12.2e\n", P, paste(sprintf("%8.5f",ev),collapse=" "),
      ns, max(abs(ev - sort(c(ns,.5,.5,1-ns))))))
}

cat("\n=== 2. the minimum over SEPARABLE pairs, with the contact condition imposed ===\n")
cat("  A linear functional is minimised on an extreme point, so separable reduces to product\n")
cat("  pure states. The step from <Q> to the production integral needs <N_+> = <N_->, so that\n")
cat("  constraint must be imposed on the competitor too. Omitting it was my first error and it\n")
cat("  gave 393.8 PeV instead of 270.3. By symmetry the minimum sits at |beta_a|^2 = |beta_b|^2\n")
cat("  = u, and the constraint becomes P u - P/2 + sqrt(P(1-P)) u(1-u) = 0.\n\n")
nsym <- function(P){
  if (P < 1e-10) return(sqrt(P)/2)          # asymptotic branch; the root degenerates as P -> 0
  cs <- sqrt(P*(1-P))
  tryCatch(uniroot(function(u) P*u - P/2 + cs*u*(1-u), c(1e-14,0.5), tol=1e-15)$root,
           error=function(e) sqrt(P)/2)
}
cat("      P        n_*        n_sep     ratio\n")
for (P in c(0.05,0.2,0.5,0.8,0.95)) {
  ns <- (1-sqrt(1-P))/2; np <- nsym(P)
  cat(sprintf("  %6.3f %10.6f %10.6f %8.3f\n", P, ns, np, np/ns))
}
cat("\n  Small-P behaviour matters for the integral: n_sep -> sqrt(P)/2 while n_* -> P/4, so the\n")
cat("  separable floor decays as exp(-x^2/2) against the adopted state's exp(-x^2).\n")

I0 <- integrate(function(x) x^2*(1-sqrt(1-exp(-x^2)))/2, 0, 40, subdivisions=8000)$value/pi^2
Is <- integrate(function(x) x^2*sapply(x, function(y) nsym(exp(-y^2))), 0, 40, subdivisions=8000)$value/pi^2
cat(sprintf("\n=== 3. the consequence ===\n  I_0   = %.12f\n  I_sep = %.12f   ratio %.6f\n", I0, Is, Is/I0))
cat(sprintf("  M_sep = 491.6 (I_0/I_sep)^(2/5) = %.4f PeV, half-mass line %.4f PeV\n",
    491.6*(I0/Is)^(2/5), 491.6*(I0/Is)^(2/5)/2))
cat("\n  READING. If every out-pair marginal were separable the abundance would cap the mass at\n")
cat("  270.30 PeV, not 491.6. A securely inferred mass above 270.30, or a primary line above\n")
cat("  135.15, therefore requires entanglement in pairs carrying nonzero weight in the integral.\n")
cat("  The dark-matter mass becomes an entanglement witness.\n")
cat("\n  WHAT IT IS NOT. It inherits every assumption of the bound it refines: free Bogoliubov\n")
cat("  evolution, the stated asymptotic particle definitions, per-block equality of in- and\n")
cat("  out-region expectations. It witnesses entanglement in the OUT-PAIR marginals of this\n")
cat("  model, not entanglement between the two sheets, and not a measured mass, since no mass\n")
cat("  has been measured. It says what a measurement would certify if one were made.\n")
