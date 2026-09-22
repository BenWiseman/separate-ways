# Does minimising marginal entropy pick out the minimum-occupation member of the band?
# state_selection.R answers yes using the BOSONIC marginal entropy of a two-mode squeezed
# pair. The paper's produced quanta are fermions, and its own entanglement passage uses the
# fermionic form. Check whether the argument survives the change of statistics.
Sb <- function(n) ifelse(n <= 0, 0, (n+1)*log(n+1) - n*log(n))   # boson, as used in state_selection.R
Sf <- function(n) ifelse(n <= 0 | n >= 1, 0, -n*log(n) - (1-n)*log(1-n))  # fermion

cat("  P      n_min      n_max=1-n_min    S_boson(min)  S_boson(max)   S_fermi(min)  S_fermi(max)\n")
for (P in c(0.05, 0.2, 0.5, 0.8, 0.95)) {
  nmin <- (1-sqrt(1-P))/2; nmax <- 1-nmin
  cat(sprintf(" %4.2f  %9.6f  %13.6f  %13.6f %13.6f  %13.6f %13.6f\n",
      P, nmin, nmax, Sb(nmin), Sb(nmax), Sf(nmin), Sf(nmax)))
}
cat("\n  Boson: S is strictly increasing, so the minimum-occupation member is the minimum-entropy\n")
cat("  member and the argument selects. Fermion: S(n) = S(1-n) identically, so the two ends of\n")
cat("  the band carry EXACTLY the same marginal entropy and the argument selects nothing.\n\n")
d <- sapply(c(0.05,0.2,0.5,0.8,0.95), function(P){n<-(1-sqrt(1-P))/2; Sf(n)-Sf(1-n)})
cat(sprintf("  max |S_fermi(n_min) - S_fermi(n_max)| over those P: %.3e\n", max(abs(d))))
cat("\n  VERDICT: for fermions the entropy argument is vacuous. What still selects the minimum\n")
cat("  is the energy prescription, and separately the divergent number density as n_max -> 1.\n")
cat("  Both are already in the paper; the entropy sentence is the one that has to go.\n")
