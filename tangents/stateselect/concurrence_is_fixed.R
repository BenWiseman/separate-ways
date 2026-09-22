# How MUCH entanglement does the contact condition require, not just whether any is required.
#
# On the even block {|00>,|11>} write the occupation per member as n = p_11 and the coherence as
# rho_01. From the reconstructed Q,
#     <N_+> = n,   <Q> = (P/2)(1-n) + ((2-P)/2) n + sqrt(P(1-P)) Re(rho_01),
# so the contact condition <N_+> = <Q> forces
#     Re(rho_01) = P (n - 1/2) / sqrt(P(1-P)),
# and the concurrence of the block, C = 2|rho_01|, therefore obeys
#     C >= 2 |n - 1/2| sqrt(P/(1-P)).
# The condition does not merely permit entanglement. It prescribes a floor for it.
Cmin <- function(n,P) 2*abs(n-0.5)*sqrt(P/(1-P))
nstar <- function(P) (1-sqrt(1-P))/2
Cpure <- function(n) 2*sqrt(n*(1-n))              # concurrence of sqrt(1-n)|00> + sqrt(n)|11>

cat("  At the adopted least-occupied state the floor is saturated, and equals sqrt(P):\n\n")
cat("      P       n_*        C_min(n_*)    C of the pure state    sqrt(P)\n")
for (P in c(0.05,0.2,0.5,0.8,0.95)) {
  ns <- nstar(P)
  cat(sprintf("  %5.2f %10.6f %13.8f %22.8f %11.6f\n", P, ns, Cmin(ns,P), Cpure(ns), sqrt(P)))
}
cat("\n  Algebra, not coincidence: n_*(1-n_*) = P/4 exactly, so C = 2 sqrt(P/4) = sqrt(P), and\n")
cat("  the floor 2|n_*-1/2| sqrt(P/(1-P)) = sqrt(1-P) sqrt(P/(1-P)) = sqrt(P) as well.\n")

cat("\n  So the concurrence of each produced pair is fixed by its momentum:\n")
cat("      C(x) = sqrt(P) = exp(-x^2/2)\n\n")
cat("        x        P = e^{-x^2}      C(x)\n")
for (x in c(0,0.5,1,1.5,2,3)) cat(sprintf("  %7.1f %17.6e %12.6f\n", x, exp(-x^2), exp(-x^2/2)))

cat("\n  The entanglement is not a free parameter of the model. It is the square root of the\n")
cat("  Landau-Zener transition probability, mode by mode, and it is maximal in the deep infrared\n")
cat("  where C -> 1 and the pair is a Bell state.\n")

W <- function(f) integrate(function(x) x^2*f(x), 0, 30, subdivisions=8000)$value
cat(sprintf("\n  Mode-weighted mean concurrence, int x^2 C dx / int x^2 dx over the support that\n"))
cat(sprintf("  carries the abundance (x < 3): %.6f\n",
    integrate(function(x) x^2*exp(-x^2/2),0,3)$value/integrate(function(x) x^2,0,3)$value))

cat("\n  WHAT THIS IS NOT. It is a property of the adopted state within the two-level reduction,\n")
cat("  with the same free-Bogoliubov and contact-condition assumptions as the bound it refines.\n")
cat("  It is the concurrence of an out-pair marginal, a particle at p with an antiparticle at\n")
cat("  -p, and not entanglement between the two sheets. Nothing here is measured.\n")

cat("\n=== the same fact, a third time ===\n")
cat("  C = 2 sqrt(n(1-n)) is symmetric under n <-> 1-n, so the concurrence is blind to\n")
cat("  particle-hole conjugation exactly as the Renyi entropies are. Check it:\n\n")
cat("      P        n_min      n_max      C(n_min)     C(n_max)     difference\n")
for (P in c(0.05,0.5,0.95)) {
  nm <- (1-sqrt(1-P))/2; nx <- 1-nm
  cat(sprintf("  %5.2f %10.6f %10.6f %12.8f %12.8f %14.2e\n",
      P, nm, nx, 2*sqrt(nm*(1-nm)), 2*sqrt(nx*(1-nx)), abs(2*sqrt(nm*(1-nm))-2*sqrt(nx*(1-nx)))))
}
cat("\n  So one structural fact, the band's endpoints being particle-hole conjugates, now has\n")
cat("  three consequences in this paper: no spectral entropy can select the state, the\n")
cat("  decoherence clock is a min-entropy and therefore blind the same way, and the pair\n")
cat("  concurrence is blind too. They are not three coincidences.\n")

cat("\n  ==========================================================================\n")
cat("  SCOPE CORRECTION, 2026-09-21, from the pre-submission review.\n")
cat("  ==========================================================================\n")
cat("  The table above evaluates the LEAST-OCCUPIED state, and the conclusion drawn from it\n")
cat("  was extended to the family. It does not extend. At n = 1/2 the floor\n")
cat("  C >= 2|n-1/2| sqrt(P/(1-P)) is empty, and the family contains admissible states there\n")
cat("  with any concurrence. Two counterexamples at P = 0.2, both verified below to satisfy\n")
cat("  the contact condition:\n\n")
I2 <- diag(2); Zz <- diag(c(1,-1)); ff <- matrix(c(0,1,0,0),2,2,byrow=TRUE)
aa <- kronecker(ff,I2); bb <- kronecker(Zz,ff)
Npp <- (t(aa)%*%aa + t(bb)%*%bb)/2
Pc <- 0.2; sc <- sqrt(Pc); ccx <- sqrt(1-Pc)
amx <- ccx*aa + sc*t(bb); bmx <- ccx*bb - sc*t(aa)
Nmm <- (t(amx)%*%amx + t(bmx)%*%bmx)/2
psi <- c(1,0,0,1i)/sqrt(2); rhoA <- psi %*% t(Conj(psi))
rhoB <- diag(c(0.5,0,0,0.5))
for (nm in c("A: (|00>+i|11>)/sqrt2", "B: mixed, equal weights")) {
  rr <- if (substr(nm,1,1)=="A") rhoA else rhoB
  np <- Re(sum(diag(rr %*% Npp))); nm2 <- Re(sum(diag(rr %*% Nmm)))
  CC <- if (substr(nm,1,1)=="A") 2*abs(psi[1]*Conj(psi[4])) else 0
  cat(sprintf("     %-24s <N_+>=%.6f  <N_->=%.6f  |diff|=%.1e   C=%.4f\n", nm, np, nm2, abs(np-nm2), CC))
}
cat(sprintf("\n     the saturating state would give C = sqrt(P) = %.4f; these give 1 and 0.\n", sqrt(Pc)))
cat("\n  So C(x) = exp(-x^2/2) is a property of the SATURATING state, which is the state the\n")
cat("  ceiling already assumes, and not of every admissible state. The paper now says so, and\n")
cat("  the abstract no longer calls the concurrence 'fixed'.\n")
