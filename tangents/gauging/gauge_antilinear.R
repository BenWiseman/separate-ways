# Ben: what happens if you gauge the ANTILINEAR map rather than the linear one? 2.4 follows
# Harlow-Numasawa and gauges P_perp, which is linear. 2.5's limit says gauging Theta is a
# different operation we do not perform. What IS that operation, and what does it force?
#
# An antiunitary involution with Theta^2 = 1 defines a REAL STRUCTURE: the invariant vectors
# form a real Hilbert space H_R with H = H_R (x) C. So gauging Theta is not "restrict to a
# subspace of the same kind", it is "make the theory real". That is a much stronger move
# than gauging a linear Z2, and it has consequences that are checkable.

set.seed(31); N <- 12
Th <- function(v) Conj(v)                           # Theta = complex conjugation in this basis
proj <- function(v) (v + Th(v))/2                   # (1 + Theta)/2

cat("=== 1. the gauged sector is a REAL Hilbert space of half the real dimension\n\n")
V <- matrix(rnorm(N*5) + 1i*rnorm(N*5), N, 5)
P <- apply(V, 2, proj)
cat(sprintf("   arbitrary complex space: complex dim 5, real dim %d\n", 2*5))
cat(sprintf("   Theta-invariant part:    measured real rank %d\n",
    qr(rbind(Re(P), Im(P)))$rank))
cat(sprintf("   and it is genuinely real: max |Im| = %.2e\n", max(abs(Im(P)))))

cat("\n=== 2. what that FORCES: every Theta-ODD observable has zero expectation\n\n")
cat("   For A with Theta A Theta^-1 = -A, antiunitarity gives <A> = -<A>*, and <A> is real\n")
cat("   for Hermitian A, so <A> = 0 exactly. Test on random Hermitian operators of each\n")
cat("   parity, in random gauged states.\n\n")
mkH <- function() { M <- matrix(rnorm(N*N)+1i*rnorm(N*N),N,N); (M+Conj(t(M)))/2 }
odd  <- function() { H <- mkH(); (H - Conj(H))/2 }   # Theta A Theta = -A  (purely imaginary)
even <- function() { H <- mkH(); (H + Conj(H))/2 }   # Theta A Theta = +A  (purely real)
psi <- function() { v <- matrix(rnorm(N),N,1); v/sqrt(sum(v^2)) }   # real = gauged
eo <- replicate(3000, { A <- odd();  p <- psi(); Re(t(Conj(p)) %*% A %*% p) })
ee <- replicate(3000, { A <- even(); p <- psi(); Re(t(Conj(p)) %*% A %*% p) })
cat(sprintf("   Theta-ODD  observables: max |<A>| over 3000 = %.2e\n", max(abs(eo))))
cat(sprintf("   Theta-EVEN observables: max |<A>| over 3000 = %.4f  (unconstrained)\n", max(abs(ee))))
cat("\n  Exactly zero for every odd observable, in every gauged state. This is not a\n")
cat("  statistical statement; it is forced by the real structure.\n")

cat("\n=== 3. which observables are Theta-odd, and why it matters\n\n")
cat("  Theta is CPT. So the gauged theory predicts <A> = 0 for every CPT-odd A, in the\n")
cat("  state of the PAIR. Baryon number is C-odd, P-even, T-even, hence CPT-ODD.\n")
cat("  So gauging Theta forces the two-sheet universe to carry zero net baryon number,\n")
cat("  with whatever our sheet holds balanced by the mirror. That is the fold's\n")
cat("  baryogenesis story, and gauging turns it from a feature into a consequence.\n")
cat("  Energy is CPT-EVEN, so it is unconstrained and the sheets do not cancel there.\n")
cat("  That asymmetry between what cancels and what does not is the content.\n")

cat("\n=== 4. flatly, what this is NOT\n\n")
cat("  It is not a derivation that the fold MUST gauge Theta. 2.4 gauges P_perp following\n")
cat("  Harlow and Numasawa, whose argument is about spacetime inversions being linear\n")
cat("  symmetries of the gravitational path integral; extending it to an antiunitary is a\n")
cat("  separate claim we have not made.\n")
cat("  Nor is the twisted-sector question settled. Gauging a linear Z2 gives sectors\n")
cat("  classified by H^1(M; Z2), which A.3 shows vanishes on the cover. For an ANTILINEAR\n")
cat("  gauging the relevant classification is not the same cohomology and A.3 does not\n")
cat("  cover it. That is a genuine open fork, not a gap we can close by citing A.3.\n")
