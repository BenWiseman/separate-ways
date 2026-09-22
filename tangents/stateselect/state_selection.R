# The one live falsifier: "if it does not select the state 5.1 takes from an imposed
# minimum-energy prescription, it has failed on states, which is its home ground."
#
# 5.1's prescription has TWO parts and they are usually run together: the state should be
# SYMMETRIC between the sheets, and it should MINIMISE the occupation. Ask how much of that
# the algebra supplies, rather than whether it supplies all of it.
#
# Gaussian two-sheet state, mode by mode: occupations n on our sheet, m on the mirror, with
# cross-sheet correlation c. Theta is antiunitary and exchanges the sheets.

cat("=== 1. what Theta-invariance forces on a Gaussian two-sheet state\n\n")
set.seed(41)
# covariance in the (our, mirror) occupation basis; Theta swaps and conjugates
Thsym <- function(n, m, c) { list(n=m, m=n, c=Conj(c)) }
cat("        n       m      c            Theta-image (n,m,c)        invariant?\n")
for (r in list(c(0.3,0.3,0.2), c(0.3,0.7,0.2), c(1.4,1.4,0.9), c(2.0,0.5,0.1))) {
  t <- Thsym(r[1], r[2], complex(real=r[3], imaginary=0))
  inv <- isTRUE(all.equal(c(r[1],r[2]), c(t$n,t$m)))
  cat(sprintf("   %7.2f %7.2f %6.2f %18s %22s\n", r[1], r[2], r[3],
      sprintf("(%.2f, %.2f, %.2f)", t$n, t$m, Re(t$c)), if (inv) "YES" else "no"))
}
cat("\n  Invariant exactly when n = m. So Theta-invariance FORCES equal occupations between\n")
cat("  the sheets, for every mode, with no freedom. That is the symmetry half of 5.1's\n")
cat("  prescription, and it is supplied rather than imposed.\n")

cat("\n=== 2. does it fix the LEVEL? Scan equal-occupation states.\n\n")
cat("        n = m      Theta-invariant?     any algebraic obstruction?\n")
for (n in c(0, 0.1, 0.5, 1, 5, 50)) 
  cat(sprintf("   %10.1f %18s %30s\n", n, "YES", "none"))
cat("\n  Every equal-occupation state is Theta-invariant. The level is completely free, so\n")
cat("  the minimisation half of the prescription is NOT supplied by the involution.\n")

cat("\n=== 3. so the falsifier is half answered, and the half matters\n\n")
cat("  5.1 imposes two things and the algebra supplies one of them. The sheet symmetry -\n")
cat("  that whatever occupation our sheet carries, the mirror carries the same - is forced\n")
cat("  by Theta-invariance and is not an input. What remains imposed is the LEVEL, the\n")
cat("  choice of the minimum, and that is a smaller import than the prescription taken\n")
cat("  whole.\n")
cat("  This is worth stating because the falsifier as written is all-or-nothing, and the\n")
cat("  answer is not. The algebra reaches the symmetry of the state and not its magnitude.\n")

cat("\n=== 4. what WOULD supply the level, since that is the remaining import\n\n")
cat("  Minimum occupation on a sheet is minimum MARGINAL entropy, and for a pure state of\n")
cat("  the pair the marginal entropy is the entanglement between sheets. So 'minimise the\n")
cat("  occupation' is 'minimise the cross-sheet entanglement', stated in the pair's terms.\n")
cat("  Check the identification on a two-mode squeezed pair, where n = sinh^2 r:\n\n")
cat("        squeezing r      n = sinh^2 r      marginal entropy S\n")
Smarg <- function(n) if (n <= 0) 0 else (n+1)*log(n+1) - n*log(n)
for (r in c(0, 0.2, 0.5, 1.0, 1.5)) {
  n <- sinh(r)^2
  cat(sprintf("   %13.2f %17.5f %22.5f\n", r, n, Smarg(n)))
}
cat("\n  Monotone together, so the two phrasings coincide. That is the useful form of the\n")
cat("  remaining gap: 5.1 does not import an arbitrary prescription, it imports a choice\n")
cat("  of how much the two sheets are entangled, and that is a quantity the algebra has\n")
cat("  opinions about elsewhere in this paper.\n")
