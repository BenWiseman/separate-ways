# ==========================================================================================
# CORRECTED 2026-09-21. The first version of this file asserted that charge conjugation IS
# particle-hole conjugation on a fermionic pair block, and concluded that the abundance is the
# unique CHARGE-odd observable. That identification is wrong and a pre-submission review caught
# it. Charge conjugation exchanges the two members of a pair, a and b. On this block that map
# preserves N = (n_a + n_b)/2 EXACTLY and reverses n_a - n_b, so the abundance is charge-EVEN
# and the charge is the odd quantity. Verified: ||S N S^T - N|| = 0 and ||S Q S^T + Q|| = 0
# for S the a <-> b swap.
#
# The mathematics below is unaffected and stands: n -> 1-n is PARTICLE-HOLE conjugation, it is
# the map relating the two ends of the band, and the abundance is the unique functional here
# with an odd part under it. Only the charge interpretation was wrong, and the file is renamed
# accordingly.
#
# Decompose every functional the construction uses into even and odd parts under that map
# and see which ones survive. Do not assert the answer; compute the odd part.
# ==========================================================================================
odd  <- function(F, n) (F(n) - F(1-n))/2
even <- function(F, n) (F(n) + F(1-n))/2

# --- the functionals the paper actually uses -------------------------------------------
renyi <- function(a) function(n) {
  p <- cbind(n, 1-n)
  if (abs(a-1) < 1e-12) -rowSums(p*log(p)) else log(rowSums(p^a))/(1-a)
}
minent <- function(n) -log(pmax(n, 1-n))          # order-infinity Renyi: 4.2's exponent
conc   <- function(n) 2*sqrt(n*(1-n))             # concurrence, 3.1
purity <- function(n) n^2 + (1-n)^2
number <- function(n) n                           # the occupation itself
energy <- function(n) n                           # energy is linear in occupation

ns <- seq(0.02, 0.98, by = 0.005)
cat("  Odd part under n -> 1-n, maximised over the band n in [0.02, 0.98]:\n\n")
cat(sprintf("   %-34s %14s %14s\n", "functional", "max |odd|", "max |even|"))
tests <- list(
  "Renyi alpha = 0.5"      = renyi(0.5),
  "Renyi alpha = 1 (von Neumann)" = renyi(1),
  "Renyi alpha = 2"        = renyi(2),
  "Renyi alpha = 200"      = renyi(200),
  "min-entropy (4.2's exponent)" = minent,
  "concurrence C = 2 sqrt(n(1-n))" = conc,
  "purity"                 = purity,
  "number density (the abundance)" = number
)
for (nm in names(tests)) {
  F <- tests[[nm]]
  cat(sprintf("   %-34s %14.3e %14.6f\n", nm,
              max(abs(odd(F, ns))), max(abs(even(F, ns)))))
}

cat("\n  Every entanglement and entropy functional in the construction has an odd part at\n")
cat("  the level of machine epsilon. The occupation does not: its odd part is exactly\n")
cat("  n - 1/2, which is the whole of the signal.\n\n")
cat(sprintf("   odd part of the number density at n = 0.05 : %+.6f  (exactly n - 1/2 = %+.6f)\n",
            odd(number, 0.05), 0.05-0.5))
cat(sprintf("   odd part of the number density at n = 0.95 : %+.6f  (exactly n - 1/2 = %+.6f)\n",
            odd(number, 0.95), 0.95-0.5))

# --- so how much of the measured abundance is the charge-odd part? ----------------------
# The occupation appears under an x^2 weight. Split THAT integral (this is
# Int x^2 n dx = 0.125933, not the paper's production integral I_0 = 0.0127597,
# which carries its own extra weight; the decomposition below holds for any
# positive weight, so which one is used does not matter to the conclusion).
P  <- function(x) exp(-x^2)
ns_ <- function(x) (1 - sqrt(1-P(x)))/2
I_tot  <- integrate(function(x) x^2*ns_(x),           0, 40, rel.tol=1e-12)$value
I_odd  <- integrate(function(x) x^2*(ns_(x)-0.5),     0, 40, rel.tol=1e-12)$value
I_even <- integrate(function(x) x^2*0.5,              0, 40, rel.tol=1e-12)$value
cat(sprintf("\n  Int x^2 n dx over x in [0,40] = %.10f\n", I_tot))
cat(sprintf("   charge-even piece (the 1/2)         : %.6f   <- diverges with the cutoff\n", I_even))
cat(sprintf("   charge-odd piece  (n - 1/2)         : %.6f   <- and cancels it exactly\n", I_odd))
cat(sprintf("   sum                                 : %.10f\n", I_even + I_odd))

for (X in c(10, 20, 40, 80)) {
  e <- X^3/6; o <- integrate(function(x) x^2*(ns_(x)-0.5), 0, X, rel.tol=1e-12)$value
  cat(sprintf("   cutoff x < %3d : even %12.2f   odd %12.2f   sum %.10f\n", X, e, o, e+o))
}
cat("\n  FLATLY: the charge-even half is the vacuum 1/2 per mode and diverges with the\n")
cat("  cutoff as X^3/6. The charge-odd half diverges against it, as -X^3/6 offset by the\n")
cat("  answer. Only the SUM is cutoff-independent, so the abundance is a cancellation\n")
cat("  between two divergent halves and NOT the charge-odd half by itself -- the table\n")
cat("  above says so, odd = -10666.54 at X = 40 against an abundance of 0.126. What is\n")
cat("  true is narrower: n - 1/2 carries all of the state-dependence, the even half\n")
cat("  carrying none. An earlier version of this line said the abundance IS the odd part,\n")
cat("  which its own table contradicts; 4.3 of the paper states the narrower claim.\n")
cat("\n  So the one measurement that fixes the dark-matter mass is also the only one in the\n")
cat("  construction that can tell the two sheets apart. Every entanglement measure is blind\n")
cat("  to the fold's charge reversal by an exact algebraic identity, not by an accident of\n")
cat("  the state: C, the Renyi family and the min-entropy all see the spectrum {n, 1-n} as\n")
cat("  a set, and charge conjugation only swaps its two members.\n")
