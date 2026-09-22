# A hostile referee's kill shot: the bang-adiabatic state is itself Theta-invariant, sitting at
# mu = pi in A.18's own family, so the fold cannot be what separates it from the adopted state
# and the abstract's "what derives it is the difference between one sheet and two" is false.
# Adjudicate it rather than argue. Three questions, each answerable by arithmetic.

gamma <- 1
cat("=== 1. is the bang-adiabatic state in the Theta-invariant family?\n\n")
p <- 0.6
E <- sqrt((gamma*0)^2 + p^2)
v <- c(p+0i, -gamma*0 - E); v <- v/sqrt(sum(Mod(v)^2))
cat(sprintf("   lower adiabatic state at eta = 0, p = %.1f : (%.4f, %.4f)\n", p, Re(v[1]), Re(v[2])))
cat(sprintf("   |psi_1| = %.6f, |psi_2| = %.6f, both 1/sqrt2 = %.6f\n",
    Mod(v[1]), Mod(v[2]), 1/sqrt(2)))
mu <- Arg(v[2]/v[1]) %% (2*pi)
cat(sprintf("   relative phase mu = %.4f, and pi = %.4f\n\n", mu, pi))
cat("   ANSWER: YES. The referee is right, and A.18 already said so in its own words:\n")
cat("   'The bang-adiabatic state sits at mu = pi and is not the minimum.' So Theta-invariance\n")
cat("   does NOT separate the two states. Anything claiming it does is wrong.\n")

cat("\n=== 2. then what does separate them?\n\n")
nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
cat("   Occupation tails: mu = pi gives n ~ gamma^2/16p^4, mu_*(p) gives n ~ exp(-x^2)/4.\n")
cat("   The separating quantity is the ENERGY integral, int x^2 * x * n dx:\n\n")
cat("      state            tail        int_1^X x^3 n dx at X = 10, 100, 1000\n")
f_pi <- function(x) x^3 * (1/(16*x^4))
f_mn <- function(x) x^3 * nmin(x)
for (X in c(10,100,1000)) {
  a <- integrate(f_pi, 1, X)$value; b <- integrate(f_mn, 1, X, subdivisions=2000)$value
  cat(sprintf("      mu = pi          p^-4     %12.4f\n", a))
  cat(sprintf("      mu = mu_*(p)     Gaussian %12.4e   (X = %d)\n", b, X))
}
cat("\n   The first grows without bound, logarithmically, and the second is converged by x ~ 3.\n")
cat("   So the separation is by HADAMARD REGULARITY through the energy integral, exactly as\n")
cat("   5.1 says when it excludes the bang-adiabatic state. Not by the fold.\n")

cat("\n=== 3. so is the one-sheet / two-sheet framing empty?\n\n")
cat("   No, but it is weaker than 'derives'. The adopted state (1, e^{i mu_*})/sqrt2 is a\n")
cat("   perfectly good initial condition on a half-line too: nothing forbids a one-sheet\n")
cat("   universe from positing it. What a one-sheet universe LACKS is a reason to. Its natural\n")
cat("   choice is the adiabatic vacuum at the bang, which is mu = pi, and that one fails on\n")
cat("   energy. The fold supplies the magnitude condition |psi_1| = |psi_2| that makes the\n")
cat("   family the right object to minimise over; it does not by itself pick the member.\n\n")
cat("   HONEST CLAIM: the fold SUPPLIES the condition, and Hadamard regularity then selects\n")
cat("   within it. 'What derives it is the difference between one sheet and two' overstates\n")
cat("   that and must be corrected. The referee's shot lands on the wording, not on the result.\n")
