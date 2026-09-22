# The companion computes the two-branch overlap by DELETING odd-parity oscillators, on the
# grounds that gauging P_perp restricts to its invariant sector. A verification pass says
# projection deletes nothing, because the Gaussian branch states are already invariant.
# Check the three steps.
a <- 1.3; b <- 0.8                       # arbitrary; nothing below depends on the values
psi <- function(q, sgn) (a/pi)^(1/4) * exp(-(a + sgn*1i*b)*q^2/2)

cat("=== 1. is the branch state even in q, hence P_perp-invariant on an odd mode? ===\n")
qs <- c(0.3, 1.1, 2.7)
cat("      q      |psi(+q) - psi(-q)|\n")
for (q in qs) cat(sprintf("  %6.2f   %18.3e\n", q, Mod(psi(q,+1) - psi(-q,+1))))
cat("  A Gaussian in q^2 cannot tell +q from -q. The projector (1 + U(P_perp))/2 therefore acts\n")
cat("  as the identity on it, whatever the parity of the spatial harmonic it belongs to.\n")

cat("\n=== 2. what each mode actually contributes to the overlap ===\n")
num <- integrate(function(q) Re(Conj(psi(q,+1))*psi(q,-1)), -Inf, Inf)$value
num <- num + 1i*integrate(function(q) Im(Conj(psi(q,+1))*psi(q,-1)), -Inf, Inf)$value
closed <- (1 + b^2/a^2)^(-1/4)
cat(sprintf("  numerical |<psi+|psi->| = %.10f\n  closed form (1+b^2/a^2)^(-1/4) = %.10f\n",
            Mod(num), closed))
cat(sprintf("  projected state is identical, so the projected overlap is the same: %.10f\n", closed))
cat("  Deleting the mode would replace this by 1, which is strictly larger.\n")

cat("\n=== 3. the effect on the threshold ===\n")
cat("  Every deleted mode multiplies the overlap by 1 instead of by a number below 1, so the\n")
cat("  overlap decays more slowly and the threshold moves LATER. That is the whole of the\n")
cat("  claimed 19 per cent delay: it is the deletion, not the gauging.\n")
r <- closed
for (k in c(1, 3, 10)) cat(sprintf("   %2d deleted modes: overlap ratio kept/deleted = %.6f\n", k, r^k))
cat("\n  VERDICT: gauging physical states does not delete these oscillators. Recovering the\n")
cat("  deleted-mode sum needs an additional restriction on the FIELD, not on the states, and\n")
cat("  the companion does not derive one. The ungauged threshold is what the stated projection\n")
cat("  gives; the 19 per cent delay is withdrawn.\n")
