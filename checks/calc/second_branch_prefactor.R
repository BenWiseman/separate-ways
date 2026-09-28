#!/usr/bin/env Rscript
# second_branch_prefactor.R -- referee item 2: which prefactor carries the OUTGOING horizon
# branch, settled from asymptotics alone so the next attempt at the second continued fraction
# starts from the right ansatz.
#
# A.14's argument needs a reflecting horizon, psi ~ e^{-i omega r*} + R e^{+i omega r*}, so it
# needs BOTH Frobenius solutions at the horizon. leaver_qnm.R has the ingoing one. The
# outgoing one has never been built, and REFEREE_OPEN records four attempts that failed by
# order one on the Regge-Wheeler residual.
#
# All four kept Leaver's coefficients a_n and varied the prefactor. That cannot work: the two
# Frobenius solutions at a regular singular point have DIFFERENT recursions, and substituting a
# different prefactor into the same series is not the second solution. What a prefactor test can
# settle, and what this file settles, is which prefactor the second recursion must be built
# under. That is forced, and it is not among the four tried.
#
# Units 2M = 1, rho = -i omega, tortoise r* = r + log(r-1).
#   outgoing at infinity: psi ~ e^{+i omega r*} = e^{-rho r*} = r^{-rho} e^{-rho r}
#   outgoing at horizon:  psi ~ (r-1)^{-rho}
# A trial (r-1)^A r^B e^{-rho(r-1)} has r-power A + B as r -> inf, so infinity forces A + B =
# -rho. Leaver takes A = +rho and therefore needs B = -2rho. The outgoing branch takes A = -rho
# and therefore needs B = 0: no r power at all.

om <- 0.4 + 0i; rho <- -1i*om
rstar <- function(r) r + log(r - 1)
P <- function(r, A, B) (r-1)^A * r^B * exp(-rho*(r-1))
far <- c(1e3, 1e5, 1e7)

cat("=== the outgoing-at-both-ends prefactor is forced by the two asymptotics ===\n\n")
cat("   trial                            psi / e^{-rho r*} at r = 1e3, 1e5, 1e7\n")
trials <- list(list(-rho, 0,       "(r-1)^-rho r^0      [outgoing branch]"),
               list(-rho, -2*rho,  "(r-1)^-rho r^-2rho"),
               list(-rho,  2*rho,  "(r-1)^-rho r^+2rho"),
               list( rho, -2*rho,  "(r-1)^+rho r^-2rho  [Leaver, ingoing]"))
spread <- c()
for (t in trials) {
  v <- sapply(far, function(r) P(r, t[[1]], t[[2]]) / exp(-rho*rstar(r)))
  spread <- c(spread, max(abs(v - v[length(v)])))
  cat(sprintf("   %-38s %s\n", t[[3]], paste(sprintf("%.4f%+.4fi", Re(v), Im(v)), collapse="  ")))
}
cat(sprintf("\n   drift across six decades of r: %s\n",
            paste(sprintf("%.1e", spread), collapse = "  ")))
# the two that satisfy A + B = -rho are flat; the other two wander by order one
# 1e-8, not 1e-12: at r = 1e7 a complex power carries that much rounding on its own
stopifnot(spread[1] < 1e-8, spread[4] < 1e-3, spread[2] > 0.5, spread[3] > 0.5)

near <- c(1+1e-4, 1+1e-6, 1+1e-8)
v <- sapply(near, function(r) P(r, -rho, 0) / (r-1)^(-rho))
cat(sprintf("\n   at the horizon, psi/(r-1)^{-rho} at r-1 = 1e-4, 1e-6, 1e-8:\n     %s\n",
            paste(sprintf("%.8f%+.8fi", Re(v), Im(v)), collapse = "  ")))
stopifnot(max(abs(v - 1)) < 1e-3)

# The check must be able to fail: a prefactor with the wrong horizon exponent does not tend to
# the outgoing behaviour there, whatever it does at infinity.
# For real omega these exponents are pure phase, so the control has to look at the phase and
# not the modulus: |(r-1)^{-i omega}| is 1 whichever sign it carries. Leaver's prefactor divided
# by the outgoing exponent leaves (r-1)^{2 rho}, which winds without bound as r -> 1.
w <- sapply(near, function(r) P(r, rho, -2*rho) / (r-1)^(-rho))
cat(sprintf("   the same ratio for Leaver's prefactor, which is INgoing there:\n     %s\n",
            paste(sprintf("%.5f%+.5fi", Re(w), Im(w)), collapse = "  ")))
cat(sprintf("   its phase moves %.2f radians over those two decades; the outgoing one moves %.1e\n",
            max(abs(diff(Arg(w)))), max(abs(diff(Arg(v))))))
stopifnot(max(abs(diff(Arg(w)))) > 1, max(abs(diff(Arg(v)))) < 1e-3)

cat("\n=== flatly ===\n\n")
cat("  The second branch's ansatz is (r-1)^{-rho} e^{-rho(r-1)} times a series in (r-1)/r,\n")
cat("  with NO power of r, and that is forced rather than chosen. What is still not built is\n")
cat("  the three-term recursion for that series, which has to come from Leaver's confluent\n")
cat("  reduction under this prefactor and not from his coefficients under another one.\n")
