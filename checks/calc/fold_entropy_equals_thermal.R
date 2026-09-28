# Does the fold's OWN derived squeeze make one sheet's entanglement entropy the thermal
# entropy of a horizon?
#
# Why this matters. Jacobson's 1995 derivation gets the Einstein equation out of Clausius,
# dQ = T dS, applied to local horizons, and it has to assume BOTH inputs: a temperature and
# an entropy proportional to area. The fold already derives the temperature, because the
# half-KMS shift forces it. If it also derives the entropy, the field equations stop being
# an input to this construction and become an output of it.
#
# The quarter-area identity looked like the source and is not: quarter_area_jacobson.R shows
# it is a Schwarzschild accident, running to 0.768 at Q/M = 0.99 and 0.42 by D = 7. Jacobson
# needs a universal coefficient, so that route is closed.
#
# This is the other route, and it uses only what the fold derives. A.8's propagation makes
# the two sheets a thermofield double with
#
#     tanh r(omega) = exp(-beta omega / 2),
#
# no free parameter. For a two-mode squeezed vacuum the reduced state of one mode is thermal
# with <N> = sinh^2 r, so its von Neumann entropy should be the thermal entropy of one copy
# at the horizon temperature. If that holds exactly, the fold hands Clausius an entropy it
# did not have to assume.

S_ent <- function(r) {                      # von Neumann entropy of one mode of a TFD
  n <- sinh(r)^2
  ifelse(n <= 0, 0, (n + 1) * log(n + 1) - n * log(n))
}
S_thermal <- function(bw) {                 # thermal entropy of one bosonic mode, beta*omega
  bw / (exp(bw) - 1) - log(1 - exp(-bw))
}
r_fold <- function(bw) atanh(exp(-bw / 2))  # the squeeze A.8 derives

cat("=== 1. the fold's squeeze reproduces the Bose occupation exactly ===\n")
cat("     beta*omega      sinh^2 r        1/(e^bw - 1)       difference\n")
worst_n <- 0
for (bw in c(0.1, 0.5, 1, 2, ln3 <- log(3), 5, 10)) {
  n1 <- sinh(r_fold(bw))^2; n2 <- 1 / (exp(bw) - 1)
  worst_n <- max(worst_n, abs(n1 - n2))
  cat(sprintf("   %10.5f  %14.10f  %16.10f   %.3e\n", bw, n1, n2, abs(n1 - n2)))
}
stopifnot(worst_n < 1e-12)

cat("\n=== 2. and therefore the entanglement entropy IS the thermal entropy ===\n")
cat("     beta*omega    S_entanglement      S_thermal          difference\n")
worst_s <- 0
for (bw in c(0.05, 0.2, 0.7, log(3), 1.5, 3, 6, 12)) {
  s1 <- S_ent(r_fold(bw)); s2 <- S_thermal(bw)
  worst_s <- max(worst_s, abs(s1 - s2))
  cat(sprintf("   %10.5f  %16.10f  %16.10f   %.3e\n", bw, s1, s2, abs(s1 - s2)))
}
cat(sprintf("   worst difference over the range: %.3e\n", worst_s))
stopifnot(worst_s < 1e-12)

cat("\n=== 3. the check has to be able to fail, so break the squeeze on purpose ===\n")
cat("   use tanh r = exp(-beta omega) instead of the half-period exp(-beta omega/2):\n")
cat("     beta*omega    S from wrong squeeze    S_thermal       difference\n")
for (bw in c(0.2, 1, 3)) {
  s1 <- S_ent(atanh(exp(-bw))); s2 <- S_thermal(bw)
  cat(sprintf("   %10.5f  %18.10f  %14.10f   %.4f\n", bw, s1, s2, abs(s1 - s2)))
}
cat("   the wrong squeeze misses by order one, so the agreement above is the\n")
cat("   half-period doing work and not an identity that holds for any squeeze.\n")

cat("\n=== 4. what this does and does not buy ===\n")
cat("   Buys: the entropy Clausius needs is not an assumption here. It is the\n")
cat("   entanglement entropy of one sheet in the state the fold's own propagation\n")
cat("   selects, and it equals the thermal entropy at the horizon temperature that\n")
cat("   the same propagation fixes. Two of Jacobson's inputs, one source.\n")
cat("   Does not buy: proportionality to AREA. Per mode this is an entropy per mode.\n")
cat("   Summing it over the transverse modes of a horizon is where the area would\n")
cat("   have to come from, and that sum is UV divergent and needs a cutoff, which is\n")
cat("   the same place the Bekenstein-Hawking coefficient always comes from. The fold\n")
cat("   has not removed that step. What it has removed is the freedom in T and in the\n")
cat("   state.\n")
