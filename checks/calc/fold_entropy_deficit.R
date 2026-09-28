# If the bang ceiling binds at a horizon, the fold predicts LESS entropy than thermal,
# by a pure number.
#
# The chain. A.8's propagation fixes the horizon squeeze at tanh r = exp(-beta omega / 2),
# with no free parameter. The crossing at the bang caps the squeeze at r <= arcsinh(1/sqrt2),
# because Theta-invariance forces half-and-half occupancy and the pair-block floor tops out
# at n_* = 1/2. Those two agree at beta omega = ln 3, which the paper reports as a
# coincidence nothing was arranged to produce.
#
# Read it as a constraint instead of a coincidence and it says something sharp: a horizon
# mode with beta omega < ln 3 would need a squeeze the crossing cannot supply. Such modes
# cannot be in the fold-selected state. The horizon's entropy is then the thermal entropy
# with its soft modes removed, and the removed fraction is a pure number, independent of
# the ultraviolet cutoff that the area law always needs.
#
# This is a deviation from Bekenstein-Hawking, so it had better be computed rather than
# asserted, and it had better be attacked.

S_th <- function(x) ifelse(x <= 0, 0, x / (exp(x) - 1) - log(1 - exp(-x)))   # x = beta*omega
LN3  <- log(3)

cat("=== 1. the ceiling, and where it bites ===\n")
r_ceiling <- asinh(1 / sqrt(2))
bw_cross  <- -2 * log(tanh(r_ceiling))
cat(sprintf("   bang ceiling r_max          = %.12f\n", r_ceiling))
cat(sprintf("   horizon squeeze meets it at beta*omega = %.12f\n", bw_cross))
cat(sprintf("   ln 3                        = %.12f   difference %.2e\n",
            LN3, abs(bw_cross - LN3)))
stopifnot(abs(bw_cross - LN3) < 1e-12)

cat("\n=== 2. the removed fraction, for a massless field in three spatial dimensions ===\n")
# entropy density in frequency: rho(omega) ~ omega^2 for a 3D gas of modes
num <- integrate(function(x) x^2 * S_th(x), 0, LN3,  rel.tol = 1e-12)$value
den <- integrate(function(x) x^2 * S_th(x), 0, Inf, rel.tol = 1e-12)$value
cat(sprintf("   entropy below ln 3 : %.10f\n", num))
cat(sprintf("   entropy, all modes : %.10f\n", den))
cat(sprintf("   REMOVED FRACTION   : %.10f   (%.4f per cent)\n", num / den, 100 * num / den))

cat("\n=== 3. is it an artefact of the mode measure? sweep the density of states ===\n")
cat("      rho(omega) ~ omega^p     removed fraction\n")
for (p in 0:4) {
  a <- integrate(function(x) x^p * S_th(x), 0, LN3,  rel.tol = 1e-11)$value
  b <- integrate(function(x) x^p * S_th(x), 0, Inf, rel.tol = 1e-11)$value
  cat(sprintf("        p = %d                 %.8f   (%.3f%%)\n", p, a / b, 100 * a / b))
}
cat("   the fraction falls as the measure weights high frequencies more, which is what a\n")
cat("   low-frequency cut must do. It is not a fixed number: it depends on p, so the\n")
cat("   prediction is only as sharp as the horizon's density of states.\n")

cat("\n=== 4. attacking it: is the cut at ln 3 even the right reading? ===\n")
cat("   The ceiling is a statement about the state the CROSSING produces at the bang.\n")
cat("   The horizon squeeze is a statement about equilibrium at a horizon. For the cut\n")
cat("   to be real, the bang ceiling has to bound the horizon state, which needs the\n")
cat("   horizon state to be the evolved bang state and nothing else to have happened in\n")
cat("   between. A.8 says exactly that the state at a horizon is whatever the selected\n")
cat("   bang state evolves into, so the reading is the paper's own, but equilibration\n")
cat("   with local matter is not modelled anywhere and would refill the soft modes.\n")
cat("   So: the cut is what the construction says IF the horizon has not equilibrated\n")
cat("   with anything since the bang. For a cosmological horizon that is arguable. For\n")
cat("   a black hole that formed from a star it is not.\n")

cat("\n=== 5. the same arithmetic, said the other way round ===\n")
s_at_cross <- S_th(LN3)
cat(sprintf("   entropy of one mode exactly at the crossing: %.10f\n", s_at_cross))
cat(sprintf("   occupation there: sinh^2(r_max) = %.10f, exactly one half\n",
            sinh(r_ceiling)^2))
cat("   So ln 3 is the frequency at which a horizon mode is half occupied, and the\n")
cat("   bang ceiling is the statement that half occupancy is the most the crossing can\n")
cat("   make. Those are the same sentence, which is why the two curves meet. That is a\n")
cat("   better description of the coincidence than 'nothing was arranged to produce it'.\n")
