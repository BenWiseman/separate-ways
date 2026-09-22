# A review points out Elbers et al. publish a profile approximation for the DR2 BAO+CMB
# neutrino-mass posterior: mu0 = -0.036 eV, sigma = 0.043 eV, Gaussian and allowed negative.
# If that is right we do NOT need their chain, and the paper's shape-dependent 6.4-8.4 per cent
# range can be replaced by a number from the survey's own profile.
# The test of whether we may use it: does it reproduce their PUBLISHED 95 per cent limit?

mu <- -0.036; sig <- 0.043; floor <- 0.0588
cat(sprintf("   profile: mu0 = %.3f eV, sigma = %.3f eV\n", mu, sig))

# physical prior: truncate at zero and renormalise
P_pos <- 1 - pnorm((0 - mu)/sig)
cat(sprintf("   P(Sigma m_nu > 0) before renormalising = %.4f\n", P_pos))
L <- uniroot(function(L) (pnorm((L-mu)/sig) - pnorm((0-mu)/sig))/P_pos - 0.95, c(0.01, 0.5))$root
cat(sprintf("   implied 95%% upper limit = %.1f meV\n", 1000*L))
cat(sprintf("   Elbers et al. published   = 64.2 meV      difference %.1f meV\n\n", 1000*L-64.2))

tail <- (1 - pnorm((floor - mu)/sig))/P_pos
cat(sprintf("   P(Sigma m_nu > 58.8 meV | physical) = %.4f, i.e. %.1f per cent\n", tail, 100*tail))
cat(sprintf("   the paper currently quotes a shape-dependent range of 6.4 to 8.4 per cent\n"))
cat(sprintf("   the profile value sits %s that range\n\n", if (tail>0.064 && tail<0.084) "INSIDE" else "OUTSIDE"))

cat("=== the three-way partition the comparison actually wants\n\n")
for (w in c(0.002, 0.005, 0.010)) {
  lo <- floor - w; hi <- floor + w
  below  <- (pnorm((lo-mu)/sig) - pnorm((0-mu)/sig))/P_pos
  within <- (pnorm((hi-mu)/sig) - pnorm((lo-mu)/sig))/P_pos
  above  <- (1 - pnorm((hi-mu)/sig))/P_pos
  cat(sprintf("   equivalence half-width %4.1f meV : below %.3f  within %.3f  above %.3f\n",
      1000*w, below, within, above))
}
cat("\n   The width must be fixed in advance by the floor's own propagated uncertainty, not\n")
cat("   chosen after seeing these. The rows are shown to make the sensitivity visible.\n")

cat("\n=== verdict\n\n")
cat("   The profile reproduces the published limit to well under a meV, so it is a faithful\n")
cat("   summary of their posterior and the chain is not needed for this diagnostic. The\n")
cat("   partition is therefore computable now, from published numbers, and the paper can\n")
cat("   replace a shape-dependent range with the survey's own profile.\n")
