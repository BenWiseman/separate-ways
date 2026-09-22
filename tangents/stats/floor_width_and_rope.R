# Ben's call: set the equivalence width ourselves rather than wait on a methodologist, then ask
# for thoughts on the choice rather than for the analysis. The width must come from the FLOOR's
# own propagated uncertainty, fixed before looking at the partition. Do it in that order.
#
# Normal ordering with m1 = 0:  Sigma m_nu = sqrt(dm21^2) + sqrt(dm31^2).
# Inputs: global-fit central values with symmetric 1-sigma errors. These exact values are now
# displayed in section 5.6 of the paper, so the paper and this script quote ONE input set. A
# verification pass found the paper carrying a second, asymmetric width (+0.24/-0.22 meV) for the
# same floor from a different fit. Corrected 2026-09-21: that comment claimed the second width
# was gone when it was still in the companion, so the claim was false for a day. Both documents
# now carry 58.78 +/- 0.32 meV from the inputs displayed above, under independent first-order
# propagation, and no other width appears anywhere.
dm21 <- 7.53e-5; s21 <- 0.18e-5
dm31 <- 2.510e-3; s31 <- 0.030e-3

m2 <- sqrt(dm21); m3 <- sqrt(dm31); S <- m2 + m3
cat(sprintf("   m2 = %.3f meV, m3 = %.3f meV\n", 1000*m2, 1000*m3))
cat(sprintf("   Sigma m_nu floor = %.3f meV   (paper quotes 58.8)\n", 1000*S))
stopifnot(abs(1000*S - 58.8) < 0.2)

# propagate: d(sqrt(x)) = dx/(2 sqrt(x)), added in quadrature
e2 <- s21/(2*m2); e3 <- s31/(2*m3); eS <- sqrt(e2^2+e3^2)
cat(sprintf("   contributions: from dm21 %.3f meV, from dm31 %.3f meV\n", 1000*e2, 1000*e3))
cat(sprintf("   PROPAGATED FLOOR UNCERTAINTY = %.3f meV, and dm31 supplies %.0f%% of it\n\n",
    1000*eS, 100*e3^2/(e2^2+e3^2)))

cat("=== the width, fixed here and not revisited\n\n")
w <- 2*eS
cat(sprintf("   Equivalence half-width taken as TWO propagated sigma: %.2f meV.\n", 1000*w))
cat("   Two sigma rather than one because the region should contain values a measurement\n")
cat("   could not honestly distinguish from the floor, and one sigma is too tight for that.\n")
cat("   Chosen before the partition below was computed, and not adjusted after.\n\n")

mu <- -0.036; sig <- 0.043                      # Elbers et al. profile, physical prior
Ppos <- 1 - pnorm((0-mu)/sig)
part <- function(lo,hi) c(below=(pnorm((lo-mu)/sig)-pnorm((0-mu)/sig))/Ppos,
                          within=(pnorm((hi-mu)/sig)-pnorm((lo-mu)/sig))/Ppos,
                          above=(1-pnorm((hi-mu)/sig))/Ppos)
p <- part(S-w, S+w)
cat(sprintf("   ROPE = [%.2f, %.2f] meV around the floor\n", 1000*(S-w), 1000*(S+w)))
cat(sprintf("   posterior mass BELOW  the region: %.4f\n", p["below"]))
cat(sprintf("   posterior mass WITHIN the region: %.4f\n", p["within"]))
cat(sprintf("   posterior mass ABOVE  the region: %.4f\n", p["above"]))

cat("\n=== how to read it, and what it is not\n\n")
cat(sprintf("   %.1f per cent of DESI's posterior lies below a region the floor cannot be\n", 100*p["below"]))
cat("   distinguished from. That is the honest form of the tension: not that the model is\n")
cat("   excluded, but that the data prefer a sum the model cannot reach, and the model sits\n")
cat(sprintf("   in the %.1f per cent that remains.\n", 100*(p["within"]+p["above"])))
cat("   NOT a hypothesis test. The partition is descriptive; no decision rule is attached and\n")
cat("   none should be read in. And the width is a judgement, stated so it can be disagreed\n")
cat("   with: anyone who prefers one sigma or three can recompute from the two numbers above.\n")
