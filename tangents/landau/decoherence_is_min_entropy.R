# The decoherence exponent has a name, and the paper does not use it.
#
# 4.2 integrates -log max(n, 1-n) over modes to get R, the normalisation of the decoherence time.
# That function is exactly the order-infinity Renyi entropy of the pair spectrum {n, 1-n}, i.e.
# the MIN-ENTROPY H_inf. Check the convergence, then note two things it buys for free.
Sr <- function(n,a) if (is.infinite(a)) -log(max(n,1-n)) else log(n^a+(1-n)^a)/(1-a)
cat("      n      alpha=10     alpha=50    alpha=200       alpha=Inf   -log max(n,1-n)\n")
ns <- c(0.05,0.2,0.35,0.5,0.7,0.95)
for (n in ns)
  cat(sprintf("  %5.2f %12.7f %12.7f %12.7f %15.7f %17.7f\n",
      n, Sr(n,10), Sr(n,50), Sr(n,200), Sr(n,Inf), -log(max(n,1-n))))
cat(sprintf("\n  max |S_200 - S_Inf| over this range: %.2e ; the family converges from above.\n",
    max(sapply(ns, function(n) abs(Sr(n,200)-Sr(n,Inf))))))

cat("\n  TWO THINGS THIS BUYS.\n")
cat("  1. The min-entropy is the one-shot quantity: -log of the largest outcome probability, the\n")
cat("     guessing entropy. So R is not an average uncertainty over many trials, it is the\n")
cat("     single-shot distinguishability of the two branches, summed over modes. That is the\n")
cat("     right object for a decoherence threshold and it was already being used as one.\n")
cat("  2. It explains the particle-hole blindness found separately. H_inf depends on the\n")
cat("     spectrum {n, 1-n} only through its maximum, so it cannot tell n from 1-n either:\n")
cat(sprintf("       n = 0.05 gives %.7f and n = 0.95 gives %.7f, difference %.1e\n",
    Sr(0.05,Inf), Sr(0.95,Inf), abs(Sr(0.05,Inf)-Sr(0.95,Inf))))
cat("     The same degeneracy that kills entropic state selection also fixes the decoherence\n")
cat("     clock's functional form. One fact, two consequences.\n")
cat("\n  WHAT IT IS NOT: a new result. It renames a function the paper already integrates, and the\n")
cat("  numbers do not move. It is worth one sentence because it tells a reader which quantity\n")
cat("  they are looking at and connects two sections that currently look unrelated.\n")
