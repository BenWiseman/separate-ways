# Third and last of the paper's observational comparisons, checked for the same defect the other
# two had: a sharp prediction set against a soft measurement as though the measurement were sharp.
# D.1 predicts sigma_n <~ 1.293e-72 cm^2 and about 2.3e-30 events for LZ's exposure. This
# file carried 1.35e-72, which is the superseded normalisation; the margin below is the same
# to the digit quoted either way, but the input should be the paper's own number.

pred <- 1.293e-72         # cm^2, the paper's contact benchmark (D.1)
lz   <- 2.2e-48           # cm^2, LZ's published spin-independent sensitivity, order of magnitude
cat(sprintf("   predicted cross-section   %.2e cm^2\n", pred))
cat(sprintf("   current best sensitivity  %.2e cm^2 (order of magnitude)\n", lz))
cat(sprintf("   margin                    %.1f orders of magnitude\n\n", log10(lz/pred)))

cat("   Ask what error bar could close that margin. Each input enters the benchmark as a\n")
cat("   power, so give every one of them a generous factor and compound them:\n\n")
cat("      input                          generous uncertainty     effect on sigma_n\n")
facs <- list(c("relic abundance", 1.1), c("contact-rate normalisation", 10),
             c("scattering model", 100), c("radiation history / g_*", 3),
             c("nuclear form factor", 10), c("detector efficiency", 2))
tot <- 1
for (f in facs) { tot <- tot*as.numeric(f[[2]])
  cat(sprintf("      %-30s %10s %25s\n", f[[1]], paste0("x", f[[2]]), paste0("x", f[[2]]))) }
cat(sprintf("\n      compounded, every one taken in the same direction: x%.0e, i.e. %.1f orders\n",
    tot, log10(tot)))
cat(sprintf("      remaining margin after that: %.1f orders of magnitude\n", log10(lz/pred)-log10(tot)))

cat("\n=== verdict, and it is the opposite of the other two\n\n")
cat("  FLAT: there is no error-bar mismatch here, and the reason is that the margin is not a\n")
cat("  few per cent but about twenty-four orders of magnitude. Compounding six generous\n")
cat("  uncertainties, every one pushed the same way, moves the prediction by five orders and\n")
cat("  leaves nineteen. No plausible uncertainty bridges that, so the comparison is robust in\n")
cat("  a way the KM3NeT and DESI comparisons were not.\n")
cat("  WHAT THAT COSTS: the benchmark is therefore not a test. It cannot fail against any\n")
cat("  foreseeable experiment, so it belongs in the paper as a statement of what the model\n")
cat("  forbids rather than as a prediction awaiting data, and 5.6 should not sit it beside the\n")
cat("  endpoint and the corner as though the three were comparable.\n")
cat("  All three observational comparisons in the paper have now been checked this way: two\n")
cat("  needed the posterior stating and one needs its status downgraded.\n")
