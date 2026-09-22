# assignment.R uses n = (k*0.5/(p-0.5))^2, the normal approximation with sigma = 1/2.
# A statistician referee will ask for the exact binomial. Compute both and see if it matters.
approx <- function(p,k) ceiling((k*0.5/(p-0.5))^2)

# exact: smallest n such that a level-alpha one-sided binomial test of H0: p=0.5 has power >= 1/2
# at the true p. Use the conventional "k sigma" reading: alpha = one-sided normal tail at k.
exact <- function(p, k, target_power = 0.5) {
  alpha <- pnorm(k, lower.tail = FALSE)
  for (n in 1:200000) {
    crit <- qbinom(alpha, n, 0.5, lower.tail = FALSE)      # reject if X > crit
    pow  <- pbinom(crit, n, p, lower.tail = FALSE)
    if (pow >= target_power) return(n)
  }
  NA
}
cat("  Galactic-only p = 0.6987, diluted p = 0.6138\n\n")
cat(sprintf("  %-28s %8s %8s %8s\n", "case", "k", "normal", "exact"))
for (nm in c("pure Galactic","with extragalactic")) {
  p <- if (nm=="pure Galactic") 0.6987 else 0.6138
  for (k in c(3,5))
    cat(sprintf("  %-28s %8.0f %8d %8d\n", nm, k, approx(p,k), exact(p,k)))
}
cat("\n  The two agree. At 3 sigma the exact requirement is 58 against the approximation's\n")
cat("  58 (assignment.R uses p = 0.6985, where the ceiling lands on 58); at 5 sigma it is\n")
cat("  157 against 159, and the diluted pair are 177 against 174 and 483 against 483.\n")
cat("  The exact figure is sometimes above and sometimes below the approximation, because\n")
cat("  the binomial is discrete and the critical value moves in steps. Spread is under two\n")
cat("  per cent either way.\n")
cat("\n  FLATLY: the normal approximation is adequate at this precision and the paper's\n")
cat("  numbers do not move. A first version of this file concluded that the exact figure is\n")
cat("  smaller in every case, which its own table contradicts in two rows of four.\n")
