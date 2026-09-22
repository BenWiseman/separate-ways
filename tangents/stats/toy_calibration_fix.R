# The paper quotes a 6.4-8.4 per cent range for the posterior mass above the 58.8 meV floor,
# from three toy shapes all said to be calibrated so that 64.2 meV is the 95th percentile.
# The uniform toy is not calibrated that way: it ENDS at 64.2, which makes 64.2 its 100th
# percentile. Recalibrate it the way the other two are and see what the range becomes.
L <- 64.2; F <- 58.8
cat(sprintf("  DESI 95%% upper limit L = %.1f meV, floor F = %.1f meV\n\n", L, F))

# exponential: 1 - exp(-lam L) = 0.95
lam <- log(20)/L
p_exp <- exp(-lam*F)
# half-Gaussian: erf(L/(s sqrt2)) = 0.95
erf  <- function(x) 2*pnorm(x*sqrt(2)) - 1
erfi <- function(p) qnorm((p+1)/2)/sqrt(2)
s <- L/(erfi(0.95)*sqrt(2))
p_hg <- 1 - erf(F/(s*sqrt(2)))
# uniform, as the paper has it: support [0, L], so L is the MAXIMUM
p_uni_bad <- (L-F)/L
# uniform, calibrated like the others: 95th percentile at L, so support [0, L/0.95]
W <- L/0.95
p_uni_ok <- (W-F)/W

cat(sprintf("  exponential   (95th pct at L) : %8.4f %%\n", 100*p_exp))
cat(sprintf("  half-Gaussian (95th pct at L) : %8.4f %%\n", 100*p_hg))
cat(sprintf("  uniform, support [0, L]       : %8.4f %%   <- L is its 100th percentile\n", 100*p_uni_bad))
cat(sprintf("  uniform, 95th pct at L        : %8.4f %%   <- calibrated like the other two\n", 100*p_uni_ok))
cat(sprintf("\n  paper's quoted range  : %.1f to %.1f per cent\n", 100*p_exp, 100*p_uni_bad))
cat(sprintf("  consistently calibrated: %.1f to %.1f per cent\n", 100*p_exp, 100*p_uni_ok))
cat("\n  The upper end of the range roughly doubles once the uniform is calibrated the way the\n")
cat("  text says all three are. A single upper quantile does not pin the tail fraction, which\n")
cat("  is the honest reading of a spread this wide.\n")
