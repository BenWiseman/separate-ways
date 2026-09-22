# Does the surviving band leave a signature beyond the abundance integral? There are two
# functionals of n in the paper, not one:
#     I = (1/pi^2) int x^2 n dx            -> the mass, M_1 ~ I^(-2/5)
#     R = int x^2 [-log(1-n)] dx / int x^2 n dx   -> enters t_dec
# 4.2 eliminates I between them to get t_dec ~ M_1^(2/3) "with no free parameter". Carrying R
# through: t_dec ~ M_1^(-1) (R I)^(-2/3) and I ~ M_1^(-5/2), so
#     t_dec ~ M_1^(2/3) R^(-2/3).
# The claim is parameter-free only if R is FIXED. R is a different functional of n from I, so a
# band excursion need not move them together. Compute how much R moves.

nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmax <- function(x) (1+sqrt(1-exp(-x^2)))/2
w    <- function(f, lo, hi) integrate(function(x) x^2*f(x), lo, hi, subdivisions=4000)$value
# CORRECTED 2026-09-20: the phase-averaged exponent is -log max(n,1-n) (Jensen), not -log(1-n).
# The two agree for n <= 1/2 and differ on the band EDGE, where n > 1/2, which reversed the sign
# of the whole effect. See R_phase_average_fix.R.
wl   <- function(f, lo, hi) integrate(function(x) x^2*(-log(pmax(f(x),1-f(x)))), lo, hi, subdivisions=4000)$value

num0 <- wl(nmin, 0, Inf); den0 <- w(nmin, 0, Inf)
cat(sprintf("   minimum state:  R = %.5f   (4.2 quotes 1.07037)   I = %.7f  (quotes 0.0127597)\n",
    num0/den0, den0/pi^2))

cat("\n   worst admissible excursion: band edge to x_c, minimum beyond\n\n")
cat("      x_c        I         I/I_0     R        R/R_0    M_1 (PeV)   t_dec factor   net t_dec\n")
R0 <- num0/den0; I0 <- den0/pi^2
for (xc in c(0.10,0.25,0.50,0.644,0.75,1.00,1.50)) {
  den <- w(nmax,0,xc) + w(nmin,xc,Inf)
  num <- wl(nmax,0,xc) + wl(nmin,xc,Inf)
  I <- den/pi^2; R <- num/den
  mfac <- (I/I0)^(-2/5)                       # M_1 relative
  tfac <- mfac^(2/3) * (R/R0)^(-2/3)          # t_dec relative, carrying R
  cat(sprintf("   %6.3f %10.6f %9.4f %9.5f %8.4f %11.1f %13.4f %11.4f\n",
      xc, I, I/I0, R, R/R0, 491.6*mfac, (R/R0)^(-2/3), tfac))
}

cat("\n=== reading\n\n")
cat("  R is NOT constant across the band, so the band is not invisible: it leaves a second\n")
cat("  signature besides the abundance. Whether that matters depends on the size.\n")
d <- sapply(c(0.25,0.50,0.619), function(xc){
  den <- w(nmax,0,xc)+w(nmin,xc,Inf); num <- wl(nmax,0,xc)+wl(nmin,xc,Inf)
  (num/den)/R0 })
cat(sprintf("  Over the range KM3NeT already allows (x_c up to 0.619), R/R_0 runs %.4f to %.4f,\n",
    min(d), max(d)))
cat(sprintf("  so R^(-2/3) moves t_dec by at most %.1f per cent.\n", 100*abs(max(d)^(-2/3)-1)))
cat("  The 2/3 POWER itself is untouched: it comes from eliminating I and does not involve R.\n")
cat("  What acquires a dependence is the PREFACTOR. So 4.2's scaling law survives and its\n")
cat("  normalisation carries a state-dependence the section does not currently flag.\n")
