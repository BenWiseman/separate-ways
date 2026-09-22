# R_under_band.R uses -log(1-n) as the phase-averaged decoherence exponent for EVERY occupation.
# That form is only the n <= 1/2 branch. The phase average of log|a1 + e^{i phi} a2| over a full
# period is log max(|a1|,|a2|) (Jensen), so with |a1|^2 = 1-n and |a2|^2 = n the exponent is
#     -<log|D|^2> = -log max(n, 1-n),
# which equals -log(1-n) only when n <= 1/2. The band EDGE has n > 1/2, which is exactly where
# the script applies the wrong branch. Check the identity first, then redo the table.

cat("=== 1. the phase-average identity, checked numerically ===\n")
ph <- seq(0, 2*pi, length.out=200001)[-200001]      # drop the duplicate endpoint
cat("     n      mean log|D|^2      -log max(n,1-n)      -log(1-n)\n")
for (n in c(0.05, 0.3, 0.5, 0.7, 0.95)) {
  a1 <- sqrt(1-n); a2 <- sqrt(n)
  emp <- mean(log(Mod(a1 + exp(1i*ph)*a2)^2))
  cat(sprintf("  %5.2f  %15.9f  %18.9f  %14.9f\n",
      n, -emp, -log(max(n,1-n)), -log(1-n)))
}

cat("\n=== 2. R with the corrected exponent ===\n")
nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmax <- function(x) (1+sqrt(1-exp(-x^2)))/2
w  <- function(f,lo,hi) integrate(function(x) x^2*f(x), lo, hi, subdivisions=4000)$value
# old: -log(1-n).  new: -log(max(n,1-n))
wl_old <- function(f,lo,hi) integrate(function(x) x^2*(-log(1-f(x))), lo,hi, subdivisions=4000)$value
wl_new <- function(f,lo,hi) integrate(function(x) x^2*(-log(pmax(f(x),1-f(x)))), lo,hi, subdivisions=4000)$value
I0 <- w(nmin,0,Inf)/pi^2
R0_old <- wl_old(nmin,0,Inf)/w(nmin,0,Inf)
R0_new <- wl_new(nmin,0,Inf)/w(nmin,0,Inf)
cat(sprintf("  minimum state: R_old = %.6f  R_new = %.6f  (identical: n_min < 1/2 everywhere)\n\n",
            R0_old, R0_new))
cat("    x_c     I/I_0     R_old   R_old/R0   R_new   R_new/R0   t_dec factor old   new\n")
for (xc in c(0.25, 0.50, 0.619, 0.644, 1.00)) {
  den <- w(nmax,0,xc) + w(nmin,xc,Inf)
  Ro <- (wl_old(nmax,0,xc) + wl_old(nmin,xc,Inf))/den
  Rn <- (wl_new(nmax,0,xc) + wl_new(nmin,xc,Inf))/den
  I <- den/pi^2; mfac <- (I/I0)^(-2/5)
  cat(sprintf("  %6.3f %9.6f %8.5f %9.5f %8.5f %9.5f %14.4f %8.4f\n",
      xc, I/I0, Ro, Ro/R0_old, Rn, Rn/R0_new,
      mfac^(2/3)*(Ro/R0_old)^(-2/3), mfac^(2/3)*(Rn/R0_new)^(-2/3)))
}
cat("\n  The sign of the R excursion REVERSES. With the correct branch R falls below its\n")
cat("  minimum-state value instead of rising above it, so the R^(-2/3) prefactor moves the\n")
cat("  other way and the net decoherence time rises rather than falls.\n")
