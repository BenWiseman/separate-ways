# SUPERSEDED 2026-09-21. Every profile this file scans, n_min + A exp(-(x-c)^2/2w^2), violates
# the band at the origin: the band closes to n=1/2 there, so any positive bump exceeds n_max(0).
# The grid starts at x = 0.005 and never checks, so the 'admissible' column is not admissibility.
# Use same_I_admissible.R, which builds inside the band by construction. The 1.4 and 2.7 per cent
# figures this file produced are withdrawn from the paper.
# Does the measured abundance determine the decoherence prefactor? M_1 ~ I^(-2/5) is the only
# measured functional, and t_dec ~ M_1^(2/3) R^(-2/3), so the question is whether I fixes R.
#
# FIRST, a constraint that decides half of it. The Theta-invariant band has FLOOR n_min, and
# every admissible excursion therefore RAISES I. So holding I fixed while reshaping the profile
# would need a dip below n_min, which leaves the band. At I = I_0 the minimum state is the
# UNIQUE admissible profile and R = R_0 is determined. The degeneracy, if any, lives at I > I_0.

nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmax <- function(x) (1+sqrt(1-exp(-x^2)))/2
g    <- function(x,c,w) exp(-(x-c)^2/(2*w^2))
W2   <- function(f) integrate(function(x) x^2*f(x), 0, 30, subdivisions=6000)$value
# CORRECTED 2026-09-21: the phase-averaged exponent is -log max(n,1-n), not -log(1-n). The two
# agree only for n <= 1/2. R_under_band.R was fixed yesterday and THIS file was missed, which is
# the same defect in a second place. See R_phase_average_fix.R.
WL   <- function(f) integrate(function(x) x^2*(-log(pmax(pmin(pmax(f(x),0),1-1e-13),
                                                        1-pmin(pmax(f(x),0),1-1e-13)))),
                              0, 30, subdivisions=6000)$value
I0 <- W2(nmin); R0 <- WL(nmin)/I0
xs <- seq(0.005, 12, length.out=6000)

cat(sprintf("  reference I = %.7f (5.1), R = %.5f (4.2)\n", I0/pi^2, R0))
cat("  Checked: a profile holding I at I_0 while reshaping must dip below n_min and leaves\n")
cat("  the band, so at the endpoint the state is unique and R is not free.\n")

for (tgt in c(1.05, 1.15, 1.30)) {
  cat(sprintf("\n=== admissible profiles at I/I_0 = %.2f, bump only, no dip\n\n", tgt))
  cat("      bump at x1    width    amplitude      I/I_0        R        R/R_0   in band?\n")
  rr <- c()
  for (x1 in c(0.4, 0.8, 1.2, 1.8, 2.4)) for (w in c(0.25, 0.5)) {
    f <- function(A) W2(function(x) nmin(x)+A*g(x,x1,w))/I0 - tgt
    if (f(0)*f(3) > 0) next
    A <- uniroot(f, c(0, 3), tol=1e-12)$root
    n <- function(x) nmin(x) + A*g(x,x1,w)
    ok <- all(n(xs) <= nmax(xs) + 1e-12)
    I <- W2(n); R <- WL(n)/I
    if (ok) rr <- c(rr, R/R0)
    cat(sprintf("   %10.2f %8.2f %12.5f %11.6f %9.5f %9.5f %10s\n",
        x1, w, A, I/I0, R, R/R0, if (ok) "yes" else "NO"))
  }
  if (length(rr) > 1)
    cat(sprintf("\n   R/R_0 over admissible profiles: %.5f to %.5f, so at FIXED M_1 the t_dec\n   prefactor R^(-2/3) spans %.1f per cent.\n",
        min(rr), max(rr), 100*abs(min(rr)^(-2/3)/max(rr)^(-2/3) - 1)))
}

cat("\n=== reading, flatly\n\n")
cat("  At the endpoint itself the abundance determines everything: the minimum profile is the\n")
cat("  unique admissible state with I = I_0, so R = R_0 and the t_dec prefactor is fixed.\n")
cat("  Away from the endpoint it does not: profiles with the same I but different shapes give\n")
cat("  different R, so a measured M_1 below the endpoint leaves the prefactor undetermined.\n")
cat("  Since 5.1 reads M_1 as an UPPER ENDPOINT rather than a measurement, the operative case\n")
cat("  is the second, and 4.2's caveat is complete as written rather than provisional.\n")
