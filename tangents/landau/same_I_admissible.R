# The fixed-abundance spread in R, over profiles that are ADMISSIBLE EVERYWHERE.
#
# same_I_different_R.R scanned n = n_min + A*exp(-(x-c)^2/2w^2). Every such profile violates the
# band at the origin: the band closes there, n_min(0) = n_max(0) = 1/2, so ANY positive bump gives
# n(0) > n_max(0). Its grid starts at x = 0.005 and never looks. Its "admissible" flags are
# therefore not admissibility, and the 1.4 and 2.7 per cent it reported are not a spread over the
# allowed family.
#
# Build inside the band by construction instead: n = n_min + Delta*w(x) with Delta = n_max - n_min
# and w in [0,1]. Then n(0) = 1/2 automatically, since Delta(0) = 0.
nmin  <- function(x) (1-sqrt(1-exp(-x^2)))/2
Delta <- function(x) sqrt(1-exp(-x^2))
W2 <- function(f) integrate(function(x) x^2*f(x), 0, 30, subdivisions=8000)$value
WL <- function(f) integrate(function(x){u<-pmin(pmax(f(x),0),1-1e-14)
                                        x^2*(-log(pmax(u,1-u)))}, 0, 30, subdivisions=8000)$value
edge <- function(xc) function(x) nmin(x) + Delta(x)*exp(-(x/xc)^40)        # near top-hat
bump <- function(A)  function(x) nmin(x) + Delta(x)*A*exp(-(x^2-2.4^2)^2/(8*2.4^2*0.5^2))
I0 <- W2(nmin); R0 <- WL(nmin)/I0
cat(sprintf("  reference: I0 = %.10f, R0 = %.8f\n\n", I0/pi^2, R0))
cat("  admissibility at the origin (must be exactly 1/2):\n")
for (f in list(edge(0.54), bump(0.003)))
  cat(sprintf("     n(0) = %.12f\n", f(1e-12)))
cat("\n   I/I_0    x_c(edge)    A(bump)      R_edge      R_bump    fixed-I time ratio - 1\n")
for (tgt in c(1.15, 1.30)) {
  xc <- uniroot(function(v) W2(edge(v))/I0 - tgt, c(0.05, 3), tol=1e-12)$root
  A  <- uniroot(function(v) W2(bump(v))/I0 - tgt, c(1e-6, 0.999), tol=1e-12)$root
  Re <- WL(edge(xc))/W2(edge(xc)); Rb <- WL(bump(A))/W2(bump(A))
  cat(sprintf("   %5.2f  %11.8f  %11.8f  %10.8f  %10.8f  %18.3f%%\n",
      tgt, xc, A, Re, Rb, 100*((Re/Rb)^(-2/3)-1)))
}
cat("\n  Both profiles sit inside the band at every momentum and return the same abundance, and\n")
cat("  their decoherence prefactors differ by far more than the 1.4 and 2.7 per cent the old\n")
cat("  scan reported. Those figures are withdrawn. What the examples establish is that the\n")
cat("  freedom EXISTS at fixed I; they do not bound its size, and neither does any finite scan.\n")
