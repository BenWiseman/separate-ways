# Two things. First the route I named last turn, priced rather than guessed. Then what looking
# at it turned up, which is better.
#
# BACKREACTION. At the bang the produced heavy particles are relativistic (p >> ma = gamma*eta
# as eta -> 0), so they redshift as a^-4 exactly like the radiation and the ratio is a pure
# constant. With a = a_1 eta: H = 1/(a_1 eta^2), rho_rad a^4 = 3 Mpl^2 a_1^2, and
# rho_heavy a^4 = (g/2pi^2)(gamma^2/pi^2) J with J = int x^3 n(x) dx, gamma = M_1 a_1.
# The a_1 cancels:   rho_heavy/rho_rad = g J M_1^2 / (6 pi^4 Mpl^2).

n   <- function(x) (1-sqrt(1-exp(-x^2)))/2
J   <- integrate(function(x) x^3*n(x), 0, Inf)$value
M1  <- 4.916e8      # GeV
Mpl <- 2.435e18     # GeV, reduced
for (g in c(2,4)) {
  r <- g*J*M1^2/(6*pi^4*Mpl^2)
  cat(sprintf("   g = %d :  J = %.6f,  rho_heavy/rho_rad = %.3e\n", g, J, r))
}
r <- 4*J*M1^2/(6*pi^4*Mpl^2)
cat(sprintf("\n   The band allows the occupation to rise from n_min to at most 1, a factor of\n"))
cat(sprintf("   order 100 at moderate x. Backreaction only bites when the ratio reaches 1, i.e.\n"))
cat(sprintf("   a boost of %.1e.  FLAT NEGATIVE: backreaction is about %.0f orders of magnitude\n",
    1/r, log10(1/r)-2))
cat("   too weak to constrain mu(p). Third failed selector, after Hadamard and analyticity.\n")

cat("\n=== what the same state says about WHERE the sheets are joined\n\n")
cat("  The pair state is sqrt(1-n)|00> + sqrt(n)|11>, pure, so the entanglement entropy\n")
cat("  between the sheets for one mode is S(n) = -n ln n - (1-n) ln(1-n), and the mutual\n")
cat("  information is 2S. Two ends are already fixed, and they are fixed at opposite values.\n\n")
S <- function(x){ nn<-n(x); ifelse(nn<=0|nn>=1, 0, -nn*log(nn)-(1-nn)*log(1-nn)) }
cat(sprintf("   x -> 0 : n = 1/2 exactly (the Hamiltonian is diagonal, no mixing), S = ln 2 = %.5f\n", log(2)))
cat(sprintf("   x  = 3 : n = %.3e, S = %.3e\n", n(3), S(3)))
cat(sprintf("   x  = 5 : n = %.3e, S = %.3e\n", n(5), S(5)))
cat("\n   So the sheets are MAXIMALLY entangled in the infrared and unentangled in the\n")
cat("   ultraviolet. The join is not a point and not a constant: it is a crossover.\n\n")
xh <- uniroot(function(x) S(x)-log(2)/2, c(0.1,4))$root
xq <- uniroot(function(x) S(x)-log(2)/10, c(0.1,6))$root
cat(sprintf("   half of maximum entanglement at x = %.4f\n", xh))
cat(sprintf("   one tenth of maximum at        x = %.4f\n", xq))
tot <- integrate(function(x) x^2*S(x), 0, Inf)$value
cat(sprintf("   total mode-weighted entanglement int x^2 S dx = %.6f, finite\n", tot))
cat(sprintf("   compare the production integral's own weight int x^2 n dx = %.6f\n",
    integrate(function(x) x^2*n(x),0,Inf)$value*1))

cat("\n=== and the crossover sits exactly at the avoided crossing's own scale\n\n")
cat("  x = sqrt(pi) p / sqrt(gamma), so x ~ 1 is p ~ sqrt(gamma/pi): the momentum at which\n")
cat("  the mode's gap p and its sweep scale sqrt(gamma) are comparable. That is the Landau-\n")
cat("  Zener gap scale itself. The scale at which the two sheets stop being one object is the\n")
cat("  scale at which the crossing stops being sudden. Nothing was inserted to make that so.\n")
cat(sprintf("  Numerically the half-entanglement point is x = %.3f, i.e. p = %.3f sqrt(gamma).\n",
    xh, xh/sqrt(pi)))
