# The companion argues an even horizon response can never mix l with l+-1, so a ringdown showing
# l = 2 <-> 3 mixing at fixed m would falsify it. That argument uses SCALAR Legendre parity.
# Gravitational Teukolsky modes are SPIN-WEIGHTED (s = -2), and spin-weighted harmonics do not
# carry the scalar parity that the argument needs. Test the decisive matrix element.
#
#   sYlm ∝ d^l_{m,-s}(theta) e^{i m phi},  so for s = -2, m = 2 we need d^l_{2,2}.
d22 <- function(l, th) {                      # Wigner small-d, m' = m = 2
  u <- cos(th); c2 <- ((1+u)/2)^2
  if (l == 2) return(c2)
  if (l == 3) return(c2*(3*u-2))
  stop("only l = 2,3 needed here")
}
chk <- function(l) abs(d22(l,0) - 1)          # d^l_{mm}(0) = 1
cat(sprintf("  normalisation d^l_22(0) = 1 :  l=2 err %.1e   l=3 err %.1e\n", chk(2), chk(3)))

orth <- integrate(function(th) d22(2,th)*d22(3,th)*sin(th), 0, pi)$value
cat(sprintf("  orthogonality <2|3> (must vanish)          : %.3e\n", orth))

# <_-2Y_22 | cos^2 theta | _-2Y_32>
rad <- integrate(function(th) d22(2,th)*d22(3,th)*cos(th)^2*sin(th), 0, pi)$value
val <- (sqrt(5*7)/2) * rad
cat(sprintf("\n  radial integral                            : %.10f   (exact 2/21 = %.10f)\n",
            rad, 2/21))
cat(sprintf("  <_-2Y_22 | cos^2 | _-2Y_32>                : %.10f\n", val))
cat(sprintf("  sqrt(35)/21                                : %.10f\n", sqrt(35)/21))
cat(sprintf("  difference                                 : %.2e\n", abs(val - sqrt(35)/21)))

cat("\n  Non-zero. An even response CAN mix l = 2 and l = 3 at fixed m = 2 for spin-weighted\n")
cat("  modes, so 'never l +- 1' is not a gravitational falsifier. For contrast, the same\n")
cat("  element built from SCALAR harmonics does vanish, which is why the scalar calculation\n")
cat("  looked conclusive:\n")
P <- function(l,x) if (l==2) (3*x^2-1)/2 else (5*x^3-3*x)/2
sc <- integrate(function(th) P(2,cos(th))*P(3,cos(th))*cos(th)^2*sin(th), 0, pi)$value
cat(sprintf("     scalar <P2|cos^2|P3> = %.3e  (vanishes by parity)\n", sc))
cat("\n  The scalar coefficients in the companion are correct; what fails is carrying their\n")
cat("  parity selection rule over to the gravitational case.\n")
