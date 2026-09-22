# two_port_noise.R -- the fold-even / fold-odd force-noise ports of Paper 2 sec 2.10 (second paragraph).
# Base R only. Units H = 1, scalar charge lambda = 1, one transverse component.
# Toy: a comoving scalar charge (action -m int ds + lambda int ds phi) in the conformal Bunch-Davies vacuum
# of dS4, and its literal antipodal image in the push-forward frame (e' = -e, the A-gluing of App. A.1).
# Conformal BD Wightman function G = H^2/(8 pi^2 (1 - Z)), Z = cosh(H dt) on the geodesic, image Z -> -Z:
#   N_0(t) =  1/(32 pi^2 sinh^4((t - i eps)/2))   local force kernel (Bunch-Davies boundary value)
#   N_I(t) = -1/(32 pi^2 cosh^4(t/2))              image kernel   = -N_0(t - i beta/2), beta = 2 pi
#   A_0 = int N_0 dt = 1/(12 pi^2);  A_I = int N_I dt = -1/(12 pi^2);  ports A_0 +/- A_I = 0, 1/(6 pi^2)
#   S_I(w) = -w (w^2+1)/(12 pi sinh(pi w));  S_0^sym(w) = w (w^2+1) coth(pi w)/(12 pi)
#   S_even = w (w^2+1) tanh(pi w/2)/(12 pi) -> w^2/24;  S_odd = w (w^2+1) coth(pi w/2)/(12 pi) -> 1/(6 pi^2)
#   spectral function rho(w) = [S_0^W(w) - S_0^W(-w)]/2 = w (w^2+1)/(12 pi)  (w^3/(12 pi) in the flat limit)
# Source of the expected values: calc/tangents/h15_two_port/two_port_psd.py and its output.txt (2026-09-14);
# An independent derivation gives A_0, A_I and the eigenvalues 0, H^3/(6 pi^2) independently.
source("helpers.R")
N_I <- function(t) -1/(32*pi^2*cosh(t/2)^4)
# local kernel on the line Im t = -eps, real part (the imaginary part is odd and integrates to zero)
N_0_re <- function(t, eps) { z <- complex(real = t/2, imaginary = -eps/2); Re(1/(32*pi^2*sinh(z)^4)) }
A_I <- integrate(N_I, -Inf, Inf, rel.tol = 1e-12)$value
A_0_kms <- integrate(function(t) 1/(32*pi^2*cosh(t/2)^4), -Inf, Inf, rel.tol = 1e-12)$value  # Im t = -pi: sinh^4 -> cosh^4
A_0_eps <- integrate(function(t) N_0_re(t, 0.3), -60, 60, rel.tol = 1e-10, subdivisions = 4000)$value
ev <- eigen(matrix(c(A_0_kms, A_I, A_I, A_0_kms), 2))$values
S_I <- function(w) -w*(w^2+1)/(12*pi*sinh(pi*w))
S_0 <- function(w)  w*(w^2+1)/(12*pi*tanh(pi*w))
S_even <- function(w) S_0(w) + S_I(w); S_odd <- function(w) S_0(w) - S_I(w)
S_I_quad <- function(w) integrate(function(t) cos(w*t)*N_I(t), -60, 60, rel.tol = 1e-10, subdivisions = 4000)$value
# Wightman power via the shifted contour: S_0^W(w) = e^{pi w} * FT[1/(32 pi^2 cosh^4(t/2))](w)
S_0W <- function(w) exp(pi*w)*integrate(function(t) cos(w*t)/(32*pi^2*cosh(t/2)^4), -60, 60, rel.tol = 1e-10, subdivisions = 4000)$value
rho <- function(w) (S_0W(w) - S_0W(-w))/2
cat(sprintf("A_I = %.10f (exact -1/(12 pi^2) = %.10f); A_0 (Im t = -pi) = %.10f; A_0 (Im t = -0.3) = %.10f\n",
            A_I, -1/(12*pi^2), A_0_kms, A_0_eps))
cat(sprintf("two-port eigenvalues: %.3e and %.10f (1/(6 pi^2) = %.10f)\n", min(abs(ev)), max(ev), 1/(6*pi^2)))
cat(sprintf("S_even(0.01)/(w^2/24) = %.6f ; S_odd(1e-3) = %.10f ; rho(1) = %.10f vs w(w^2+1)/(12 pi) = %.10f\n",
            S_even(0.01)/(0.01^2/24), S_odd(1e-3), rho(1), 2/(12*pi)))
cat("\nchecks:\n")
report("A_I = -H^3/(12 pi^2) (sec 2.10)", expected = -1/(12*pi^2), reproduced = A_I)
report("A_0 = +H^3/(12 pi^2), contour Im t = -pi (sec 2.10)", expected = 1/(12*pi^2), reproduced = A_0_kms)
report("A_0 on Im t = -0.3 (any 0 < eps < 2 pi gives the same)", expected = 1/(12*pi^2), reproduced = A_0_eps)
report("even port: 1 + (A_0 + A_I)/(A_0 - A_I) = 1 (sec 2.10; eigenvalue exactly zero)", expected = 1,
       reproduced = 1 + (A_0_kms + A_I)/(A_0_kms - A_I), tol = 1e-9, note = "the image integral is minus the local one")
report("odd-port eigenvalue A_0 - A_I = H^3/(6 pi^2) (sec 2.10)", expected = 1/(6*pi^2), reproduced = max(ev))
report("S_I(w=1) closed form vs quadrature", expected = S_I(1), reproduced = S_I_quad(1))
report("S_I(w=0.25) closed form vs quadrature", expected = S_I(0.25), reproduced = S_I_quad(0.25))
report("S_even -> H w^2/24 at small w (w = 0.01)", expected = 1, reproduced = S_even(0.01)/(0.01^2/24))
report("S_odd(w -> 0) = H^3/(6 pi^2) (w = 1e-3)", expected = 1/(6*pi^2), reproduced = S_odd(1e-3))
report("spectral function rho(1) = w(w^2+1)/(12 pi)", expected = 2/(12*pi), reproduced = rho(1))
report("spectral function rho(2) = w(w^2+1)/(12 pi)", expected = 10/(12*pi), reproduced = rho(2))
report("spectral function rho(3) = w(w^2+1)/(12 pi)", expected = 30/(12*pi), reproduced = rho(3))
# flat-space kernel 1/(2 pi^2 (t - i eps)^4): the same contour construction must give rho_flat(w) = w^3/(12 pi),
# the coefficient of the scalar radiation-reaction force (lambda^2/12 pi) d^3x/dt^3
F_flat <- function(t, eps) 1/(2*pi^2*complex(real = t, imaginary = -eps)^4)
SW_flat <- function(w, eps = 1) {
  re <- integrate(function(t) Re(exp(1i*w*t)*F_flat(t, eps)), -200, 200, rel.tol = 1e-10, subdivisions = 8000)$value
  exp(w*eps)*re
}
rho_flat <- function(w) (SW_flat(w) - SW_flat(-w))/2
report("flat-space rho(2) = w^3/(12 pi) (radiation-reaction coefficient)", expected = 8/(12*pi), reproduced = rho_flat(2))
