# graviton_whichpath.R -- standard which-path decoherence from graviton emission during a
# split/recombine of a laboratory mass (Paper 2 sec 5, laboratory paragraph). Base R only.
#
# Formula (linearised gravity, quadrupole order): the decoherence exponent between two
# branches equals half the mean number of gravitons radiated by the branch DIFFERENCE of
# the quadrupole (Matsui, arXiv:2607.20867, abstract: "the decoherence exponent equals half
# the mean number of gravitons").  With Q_ij = q(t) diag(2/3,-1/3,-1/3) (motion along x),
#   E_rad = (G/5c^5) int Qddd_ij Qddd_ij dt = (G/5c^5)(2/3)(1/pi) int_0^inf w^6 |q^(w)|^2 dw
#   N     = int (dE/dw)/(hbar w)  = (2G/(15 pi hbar c^5)) int_0^inf w^5 |q^(w)|^2 dw
#   Gamma = N/2.
# Two equal masses m at separation D along x, reduced mass mu = m/2, relative coordinate
# r = D + Delta_x g(t/T) on the mixed branch pair (each mass displaced by +/- Delta_x/2):
#   Q_ij = mu (r_i r_j - delta_ij r^2/3) = mu r^2 diag(2/3,-1/3,-1/3),  so  q = mu r^2  and
#   Delta q = mu[(D + Delta_x g)^2 - D^2] = m D Delta_x g + (m/2) Delta_x^2 g^2   (exactly).
# Leading order in Delta_x/D:  q(t) = m D Delta_x g(t/T).
# CORRECTION 2026-09-12 09:4x: the first version used (4/3) m D Delta_x, mixing the xx tensor
# component with the STF scalar; the exponent was high by 16/9. Caught by an independent check
# (2026-09-12) and reproduced here.
# Profile g = sin^4(pi t) on [0,1] (smooth to third derivative at both ends).
# Spectral constants (A,B,C) = int_0^inf u^5 (|g^|^2, |g2^|^2, Re g^ g2^*) du with g2 = g^2 are
# the independently derived Gram constants (11656.84357309, 42092.91505331, 16883.19030278).
source("helpers.R")
G <- 6.67430e-11; hbar <- 1.054571817e-34; c <- 2.99792458e8
# --- spectral constant by FFT on a zero-padded grid (window [0,1] inside [0,L]) ---
L <- 64; N <- 2^20; dt <- L/N; t <- (0:(N-1))*dt
g <- ifelse(t <= 1, sin(pi*pmin(t,1))^4, 0)
g2 <- g^2
gh <- fft(g)*dt; g2h <- fft(g2)*dt               # continuous FT approximations
w <- 2*pi*(c(0:(N/2-1), -(N/2):-1))/L            # angular frequency grid
dw <- w[2]-w[1]
A <- sum(abs(w)^5*Mod(gh)^2)*dw/2                # = int_0^inf w^5 |g^|^2 dw (one-sided)
B <- sum(abs(w)^5*Mod(g2h)^2)*dw/2
C <- sum(abs(w)^5*Re(gh*Conj(g2h)))*dw/2
cat(sprintf("(A,B,C) for g = sin^4(pi t) on [0,1]: %.8f %.8f %.8f\n", A, B, C))
# --- Bose et al. design (arXiv:1707.06050): m = 1e-14 kg, D = 450 um, Delta_x = 250 um; T = 1 s ---
m <- 1e-14; D <- 450e-6; Dx <- 250e-6; T <- 1
pref <- (2*G/(15*pi*hbar*c^5)) / 2                 # Gamma = N/2
Gam_lead <- pref * (m*D*Dx)^2 * A / T^4            # leading order in Delta_x/D (closed-source estimate)
Gam_full <- pref * m^2 * (D^2*Dx^2*A + D*Dx^3*C + Dx^4*B/4) / T^4   # ideal closed two-body continuation
cat(sprintf("Gamma leading = %.7e ; Gamma with the g^2 terms (ideal two-body, actuator-dependent) = %.7e\n", Gam_lead, Gam_full))
cat(sprintf("Delta_x/D = %.3f, so the nonlinear terms are not small; a separate analysis shows they depend on the actuator\n", Dx/D))
cat("\nchecks:\n")
report("spectral constant A (independently derived)", expected = 11656.84357309, reproduced = A, tol = 1e-4, mode = "rel",
       note = "2026-09-12 Gram constant A, sin^4 profile, T = 1")
report("spectral constant B (independently derived)", expected = 42092.91505331, reproduced = B, tol = 1e-4, mode = "rel")
report("spectral constant C (independently derived)", expected = 16883.19030278, reproduced = C, tol = 1e-4, mode = "rel")
report("Gamma leading, Bose design, T = 1 s (paper sec 5)", expected = 8.2e-59, reproduced = Gam_lead, tol = 0.01, mode = "rel",
       note = "paper quotes 8e-59 (one figure); independent value 8.1822158e-59")
report("Gamma ideal two-body", expected = 1.7045731e-58, reproduced = Gam_full, tol = 1e-4, mode = "rel",
       note = "actuator-dependent; not quoted in the paper as a value")
