# dp_design.R -- sphere-smeared Diosi-Penrose collapse time for the Bose et al. design, at the
# 100 nm benchmark and at the design separation (Paper 2 sec 5, laboratory paragraph). Base R.
# tau_DP = hbar / Delta E_G,  Delta E_G = 2[U(0) - U(d)],  two identical uniform spheres (m, R):
#   U(d) = -(G m^2/R)[6/5 - x^2/2 + 3x^3/16 - x^5/160], x = d/R  (d <= 2R);  U(d) = -G m^2/d (d >= 2R).
# Same formula as the record's beyond/decoherence_floor.py (checked continuous at d = 2R there).
source("helpers.R")
G <- 6.674e-11; hbar <- 1.0546e-34
U <- function(d, m, R) { x <- d/R; if (d <= 2*R) -(G*m^2/R)*(6/5 - x^2/2 + 3*x^3/16 - x^5/160) else -G*m^2/d }
tau <- function(m, R, d) hbar / (2*(U(d, m, R) - U(0, m, R)))
m <- 1e-14
R_si <- (3*m/(4*pi*2200))^(1/3)      # silica, the record's benchmark density
R_di <- (3*m/(4*pi*3500))^(1/3)      # diamond, Bose et al. 1707.06050 microcrystals
t1 <- tau(m, R_si, 100e-9); t2 <- tau(m, R_di, 250e-6)
cat(sprintf("silica  R = %.2f um: tau_DP(100 nm) = %.3f s\n", R_si*1e6, t1))
cat(sprintf("diamond R = %.2f um: tau_DP(250 um, design separation) = %.4f s ; exponent per second = %.0f\n", R_di*1e6, t2, 1/t2))
cat(sprintf("fold's leading which-path exponent for the same protocol (graviton_whichpath.R): 8.2e-59 ; log10 ratio = %.1f\n", log10((1/t2)/8.2e-59)))
cat("\nchecks:\n")
report("tau_DP, 10^-14 kg silica, 100 nm [s]", expected = 1.8, reproduced = t1, tol = 0.02, mode = "rel",
       note = "beyond/DECOHERENCE_FLOOR.md table: 1.8 s; paper sec 5: 1.8 s")
report("tau_DP, 10^-14 kg diamond, 250 um [s]", expected = 6e-3, reproduced = t2, tol = 0.05, mode = "rel",
       note = "paper sec 5 quotes 6 ms (one figure)")
