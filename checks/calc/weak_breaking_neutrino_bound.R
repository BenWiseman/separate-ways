# Section 7: the weakly broken implementation, its Yukawa and the mass it adds.
#
# The section quotes three numbers with no arithmetic shown: M_1 = 4.916e8 GeV (carried in
# from the cosmology paper's leptogenesis scale), ||y|| = 1.8e-30 and a bound of about
# 2.1e-55 eV. The first is an input. This file checks that the other two follow from it and
# that the displayed formula is consistent with the two definitions above it, which is the
# part a referee can actually catch.
#
#   Gamma_0        = q M_1 / (8 pi),        q = y^dagger y = ||y||^2
#   delta m_nu     = -v^2 y y^T / (2 M_1)   so  ||delta m_nu||_2 = v^2 q / (2 M_1)
#   tau            = hbar / Gamma_0
# Eliminating q:
#   ||delta m_nu||_2 = 4 pi v^2 hbar / (M_1^2 tau)      <- the displayed bound

hbar_GeV_s <- 6.582119569e-25    # CODATA, exact by definition of hbar and eV
v          <- 246.22             # GeV, the Higgs vev as Section 7 states it
M1         <- 4.916e8            # GeV, carried in from the cosmology paper
tau        <- 1e28               # s, the example lifetime, inside [22]'s neutrino-telescope reach

cat("=== Section 7: weakly broken implementation ===\n\n")

q    <- 8 * pi * hbar_GeV_s / (tau * M1)
ynorm<- sqrt(q)
cat(sprintf("  q = 8 pi hbar / (tau M_1)        = %.6e\n", q))
cat(sprintf("  ||y|| = sqrt(q)                  = %.6e     (paper: 1.8e-30)\n", ynorm))

# Two routes to the bound, which is the consistency check worth having.
b_direct  <- v^2 * q / (2 * M1)                        # from ||delta m_nu||_2 = v^2 q / (2 M_1)
b_display <- 4 * pi * v^2 * hbar_GeV_s / (M1^2 * tau)  # from the displayed formula
cat(sprintf("\n  bound from v^2 q / (2 M_1)       = %.6e GeV = %.4e eV\n", b_direct,  b_direct*1e9))
cat(sprintf("  bound from 4 pi v^2 hbar/(M_1^2 tau) = %.6e GeV = %.4e eV   (paper: 2.1e-55 eV)\n",
            b_display, b_display*1e9))
stopifnot(abs(b_direct/b_display - 1) < 1e-12)

cat(sprintf("\n  rounded: ||y|| = %.1e, bound = %.1e eV\n", ynorm, b_display*1e9))

# Scale: how far below the measured splittings this sits. Sum of masses 58.78 meV.
cat(sprintf("\n  for scale, the bound is %.1e of the 58.78 meV mass sum, so the massless\n",
            b_display*1e9 / 58.78e-3))
cat("  state stays massless for any purpose an experiment can reach.\n")

stopifnot(abs(ynorm - 1.8e-30) < 0.05e-30)
stopifnot(abs(b_display*1e9 - 2.1e-55) < 0.05e-55)
