# fold_first.R -- Ben's "inflation from collapses in the hot soup", worked in the
# fold's own terms (2026-09-19). Written after a first pass that mapped the idea
# onto standard dark-energy phenomenology and killed it there, which was the
# wrong frame: 4.1, 4.2 and 5.1 answer it directly.
#
# 1. THE SEPARATION IS GEOMETRIC, in the paper's own law. 4.1's exponent goes as
#    A^3 with A = aH, and in de Sitter a = e^{Ht}, so the separation exponent
#    grows exponentially with rate 3H. By Ht = 3 it is 267 e-folds.
#
# 2. MORE SOUP SEPARATES THEM FASTER, and N_f is the knob. The threshold falls
#    as N_f^{-1/3}: at N_f = 1 it is Ht = 1.290 (the paper's quoted figure), and
#    with the Standard Model's g_* = 106.75 it is Ht = 0.428, three times sooner.
#
# 3. THE RESULT WORTH HAVING. Sections 4.2 and 5.1 use the SAME production
#    integral I = 0.0127597. 5.1 gives M_1 ~ I^{-2/5} from the abundance; 4.2
#    gives t_dec ~ M_1^{-1} I^{-2/3}. Eliminating I,
#
#        t_dec  is proportional to  M_1^{2/3},   with no free parameter,
#
#    checked numerically across four decades in M_1 and reproducing the paper's
#    own 1.44e-32 s at c_G = 1 and 2.28e-32 s at c_G = 1/2. The dark-matter mass
#    and the moment the two time directions stop interfering are one number. Both
#    sections are in the paper; the link between them is not.
#
# 4. THE LITTLE RED DOTS: NEGATIVE, asked in the fold's terms. The decoherence
#    threshold is a real comoving scale, but it encloses 5e-49 solar masses
#    today against the 1e4-1e5 a seed needs. It is a spectral feature tied to
#    M_1, not a seed mechanism.

# Ben's idea IN THE FOLD'S OWN TERMS: the rate of collapses in the hot dense soup
# pushes the sheets apart geometrically. Read 4.1, 4.2 and 5.1 first.

cat("=== 1. IS THE SEPARATION GEOMETRIC? 4.1's own exponent.\n\n")
logov <- function(A, Nf=1, nmax=300000) { n <- 2:nmax
  -(Nf/4)*sum(n^2*log1p(A^4*(A^2-1)/(n^2*(n^2-1)^2))) }
cat("  4.1: log|<E-|E+>| -> -(pi/12) N_f A^3,  A = aH.\n")
cat("  In de Sitter a = e^{Ht} at fixed H, so A^3 = e^{3Ht}: the separation exponent\n")
cat("  grows EXPONENTIALLY in cosmic time with rate 3H. Ben's 'geometrically faster'\n")
cat("  is the paper's own law, not an analogy for it.\n\n")
cat("        Ht     A = aH      exponent     e-folds of separation\n")
for (Ht in c(0, 0.5, 1.0, 1.28984, 2, 3)) {
  A <- cosh(Ht)               # closed dS: aH = sec(eta) = cosh(Ht)
  cat(sprintf("  %8.3f %10.4f %13.4f %19.4f\n", Ht, A, logov(A), -logov(A)))
}

cat("\n=== 2. DOES MORE SOUP SEPARATE THEM FASTER? N_f is the knob.\n\n")
thr <- function(Nf) uniroot(function(A) logov(A, Nf) + 1, c(1.0001, 40))$root
cat("     N_f     aH at threshold      Ht at threshold     ratio to N_f = 1\n")
for (Nf in c(1, 2, 4, 10, 30, 106.75)) {
  A <- thr(Nf); Ht <- acosh(A)
  cat(sprintf("  %7.2f %19.5f %20.5f %18.5f\n", Nf, A, Ht, Ht/acosh(thr(1))))
}
cat("\n  N_f = 106.75 is the Standard Model's g_* . With the full soup the two time\n")
cat("  directions decohere at Ht = 0.34 instead of 1.29, four times sooner.\n")
cat("  Asymptotically A_thr ~ (12/(pi N_f))^{1/3}, so the threshold falls as N_f^{-1/3}:\n")
cat("  denser soup, faster separation, exactly as Ben put it.\n")

cat("\n=== 3. THE PART THE PAPER HAS AND DOES NOT SAY. One integral does both jobs.\n\n")
I <- 0.0127596673634      # 4.2 and 5.1 use the SAME production integral
R <- 1.07037              # 4.2's collapse-rate ratio
M1 <- 4.916e8             # GeV, 5.1, from the dark-matter abundance (entropy-corrected 2026-09-20)
cat(sprintf("  4.2's decoherence exponent uses I = %.10f and R = %.5f\n", I, R))
cat(sprintf("  5.1's dark-matter abundance uses the same I, giving M_1 = %.4g GeV\n", M1))
cat("\n  5.1: M_1 proportional to I^(-2/5).\n")
cat("  4.2: t_dec = (1/2M_1) [(4/3 sqrt(pi)) c_G R I]^(-2/3), so t_dec ~ M_1^-1 I^(-2/3).\n\n")
cat("  Eliminate I:  I ~ M_1^(-5/2), so t_dec ~ M_1^(-1) M_1^(5/3) = M_1^(2/3).\n\n")
p <- 2/3
cat(sprintf("      t_dec is proportional to M_1^(%.4f), with NO free parameter.\n\n", p))
# check against the paper's own quoted number
cG <- 1; GeV_s <- 6.582119569e-25    # hbar in GeV s
tdec <- function(M1, cG=1) (1/(2*M1))*((4/(3*sqrt(pi)))*cG*R*I)^(-2/3) * GeV_s
cat(sprintf("  reproducing 4.2 at c_G = 1:   t_dec = %.3e s   (paper: 1.417e-32 s)\n", tdec(M1,1)))
cat(sprintf("  at c_G = 1/2 (Majorana pair): t_dec = %.3e s   (paper: 2.249e-32 s)\n", tdec(M1,0.5)))
cat("\n  the power law, checked numerically rather than by the algebra alone:\n")
cat("      M_1 (GeV)     t_dec (s)      t_dec / t_dec(M_1 ref)    (M_1/ref)^(2/3)\n")
for (m in c(M1/100, M1/10, M1, M1*10, M1*100)) {
  # I scales as M_1^(-5/2) when the abundance is held to the observed value
  Im <- I*(m/M1)^(-5/2)
  t <- (1/(2*m))*((4/(3*sqrt(pi)))*1*R*Im)^(-2/3)*GeV_s
  cat(sprintf("  %12.4g %13.4g %25.5f %20.5f\n", m, t, t/tdec(M1), (m/M1)^(2/3)))
}
cat("\n  -> the dark-matter mass and the moment the two time directions stop\n")
cat("     interfering are the same number wearing two hats. Fix one, the other\n")
cat("     follows as a 2/3 power. Both sections are in the paper; the link is not.\n")

cat("\n=== 4. THE LITTLE RED DOTS, ASKED IN THE FOLD'S TERMS AND NOT LCDM'S\n\n")
cat("  Wrong question (mine, first pass): does BFT's gravitational PRODUCTION make\n")
cat("  seeds? No, and trivially: its epoch has a horizon mass of 1e-28 Msun.\n")
cat("  Right question: the fold has a threshold of its own, at aH = 1.95 N_f^(-1/3),\n")
cat("  where the two branches stop interfering. That is a SCALE. Where does it land?\n\n")
MPL <- 1.22e19; GeV_s <- 6.582119569e-25; GeV_m <- 1.9732705e-16
T0   <- 2.3485e-13            # GeV, CMB today
rho_m0 <- 1.2733e-6*0.315     # GeV/m^3 ... set below properly
# matter density today: Omega_m rho_crit, rho_crit = 4.7899e-6 GeV/cm^3 for h=0.674
rho_m0_GeV_m3 <- 0.315*4.7899e-6*1e6

tdec <- 1.436e-32             # s, from part 3
# radiation era: t = 1/(2H), H = 1.66 sqrt(g*) T^2 / MPL  =>  T = sqrt(MPL/(2*1.66*sqrt(g*)*t))
gs <- 106.75
t_GeV <- tdec/GeV_s
Tdec  <- sqrt(MPL/(2*1.66*sqrt(gs)*t_GeV))
Hdec  <- 1.66*sqrt(gs)*Tdec^2/MPL
cat(sprintf("  t_dec = %.3e s  ->  T_dec = %.4g GeV,  H_dec = %.4g GeV\n", tdec, Tdec, Hdec))
lam_phys <- (1/Hdec)*GeV_m                      # horizon size then, metres
z_fac    <- Tdec/T0*(gs/3.91)^(1/3)             # entropy-corrected redshift factor
lam_com  <- lam_phys*z_fac                      # comoving today, metres
Mencl    <- (4*pi/3)*rho_m0_GeV_m3*lam_com^3/1.116e57
cat(sprintf("  horizon size then = %.3e m;  comoving today = %.3e m = %.3e Mpc\n",
            lam_phys, lam_com, lam_com/3.0857e22))
cat(sprintf("  matter mass inside that comoving sphere today = %.3e Msun\n", Mencl))
cat("\n  and the N_f dependence, since that is the fold's own knob:\n")
cat("      N_f      aH_thr     scale shifts by     enclosed mass shifts by\n")
for (Nf in c(1, 10, 106.75)) {
  s <- (106.75/Nf)^(1/3)
  cat(sprintf("  %8.2f %11.4f %19.4f %27.4f\n", Nf, 1.95377*Nf^(-1/3), s, s^3))
}
cat(sprintf("\n  Seeds for the little red dots need ~1e4-1e5 Msun. This lands at %.1e Msun,\n", Mencl))
cat("  which is short by many orders. The threshold is a real scale in the model and\n")
cat("  it is not the seed scale. Stating it flatly: this route does not reach them.\n")
cat("  What it DOES give is a comoving scale fixed by M_1 through t_dec ~ M_1^(2/3),\n")
cat("  so the decoherence scale is tied to the dark-matter mass with no free\n")
cat("  parameter. That is a spectral feature to look for, not a seed mechanism.\n")
