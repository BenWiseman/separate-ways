# ==========================================================================================
# Results the paper's own machinery already implies and does not state.  Cheap arithmetic only,
# using objects the paper has: n(x), C(x) = e^{-x^2/2}, the min-entropy decoherence exponent,
# the operator bound Q >= n_* 1, the separable bound, the closed budget.
# ==========================================================================================
n   <- function(x) (1-sqrt(1-exp(-x^2)))/2
Cx  <- function(x) exp(-x^2/2)
h   <- function(x) -log(1-n(x))                 # min-entropy decoherence exponent per mode
W   <- function(f,a=0,b=40) integrate(function(x) x^2*f(x), a, b, subdivisions=20000, rel.tol=1e-12)$value
I2  <- W(n); Rh <- W(h)/I2
cat(sprintf("  baseline checks:  I = %.10f (paper 0.0127597),  R = %.6f (paper 1.07037)\n\n",
            I2/pi^2, Rh))

cat("==========================================================================\n")
cat(" A.  THE DECOHERENCE EXPONENT IS A FUNCTION OF THE CONCURRENCE, AND IT IS EXACT\n")
cat("==========================================================================\n\n")
cat("   n = (1 - sqrt(1-C^2))/2  =>  1-n = (1 + sqrt(1-C^2))/2, so\n")
cat("        h(C) = log 2 - log(1 + sqrt(1 - C^2)).\n\n")
hC <- function(C) log(2) - log(1+sqrt(1-C^2))
xs <- c(0,0.25,0.5,1,1.5,2,3)
cat("        x        C(x)          h from n          h from C          difference\n")
for (x in xs) cat(sprintf("   %6.2f %13.8f %17.10f %17.10f %18.2e\n",
                          x, Cx(x), h(x), hC(Cx(x)), abs(h(x)-hC(Cx(x)))))
cat(sprintf("\n   At x = 0 the pair is a Bell state (C = 1) and h = log 2 = %.10f EXACTLY:\n", log(2)))
cat("   each maximally entangled pair the bang makes contributes exactly ONE BIT to the\n")
cat("   exponent that separates the two time orientations.\n")
cat(sprintf("   The threshold Gammabar = 1 nat is therefore %.4f Bell pairs' worth of entanglement\n",
            1/log(2)))
cat("   inside a Hubble volume.  Not a new calculation: a reading of the one in 3.5.\n")

cat("\n==========================================================================\n")
cat(" B.  THE SAME OPERATOR INEQUALITY THAT CAPS THE MASS CAPS THE DECOHERENCE TIME\n")
cat("==========================================================================\n\n")
cat("   Q >= n_* 1 gives <N_k> >= n_*(p_k) for every admissible state.  h = -log(1-n) is\n")
cat("   strictly increasing in n, so Gammabar = sum_k h(<N_k>) >= sum_k h(n_*(p_k)).\n")
cat("   The branch-distinguishability exponent inherits the SAME bound, in the same direction,\n")
cat("   and t_dec ~ Gammabar^{-2/3}, so at fixed mass and Hubble rate\n\n")
cat("        t_dec  <=  t_dec(least-occupied state).\n\n")
cat("   Check monotonicity of h on the admissible band [n_min, n_max = 1-n_min]:\n\n")
cat("        n        h(n) = -log(1-n)     dh/dn\n")
for (nn in c(0.01,0.1,0.25,0.4,0.5)) cat(sprintf("   %7.3f %18.8f %12.4f\n", nn, -log(1-nn), 1/(1-nn)))
cat("\n   Strictly positive throughout, so the inequality is one-sided with no tolerance needed.\n")
cat(sprintf("   The paper's %s s is therefore a CEILING on the decoherence time, not a\n","1.417e-32"))
cat("   central value: the two time directions are separated by then at the latest.\n")
cat("   (The paper's own band scan gives +12 per cent when M_1 is allowed to move with the\n")
cat("    state; the statement here is the fixed-mass one, and it is an inequality, not a scan.)\n")

cat("\n==========================================================================\n")
cat(" C.  HOW MUCH ENTANGLEMENT IS IN THE DARK MATTER, IN UNITS A PERSON CAN PICTURE\n")
cat("==========================================================================\n\n")
S <- function(x){ nn<-n(x); ifelse(nn<=0|nn>=1, 0, -nn*log(nn)-(1-nn)*log(1-nn)) }
Sint <- W(S)
cat(sprintf("   Int x^2 S dx = %.6f (paper 0.4509);  Int x^2 n dx = %.6f\n", Sint, I2))
cat(sprintf("   ratio Int x^2 S / Int x^2 n = %.4f nats per unit of the production integral\n", Sint/I2))
cat("   n is the occupation PER PAIR MEMBER (3.1), so a pair carries 2n particles and S(n)\n")
cat(sprintf("   nats: entanglement per relic particle = %.4f nats = %.4f bits\n",
            Sint/I2/2, Sint/I2/2/log(2)))
eV_kg <- 1.78266192e-36; M1kg <- 491.6e15*eV_kg
rho <- 0.12*1.878e-26                      # kg m^-3
ndens <- rho/M1kg
cat(sprintf("   number density %.3e m^-3 = %.2f particles per km^3 (paper: 2.6)\n", ndens, ndens*1e9))
cat(sprintf("   so the dark matter carries %.2f bits of pair entanglement per cubic kilometre,\n",
            ndens*1e9*Sint/I2/2/log(2)))
cat(sprintf("   and %.3e bits inside the Earth's volume (1.08e21 m^3).\n", ndens*1.08e21*Sint/I2/2/log(2)))
cat("   Every one of those particles is one half of a pair whose partner has exactly opposite\n")
cat("   momentum and whose concurrence is fixed, not fitted: C(x) = exp(-x^2/2).\n")
frac <- W(function(x) n(x)*(Cx(x)>0.5))/I2
cat(sprintf("\n   Fraction of the production integral carried by modes still more than half\n"))
cat(sprintf("   entangled (C > 1/2, i.e. x < sqrt(2 log 2) = %.4f): %.4f\n", sqrt(2*log(2)), frac))

cat("\n==========================================================================\n")
cat(" D.  THE LINE AND THE NEUTRINO-MASS FLOOR CANNOT BE TRADED AGAINST EACH OTHER\n")
cat("==========================================================================\n\n")
cat("   The rule must be weakly broken for the line to exist.  The SAME Yukawa column that\n")
cat("   breaks it both sets the lifetime and lifts the massless neutrino.  Eliminate it.\n\n")
# CORRECTED 2026-09-21. A pre-submission review found two errors here that partially
# cancelled and left the number a factor of two too large.
#   - Gamma(N -> h nu) = y^2 M_1/(32 pi) is ONE channel. The 1:1:2 structure has four, so
#     the total is y^2 M_1/(8 pi) and a given lifetime constrains y^2 four times harder.
#   - m = y^2 v^2/(2M) is the v = 246 GeV convention. With <H0> = 174 the seesaw mass is
#     y^2 <H0>^2/M, with no extra factor of two.
cat("     Gamma_total = y^2 M_1/(8 pi) over the four channels, and m_1 = y^2 <H0>^2/M_1\n")
cat("     with <H0> = 174 GeV, so   m_1 <= 8 pi <H0>^2 / (tau M_1^2).\n")
cat("     An INEQUALITY, not an equality: a Yukawa column lying in the existing rank-two\n")
cat("     span leaves the lightest mass zero despite a nonzero width.\n\n")
v <- 174; M1g <- 4.916e8                         # GeV
GeV_s <- 6.582119569e-25                         # hbar in GeV s
m1_of_tau <- function(tau_s) 8*pi*v^2/((tau_s/GeV_s)*M1g^2)    # GeV, upper bound
cat("        lifetime tau (s)      induced m_1 (eV)        as a fraction of 58.78 meV\n")
for (tau in c(1e26,1e27,1e28,1e29,1e30)) {
  m1 <- m1_of_tau(tau)*1e9
  cat(sprintf("   %18.0e %22.3e %28.3e\n", tau, m1, m1/58.78e-3))
}
tau_bite <- 8*pi*v^2/((1e-3*1e-9)*M1g^2)*GeV_s
cat(sprintf("\n   To lift m_1 to even 1 meV the lifetime would have to be tau = %.2e s = %.2e yr.\n",
            tau_bite, tau_bite/3.156e7))
cat("   The particle would decay in far less than one expansion time and would not be the dark\n")
cat("   matter at all.  So over the WHOLE range in which the particle is dark matter, breaking\n")
cat("   the rule hard enough to produce a line leaves the 58.78 meV floor untouched to fifty\n")
cat("   orders of magnitude.  The two commitments are independent, and the DESI squeeze of 3.3\n")
cat("   cannot be relieved by turning up the decay.\n")

cat("\n==========================================================================\n")
cat(" E.  KM3-230213A ALREADY SITS ABOVE THE SEPARABLE CEILING\n")
cat("==========================================================================\n\n")
m<-220; lo<-110; hi<-790
s1<-log(m/lo); s2<-log(hi/m)
Psplit <- function(E) if (E>=m) 1-pnorm(log(E/m)/s2) else 1-pnorm(log(E/m)/s1)
for (E in c(245.8, 135.15)) {
  cat(sprintf("   P(E_nu > %7.2f PeV) = %.3f  (split lognormal)  |  %.3f (upper-half sigma) | %.3f (lower)\n",
      E, Psplit(E), 1-pnorm(log(E/m)/s2), 1-pnorm(log(E/m)/s1)))
}
cat("\n   135.15 PeV is the half-mass line of the separable ceiling M_1 <= 270.30 PeV.  If the\n")
cat("   event is the decay line at all, the probability that it lies above the energy a\n")
cat("   separable out-pair marginal could reach is 0.65 to 0.76 across the parametrisations,\n")
cat("   against a near-even 0.44 to 0.47 for the entangled ceiling.  One posterior, two\n")
cat("   thresholds, and the separable one is roughly two-to-one decided while the other is not.\n")

cat("\n==========================================================================\n")
cat(" F.  THE SEPARABLE COMPETITOR HAS ITS OWN CLOCK, AND IT IS NOT THE SAME ONE\n")
cat("==========================================================================\n\n")
nsym <- function(P){ if (P < 1e-10) return(sqrt(P)/2)
  cs <- sqrt(P*(1-P))
  tryCatch(uniroot(function(u) P*u - P/2 + cs*u*(1-u), c(1e-14,0.5), tol=1e-15)$root,
           error=function(e) sqrt(P)/2) }
nsep <- function(x) sapply(x, function(y) nsym(exp(-y^2)))
Wq <- function(f) integrate(function(x) x^2*f(x), 0, 40, subdivisions=8000)$value
Isep <- Wq(nsep); Rsep <- Wq(function(x) -log(1-nsep(x)))/Isep
cat(sprintf("   I_sep/I_0 = %.5f (paper 4.4610),  R_sep = %.5f against R_0 = %.5f\n",
            Isep/I2, Rsep, Rh))
cat(sprintf("   mass ratio (I_0/I_sep)^{2/5} = %.5f -> M_1 <= %.2f PeV (paper 270.30)\n",
            (I2/Isep)^0.4, 491.6*(I2/Isep)^0.4))
cat(sprintf("   t_dec ~ M_1^{-1}(R I)^{-2/3}: the separable world decoheres %.3f times as fast,\n",
            ( (Rsep*Isep)/(Rh*I2) )^(2/3) * (491.6*(I2/Isep)^0.4)/491.6 ))
cat("   so the same measurement that would witness entanglement also fixes which clock ran.\n")

cat("\n==========================================================================\n")
cat(" G.  ONLY ONE OBSERVABLE IN THE WHOLE CONSTRUCTION IS PARTICLE-HOLE ODD\n")
cat("==========================================================================\n\n")
cat("   Under n -> 1-n:  every Renyi entropy fixed (paper), the min-entropy fixed (paper),\n")
cat("   the concurrence C = 2 sqrt(n(1-n)) fixed (paper).  Check what is NOT:\n\n")
for (nn in c(0.05,0.2,0.4)) {
  cat(sprintf("   n = %.2f vs %.2f:  number %8.4f vs %8.4f   |   h = -log(1-n) %8.4f vs %8.4f\n",
      nn, 1-nn, nn, 1-nn, -log(1-nn), -log(nn)))
}
cat("\n   The number density (and the energy, linear in it) is the only functional in the paper\n")
cat("   that distinguishes the band's two ends, and the decoherence exponent -log(1-n) is odd\n")
cat("   too once the max is dropped.  So the MEASURED DARK-MATTER ABUNDANCE is the unique\n")
cat("   thing that breaks the particle-hole degeneracy every entanglement measure is blind to.\n")
