# ==========================================================================================
# An editorial pass objected that the absolute 1.417e-32 s "depends on gamma = M_1 a_1, which
# the paper does not name", and asked for the power law alone. Check whether gamma enters.
#
# It does not, and the reason is that the threshold ELIMINATES H. From 3.5,
#     Gammabar_H = (4/(3 sqrt(pi))) c_G R I (M_1/H)^{3/2}
# and setting Gammabar_H = 1 solves for H rather than requiring it as an input:
#     M_1/H = [ (4/(3 sqrt(pi))) c_G R I ]^{-2/3}
# In radiation domination t = 1/(2H), so
#     t_dec = (1/(2 M_1)) [ (4/(3 sqrt(pi))) c_G R I ]^{-2/3},
# which contains M_1, c_G, R and I and nothing else. All four are quoted in the paper.
hbar_GeV_s <- 6.582119569e-25
M1  <- 491.6e6      # GeV, from 2.3
R   <- 1.07037      # from 4.2
I   <- 0.0127597    # from 2.2
cG  <- 1            # one channel per (p,h) mode

pref <- 4/(3*sqrt(pi))
inner <- pref*cG*R*I
t_dec <- 0.5 * (hbar_GeV_s/M1) * inner^(-2/3)

cat(sprintf("  4/(3 sqrt pi)          = %.8f\n", pref))
cat(sprintf("  (4/(3 sqrt pi)) c_G R I= %.8f\n", inner))
cat(sprintf("  [ ... ]^(-2/3)         = %.6f\n", inner^(-2/3)))
cat(sprintf("  hbar / M_1             = %.6e s\n", hbar_GeV_s/M1))
cat(sprintf("\n  t_dec = %.6e s   (paper quotes 1.417e-32 s)\n", t_dec))
cat(sprintf("  agreement: %.3f per cent\n", 100*abs(t_dec-1.417e-32)/1.417e-32))
cat(sprintf("\n  and with c_G = 1/2, a Majorana pair counted once: %.6e s (paper 2.249e-32)\n",
            0.5*(hbar_GeV_s/M1)*(pref*0.5*R*I)^(-2/3)))
cat("\n  FLATLY: gamma = M_1 a_1 does not appear. The threshold condition Gammabar_H = 1 is\n")
cat("  what determines H, so H is an OUTPUT and not an input, and the absolute time follows\n")
cat("  from M_1, c_G, R and I alone. The objection mistook a quantity that was solved for\n")
cat("  for one that had to be supplied.\n")
cat("\n  What IS a convention, and 3.5 already says so, is the order-unity threshold itself\n")
cat("  and c_G. Those scale the answer as c_G^(-2/3) and are declared.\n")
