rho_DM <- 3.3406e10          # Msun / Mpc^3, comoving
n_lrd  <- 1e-4; M_seed <- 1e5
f_lrd  <- n_lrd*M_seed/rho_DM
cat(sprintf("f at the observed LRD seed budget      : %.3e\n", f_lrd))
cat(sprintf("margin to the TOTAL budget   (f=1)     : %.2f orders\n", log10(1/f_lrd)))
cat(sprintf("margin to the SHARP ceiling  (f=0.0101): %.2f orders\n", log10(0.0101/f_lrd)))
cat(sprintf("tightest PBH accretion limit 3e-9      : %.2f orders ABOVE the LRD seed budget\n",
            log10(3e-9/f_lrd)))
cat("\n  So observation does NOT exclude PBH seeds at the observed LRD abundance:\n")
cat(sprintf("  they need f = %.1e and the tightest limit allows 3e-9, one order higher.\n", f_lrd))
cat("  The defensible claim is that PBHs cannot MOVE THE ENDPOINT, not that the seeds\n")
cat("  are astrophysical.\n")
