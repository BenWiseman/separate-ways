# Recompute every figure the 2026-09-20 verification pass flagged as stale or wrong.
# The entropy correction moved the endpoint 484.8 -> 491.6 PeV. Anything derived from
# the endpoint and not recomputed is stale. This does the recomputation independently
# of the reviewer's report, so the reviewer is checked rather than copied.
cat("=== derived from the endpoint move ===\n")
M_old <- 484.8; M_new <- 491.6; k <- M_new/M_old
cat(sprintf("scale factor M_new/M_old = %.8f\n\n", k))

# Decoherence times: at fixed I and R the time goes as 1/M.
for (t in c(1.436, 2.280))
  cat(sprintf("time  %.4f -> %.5f e-32 s\n", t, t/k))

# Mass at fixed mixing fraction f: strictly proportional.
cat("\n")
for (m in c(464.8, 367.4, 193.0))
  cat(sprintf("f-mass %6.1f -> %7.2f PeV   (half %7.2f)\n", m, m*k, m*k/2))

# Particle mass in kg.
kg_per_eV <- 1.78266192e-36
cat(sprintf("\nparticle mass: old %.3e kg, new %.3e kg\n",
            M_old*1e15*kg_per_eV, M_new*1e15*kg_per_eV))

cat("\n=== claims that are wrong independently of the entropy fix ===\n")
# "seven orders of magnitude" between 0.115 and 4e-17
cat(sprintf("0.115 vs 4e-17 spans %.2f orders (paper says seven)\n", log10(0.115/4e-17)))

# pi G M / c^3
GMsun_c3 <- 4.925490947e-6   # s, IAU nominal
for (Ms in c(3, 6.5e9)) {
  s <- pi*GMsun_c3*Ms
  cat(sprintf("pi G M/c^3 at %-8.1e Msun = %.4g s = %.3f hours = %.1f us\n",
              Ms, s, s/3600, s*1e6))
}

# CMB temperature contribution to the mass error
Tt <- 2.7255; sT <- 0.0006
rel_s <- 3*sT/Tt          # s ~ T^3
rel_M <- (2/5)*rel_s      # M ~ s^(-2/5)
cat(sprintf("\nT = %.4f +/- %.4f K -> sigma_s/s = %.4f%%, sigma_M/M = %.4f%% (paper says 0.008%%)\n",
            Tt, sT, 100*rel_s, 100*rel_M))
cat(sprintf("  that is %.4f PeV on the endpoint\n", rel_M*M_new))
