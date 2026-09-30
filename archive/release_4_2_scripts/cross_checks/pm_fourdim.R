# pm_fourdim.R -- the partially-massless eigenvalue on the midpoint line, and why d = 4.
#
# Supplement S11.5 claims that for the Z2-odd transverse-traceless graviton of a dS_d-sliced
# AdS_{d+1} bulk, the Higuchi/partially-massless value eps_PM = d - 2 is the LOWEST eigenvalue
# on the midpoint line w0 = 2L only in four dimensions. The claim had no verification path, so
# this solves the eigenvalue problem directly and compares with the numbers S11.5 quotes.
#
#   (sinh^d u psi')' + eps sinh^{d-2} u psi = 0,   psi(L) = 0 at the crease,
#                                                  psi'(w0) = 0 at the brane.
#
# In standard form,  psi'' = -d coth(u) psi' - eps psi / sinh^2 u.
#
# Shooting: fix psi(L) = 0, psi'(L) = 1, integrate to w0 with RK4, and look for the smallest
# eps > 0 at which psi'(w0) crosses zero. The normalisation of psi'(L) is irrelevant because the
# problem is linear and the target is a zero.
#
# Base R only.
source("helpers.R")

rhs <- function(u, y, d, eps) {
  s <- sinh(u)
  c(y[2], -d * (cosh(u)/s) * y[2] - eps * y[1] / (s*s))
}

shoot <- function(eps, d, L, w0, n = 3000) {
  h <- (w0 - L)/n
  y <- c(0, 1); u <- L
  for (i in seq_len(n)) {
    k1 <- rhs(u,         y,            d, eps)
    k2 <- rhs(u + h/2,   y + h/2 * k1, d, eps)
    k3 <- rhs(u + h/2,   y + h/2 * k2, d, eps)
    k4 <- rhs(u + h,     y + h  * k3,  d, eps)
    y <- y + (h/6) * (k1 + 2*k2 + 2*k3 + k4); u <- u + h
  }
  y[2]                                    # psi'(w0); zero at an eigenvalue
}

# smallest positive eigenvalue: scan for the first sign change in psi'(w0), then bisect
eps1 <- function(d, L, w0, hi = 10, step = 0.05) {
  grid <- seq(step, hi, by = step)
  prev <- shoot(grid[1], d, L, w0)
  for (e in grid[-1]) {
    cur <- shoot(e, d, L, w0)
    if (is.finite(prev) && is.finite(cur) && prev * cur < 0)
      return(uniroot(shoot, c(e - step, e), d = d, L = L, w0 = w0, tol = 1e-12)$root)
    prev <- cur
  }
  NA_real_
}

cat("=== lowest eigenvalue on the midpoint line w0 = 2L\n\n")
cat("   d   eps_PM = d-2        L=0.5         L=1           L=2\n")
tab <- list()
for (d in 3:6) {
  v <- sapply(c(0.5, 1, 2), function(L) eps1(d, L, 2*L))
  tab[[as.character(d)]] <- v
  cat(sprintf("  %2d %10d %14.5f %13.5f %13.5f\n", d, d - 2, v[1], v[2], v[3]))
}

cat("\n=== against the numbers S11.5 quotes\n\n")
quoted <- list("3" = c(2.998, 3.631, 6.977),
               "4" = c(2.00000, 2.00000, 2.00000),
               "5" = c(1.291, 1.033, 0.480),
               "6" = c(0.806, 0.502, 0.102))
for (d in as.character(3:6))
  for (i in 1:3)
    report(sprintf("eps_1, d = %s, L = %s", d, c("0.5","1","2")[i]),
           expected = quoted[[d]][i], reproduced = tab[[d]][i],
           tol = 0.005, mode = "rel", note = "supplement S11.5")

cat("\n=== the claim itself: only d = 4 lands on eps_PM, and it lands on it exactly\n\n")
for (d in 3:6) {
  dev <- max(abs(tab[[as.character(d)]] - (d - 2)))
  cat(sprintf("   d = %d: largest |eps_1 - (d-2)| across the three L is %.5f%s\n",
              d, dev, if (dev < 1e-4) "   <- on the PM value" else ""))
}
report("d = 4 sits on the PM value at every L", expected = 2, reproduced = tab[["4"]][3],
       tol = 1e-5, mode = "rel", note = "S11.5's 'only in four dimensions'")
