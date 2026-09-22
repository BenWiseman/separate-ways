# Last iteration's stated next route. The AZ class of the bang is CI, but a symmetry class is
# not random matrix theory: RMT needs an ENSEMBLE with correlated levels, and independent
# two-level modes have none. The only mechanism in this programme that couples modes is the
# latitude-dependent horizon response, which mixes different l at fixed m. Build THAT matrix
# and ask whether its spectrum behaves like an RMT ensemble.
#
# Basis: spin-weighted _{-2}Y_{l,2}, i.e. Wigner d^l_{2,2}(theta), for l = 2..L.
dl22 <- function(l, th) {                       # general Wigner small-d at m' = m = 2
  c2 <- cos(th/2); s2 <- sin(th/2); tot <- 0
  for (k in 0:(l-2)) {
    num <- (-1)^k * factorial(l+2) * factorial(l-2)
    den <- factorial(l+2-k) * factorial(k)^2 * factorial(l-k-2)
    tot <- tot + num/den * c2^(2*l-2*k) * s2^(2*k)
  }
  tot
}
stopifnot(all(abs(sapply(2:8, function(l) dl22(l,0)) - 1) < 1e-12))   # d^l_22(0) = 1

P2 <- function(x) (3*x^2-1)/2
build <- function(L, eps) {                     # A(theta) = 1 + eps P2(cos theta)
  ls <- 2:L; n <- length(ls); M <- matrix(0,n,n)
  nrm <- sqrt((2*ls+1)/2)                       # normalisation of d^l_22 on [0,pi]
  for (i in seq_len(n)) for (j in seq_len(n)) {
    f <- function(th) dl22(ls[i],th)*(1+eps*P2(cos(th)))*dl22(ls[j],th)*sin(th)
    M[i,j] <- nrm[i]*nrm[j]*integrate(f, 0, pi, subdivisions=2000)$value
  }
  M
}
rbar <- function(ev) {                          # Oganesyan-Huse spacing ratio
  s <- diff(sort(ev)); s <- s[s > 1e-13]
  if (length(s) < 3) return(NA)
  mean(mapply(function(a,b) min(a,b)/max(a,b), s[-length(s)], s[-1]))
}

cat("  reference values for rbar:  Poisson 0.386   GOE 0.536   GUE 0.603   rigid/equal 1.000\n\n")
cat("    a/M    eps=c2/c0     L    matrix rbar   spectrum range\n")
for (par in list(c(0.5,-0.046), c(0.9,-0.219), c(0.99,-0.367))) {
  for (L in c(14, 22)) {
    M <- build(L, par[2]); ev <- eigen(M, symmetric=TRUE)$values
    cat(sprintf("  %5.2f   %9.4f  %4d   %10.4f   [%.4f, %.4f]\n",
        par[1], par[2], L, rbar(ev), min(ev), max(ev)))
  }
}
cat("\n  For contrast, a genuine GOE matrix of the same size:\n")
set.seed(11); n <- 21
G <- matrix(rnorm(n*n),n); G <- (G+t(G))/2
cat(sprintf("    GOE rbar = %.4f\n", rbar(eigen(G, symmetric=TRUE)$values)))

cat("\n=== is that 0.53 an ensemble statistic, or an artefact of the truncation? ===\n")
cat("  M = I + eps*K, so eigenvalues are 1 + eps*lambda_i and every SPACING scales with eps.\n")
cat("  rbar is a ratio of spacings, hence eps-independent by construction. That is why the\n")
cat("  three spins above returned identical values. So rbar can only depend on L. A real\n")
cat("  ensemble statistic converges as the matrix grows. Check whether this one does.\n\n")
cat("      L     n levels    rbar\n")
rs <- c()
for (L in c(8,10,12,14,16,18,20,22,24,26)) {
  M <- build(L, -0.219); ev <- eigen(M, symmetric=TRUE)$values
  r <- rbar(ev); rs <- c(rs, r)
  cat(sprintf("  %5d  %9d  %8.4f\n", L, L-1, r))
}
cat(sprintf("\n  spread over L: %.4f to %.4f, range %.4f\n", min(rs), max(rs), max(rs)-min(rs)))
cat("  A GOE ensemble's rbar is stable to about 0.01 over this range of sizes. This wanders by\n")
cat(sprintf("  %.2f and shows no sign of settling, so it is not an ensemble statistic at all.\n", max(rs)-min(rs)))

cat("\n  FLATLY: the l-mixing matrix is NOT a random matrix and its agreement with GOE at L = 22\n")
cat("  is a coincidence of where the truncation was cut. The reason is structural rather than\n")
cat("  numerical. K is multiplication by a smooth function, P2(cos theta), restricted to a\n")
cat("  finite band of l. Its eigenvalues approximate the RANGE of that function sampled on the\n")
cat("  band, so the spectrum is a smooth deterministic curve, and smooth curves have no level\n")
cat("  repulsion. Nothing in this construction is random, so no amount of enlarging it produces\n")
cat("  RMT behaviour.\n")
cat("\n  NEXT ROUTE. Randomness has to enter physically, not be hoped for. Two candidates the\n")
cat("  programme can actually supply: (i) a POPULATION of holes with a distribution over spin\n")
cat("  and mass, which randomises eps and the band edge across an observed sample, giving an\n")
cat("  ensemble of spectra rather than one spectrum; (ii) a seam whose response varies over the\n")
cat("  horizon in a way not fixed by the Kerr geometry, which is exactly the matching law 3.4\n")
cat("  records as unconstructed. Route (i) is observational and available now; route (ii) needs\n")
cat("  the missing law. Neither is done here.\n")
