# B.1 says the bulk divergence of the bang-adiabatic state has log coefficient gamma^2/16 pi^2
# and that it "is reproduced at eta = 0.5, 1, 2, 5 as 0.006330". OPEN_FOR_BEN_20260922.md
# section 3 filed that number as having no script anywhere in the repo. This is the script.
#
# The half sweep leaves n(p) -> gamma^2/(16 p^4), so in the measure the paper uses throughout,
# pi^-2 int p^2 (...) dp, the energy integral pi^-2 int p^3 n dp has integrand -> A/p with
#
#     A = lim_{p->inf} p^4 n(p) / pi^2 = gamma^2 / (16 pi^2)
#
# and therefore diverges logarithmically with coefficient A. Computing A directly is the same
# check as fitting the log and does not need a cutoff. gamma = 1 throughout, as elsewhere.

gamma <- 1

rk4 <- function(p, eta0, eta1, psi0, nstep) {
  h <- (eta1 - eta0) / nstep
  f <- function(eta, psi)
    rbind(-1i*(gamma*eta*psi[1,] + p*psi[2,]),
          -1i*(p*psi[1,] - gamma*eta*psi[2,]))
  psi <- psi0; eta <- eta0
  for (s in seq_len(nstep)) {
    k1 <- f(eta,     psi); k2 <- f(eta + h/2, psi + h/2*k1)
    k3 <- f(eta+h/2, psi + h/2*k2); k4 <- f(eta + h, psi + h*k3)
    psi <- psi + h/6*(k1 + 2*k2 + 2*k3 + k4); eta <- eta + h
  }
  psi
}

# instantaneous eigenvectors of H(eta) = [[g eta, p], [p, -g eta]]
eigvec <- function(eta, p, upper) {
  E <- sqrt((gamma*eta)^2 + p^2) * if (upper) 1 else -1
  v <- rbind(p + 0i, E - gamma*eta)
  sweep(v, 2, sqrt(colSums(Mod(v)^2)), "/")
}

# occupation of a mode started at eta = 0 in the lower adiabatic state and evolved to eta
n_of <- function(p, eta, per_osc = 240) {
  psi <- eigvec(0, p, FALSE)
  nstep <- max(2000, ceiling(per_osc * eta * sqrt(p^2 + (gamma*eta)^2)))
  psi <- rk4(p, 0, eta, psi, nstep)
  up <- eigvec(eta, p, TRUE)
  Mod(colSums(Conj(up) * psi))^2
}

pred <- gamma^2 / (16*pi^2)
cat(sprintf("=== predicted log coefficient  gamma^2/(16 pi^2) = %.7f\n\n", pred))
cat("     eta        p=8        p=12       p=16       p=20     mean(12,16,20)\n")
best <- c()
for (eta in c(0.5, 1, 2, 5)) {
  ps <- c(8, 12, 16, 20)
  A <- sapply(ps, function(p) p^4 * n_of(p, eta) / pi^2)
  m <- mean(A[2:4]); best <- c(best, m)
  cat(sprintf("  %6.1f %10.6f %10.6f %10.6f %10.6f %12.6f\n", eta, A[1], A[2], A[3], A[4], m))
}

cat(sprintf("\n  across the four epochs: %.6f to %.6f, against the predicted %.7f\n",
            min(best), max(best), pred))
cat("\n=== flatly, and this is a NEGATIVE result\n\n")
cat("  The finite-eta reading does NOT converge to gamma^2/(16 pi^2). It scatters over a\n")
cat("  factor of ten, and it should: at eta = 0.5 to 5 a mode of momentum 8 to 20 is still\n")
cat("  inside its own crossing, whose width in conformal time is eta ~ p/gamma, so the\n")
cat("  adiabatic occupation there oscillates and has no settled value to read off.\n\n")
cat("  So B.1's 'reproduced at eta = 0.5, 1, 2, 5 as 0.006330' could not be reproduced here,\n")
cat("  and the sentence has been rewritten rather than carried. What IS established, and by\n")
cat("  a computation already in the paper, is the asymptotic tail: lz_bang.R gives\n")
cat("  n_B/(gamma^2/16 p^4) falling from 1.114 to 1.0010 across p/sqrt(gamma) = 2 to 6, so\n")
cat("  the coefficient approaches gamma^2/(16 pi^2) = 0.0063326 from above. B.1 now says that.\n\n")
cat("  Kept rather than deleted because a failed reproduction is the reason the text changed,\n")
cat("  and a reader who wonders where 0.006330 went should find this file.\n")
