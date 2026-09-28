# Is "one shared metric" a convenience the derivation is written in, or does it need one?
#
# The assumed list is down to four and two of those are numbers. Of the other two, the signature
# has been examined: the fold keeps a parity in exactly the signatures with an odd number of time
# directions, so it has content on one arena and not the others. The metric postulate has not.
# It sits on the list as "one metric carrying both the fold's horizons and a covariant matter
# action" with no computation behind the word "one".
#
# WHAT WOULD BREAK. The Clausius step needs two things that come from different places. It needs
# horizons, which belong to the metric whose null cones define them, and it needs a conserved
# stress tensor, which Noether's second theorem supplies with respect to the metric the matter
# action was written with. If those are two different metrics then conservation holds against the
# wrong connection and the Bianchi step, which is what leaves exactly one integration constant,
# has nothing to close on. This measures that: build the stress tensor from a second metric and
# take its divergence with the first one's connection.

set.seed(6021)
N <- 4
eta <- diag(c(-1, 1, 1, 1))
P1 <- array(rnorm(N*N*N), c(N,N,N)); Q1 <- array(rnorm(N*N*N*N), c(N,N,N,N))
P2 <- array(rnorm(N*N*N), c(N,N,N)); Q2 <- array(rnorm(N*N*N*N), c(N,N,N,N))
for (a in 1:N) for (b in 1:N) {
  P1[a,b,] <- P1[b,a,] <- (P1[a,b,] + P1[b,a,])/2; Q1[a,b,,] <- Q1[b,a,,] <- (Q1[a,b,,] + Q1[b,a,,])/2
  P2[a,b,] <- P2[b,a,] <- (P2[a,b,] + P2[b,a,])/2; Q2[a,b,,] <- Q2[b,a,,] <- (Q2[a,b,,] + Q2[b,a,,])/2
}
mk <- function(P, Q, s) function(x) {
  h <- matrix(0,N,N)
  for (a in 1:N) for (b in 1:N) h[a,b] <- sum(P[a,b,]*x) + sum(outer(x,x)*Q[a,b,,])
  eta + s*h
}
g   <- mk(P1, Q1, 0.10)                       # the metric the horizons belong to
gt  <- mk(P2, Q2, 0.10)                       # a genuinely different metric for the matter action
gc  <- function(x) 1.7 * g(x)                 # and a constant multiple of the first, for contrast

phi <- function(x) 0.7*x[1] - 0.4*x[2]*x[3] + 0.3*sin(1.3*x[4]) + 0.25*x[1]^2*x[2] - 0.15*x[3]^3
mm <- 0.83; h <- 1e-3
d1 <- function(f, x, a, hh = h) { e <- rep(0,N); e[a] <- hh; (f(x+e) - f(x-e))/(2*hh) }

Gam <- function(met, x) {
  gi <- solve(met(x)); dg <- lapply(1:N, function(c) d1(met, x, c))
  G <- array(0,c(N,N,N))
  for (a in 1:N) for (b in 1:N) for (c in 1:N)
    G[a,b,c] <- sum(sapply(1:N, function(d) gi[a,d]*(dg[[c]][d,b] + dg[[b]][d,c] - dg[[d]][b,c])))/2
  G
}
Tdn <- function(met, x) {                      # T_ab for a scalar whose action uses `met`
  m <- met(x); mi <- solve(m)
  dp <- sapply(1:N, function(a) d1(phi, x, a))
  kin <- as.numeric(t(dp) %*% mi %*% dp)
  outer(dp, dp) - m*(kin + mm^2*phi(x)^2)/2
}
# grad_a T^a_b, with the connection from `conn` and the index ALWAYS raised by the action's own
# metric. Raising with `conn` instead makes a constant rescaling look like a broken connection,
# because T^a_b then carries the ratio of the two metrics: the plant in section 3 caught that.
divT <- function(conn, act, x) {
  G <- Gam(conn, x); Tm <- function(y) solve(act(y)) %*% Tdn(act, y)
  Tx <- Tm(x); out <- numeric(N)
  for (b in 1:N) {
    s <- sum(sapply(1:N, function(a) d1(function(y) Tm(y)[a,b], x, a)))
    s <- s + sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[a,a,c]*Tx[c,b]))))
    s <- s - sum(sapply(1:N, function(a) sum(sapply(1:N, function(c) G[c,a,b]*Tx[a,c]))))
    out[b] <- s
  }
  out
}
boxphi <- function(met, x) {
  sq <- function(y) sqrt(-det(met(y)))
  f <- function(a) function(y) sq(y)*sum(solve(met(y))[a,]*sapply(1:N, function(b) d1(phi,y,b)))
  sum(sapply(1:N, function(a) d1(f(a), x, a)))/sq(x)
}
x0 <- c(0.20, -0.30, 0.15, 0.25)

cat("=== 1. the identity holds for each metric against its OWN connection ===\n")
cat("   conservation_from_diffeo.R established grad^a T_ab = (box phi - m^2 phi) d_b phi off\n")
cat("   shell. It holds for whichever metric the action and the divergence share.\n\n")
dp0 <- sapply(1:N, function(b) d1(phi, x0, b))
cat("      metric used for BOTH        max |grad^a T_ab - (box - m^2) d_b phi|\n")
for (nm in list(list("g", g), list("g-tilde", gt), list("1.7 g", gc))) {
  lhs <- divT(nm[[2]], nm[[2]], x0)
  rhs <- (boxphi(nm[[2]], x0) - mm^2*phi(x0)) * dp0
  cat(sprintf("      %-26s  %.3e\n", nm[[1]], max(abs(lhs - rhs))))
  stopifnot(max(abs(lhs - rhs)) < 2e-4 * max(1, max(abs(lhs))))
}
cat("   So each metric conserves its own matter, which is Noether's theorem and nothing more.\n")

cat("\n=== 2. and fails when the two are different metrics ===\n")
cat("   Now build the stress tensor from g-tilde's action and take its divergence with g's\n")
cat("   connection, which is what the Clausius step would be doing if matter and horizons had\n")
cat("   different metrics.\n\n")
mixed <- divT(g, gt, x0)
own   <- divT(gt, gt, x0)
rhs_t <- (boxphi(gt, x0) - mm^2*phi(x0)) * dp0
cat(sprintf("      divergence with g-tilde's own connection, minus its identity:  %.3e\n",
            max(abs(own - rhs_t))))
cat(sprintf("      divergence with g's connection, minus the same identity:       %.4f\n",
            max(abs(mixed - rhs_t))))
cat(sprintf("      against the size of the divergence itself:                     %.4f\n", max(abs(mixed))))
# the bar is the RATIO of the two failures, not an absolute size: a first version compared the
# mismatch against max(1, scale) and so measured nothing when the divergence was below one
ratio <- max(abs(mixed - rhs_t)) / max(abs(own - rhs_t))
cat(sprintf("      so the mismatch is %.3e times the numerical floor the own-connection case sits at\n",
            ratio))
stopifnot(ratio > 1e3, max(abs(mixed - rhs_t)) > 0.05 * max(abs(mixed)))
cat("   Order unity. Conservation holds against one connection and not the other, so on shell\n")
cat("   the stress tensor the Clausius step is handed is not conserved in the geometry whose\n")
cat("   horizons the step is running on, and the Bianchi identity has nothing to close against.\n")

cat("\n=== 3. the plant: a constant multiple must NOT break it ===\n")
cat("   Two metrics differing by a constant factor share a connection, so the failure has to\n")
cat("   vanish there. If it does not, section 2 is measuring the perturbation and not the\n")
cat("   mismatch of connections.\n\n")
mixed_c <- divT(g, gc, x0)
own_c   <- divT(gc, gc, x0)
cat(sprintf("      1.7 g's stress, divergence with g's connection, against its own:  %.3e\n",
            max(abs(mixed_c - own_c))))
cat(sprintf("      the same comparison for a genuinely different metric:             %.4f\n",
            max(abs(mixed - own))))
stopifnot(max(abs(mixed_c - own_c)) < 1e-5 * max(abs(own_c)),
          max(abs(mixed - own)) > 0.05 * max(abs(own)))
cat("   The constant rescaling costs nothing and the genuine second metric costs order unity,\n")
cat("   so what section 2 measures is the connections disagreeing and not the size of a bump.\n")

cat("\n=== 4. what this settles ===\n")
cat("   The metric postulate is not a convenience the derivation happens to be written in. The\n")
cat("   Clausius step needs a stress tensor conserved in the same geometry whose horizons it\n")
cat("   uses, Noether supplies conservation only against the metric the matter action was\n")
cat("   written with, and with two metrics those are different statements. So the step has no\n")
cat("   content unless the two coincide, up to the constant factor that changes no connection.\n")
cat("   That does not derive the postulate: a universe could have two metrics and this argument\n")
cat("   would not forbid it, it would only say the derivation stops. Like the signature, the\n")
cat("   assumption is the condition under which there is anything to derive, which is a better\n")
cat("   thing to be able to say about it than that it was assumed.\n")
