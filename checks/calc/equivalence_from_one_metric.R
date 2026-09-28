# Is the equivalence principle an assumption of this construction, or a consequence of one?
#
# The relativity ledger lists "matter following the metric's geodesics" as assumed, which reads as
# a physical postulate about gravity standing alongside the fold. Matter conservation read that
# way too until Noether's second theorem was applied to it, and the same question is worth asking
# here, because what the Clausius step actually uses is narrower than the ledger's phrasing.
#
# WHAT THE CLAUSIUS STEP USES. It needs a stress tensor T_ab to put in the flux, and T_ab means
# the variational derivative of the matter action with respect to the metric. So the input is that
# the matter action is a functional of the SAME metric whose horizons the argument is built on.
# That is the metric postulate. This asks what follows from it.
#
# WHAT FOLLOWS. Write a matter field as A exp(i S / hbar) and let hbar go to zero. The leading
# eikonal equation is g^{ab} d_a S d_b S + m^2 = 0, whose characteristics are the integral curves
# of dx^a/dlambda = g^{ab} p_b with dp_a/dlambda = -(1/2) d_a g^{bc} p_b p_c. Those are Hamilton's
# equations for H = (g^{ab} p_a p_b + m^2)/2, and the claim is that they are the geodesic equation
# of g written another way. If so, the weak equivalence principle is what the metric postulate
# gives once the matter is on shell, and it is not a separate input. Checked below, with the mass
# and the curvature coupling varied to test composition independence, and with the matter action
# built from a DIFFERENT metric as the plant.

set.seed(5501)
N <- 4
eta <- diag(c(-1, 1, 1, 1))
Pc <- array(rnorm(N*N*N), c(N,N,N)); Qc <- array(rnorm(N*N*N*N), c(N,N,N,N))
for (a in 1:N) for (b in 1:N) { Pc[a,b,] <- Pc[b,a,] <- (Pc[a,b,] + Pc[b,a,])/2
                                Qc[a,b,,] <- Qc[b,a,,] <- (Qc[a,b,,] + Qc[b,a,,])/2 }
gmet <- function(x, s = 0.10) {
  h <- matrix(0,N,N)
  for (a in 1:N) for (b in 1:N) h[a,b] <- sum(Pc[a,b,]*x) + sum(outer(x,x)*Qc[a,b,,])
  eta + s*h
}
ginv <- function(x, s = 0.10) solve(gmet(x, s))
h <- 1e-5
d_ginv <- function(x, a, s = 0.10) { e <- rep(0,N); e[a] <- h; (ginv(x+e,s) - ginv(x-e,s))/(2*h) }
d_g     <- function(x, a, s = 0.10) { e <- rep(0,N); e[a] <- h; (gmet(x+e,s) - gmet(x-e,s))/(2*h) }
Gam <- function(x, s = 0.10) {
  gi <- ginv(x,s); dg <- lapply(1:N, function(c) d_g(x,c,s)); G <- array(0,c(N,N,N))
  for (a in 1:N) for (b in 1:N) for (c in 1:N)
    G[a,b,c] <- sum(sapply(1:N, function(d) gi[a,d]*(dg[[c]][d,b] + dg[[b]][d,c] - dg[[d]][b,c])))/2
  G
}

# --- the eikonal characteristics, which is what the WAVE equation gives
hj_rhs <- function(y, s = 0.10) {
  x <- y[1:N]; p <- y[(N+1):(2*N)]
  dx <- as.numeric(ginv(x,s) %*% p)
  dp <- -0.5 * sapply(1:N, function(a) as.numeric(t(p) %*% d_ginv(x,a,s) %*% p))
  c(dx, dp)
}
# --- the geodesic equation, which is what a TEST BODY is supposed to do
geo_rhs <- function(y, s = 0.10) {
  x <- y[1:N]; v <- y[(N+1):(2*N)]; G <- Gam(x,s)
  c(v, -sapply(1:N, function(a) sum(outer(v,v) * G[a,,])))
}
rk4 <- function(rhs, y0, lam, n, s = 0.10) {
  y <- y0; hh <- lam/n
  for (i in 1:n) { k1 <- rhs(y,s); k2 <- rhs(y + hh/2*k1,s); k3 <- rhs(y + hh/2*k2,s); k4 <- rhs(y + hh*k3,s)
    y <- y + hh/6*(k1 + 2*k2 + 2*k3 + k4) }
  y
}

x0 <- c(0.10, -0.20, 0.15, 0.05)
v0 <- c(1.00,  0.30, -0.20, 0.10)                     # a 4-velocity
p0 <- as.numeric(gmet(x0) %*% v0)                     # the matching momentum, p_a = g_ab v^b

cat("=== 1. the wave's characteristics against the test body's path ===\n")
cat("   Same start, same initial direction, one integrated from the eikonal equation of the\n")
cat("   field and one from the geodesic equation of the metric.\n\n")
cat("      lambda      max |x_eikonal - x_geodesic|\n")
for (L in c(0.2, 0.5, 1.0)) {
  ye <- rk4(hj_rhs, c(x0, p0), L, 1500)
  yg <- rk4(geo_rhs, c(x0, v0), L, 1500)
  cat(sprintf("      %6.2f        %.3e\n", L, max(abs(ye[1:N] - yg[1:N]))))
}
ye <- rk4(hj_rhs, c(x0,p0), 1.0, 1500); yg <- rk4(geo_rhs, c(x0,v0), 1.0, 1500)
stopifnot(max(abs(ye[1:N] - yg[1:N])) < 1e-8)
cat("   They are the same curve to the integrator's accuracy. What the field's short-wavelength\n")
cat("   limit follows is the geodesic of the metric its action was built from.\n")

cat("\n=== 2. composition independence, which is the content of the weak principle ===\n")
cat("   The mass sits in the eikonal equation as g^{ab} d_a S d_b S + m^2 = 0, so it fixes the\n")
cat("   normalisation of p and nothing else. Two fields of different mass released with the same\n")
cat("   4-velocity must trace the same curve.\n\n")
cat("      m        max |x(m) - x(m = 1)|     (after rescaling to the same 4-velocity)\n")
base <- rk4(geo_rhs, c(x0, v0), 1.0, 1500)[1:N]
for (m in c(0.5, 1, 2, 50)) {
  pm <- m * as.numeric(gmet(x0) %*% v0)               # p = m g v, the physical momentum
  yy <- rk4(function(y,s) { r <- hj_rhs(y,s); c(r[1:N]/m, r[(N+1):(2*N)]/m) }, c(x0, pm), 1.0, 1500)
  cat(sprintf("   %7.1f        %.3e\n", m, max(abs(yy[1:N] - base))))
  stopifnot(max(abs(yy[1:N] - base)) < 1e-8)
}
cat("   Independent of the mass, and the reason should be said plainly rather than dressed as a\n")
cat("   discovery: with p = m g v, dividing the characteristic equations by m returns dx/dlambda\n")
cat("   = v exactly, so the cancellation is algebraic and the identical figures above confirm\n")
cat("   the implementation rather than establish anything new. What is NOT algebraic is section\n")
cat("   1, where a first-order Hamiltonian system in (x, p) built from the inverse metric is\n")
cat("   compared against a second-order system in x built from the Christoffels, and the two\n")
cat("   agree. The content of the weak principle is that the mass appears in the first system\n")
cat("   and in neither the second nor the path.\n")
cat("   The curvature coupling xi sits one order down in hbar from the eikonal term, entering as\n")
cat("   xi R phi, so it is absent from the characteristics by inspection of the expansion and no\n")
cat("   two couplings can separate at this order. That is stated, not computed here.\n")

cat("\n=== 3. the plant: build the matter action from a DIFFERENT metric ===\n")
cat("   If the derivation is using the shared metric and nothing else, then a matter action\n")
cat("   built from g-tilde must send the rays along g-tilde's geodesics and away from g's.\n\n")
st <- 0.16                                            # a different metric of the same family
p0t <- as.numeric(gmet(x0, st) %*% v0)
yet <- rk4(function(y, s) hj_rhs(y, st), c(x0, p0t), 1.0, 1500)
ygt <- rk4(function(y, s) geo_rhs(y, st), c(x0, v0), 1.0, 1500)
cat(sprintf("   rays of the tilde action against tilde geodesics:  %.3e  (they agree)\n",
            max(abs(yet[1:N] - ygt[1:N]))))
cat(sprintf("   rays of the tilde action against g's geodesics:    %.3e  (they do not)\n",
            max(abs(yet[1:N] - base))))
stopifnot(max(abs(yet[1:N] - ygt[1:N])) < 1e-8, max(abs(yet[1:N] - base)) > 1e-3)
cat("   The plant fires. The geodesics matter follows are the geodesics of whichever metric its\n")
cat("   action was written with, so the content of the assumption is that there is ONE metric\n")
cat("   and the matter action uses it.\n")

cat("\n=== 4. what the ledger should say ===\n")
cat("   The equivalence principle is not a separate postulate here. What is assumed is the\n")
cat("   metric postulate, that a single metric carries both the horizons the fold acts on and\n")
cat("   the matter action whose variation defines T_ab, and the Clausius step already needs that\n")
cat("   because T_ab has no other meaning. Matter following the metric's geodesics is what the\n")
cat("   short-wavelength limit of any field on that metric does, independent of its mass and of\n")
cat("   its curvature coupling. That is the weak principle, and it is a consequence.\n")
cat("   The ledger's assumed list keeps the same number of entries and gains accuracy: the entry\n")
cat("   is one shared metric, which is what 'metric theory' means, and not a physical principle\n")
cat("   about gravity standing beside the fold.\n")
