# Universality: is there a fold at every point, or only on the fold's own fixed locus?
#
# Section 3.6 supplies Jacobson's temperature and his entropy law but assumes his universality:
# a local Rindler horizon with the same eta at EVERY point in EVERY null direction. The fold's own
# fixed surfaces are particular loci -- a bifurcation surface, a conjugate locus -- and a generic
# point is not on one, so the assumption is not free. This settles it.
#
# THE OBJECT. The local fold at p is the geodesic symmetry I_p: exp_p(v) -> exp_p(-v). Its
# differential is -Id, which is the frame-independent object section 3.6 already relies on. The
# question is whether I_p is an ISOMETRY, since the Bisognano-Wichmann half-period argument needs
# one. Cartan's theorem says it is an isometry at every point exactly when grad R = 0, so in a
# generic spacetime it is not. What decides the matter is the ORDER at which it fails, set against
# the order the Clausius balance is computed at.
#
# THE PREDICTION, which this measures rather than quotes. In normal coordinates
#   g_mn(x) = d_mn - (1/3) R_manb x^a x^b - (1/6) grad_c R_manb x^a x^b x^c + O(x^4),
# even at second order and odd at third, so
#   G_mn(v) - G_mn(-v) = (1/3) grad_c R_manb v^a v^b v^c + O(v^5),
# cubic in |v| with a coefficient that is grad R at p and nothing else. In two dimensions
# R_manb = K(d_mn d_ab - d_mb d_an), which reduces the prediction to the closed form
#   G(v) - G(-v) = (1/3) (grad K . v) ( |v|^2 d_mn - v_m v_n ),
# checked below to four figures. Two dimensions suffice for the exponent and the coefficient:
# the normal-coordinate expansion does not know the dimension. Part 2 does the tensor parities
# in four, where the physics lives.

# ---------------------------------------------------------------- conformally flat 2D machinery
# g = exp(2 phi) delta. G^x_xx = phi_x, G^x_xy = phi_y, G^x_yy = -phi_x,
#                      G^y_xx = -phi_y, G^y_xy = phi_x, G^y_yy = phi_y.
geo_rhs <- function(s, st) {
  d <- s$dphi(st[1], st[2]); px <- d[1]; py <- d[2]; vx <- st[3]; vy <- st[4]
  c(vx, vy,
    -(px*vx*vx + 2*py*vx*vy - px*vy*vy),
    -(-py*vx*vx + 2*px*vx*vy + py*vy*vy))
}
expmap <- function(s, v, n = 3000) {              # exp_p(v) from the origin, unit affine parameter
  st <- c(0, 0, v[1], v[2]); h <- 1/n
  for (i in 1:n) {
    k1 <- geo_rhs(s, st); k2 <- geo_rhs(s, st + h/2*k1)
    k3 <- geo_rhs(s, st + h/2*k2); k4 <- geo_rhs(s, st + h*k3)
    st <- st + h/6*(k1 + 2*k2 + 2*k3 + k4)
  }
  st[1:2]
}
straight <- function(s, v, n = 0) v               # the deliberately wrong map, for section 5

chart_metric <- function(s, v, mp = expmap, del = 1e-6, n = 3000) {
  x <- mp(s, v, n); J <- matrix(0, 2, 2)
  for (i in 1:2) { e <- c(0,0); e[i] <- del; J[, i] <- (mp(s, v+e, n) - mp(s, v-e, n))/(2*del) }
  exp(2*s$phi(x[1], x[2])) * (t(J) %*% J)
}
residual <- function(s, v, ...) max(abs(chart_metric(s, v, ...) - chart_metric(s, -v, ...)))

# ---------------------------------------------------------------- the surfaces
# A: generic. Quadratic part sets K(p), cubic part sets grad K(p).
A <- list(phi  = function(x,y) 0.15*(x^2+y^2) + 0.11*x^3 + 0.07*x^2*y - 0.05*x*y^2 + 0.09*y^3,
          dphi = function(x,y) c(0.30*x + 0.33*x^2 + 0.14*x*y - 0.05*y^2,
                                 0.30*y + 0.07*x^2 - 0.10*x*y + 0.27*y^2),
          gradK = c(-0.56, -0.68))                # K = -e^{-2phi}(0.60 + 0.56x + 0.68y), grad at p
# B: the round sphere. K = 1 everywhere, so grad R = 0 and Cartan says exact.
B <- list(phi  = function(x,y) -log(1 + (x^2+y^2)/4),
          dphi = function(x,y) c(-x/2, -y/2)/(1 + (x^2+y^2)/4),
          gradK = c(0, 0))
# C: the sharp case. Its cubic part x^3 - 3xy^2 is HARMONIC, so grad K(p) = 0 while K is not
# constant and the surface is not symmetric. The chart still has a visible cubic. If the residual
# tracks grad R at p rather than the chart, C's cubic must vanish and the wrong map's must not.
C <- list(phi  = function(x,y) 0.15*(x^2+y^2) + 0.11*(x^3 - 3*x*y^2) + 0.06*x^4 + 0.04*x^3*y - 0.03*y^4,
          dphi = function(x,y) c(0.30*x + 0.33*x^2 - 0.33*y^2 + 0.24*x^3 + 0.12*x^2*y,
                                 0.30*y - 0.66*x*y + 0.04*x^3 - 0.12*y^3),
          gradK = c(0, 0))
curv <- function(s, x=0, y=0, h=1e-4) {
  lap <- (s$phi(x+h,y)+s$phi(x-h,y)+s$phi(x,y+h)+s$phi(x,y-h)-4*s$phi(x,y))/h^2
  -exp(-2*s$phi(x,y))*lap
}

cat("=== 1. the three surfaces: K at p, K elsewhere, and grad K at p ===\n")
cat("                        K(p)     K(0.3,0.2)   K constant?   |grad K(p)|\n")
for (nm in c("A generic","B round sphere","C harmonic cubic")) {
  s <- get(substr(nm,1,1)); k0 <- curv(s); k1 <- curv(s,0.3,0.2)
  cat(sprintf("   %-18s %8.5f  %11.5f   %-11s  %.5f\n", nm, k0, k1,
              ifelse(abs(k1-k0)<1e-6,"yes","no"), sqrt(sum(s$gradK^2))))
}
cat("   C is the one that matters: K varies, so C is not a symmetric space, but grad K\n")
cat("   vanishes AT p. Chart-level cubic present, geometric cubic absent.\n")

cat("\n=== 2. the closed-form prediction for the cubic, before measuring it ===\n")
th <- 0.7; dir <- c(cos(th), sin(th))
predict_c3 <- function(s, dir) {
  Mx <- diag(2) - outer(dir, dir)
  abs(sum(s$gradK*dir)) * max(abs(Mx)) / 3
}
cat(sprintf("   direction (cos %.1f, sin %.1f); grad K . dir = %.6f for A\n", th, th, sum(A$gradK*dir)))
cat(sprintf("   predicted residual / |v|^3 = (1/3)|grad K . dir| max|d - dd| = %.6f\n", predict_c3(A,dir)))

cat("\n=== 3. measured ===\n")
report <- function(s, nm, vs, mp = expmap) {
  r <- sapply(vs, function(t) residual(s, t*dir, mp = mp))
  cat(sprintf("   %s\n        |v|        residual       residual/|v|^3\n", nm))
  for (i in seq_along(vs)) cat(sprintf("     %8.4f   %13.3e   %14.6f\n", vs[i], r[i], r[i]/vs[i]^3))
  ok <- r > 1e-9
  if (sum(ok) >= 3) { e <- unname(coef(lm(log(r[ok]) ~ log(vs[ok])))[2])
    cat(sprintf("     fitted exponent %.4f\n\n", e)); list(e = e, c3 = r[length(r)]/vs[length(vs)]^3)
  } else { cat("     residual at numerical noise throughout\n\n"); list(e = NA, c3 = 0) }
}
vs <- c(0.24, 0.17, 0.12, 0.085, 0.06, 0.0425, 0.03)
rA <- report(A, "A generic", vs)
rB <- report(B, "B round sphere", vs)
rC <- report(C, "C harmonic cubic, grad K(p) = 0", c(0.30, 0.24, 0.19, 0.15, 0.12))

cat("=== 4. reading it ===\n")
cat(sprintf("   A: exponent %.4f against 3, coefficient %.6f against the predicted %.6f,\n",
            rA$e, rA$c3, predict_c3(A, dir)))
cat(sprintf("      agreeing to %.2e. The cubic term IS grad R at p, measured not assumed.\n",
            abs(rA$c3 - predict_c3(A,dir))))
cat("   B: no residual at any order. Cartan's theorem in the direction that is easy to check.\n")
cat(sprintf("   C: exponent %.4f, so the cubic is gone even though the chart has one. K varies,\n", rC$e))
cat("      so C is not a symmetric space; what the residual tracks is grad R at p alone.\n")
stopifnot(abs(rA$e - 3) < 0.05,
          abs(rA$c3 - predict_c3(A, dir)) < 1e-3,
          is.na(rB$e),
          rC$e > 4.5)

cat("\n=== 5. the plant: hand the same surface a map that is not the exponential ===\n")
cat("   Replace exp_p(v) by the straight line v in the conformal chart. C's chart has a visible\n")
cat("   cubic, so a test that were reading the chart rather than the geometry must now find one.\n")
rP <- report(C, "C, straight-line map", c(0.30,0.24,0.19,0.15,0.12), mp = straight)
cat(sprintf("   with the wrong map the exponent is %.4f, not %.4f. The plant fires: the geodesic\n",
            rP$e, rC$e))
cat("   symmetry and the chart parity are different maps and the test tells them apart.\n")
stopifnot(abs(rP$e - 3) < 0.3)

cat("\n=== 6. and the finite differencing is not what is being measured ===\n")
for (d in c(1e-5, 1e-6, 1e-7)) cat(sprintf("      step %8.0e  ->  residual %.8e\n", d, residual(A, 0.12*dir, del = d)))
cat("      stable across two decades of step size.\n")

cat("\n=== 7. what this buys ===\n")
cat("   At every point of every spacetime the local fold is an isometry through second order in\n")
cat("   the normal coordinates, and its first failure is third order with coefficient grad R.\n")
cat("   The Clausius balance at p is built from R_kk, T_kk and the transverse area element, all\n")
cat("   of them second-order data. So the fold acts, as a symmetry, to exactly the order the\n")
cat("   balance is computed at, at every point and in every null direction. The failure enters\n")
cat("   one order further out, where the curvature-squared corrections live.\n")

cat("\n=== 8. is second order good enough? the balance is taken in a limit ===\n")
cat("   Jacobson's balance on a patch of size L has both sides of order L^2 with coefficients\n")
cat("   R_kk and T_kk. The fold's error is the cubic, order L^3. What decides the argument is not\n")
cat("   that the cubic exists but its RATIO to the term the balance uses, which must vanish as\n")
cat("   the patch shrinks. Both coefficients have closed forms in two dimensions:\n")
cat("      curvature term  G(v) - d   = -(K/3)( |v|^2 d_mn - v_m v_n ),  so c2 = |K|/3 max|d - dd|\n")
cat("      fold error   G(v) - G(-v)  = (1/3)(grad K . v)( |v|^2 d_mn - v_m v_n ),  c3 as in [2]\n")
Mx <- max(abs(diag(2) - outer(dir, dir)))
c2p <- abs(curv(A))/3 * Mx
c3p <- predict_c3(A, dir)
cat(sprintf("\n   predicted c2 = %.6f, c3 = %.6f, so the ratio must approach %.6f x |v|\n", c2p, c3p, c3p/c2p))
cat("\n        |v|     curvature term   c2 measured    fold error    c3 measured     ratio/|v|\n")
Ls <- c(0.12, 0.085, 0.06, 0.0425, 0.03, 0.02)
m2 <- m3 <- rr <- numeric(length(Ls))
for (i in seq_along(Ls)) {
  v <- Ls[i]*dir
  g2 <- max(abs(chart_metric(A, v) - diag(2))); g3 <- residual(A, v)
  m2[i] <- g2/Ls[i]^2; m3[i] <- g3/Ls[i]^3; rr[i] <- (g3/g2)/Ls[i]
  cat(sprintf("     %8.4f   %14.3e  %11.6f  %12.3e  %11.6f   %10.6f\n", Ls[i], g2, m2[i], g3, m3[i], rr[i]))
}
# c2 measured carries its own O(|v|) contamination from the cubic, so it is extrapolated
# rather than read off the smallest point. c3 does not need it: the next correction is O(|v|^2).
rich <- function(f, h) { n <- length(h); f[n] - (f[n-1] - f[n])/(h[n-1] - h[n]) * h[n] }
c2x <- rich(m2, Ls); c3x <- rich(m3, Ls); rrx <- rich(rr, Ls)
cat(sprintf("\n   extrapolated to |v| = 0:\n"))
cat(sprintf("      c2         %.6f   against the predicted %.6f   (%.1e)\n", c2x, c2p, abs(c2x-c2p)))
cat(sprintf("      c3         %.6f   against the predicted %.6f   (%.1e)\n", c3x, c3p, abs(c3x-c3p)))
cat(sprintf("      ratio/|v|  %.6f   against the predicted %.6f   (%.1e)\n", rrx, c3p/c2p, abs(rrx-c3p/c2p)))
cat("   Both coefficients are finite and the ratio is linear in |v| by construction, so the\n")
cat("   fold's error against the term the balance uses vanishes as the patch shrinks. Jacobson's\n")
cat("   field equation is that limit, and in it the fold is an exact symmetry of what the balance\n")
cat("   sees. Universality holds.\n")
stopifnot(abs(c2x - c2p) < 1e-4, abs(c3x - c3p) < 1e-4, abs(rrx - c3p/c2p) < 1e-3)

cat("\n=== 9. the four-dimensional parities, since -Id acts on tensors by rank ===\n")
cat("   dTheta = -Id at the fixed point, so a rank-n tensor at p picks up (-1)^n. The balance\n")
cat("   is fold-invariant if and only if every object in it has even rank. Checked on a random\n")
cat("   Riemann tensor carrying all its symmetries, a random stress tensor and a random null k.\n")
set.seed(90261)
bi <- rbind(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))       # bivector index pairs, 4D
Q <- matrix(rnorm(36), 6, 6); Q <- (Q + t(Q))/2               # 21 free, pair-symmetric
Riem <- array(0, c(4,4,4,4))
for (I in 1:6) for (J in 1:6) {
  a <- bi[I,1]; b <- bi[I,2]; c1 <- bi[J,1]; d <- bi[J,2]; q <- Q[I,J]
  Riem[a,b,c1,d] <-  q; Riem[b,a,c1,d] <- -q
  Riem[a,b,d,c1] <- -q; Riem[b,a,d,c1] <-  q
}
cyc <- Riem[1,2,3,4] + Riem[1,3,4,2] + Riem[1,4,2,3]          # first Bianchi, one condition in 4D
lev <- array(0, c(4,4,4,4))
pm <- function(p) { s <- 1; for (i in 1:3) for (j in (i+1):4) if (p[i] > p[j]) s <- -s; s }
for (p in as.matrix(expand.grid(1:4,1:4,1:4,1:4))) NULL
idx <- as.matrix(expand.grid(a=1:4,b=1:4,c=1:4,d=1:4))
for (r in 1:nrow(idx)) { p <- as.integer(idx[r,]); if (length(unique(p))==4) lev[p[1],p[2],p[3],p[4]] <- pm(p) }
Riem <- Riem - (cyc/3)*lev
cat(sprintf("      first Bianchi residual after projection: %.2e\n",
            abs(Riem[1,2,3,4] + Riem[1,3,4,2] + Riem[1,4,2,3])))
eta4 <- diag(c(-1,1,1,1))
Ric <- matrix(0,4,4)
for (m in 1:4) for (n in 1:4) Ric[m,n] <- sum(sapply(1:4, function(r) sum(sapply(1:4, function(s) eta4[r,s]*Riem[r,m,s,n]))))
Tab <- matrix(rnorm(16),4,4); Tab <- (Tab + t(Tab))/2
S <- -diag(4)                                                  # dTheta = -Id
push2 <- function(M) t(S) %*% M %*% S
kk <- c(1,1,0,0)/sqrt(2)                                       # null: -1+1 = 0
cat(sprintf("      k is null: k.k = %.2e\n", as.numeric(t(kk) %*% eta4 %*% kk)))
cat(sprintf("      R_kk  %+.8f -> %+.8f   invariant\n", t(kk)%*%Ric%*%kk, t(kk)%*%push2(Ric)%*%kk))
cat(sprintf("      T_kk  %+.8f -> %+.8f   invariant\n", t(kk)%*%Tab%*%kk, t(kk)%*%push2(Tab)%*%kk))
area <- diag(4)[3:4,3:4]
cat(sprintf("      transverse area element  det %+.6f -> %+.6f   invariant\n",
            det(area), det((-diag(2)) %*% area %*% (-diag(2)))))
chi <- c(0.3, 0.9, 0, 0)                                       # the local boost field, a vector
cat(sprintf("      boost field chi   (%.1f,%.1f) -> (%.1f,%.1f)   FLIPS, rank 1\n",
            chi[1], chi[2], (S%*%chi)[1], (S%*%chi)[2]))
cat("      so kappa flips with it and beta = 2 pi / |kappa| does not move. The flip is the\n")
cat("      wedge swap itself: the modular generator is H_R - H_L, odd by construction.\n")
cat(sprintf("      grad R, rank 5, picks up %+d: an odd-rank object DOES flip, so the\n", (-1)^5))
cat("      invariance above is a fact about which tensors the balance contains, not a triviality\n")
cat("      of the map. The balance holds R_kk, T_kk and an area: all even.\n")
stopifnot(abs(t(kk)%*%Ric%*%kk - t(kk)%*%push2(Ric)%*%kk) < 1e-12,
          abs(t(kk)%*%Tab%*%kk - t(kk)%*%push2(Tab)%*%kk) < 1e-12,
          max(abs(S %*% chi + chi)) < 1e-12)

cat("\n=== 10. the verdict ===\n")
cat("   Universality was assumed. It is now derived, with one order to spare. At every point of\n")
cat("   every spacetime the geodesic symmetry is a fold: differential -Id, an isometry through\n")
cat("   second order, failing at third with coefficient grad R. The Clausius balance uses only\n")
cat("   second-order, even-rank data, and its fractional sensitivity to the fold's failure goes\n")
cat("   to zero linearly in the patch size, which is the limit the balance is taken in.\n")
cat("   What remains assumed is not universality but the two structural facts the construction\n")
cat("   stands on: a Lorentzian metric and local orthonormal frames.\n")
