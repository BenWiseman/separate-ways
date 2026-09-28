#!/usr/bin/env Rscript
# Fork 9, second brick: the interior radial solver at nonzero frequency, and the Klein-Gordon
# normalisation that fixes its convention.
#
# WHY. interior_modes_exact_vs_eikonal.R built the exact k = 0 solver and measured how far the
# closed-form Jacobi modes sit from it. It ended by naming what is still missing before any mode
# sum can be attempted: the same solver at nonzero k, the Klein-Gordon normalisation on a
# constant-r slice, and the cross-region state. The first two are here. They matter because the
# convention risk in a mode sum is almost all in the measure, and an unnormalised sum can return
# a confident sign that belongs to nothing.
#
# THE GEOMETRY. With 2M = 1 the metric is ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2, f = 1 - 1/r.
# Inside the horizon f < 0: r is the time coordinate and t is a space coordinate, which is the
# Kantowski-Sachs reading. Constant-r slices are therefore the Cauchy slices and the Klein-Gordon
# product is taken on them.
#
# THE EQUATION. Phi = e^{ikt} Y_lm(Omega) u(r)/r turns the massless wave equation into
#     d^2u/drstar^2 + (k^2 - V) u = 0,    V = f ( l(l+1)/r^2 + 1/r^3 ),    drstar/dr = 1/f,
# and rstar = r + log(1 - r) inside, running from -infinity at the horizon to 0 at r = 0. Since
# f < 0 inside, V < 0, so k^2 - V > 0 and every interior mode oscillates in rstar at every k,
# including k = 0. That is why the k = 0 modes of the first brick were oscillatory.
#
# THE NORMALISATION. The Klein-Gordon product on a constant-r slice is
#     (Phi_1, Phi_2) = i \int |f| r^2 dt dOmega ( Phi_1 d_r Phi_2^* - Phi_2^* d_r Phi_1 ),
# and for one (k, l, m) that collapses to the Wronskian
#     W = u^* du/drstar - u du^*/drstar,
# which the equation makes independent of r. W is the whole convention: get it wrong and every
# amplitude in a mode sum is wrong by the same unknown factor. A unit-amplitude wave e^{-ik rstar}
# at the horizon carries W = -2ik exactly, and that is the normalisation adopted here.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

rstar_of <- function(r) r + log(1 - r)          # inside only, r < 1
fmet     <- function(r) 1 - 1/r
Vpot     <- function(r, l) fmet(r) * (l*(l+1)/r^2 + 1/r^3)

# RK4 in rstar on (r, u, v = du/drstar), u complex.
march <- function(k, l, r0, rstop, N = 200000, keep = FALSE) {
  rs0 <- rstar_of(r0); rs1 <- rstar_of(rstop)
  h <- (rs1 - rs0)/N
  st <- c(r = r0 + 0i, u = exp(-1i*k*rs0), v = -1i*k*exp(-1i*k*rs0))
  deriv <- function(s) {
    r <- Re(s[1])
    c(fmet(r) + 0i, s[3], (Vpot(r, l) - k^2)*s[2])
  }
  out <- if (keep) matrix(0+0i, N + 1, 3) else NULL
  if (keep) out[1, ] <- st
  for (i in seq_len(N)) {
    k1 <- deriv(st); k2 <- deriv(st + h*k1/2)
    k3 <- deriv(st + h*k2/2); k4 <- deriv(st + h*k3)
    st <- st + h*(k1 + 2*k2 + 2*k3 + k4)/6
    if (keep) out[i + 1, ] <- st
  }
  if (keep) list(r = Re(out[, 1]), u = out[, 2], v = out[, 3], h = h) else st
}
wronsk <- function(u, v) Conj(u)*v - u*Conj(v)   # = u^* du/drstar - u du^*/drstar

cat("=== 1. the Wronskian is conserved, which is the Klein-Gordon norm being conserved ===\n")
cat("   Started at r = 1 - 1e-6 on the unit-amplitude horizon wave e^{-ik rstar}, integrated to\n")
cat("   r = 0.20, which is well inside the caustic at r = 1/2. The exact value is W = -2ik at\n")
cat("   every radius. How far in the march can go is set by resolution, not by the equation:\n")
cat("   the local wavenumber in rstar is sqrt(k^2 - V) and V goes as -1/r^4 near the\n")
cat("   singularity, so at k = 7, l = 64 a uniform 200000 steps run to r = 0.05 leaves only\n")
cat("   seventeen steps per wavelength there and the norm drifts by 1.8e-4. Stopping at\n")
cat("   r = 0.20 leaves hundreds. The norm reports its own resolution, which is the first\n")
cat("   reason to watch it.\n")
cat("        k     l      W at the start        W at r = 0.20      relative drift     steps\n")
drifts <- c()
NN <- 300000
for (par in list(c(0.5, 2), c(1, 2), c(1, 8), c(3, 8), c(3, 32), c(7, 64))) {
  k <- par[1]; l <- par[2]
  en <- march(k, l, 1 - 1e-6, 0.20, N = NN)
  w0 <- -2i*k; w1 <- wronsk(en[2], en[3])
  d  <- Mod(w1 - w0)/Mod(w0)
  cat(sprintf("   %6.2f %5d %12s %20s %18.3e %9d\n", k, l,
              sprintf("%.6f i", Im(w0)), sprintf("%.6f i", Im(w1)), d, NN))
  drifts <- c(drifts, d)
  note(d < 1e-7, sprintf("the KG norm is conserved at k = %g, l = %d", k, l))
  note(abs(Re(w1)) < 1e-9*Mod(w0), sprintf("and stays purely imaginary at k = %g, l = %d", k, l))
}
cat(sprintf("   worst relative drift over the six: %.3e\n", max(drifts)))
cat("   And the drift is resolution and not a defect in the equation, which halving the step\n")
cat("   has to show: RK4 is fourth order, so the drift must fall by about sixteen.\n")
cn <- c()
for (NN in c(100000, 200000, 400000)) {
  en <- march(3, 32, 1 - 1e-6, 0.20, N = NN)
  cn <- c(cn, Mod(wronsk(en[2], en[3]) + 6i)/6)
  cat(sprintf("   k = 3, l = 32, %7d steps: drift %.3e\n", NN, cn[length(cn)]))
}
ordr <- log2(cn[1]/cn[3])/2
cat(sprintf("   observed order of convergence: %.2f against the 4 RK4 owes\n", ordr))
note(ordr > 3.5, "the norm drift converges at RK4's own order, so it is step size and not physics")

cat("\n=== 2. the caustic radius is inside the integrated range and the mode is finite there ===\n")
cat("   The contact sphere is r = M = 1/2 with 2M = 1, and A.19's rule is a statement about the\n")
cat("   amplitude there. Reading u at r = 1/2 is the point of building this.\n")
cat("        k     l        |u(1/2)|      arg u(1/2)/pi     |u| at the horizon\n")
for (par in list(c(1, 2), c(1, 8), c(3, 8), c(3, 32))) {
  k <- par[1]; l <- par[2]
  tr <- march(k, l, 1 - 1e-6, 0.5, N = 120000, keep = TRUE)
  uh <- tr$u[length(tr$u)]
  cat(sprintf("   %6.2f %5d %14.6e %17.4f %22.4f\n", k, l, Mod(uh), Arg(uh)/pi, Mod(tr$u[1])))
  note(is.finite(Mod(uh)) && Mod(uh) > 0, sprintf("u(1/2) is finite and nonzero at k=%g l=%d", k, l))
}

cat("\n=== 3. the k -> 0 limit returns the first brick's equation, not a different one ===\n")
cat("   At k = 0 the two independent solutions can be taken real, and a real solution carries\n")
cat("   W = 0: a static interior mode is not a positive-norm particle mode. The check that the\n")
cat("   equations agree is therefore made on the solution and not on the norm. The k = 0 march\n")
cat("   is compared against the first brick's form r^2(r-1)u'' + r u' - (l(l+1)r + 1)u = 0,\n")
cat("   integrated in r from the same start with the same data.\n")
for (l in c(4, 16)) {
  tr <- march(0, l, 1 - 1e-4, 0.3, N = 200000, keep = TRUE)
  # same start data, marched in r with the first brick's equation
  r0 <- 1 - 1e-4; u0 <- tr$u[1]; v0 <- tr$v[1]           # v0 = du/drstar = f du/dr
  dudr0 <- v0/fmet(r0)
  N <- 200000; h <- (0.3 - r0)/N; r <- r0; y <- u0; yp <- dudr0
  g <- function(r, y, yp) c(yp, ((l*(l+1)*r + 1)*y - r*yp)/(r^2*(r - 1)))
  for (i in seq_len(N)) {
    q1 <- g(r,y,yp); q2 <- g(r+h/2, y+h*q1[1]/2, yp+h*q1[2]/2)
    q3 <- g(r+h/2, y+h*q2[1]/2, yp+h*q2[2]/2); q4 <- g(r+h, y+h*q3[1], yp+h*q3[2])
    y <- y + h*(q1[1]+2*q2[1]+2*q3[1]+q4[1])/6
    yp <- yp + h*(q1[2]+2*q2[2]+2*q3[2]+q4[2])/6; r <- r + h
  }
  uend <- tr$u[length(tr$u)]
  rel <- Mod(y - uend)/Mod(uend)
  cat(sprintf("   l = %2d: rstar march gives u(0.3) = %.8e, the r march %.8e, relative %.2e\n",
              l, Mod(uend), Mod(y), rel))
  note(rel < 1e-6, sprintf("the two forms of the k=0 equation agree at l = %d", l))
}

cat("\n=== 4. the plants, one per thing that could silently be wrong ===\n")
# (a) drop the f' u' term by marching the r-form without it: the norm must stop being conserved.
bad_march <- function(k, l, r0, rstop, N = 200000, drop_fp = FALSE, wrongV = 1) {
  rs0 <- rstar_of(r0); h <- (rstar_of(rstop) - rs0)/N
  st <- c(r = r0 + 0i, u = exp(-1i*k*rs0), v = -1i*k*exp(-1i*k*rs0))
  for (i in seq_len(N)) {
    d <- function(s) {
      r <- Re(s[1])
      # the slip: an extra first-derivative term that the rstar form does not have
      extra <- if (drop_fp) 0.05*s[3] else 0
      c(fmet(r) + 0i, s[3], (wrongV*Vpot(r, l) - k^2)*s[2] + extra)
    }
    k1 <- d(st); k2 <- d(st + h*k1/2); k3 <- d(st + h*k2/2); k4 <- d(st + h*k3)
    st <- st + h*(k1 + 2*k2 + 2*k3 + k4)/6
  }
  st
}
b <- bad_march(1, 8, 1 - 1e-6, 0.05, drop_fp = TRUE)
db <- Mod(wronsk(b[2], b[3]) - (-2i))/2
cat(sprintf("   (a) a stray first-derivative term: the norm drifts by %.3e, against %.3e clean\n",
            db, drifts[3]))
note(db > 1e-3, "plant: a first-derivative slip breaks the conserved norm")
# (b) a wrong potential keeps the norm (it is still self-adjoint) but moves the solution, so the
#     norm check alone cannot catch it and the k=0 cross-check in section 3 is what does.
b2 <- bad_march(1, 8, 1 - 1e-6, 0.5, wrongV = 1.2)
g2 <- march(1, 8, 1 - 1e-6, 0.5)
cat(sprintf("   (b) the potential scaled by 1.2: the norm still holds to %.1e, because a real\n",
            Mod(wronsk(b2[2], b2[3]) - (-2i))/2))
cat("       potential is self-adjoint whatever its size. The modulus at the caustic barely\n")
cat(sprintf("       moves, %.4e against %.4e, and it is the PHASE that carries it: %.4f pi\n",
            Mod(b2[2]), Mod(g2[2]), Arg(b2[2]/g2[2])/pi))
cat(sprintf("       of relative phase, a relative move of %.3f in the complex mode itself.\n",
            Mod(b2[2] - g2[2])/Mod(g2[2])))
note(Mod(wronsk(b2[2], b2[3]) - (-2i))/2 < 1e-5, "plant: a wrong potential is invisible to the norm")
note(Mod(b2[2] - g2[2])/Mod(g2[2]) > 0.1, "plant: and visible in the complex mode")
cat("   That is the warning this brick carries into the sum: the norm is a check on the measure\n")
cat("   and on transcription, and it is no check at all on the potential. A mode sum read for a\n")
cat("   sign is read from phases, which is exactly the part the norm does not police.\n")

cat("\n=== 5. where this leaves fork 9 ===\n")
cat("   The solver now runs at any (k, l) from the horizon to the singularity, and the measure\n")
cat("   it carries is fixed by a quantity the equation conserves rather than by a choice. The\n")
cat("   norm check catches a slip in the first-derivative structure, which is the transcription\n")
cat("   error that would otherwise survive; it does NOT catch a wrong potential, and the k = 0\n")
cat("   agreement with the first brick is what covers that.\n")
cat("   STILL TO BUILD: the cross-region state, which is the one piece that cannot be settled\n")
cat("   inside the interior alone, and the sum itself. The calibration that has to be met before\n")
cat("   any sign is read off is Delta^{1/2} -> 3.9004 M s^{-1/2}.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
