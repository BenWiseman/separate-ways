## verify_v6_modes.R
## Mode equation near a radiation bang, in the paper's variables:
##   i d(psi)/dt = (t sigma_z + x sigma_x) psi,   t = sqrt(M a1) eta,   x = p / sqrt(M a1).
## Fourth-order Runge-Kutta, vectorised over modes, with steps h = c / sqrt(t^2 + xmax^2) so that
## the fastest phase turns by about 2c per step.

rk4_modes <- function(x, p1, p2, t0, t1, c = 0.02) {
  xm <- max(x, 1); t <- t0; dir <- sign(t1 - t0); nstep <- 0L
  while (dir * (t1 - t) > 1e-12) {
    h <- dir * min(c / sqrt(t^2 + xm^2), abs(t1 - t))
    th <- t + h / 2; tt <- t + h
    k1a <- -1i * (t * p1 + x * p2);   k1b <- -1i * (x * p1 - t * p2)
    a <- p1 + (h / 2) * k1a;          b <- p2 + (h / 2) * k1b
    k2a <- -1i * (th * a + x * b);    k2b <- -1i * (x * a - th * b)
    a <- p1 + (h / 2) * k2a;          b <- p2 + (h / 2) * k2b
    k3a <- -1i * (th * a + x * b);    k3b <- -1i * (x * a - th * b)
    a <- p1 + h * k3a;                b <- p2 + h * k3b
    k4a <- -1i * (tt * a + x * b);    k4b <- -1i * (x * a - tt * b)
    p1 <- p1 + (h / 6) * (k1a + 2 * k2a + 2 * k3a + k4a)
    p2 <- p2 + (h / 6) * (k1b + 2 * k2b + 2 * k3b + k4b)
    t <- tt; nstep <- nstep + 1L
  }
  list(p1 = p1, p2 = p2, t = t, nstep = nstep)
}

## Instantaneous eigenvectors of t sigma_z + x sigma_x (eigenvalues +E excited, -E ground).
excited <- function(t, x) { th <- atan2(x, t); cbind(cos(th / 2), sin(th / 2)) }
ground  <- function(t, x) { th <- atan2(x, t); cbind(-sin(th / 2), cos(th / 2)) }

## Amplitude in the excited level at time t, with the first-order adiabatic correction
## <e~| = <e| - i (x / 4E^3) <g|, which removes the 1/t^3 ringing of a finite end time.
excited_amp <- function(t, x, a, b) {
  e <- excited(t, x); g <- ground(t, x); E <- sqrt(t^2 + x^2)
  (e[, 1] * a + e[, 2] * b) - 1i * (x / (4 * E^3)) * (g[, 1] * a + g[, 2] * b)
}
## Ground level at time t to first order, |g~> = |g> + i (x / 4E^3) |e>, normalised.
ground_state <- function(t, x) {
  e <- excited(t, x); g <- ground(t, x); E <- sqrt(t^2 + x^2); d <- x / (4 * E^3)
  a <- g[, 1] + 1i * d * e[, 1]; b <- g[, 2] + 1i * d * e[, 2]; nn <- sqrt(Mod(a)^2 + Mod(b)^2)
  list(a = a / nn, b = b / nn)
}

## Half crossing from the bang: evolve the two basis states (1,0) and (0,1) from t = 0 to t = T and
## project onto the late excited level. For psi(0) = (1, e^{i phi}) / sqrt(2) the occupation is
## n(phi) = |u1 + e^{i phi} u2|^2 / 2, whose least value over phi is (|u1| - |u2|)^2 / 2, whose
## largest is (|u1| + |u2|)^2 / 2, and whose value at phi = pi (the massless vacuum) is |u1 - u2|^2 / 2.
half_crossing <- function(x, T = 20, c = 0.02) {
  n <- length(x)
  r <- rk4_modes(c(x, x), c(rep(1 + 0i, n), rep(0 + 0i, n)), c(rep(0 + 0i, n), rep(1 + 0i, n)), 0, T, c)
  ea <- excited_amp(T, c(x, x), r$p1, r$p2)
  u1 <- ea[seq_len(n)]; u2 <- ea[n + seq_len(n)]
  list(x = x, u1 = u1, u2 = u2, nmin = 0.5 * (Mod(u1) - Mod(u2))^2, nmax = 0.5 * (Mod(u1) + Mod(u2))^2,
       nvac = 0.5 * Mod(u1 - u2)^2, nstep = r$nstep)
}

## Full sweep from t = -T to +T, started in the ground level: excitation probability.
full_sweep <- function(x, T = 20, c = 0.03) {
  g <- ground_state(-T, x)
  r <- rk4_modes(x, g$a, g$b, -T, T, c)
  list(x = x, P = Mod(excited_amp(T, x, r$p1, r$p2))^2, nstep = r$nstep)
}

## ---- closed forms ----------------------------------------------------------------------------

## half-angle state, n_* = (1 - sqrt(1 - P)) / 2 written without cancellation
n_star <- function(x) { u <- exp(-pi * x^2); u / (2 * (1 + sqrt(1 - u))) }

## complex log-gamma for Re z > 0: Stirling series after shifting z up by 12
lgamma_c <- function(z) {
  N <- 12; w <- z + N
  s <- (w - 0.5) * log(w) - w + 0.5 * log(2 * pi) + 1 / (12 * w) - 1 / (360 * w^3) +
       1 / (1260 * w^5) - 1 / (1680 * w^7) + 1 / (1188 * w^9)
  for (k in 0:(N - 1)) s <- s - log(z + k)
  s
}
## massless vacuum at the bang, the Gamma-function form the paper quotes for this half crossing:
## n = [1 - sqrt(1 - P) cos(vartheta)] / 2, vartheta = pi/4 + arg Gamma(1/2 + i x^2/4) - arg Gamma(1 + i x^2/4)
n_vac_cf <- function(x) {
  y <- x^2 / 4
  vt <- pi / 4 + Im(lgamma_c(0.5 + 1i * y)) - Im(lgamma_c(1 + 1i * y))
  0.5 * (1 - sqrt(1 - exp(-pi * x^2)) * cos(vt))
}

## ---- production integral to mass (Section 3.2) ------------------------------------------------
## M^(5/2) = (rho_DM / s0) (2 pi^2 / 45) g*^(1/4) 15^(3/4) mu^(3/2) / I
mass_from_I <- function(I, rho = K("rho_DM"), s0 = K("s0_cm3") * K("hbarc")^3, g = K("gstar"), mu = MU)
  ((rho / s0) * (2 * pi^2 / 45) * g^0.25 * 15^0.75 * mu^1.5 / I)^0.4
