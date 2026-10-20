## verify_v6_amplitudes.R
## s-channel graviton exchange, L = -h_{mu nu} T^{mu nu} / M_P, de Donder propagator:
##   M = [T_i^{mu nu} T_f,mu nu - (1/2) T_i T_f] / (M_P^2 s)
## with stress-tensor matrix elements built from explicit helicity spinors and polarisation
## vectors. Amplitudes are in units of 1/M_P^2 and spin sums in units of s^2/M_P^4 (s = 1).
## Metric (+,-,-,-); Dirac matrices in the chiral basis, psi = (psi_L, psi_R).

ETA <- diag(c(1, -1, -1, -1))
I2 <- diag(2) + 0i
Z2 <- matrix(0 + 0i, 2, 2)
SIGMA <- list(matrix(c(0, 1, 1, 0), 2) + 0i,
              matrix(c(0, 1i, -1i, 0), 2),
              matrix(c(1, 0, 0, -1), 2) + 0i)
blk <- function(A, B, C, D) rbind(cbind(A, B), cbind(C, D))
GAM <- c(list(blk(Z2, I2, I2, Z2)), lapply(SIGMA, function(s) blk(Z2, s, -s, Z2)))  # gamma^mu
G5 <- 1i * GAM[[1]] %*% GAM[[2]] %*% GAM[[3]] %*% GAM[[4]]

mdot <- function(a, b) sum(a * (ETA %*% b))                  # a.b with the metric (no conjugation)
slash <- function(p) { pl <- as.vector(ETA %*% p); Reduce(`+`, Map(`*`, GAM, pl)) }
dbar <- function(w) Conj(w) %*% GAM[[1]]                      # row vector w^dagger gamma^0

## direction n = (sin th cos ph, sin th sin ph, cos th)
nvec <- function(th, ph) c(sin(th) * cos(ph), sin(th) * sin(ph), cos(th))
mom <- function(E, th, ph) E * c(1, nvec(th, ph))

## two-component helicity spinors, sigma.n xi_pm = +/- xi_pm
xi_p <- function(th, ph) c(cos(th / 2), exp(1i * ph) * sin(th / 2))
xi_m <- function(th, ph) c(-exp(-1i * ph) * sin(th / 2), cos(th / 2))

## massless Dirac spinors of helicity h for momentum E(1, n(th, ph))
u_sp <- function(E, th, ph, h)
  sqrt(2 * E) * (if (h > 0) c(0, 0, xi_p(th, ph)) else c(xi_m(th, ph), 0, 0))
v_sp <- function(E, th, ph, h)
  sqrt(2 * E) * (if (h > 0) c(xi_m(th, ph), 0, 0) else c(0, 0, -xi_p(th, ph)))

## helicity polarisation vectors of a massless vector boson moving along n(th, ph)
eps_h <- function(th, ph, lam) {
  e_th <- c(cos(th) * cos(ph), cos(th) * sin(ph), -sin(th))
  e_ph <- c(-sin(ph), cos(ph), 0)
  c(0, (e_th + 1i * lam * e_ph) / sqrt(2))
}

## ---- stress-tensor matrix elements, upper indices -------------------------------------------

## two real scalars with momenta k1, k2 (both incoming or both outgoing); xi is the coupling to
## curvature, whose improvement term is -2 xi (P P - eta P^2), P = k1 + k2
T_scalar <- function(k1, k2, xi = 0) {
  P <- k1 + k2
  outer(k1, k2) + outer(k2, k1) - ETA * mdot(k1, k2) - 2 * xi * (outer(P, P) - ETA * mdot(P, P)) + 0i
}

## fermion (p1, helicity h1) and antifermion (p2, h2) annihilated:
##   (1/4) vbar(p2) [gamma^mu k^nu + gamma^nu k^mu] u(p1), k = p1 - p2
T_dirac_in <- function(E, th1, ph1, h1, th2, ph2, h2) {
  p1 <- mom(E, th1, ph1); p2 <- mom(E, th2, ph2); k <- p1 - p2
  vb <- dbar(v_sp(E, th2, ph2, h2)); u <- u_sp(E, th1, ph1, h1)
  J <- vapply(GAM, function(g) as.complex(vb %*% g %*% u), 0i)
  0.25 * (outer(J, k) + outer(k, J))
}
## fermion (k3, h3) and antifermion (k4, h4) created:
##   (1/4) ubar(k3) [gamma^mu k^nu + gamma^nu k^mu] v(k4), k = k3 - k4
T_dirac_out <- function(E, th3, ph3, h3, th4, ph4, h4) {
  k3 <- mom(E, th3, ph3); k4 <- mom(E, th4, ph4); k <- k3 - k4
  ub <- dbar(u_sp(E, th3, ph3, h3)); v <- v_sp(E, th4, ph4, h4)
  J <- vapply(GAM, function(g) as.complex(ub %*% g %*% v), 0i)
  0.25 * (outer(J, k) + outer(k, J))
}

## two vector bosons, from the field strengths f^{mu a} = p^mu e^a - p^a e^mu:
##   T = f1^{mu a} f2^nu_a + f2^{mu a} f1^nu_a - (1/2) eta^{mu nu} f1^{ab} f2_ab
fstr <- function(p, e) outer(p, e) - outer(e, p)
T_vector_fs <- function(f1, f2) {
  f1 %*% ETA %*% t(f2) + f2 %*% ETA %*% t(f1) - 0.5 * ETA * sum(f1 * (ETA %*% f2 %*% ETA))
}
T_vector <- function(p1, e1, p2, e2) T_vector_fs(fstr(p1, e1), fstr(p2, e2))

## de Donder contraction, in units of 1/M_P^2
trT <- function(T) sum(diag(ETA %*% T))
amp <- function(Ti, Tf, s = 1) (sum(Ti * (ETA %*% Tf %*% ETA)) - 0.5 * trT(Ti) * trT(Tf)) / s

## ---- kinematics and spin sums ---------------------------------------------------------------
## Centre of mass, s = 1: incoming p1 along +z, p2 along -z; outgoing k3 along n(theta, phi),
## k4 opposite. c = cos(theta).

initial_T <- function(ch, xi = 0) {
  E <- 0.5
  p1 <- mom(E, 0, 0); p2 <- mom(E, pi, 0)
  switch(ch,
    scalar = list(T_scalar(p1, p2, xi)),
    dirac  = { L <- list(); for (h1 in c(1, -1)) for (h2 in c(1, -1))
                 L[[length(L) + 1]] <- T_dirac_in(E, 0, 0, h1, pi, 0, h2); L },
    vector = { L <- list(); for (l1 in c(1, -1)) for (l2 in c(1, -1))
                 L[[length(L) + 1]] <- T_vector(p1, eps_h(0, 0, l1), p2, eps_h(pi, 0, l2)); L })
}

final_T <- function(dm, c, phi = 0.3, xi = 0) {
  E <- 0.5; th <- acos(c)
  k3 <- mom(E, th, phi); k4 <- mom(E, pi - th, phi + pi)
  switch(dm,
    scalar = list(T_scalar(k3, k4, xi)),
    dirac  = { L <- list(); for (h3 in c(1, -1)) for (h4 in c(1, -1))
                 L[[length(L) + 1]] <- T_dirac_out(E, th, phi, h3, pi - th, phi + pi, h4); L },
    vector = { L <- list(); for (l3 in c(1, -1)) for (l4 in c(1, -1))
                 L[[length(L) + 1]] <- T_vector(k3, Conj(eps_h(th, phi, l3)), k4,
                                                Conj(eps_h(pi - th, phi + pi, l4))); L })
}

## sum over all spins of |M|^2, units s^2/M_P^4
spin_sum <- function(ch, dm = "dirac", c, phi = 0.3, xi_in = 0, xi_out = 0) {
  Ti <- initial_T(ch, xi_in); Tf <- final_T(dm, c, phi, xi_out)
  tot <- 0
  for (a in Ti) for (b in Tf) tot <- tot + Mod(amp(a, b))^2
  tot
}

## integral over the solid angle of the final pair, by Gauss-Legendre in cos(theta) (the spin
## sums are polynomials of degree 4 in cos(theta), so 8 nodes are exact) times 2 pi
GL8 <- gauss_legendre(8)
ang_int <- function(f) 2 * pi * sum(GL8$w * vapply(GL8$x, f, 0))

## ---- Einstein-Cartan contact term -----------------------------------------------------------
## L = -(3/16) A.A / M_P^2 with A the axial current of all fermions; the cross term between a
## Standard Model Dirac fermion and N gives the amplitude (3/8) J5_in . J5_out (units 1/M_P^2).
contact_and_graviton <- function(c, phi = 0.3, coupling = 3 / 8) {
  E <- 0.5; th <- acos(c)
  out <- list()
  for (h1 in c(1, -1)) for (h2 in c(1, -1)) {
    vb <- dbar(v_sp(E, pi, 0, h2)); u <- u_sp(E, 0, 0, h1)
    J5i <- vapply(GAM, function(g) as.complex(vb %*% g %*% G5 %*% u), 0i)
    Ti <- T_dirac_in(E, 0, 0, h1, pi, 0, h2)
    for (h3 in c(1, -1)) for (h4 in c(1, -1)) {
      ub <- dbar(u_sp(E, th, phi, h3)); v <- v_sp(E, pi - th, phi + pi, h4)
      J5f <- vapply(GAM, function(g) as.complex(ub %*% g %*% G5 %*% v), 0i)
      Tf <- T_dirac_out(E, th, phi, h3, pi - th, phi + pi, h4)
      out[[length(out) + 1]] <- c(contact = coupling * mdot(J5i, J5f), graviton = amp(Ti, Tf))
    }
  }
  do.call(rbind, out)
}
