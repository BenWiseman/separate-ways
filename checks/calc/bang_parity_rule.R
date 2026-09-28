# What the fold permits at the bang, and whether Lambda can be fixed there.
#
# The de Sitter budget shows the fold's image stress is zero throughout the observable universe,
# so no part of Lambda is sourced there. The bang is the one place left: it is where the fold acts
# and it is the boundary whose datum Lambda is. This asks what the fold actually constrains there.
#
# THE CONDITION. In conformal time a flat cosmology is ds^2 = a(eta)^2(-deta^2 + dx^2) and the
# cosmological fold is Theta: eta -> -eta. It is an isometry exactly when a(-eta)^2 = a(eta)^2, so
# a is either even or odd. A bang has a(0) = 0, and an even a smooth at the origin also has
# a'(0) = 0, which is a bounce and not a radiation bang. So the fold requires a ODD in eta, which
# is Boyle, Finn and Turok's CPT-symmetric universe reached as a condition rather than a choice.
#
# WHAT THAT FORBIDS. In conformal time the Friedmann equation is a'^2 = (8 pi G/3) rho a^4, and a
# fluid of equation of state w contributes rho_i a^4 proportional to a^(1-3w). Differentiating,
#     a'' = (K/2) sum_i c_i (1 - 3 w_i) a^(-3 w_i).
# If a is odd then a'' is odd, so every term must be an odd power of a. The rule is therefore that
# a fluid is admissible at the bang exactly when -3w is an odd integer, with radiation exempt
# because its coefficient 1 - 3w vanishes. Dust has -3w = 0 and is forbidden. Everything below
# checks that, numerically, on solutions rather than on the series.

G <- 1; K <- 1                                   # units: K = (8 pi G/3) rho_r0 a0^4 set to one
# a'' = (K/2) sum c_i (1 - 3 w_i) a^{-3 w_i},  with the radiation term absent identically
accel <- function(a, comps) {
  s <- 0
  for (cc in comps) s <- s + cc$c * (1 - 3*cc$w) * a^(-3*cc$w)
  (K/2) * s
}
grow <- function(comps, emax = 1.2, n = 400000) { # integrate from the bang both ways
  h <- emax/n; e <- 0; a <- 0; v <- sqrt(K)      # a' -> sqrt(K) as a -> 0, the radiation bang
  fw <- numeric(n+1); fw[1] <- 0
  for (i in 1:n) {
    k1 <- v;                 m1 <- accel(max(a, 1e-300), comps)
    k2 <- v + h/2*m1;        m2 <- accel(max(a + h/2*k1, 1e-300), comps)
    k3 <- v + h/2*m2;        m3 <- accel(max(a + h/2*k2, 1e-300), comps)
    k4 <- v + h*m3;          m4 <- accel(max(a + h*k3, 1e-300), comps)
    a <- a + h/6*(k1 + 2*k2 + 2*k3 + k4); v <- v + h/6*(m1 + 2*m2 + 2*m3 + m4)
    fw[i+1] <- a
  }
  list(eta = seq(0, emax, length.out = n+1), a = fw)
}
# the continuation to eta < 0 is the solution of the same equation run backwards from the bang,
# which for the fold to be an isometry must satisfy a(-eta)^2 = a(eta)^2.
grow_back <- function(comps, emax = 1.2, n = 400000) {
  h <- -emax/n; e <- 0; a <- 0; v <- sqrt(K)
  bw <- numeric(n+1); bw[1] <- 0
  for (i in 1:n) {
    k1 <- v;                 m1 <- accel_signed(a, comps)
    k2 <- v + h/2*m1;        m2 <- accel_signed(a + h/2*k1, comps)
    k3 <- v + h/2*m2;        m3 <- accel_signed(a + h/2*k2, comps)
    k4 <- v + h*m3;          m4 <- accel_signed(a + h*k3, comps)
    a <- a + h/6*(k1 + 2*k2 + 2*k3 + k4); v <- v + h/6*(m1 + 2*m2 + 2*m3 + m4)
    bw[i+1] <- a
  }
  list(eta = seq(0, -emax, length.out = n+1), a = bw)
}
# on the far sheet a is negative, so a^(-3w) is continued as sign(a)^(-3w)|a|^(-3w): the odd/even
# distinction IS the continuation, and this is where the rule bites.
accel_signed <- function(a, comps) {
  s <- 0
  for (cc in comps) {
    p <- -3*cc$w
    stopifnot(abs(p - round(p)) < 1e-12)          # the parity rule is about integer powers only
    s <- s + cc$c * (1 - 3*cc$w) * sign(a)^(round(p) %% 2) * abs(a)^p
  }
  (K/2) * s
}

cat("=== 1. the fold's condition on the scale factor ===\n")
cat("   Theta: eta -> -eta is an isometry of a(eta)^2(-deta^2 + dx^2) iff a(-eta)^2 = a(eta)^2.\n")
cat("   With a(0) = 0 and a'(0) finite and nonzero, that forces a odd. The even branch needs\n")
cat("   a'(0) = 0, which is a bounce and has no radiation bang.\n")

cat("\n=== 2. which fluids keep a odd ===\n")
cat("      fluid              w        -3w    odd power?   admissible\n")
fl <- list(list("dust", 0), list("radiation", 1/3), list("Lambda", -1),
           list("spatial curvature", -1/3), list("stiff", 1), list("w = 2/3", 2/3))
for (f in fl) {
  w <- f[[2]]; p <- -3*w
  exempt <- abs(1 - 3*w) < 1e-12
  odd <- abs(p - round(p)) < 1e-12 && (round(p) %% 2 != 0)
  cat(sprintf("      %-18s %6.3f  %7.3f    %-10s   %s\n", f[[1]], w, p,
              ifelse(odd, "yes", "no"),
              ifelse(exempt, "yes, coefficient 1-3w vanishes", ifelse(odd, "yes", "NO"))))
}
cat("   Radiation and Lambda are both admissible. Dust is not.\n")

cat("\n=== 3. checked on solutions, not on the series ===\n")
cat("   Integrate a'' from the bang in both directions and measure the fold's violation,\n")
cat("   max | a(-eta)^2 - a(eta)^2 | over the run.\n\n")
cat("      content                                  violation      a''(0)\n")
cases <- list(
  list("radiation only",            list()),
  list("radiation + Lambda",        list(list(c = 0.30, w = -1))),
  list("radiation + curvature",     list(list(c = 0.30, w = -1/3))),
  list("radiation + Lambda + dust", list(list(c = 0.30, w = -1), list(c = 0.30, w = 0))),
  list("radiation + dust",          list(list(c = 0.30, w = 0))))
viol <- numeric(length(cases))
for (i in seq_along(cases)) {
  cmp <- cases[[i]][[2]]
  f <- grow(cmp); b <- grow_back(cmp)
  viol[i] <- max(abs(b$a^2 - f$a^2))
  cat(sprintf("      %-38s %.3e      %+.4f\n", cases[[i]][[1]], viol[i],
              accel_signed(1e-9, cmp)))
}
cat("\n   The three passing cases give exactly zero, and the reason is worth stating rather\n")
cat("   than letting it pass as numerical luck: with only odd powers on the right the equation\n")
cat("   itself is invariant under (eta, a) -> (-eta, -a), so the backward solution IS minus the\n")
cat("   forward one, to the last bit. That invariance is the rule. Dust breaks it at order\n")
cat("   eta^2, by a'' not vanishing at the bang, and adding Lambda alongside does not repair it.\n")
cat("   The integrator is the same one in all five rows, so the difference is the physics.\n")
stopifnot(viol[1] < 1e-9, viol[2] < 1e-9, viol[3] < 1e-9, viol[4] > 1e-3, viol[5] > 1e-3)

cat("\n=== 4. the plant: the test must not pass everything ===\n")
cat("   Section 3 already contains its own plant, since two of the five cases fail and three\n")
cat("   pass. One more, to show the measure is not saturating: scale the dust down and the\n")
cat("   violation must follow it linearly.\n\n")
cat("      dust coefficient      violation     ratio to the previous\n")
prev <- NA
for (cc in c(0.30, 0.15, 0.075, 0.0375)) {
  cmp <- list(list(c = cc, w = 0))
  f <- grow(cmp); b <- grow_back(cmp)
  v <- max(abs(b$a^2 - f$a^2))
  cat(sprintf("      %14.4f      %.4e     %s\n", cc, v,
              ifelse(is.na(prev), "", sprintf("%.4f", v/prev))))
  prev <- v
}
cat("   The violation is linear in the dust, so it measures the dust and not the integrator.\n")

cat("\n=== 5. what it means, and what it does not ===\n")
cat("   Dust is not a fundamental fluid. A species behaves as dust once it is non-relativistic,\n")
cat("   and at the bang every species is relativistic, so the forbidden component is absent for\n")
cat("   a physical reason and not by fiat. Read forwards the rule says: the fold REQUIRES the\n")
cat("   bang to be hot. A cold component at eta = 0 would break the parity the whole\n")
cat("   construction is built on. The radiation bang of section 2.2 is therefore forced rather\n")
cat("   than assumed.\n")
cat("\n   For Lambda the rule is permissive, and that is the answer. w = -1 gives -3w = 3, an odd\n")
cat("   power, so the fold admits any value of Lambda and constrains none. Spatial curvature,\n")
cat("   w = -1/3, is admitted on the same footing. The bang does not fix Lambda.\n")

cat("\n=== 6. what is NOT closed, stated precisely ===\n")
cat("   It is tempting to finish by saying the fold switches off its own image stress at the\n")
cat("   bang, since that stress vanishes for conformally invariant matter and the bang has just\n")
cat("   been forced to be radiation. That is too strong and image_stress_conformal.R says why.\n")
cat("   The breaking sits in the Hadamard coefficient V_0 = Delta^(1/2)[m^2 + (xi - 1/6)R]/2,\n")
cat("   and m^2 occupies the same slot as (xi - 1/6)R. A massive field therefore feels the image\n")
cat("   term whatever its coupling, and being relativistic at the bang makes m/T small without\n")
cat("   making m zero. So the bang carries an image stress of order m^2, which for the adopted\n")
cat("   content is of order M_1^2 and is the gamma^2 of section 2.2 rather than nothing.\n")
cat("   Whether it can look like Lambda turns on how it scales with a, since only a term\n")
cat("   independent of a is a cosmological constant. That scaling is not computed here and it\n")
cat("   is the one route to the value of Lambda this construction has not yet closed.\n")
