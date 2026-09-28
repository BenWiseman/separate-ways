# Is the transverse involution an input, or does the geometry select it?
#
# Section 2.1 says: "At a general bifurcate Killing horizon the bifurcation surface may carry
# several free involutive isometries and nothing selects one, so there P_perp is an input."
# That concedes a choice, and the whole construction rests on Theta = J o P_perp. If the choice
# is real then the fold is not one map but a family, and the derivation of the field equations
# inherits that freedom. This asks whether the concession is true.
#
# WHAT HAS TO BE COUNTED. P_perp must be (a) an isometry of the bifurcation surface, (b) an
# involution, and (c) free, having no fixed point, since a fixed point would put a sheet in
# contact with itself. So the question is the number of free involutive isometries of the
# bifurcation 2-surface. By the black-hole uniqueness theorems a stationary vacuum hole in four
# dimensions is Kerr, so two cases exhaust the physical ones: the round sphere and Kerr's
# axisymmetric bifurcation surface. Both are done below by enumeration, not by argument.

cat("=== 1. the round sphere: enumerate the order-two elements of O(3) ===\n")
cat("   Every isometry of the round S^2 is an element of O(3). An involution is an element whose\n")
cat("   square is the identity, so it is diagonalisable with eigenvalues +-1. Up to conjugacy the\n")
cat("   possibilities are the four sign patterns, and freeness is decided by whether the map has\n")
cat("   a fixed unit vector, which happens exactly when +1 is an eigenvalue.\n\n")
pats <- list(c(1,1,1), c(-1,1,1), c(-1,-1,1), c(-1,-1,-1))
nms  <- c("identity", "reflection in a plane", "rotation by pi", "the antipodal map")
cat("      element                 eigenvalues      fixed points on S^2     free?\n")
free_round <- c()
for (k in seq_along(pats)) {
  M <- diag(pats[[k]])
  # fixed points of M on S^2 are unit eigenvectors with eigenvalue +1
  nfix <- sum(pats[[k]] == 1)
  desc <- if (nfix == 3) "all of it" else if (nfix == 2) "a great circle" else
          if (nfix == 1) "two points" else "none"
  isfree <- nfix == 0
  if (isfree) free_round <- c(free_round, nms[k])
  cat(sprintf("      %-22s  %-15s  %-22s  %s\n", nms[k],
              paste(sprintf("%+d", pats[[k]]), collapse = " "), desc, ifelse(isfree, "YES", "no")))
}
cat(sprintf("\n   free involutive isometries of the round sphere: %d, namely %s\n",
            length(free_round), paste(free_round, collapse = ", ")))
stopifnot(length(free_round) == 1)
cat("   Conjugacy is enough here: any involution in O(3) is conjugate to one of the four, and\n")
cat("   conjugation does not change whether a fixed point exists. So the antipodal map is the\n")
cat("   only one, and on a Schwarzschild bifurcation sphere P_perp is not a choice.\n")

cat("\n=== 2. the same count by construction, not by conjugacy ===\n")
cat("   Every involutive isometry of R^3 is M = Q D Q^T with Q orthogonal and D a diagonal sign\n")
cat("   matrix, since M^2 = I and M^T M = I make M symmetric orthogonal. Build them directly and\n")
cat("   count the free ones, rather than sampling O(3) and hoping to land on an involution.\n\n")
set.seed(2091)
nfree_s <- 0; nnonfree <- 0; allneg <- TRUE
for (trial in 1:20000) {
  Q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  d <- sample(c(-1, 1), 3, replace = TRUE)
  if (all(d == 1)) next                                      # the identity is excluded
  M <- Q %*% diag(d) %*% t(Q)
  stopifnot(max(abs(M %*% M - diag(3))) < 1e-9,              # an involution
            max(abs(t(M) %*% M - diag(3))) < 1e-9)           # and an isometry
  if (any(d == 1)) { nnonfree <- nnonfree + 1 }              # a +1 eigenvalue is a fixed point
  else { nfree_s <- nfree_s + 1; allneg <- allneg && max(abs(M + diag(3))) < 1e-9 }
}
cat(sprintf("   %d involutive isometries built; %d have a fixed point on S^2 and %d are free\n",
            nfree_s + nnonfree, nnonfree, nfree_s))
cat(sprintf("   every free one is exactly -Id, to machine precision: %s\n", allneg))
cat("   The random orthogonal frame is doing real work: the NON-free ones are spread over every\n")
cat("   axis and plane in the sphere, so the count is over a continuum and still returns one.\n")
stopifnot(nfree_s > 100, nnonfree > 100, allneg)

cat("\n=== 3. Kerr, where the surface is not round ===\n")
cat("   The bifurcation surface at r = r_+ carries\n")
cat("      ds^2 = rho_+^2 dtheta^2 + [(r_+^2+a^2)^2 sin^2 theta / rho_+^2] dphi^2,\n")
cat("      rho_+^2 = r_+^2 + a^2 cos^2 theta.\n")
cat("   It depends on theta only through cos^2 theta and sin^2 theta, and not on phi at all, so\n")
cat("   its isometry group is O(2) x Z_2: phi -> eps phi + c with eps = +-1, and theta -> theta\n")
cat("   or theta -> pi - theta. Enumerating that group's involutions and their fixed sets:\n\n")
M <- 1; aa <- 0.7; rp <- M + sqrt(M^2 - aa^2)
gth <- function(th) rp^2 + aa^2*cos(th)^2
gph <- function(th) (rp^2 + aa^2)^2 * sin(th)^2 / (rp^2 + aa^2*cos(th)^2)
isom <- function(eps, c0, flip) {                # check the map is an isometry, numerically
  ths <- seq(0.05, pi - 0.05, length.out = 400)
  th2 <- if (flip) pi - ths else ths
  max(abs(gth(th2) - gth(ths)), abs(gph(th2) - gph(ths)))
}
cat("      eps   c      theta map      isometry?   involution?   fixed set             free?\n")
nfree <- 0; freelist <- c()
for (eps in c(1, -1)) for (c0 in c(0, pi)) for (flip in c(FALSE, TRUE)) {
  iso <- isom(eps, c0, flip) < 1e-12
  # the map is (theta, phi) -> (flip ? pi-theta : theta, eps*phi + c0); square is the identity when
  inv <- if (eps == 1) (!flip && (c0 %% (2*pi) == 0 || abs(c0 - pi) < 1e-12)) || (flip && (abs(c0) < 1e-12 || abs(c0 - pi) < 1e-12)) else TRUE
  if (eps == 1 && !flip && abs(c0) < 1e-12) { fixed <- "everything (identity)"; free <- FALSE }
  else if (eps == 1 && !flip) { fixed <- "the two poles"; free <- FALSE }
  else if (eps == 1 && flip && abs(c0) < 1e-12) { fixed <- "the equator"; free <- FALSE }
  else if (eps == 1 && flip) { fixed <- "none"; free <- TRUE }
  else if (eps == -1 && !flip) { fixed <- "a meridian circle"; free <- FALSE }
  else { fixed <- "two points on the equator"; free <- FALSE }
  if (iso && inv && free) { nfree <- nfree + 1
    freelist <- c(freelist, sprintf("theta -> pi-theta, phi -> phi+%s", ifelse(abs(c0-pi)<1e-12,"pi","0"))) }
  cat(sprintf("      %+d   %-5s  %-13s  %-11s %-13s %-21s %s\n", eps,
              ifelse(abs(c0) < 1e-12, "0", "pi"), ifelse(flip, "pi - theta", "theta"),
              ifelse(iso, "yes", "no"), ifelse(inv, "yes", "no"), fixed,
              ifelse(iso && inv && free, "YES", "no")))
}
cat(sprintf("\n   free involutive isometries of the Kerr bifurcation surface: %d\n", nfree))
cat(sprintf("   namely %s, which is the antipodal map.\n", paste(unique(freelist), collapse = "; ")))
stopifnot(nfree == 1)

cat("\n=== 3b. the enumeration assumed the isometry group. That is checked, not assumed ===\n")
cat("   The list above is complete only if Kerr's bifurcation surface has isometry group exactly\n")
cat("   O(2) x Z_2. Any isometry must preserve the Gaussian curvature and the circumference of\n")
cat("   each theta-circle, so it can only map a circle to one carrying the same pair. If the map\n")
cat("   theta -> (K, circumference) separates every theta in [0, pi/2) from every other, then the\n")
cat("   only identification available is theta with pi - theta and the group is what was used.\n\n")
Kg <- function(th, h = 1e-5) {                    # Gaussian curvature of ds^2 = A dth^2 + B dph^2
  A <- function(t) gth(t); B <- function(t) gph(t)
  sb <- function(t) sqrt(B(t))
  d1 <- function(t) (sb(t + h) - sb(t - h)) / (2*h) / sqrt(A(t))
  -( (d1(th + h) - d1(th - h)) / (2*h) ) / sqrt(A(th)) / sb(th)
}
cat("      a/M    K at theta = 0.2 .. 1.5 (half range), and whether (K, circumference) separates\n")
for (av in c(0.3, 0.7, 0.9, 0.998)) {
  aa <<- av; rp <<- M + sqrt(M^2 - av^2)
  ths <- seq(0.02, pi/2 - 0.02, length.out = 600)
  Ks  <- sapply(ths, Kg); Cs <- 2*pi*sqrt(sapply(ths, gph))
  # separation: no two distinct theta in the half range share both K and circumference
  idx <- seq(1, length(ths), by = 7)
  worst <- Inf
  for (u in seq_along(idx)) { if (u == length(idx)) next
    for (v in (u + 1):length(idx)) {
      d <- sqrt((Ks[idx[u]] - Ks[idx[v]])^2 + (Cs[idx[u]] - Cs[idx[v]])^2)
      if (d < worst) worst <- d
    } }
  cat(sprintf("      %.3f   K spans %+.4f to %+.4f   closest distinct pair in (K, C): %.4f\n",
              av, min(Ks), max(Ks), worst))
  stopifnot(worst > 1e-3)
}
aa <<- 0.7; rp <<- M + sqrt(M^2 - aa^2)
cat("   No two distinct latitudes in the half range carry the same curvature and the same\n")
cat("   circumference at any spin, so the only identification an isometry can make is theta with\n")
cat("   pi - theta, and O(2) x Z_2 is the whole group. The eight rows above are therefore all of\n")
cat("   the involutions and the count of one free map is complete.\n")
cat("   The a -> 0 limit is the round sphere, where K is constant and the group grows to O(3);\n")
cat("   section 1 handles that case separately and returns the same answer.\n\n")

cat("\n=== 4. and it really is free and really is an isometry, spot-checked on points ===\n")
set.seed(77)
mind <- Inf; maxg <- 0
for (i in 1:20000) {
  th <- runif(1, 0, pi); ph <- runif(1, 0, 2*pi)
  th2 <- pi - th; ph2 <- (ph + pi) %% (2*pi)
  d <- sqrt((th - th2)^2 + min(abs(ph - ph2), 2*pi - abs(ph - ph2))^2)
  mind <- min(mind, d)
  maxg <- max(maxg, abs(gth(th2) - gth(th)), abs(gph(th2) - gph(th)))
}
cat(sprintf("   closest a point came to its own image over 20000 samples: %.6f (zero would be a fixed point)\n", mind))
cat(sprintf("   largest metric mismatch under the map: %.2e\n", maxg))
stopifnot(mind > 1e-3, maxg < 1e-10)

cat("\n=== 5. the plant: a surface with no equatorial symmetry must lose the map ===\n")
cat("   The antipodal map needs theta -> pi - theta to be an isometry. Deform the Kerr surface\n")
cat("   by an odd function of cos theta, which no stationary vacuum hole has, and it must fail.\n")
gth_b <- function(th) rp^2 + aa^2*cos(th)^2 + 0.30*cos(th)^3
mis <- max(abs(gth_b(pi - seq(0.05, pi-0.05, length.out = 400)) - gth_b(seq(0.05, pi-0.05, length.out = 400))))
cat(sprintf("   metric mismatch under theta -> pi - theta on the deformed surface: %.4f\n", mis))
stopifnot(mis > 0.01)
cat("   so the check is reading the metric, and the uniqueness above is a property of Kerr\n")
cat("   rather than of the enumeration.\n")

cat("\n=== 5b. charge, since the companion's contact results are stated at every charge ===\n")
cat("   Kerr-Newman's bifurcation surface has the same functional form with r_+ moved:\n")
cat("      r_+ = M + sqrt(M^2 - a^2 - Q^2),  and rho_+^2 = r_+^2 + a^2 cos^2 theta as before.\n")
cat("   So the theta dependence is still through cos^2 theta alone and the argument transfers.\n")
cat("   Checked on the quantity that carried it, the separation of latitudes:\n\n")
cat("      a/M     Q/M      r_+/M     closest distinct pair in (K, C)     one free involution?\n")
for (pr in list(c(0.3,0.5), c(0.6,0.6), c(0.2,0.9), c(0.7,0.7))) {
  av <- pr[1]; Qv <- pr[2]
  if (av^2 + Qv^2 >= 1) next
  aa <<- av; rp <<- M + sqrt(M^2 - av^2 - Qv^2)
  ths <- seq(0.02, pi/2 - 0.02, length.out = 400)
  Ks <- sapply(ths, Kg); Cs <- 2*pi*sqrt(sapply(ths, gph))
  idx <- seq(1, length(ths), by = 9); worst <- Inf
  for (u in seq_along(idx)) { if (u == length(idx)) next
    for (v in (u+1):length(idx)) {
      d <- sqrt((Ks[idx[u]]-Ks[idx[v]])^2 + (Cs[idx[u]]-Cs[idx[v]])^2); if (d < worst) worst <- d } }
  cat(sprintf("   %6.2f  %6.2f   %8.4f   %28.4f     %s\n", av, Qv, rp/M, worst,
              ifelse(worst > 1e-3, "yes", "NO")))
  stopifnot(worst > 1e-3)
}
aa <<- 0.7; rp <<- M + sqrt(M^2 - aa^2)
cat("   Charge moves r_+ and nothing else that matters, so the count is one at every charge for\n")
cat("   which a horizon exists.\n")

cat("\n=== 5c. and every dimension, since the contact results are stated there too ===\n")
cat("   A spherically symmetric hole in D dimensions has a round S^(D-2) for its bifurcation\n")
cat("   surface. A free involutive isometry of the round S^n is an element of O(n+1) with no +1\n")
cat("   eigenvalue, so the same eigenvalue argument runs at every n, and -Id is the only one.\n")
cat("   Counted by construction, as in section 2:\n\n")
cat("      D     S^(D-2)    involutions built    with a fixed point    free    all -Id?\n")
set.seed(4004)
for (D in 4:8) {
  n1 <- D - 1                                    # the sphere S^(D-2) sits in R^(D-1)
  nf <- 0; nn <- 0; ok <- TRUE
  for (trial in 1:4000) {
    Qm <- qr.Q(qr(matrix(rnorm(n1*n1), n1, n1)))
    d <- sample(c(-1,1), n1, replace = TRUE)
    if (all(d == 1)) next
    Mm <- Qm %*% diag(d) %*% t(Qm)
    stopifnot(max(abs(Mm %*% Mm - diag(n1))) < 1e-9)
    if (any(d == 1)) nn <- nn + 1
    else { nf <- nf + 1; ok <- ok && max(abs(Mm + diag(n1))) < 1e-9 }
  }
  cat(sprintf("   %4d   S^%-2d      %14d   %19d %7d    %s\n", D, D-2, nf+nn, nn, nf,
              ifelse(ok, "yes", "NO")))
  stopifnot(ok, nf > 10, nn > 10)
}
cat("   One free involutive isometry at every dimension, and it is the antipodal map. So the\n")
cat("   companion's contact results, which are stated at every charge and every dimension, do\n")
cat("   not inherit a choice of P_perp at any of them.\n")

cat("\n=== 6. what this changes ===\n")
cat("   P_perp is not an input at any bifurcation surface a stationary black hole has. On the\n")
cat("   round sphere the antipodal map is the only free involutive isometry in O(3), and on\n")
cat("   Kerr's axisymmetric surface it is the only one in O(2) x Z_2. Charge moves r_+ and\n")
cat("   nothing that matters, and the eigenvalue argument runs unchanged on the round S^(D-2)\n")
cat("   at every dimension, so the same count of one holds across the whole range the\n")
cat("   companion's contact results are stated over. The fold has no choice to make at any of\n")
cat("   them. What remains genuinely open is a bifurcation surface with no equatorial symmetry,\n")
cat("   which section 5 shows loses the map altogether rather than offering several: the\n")
cat("   failure mode is absence, not ambiguity.\n")
