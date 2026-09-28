#!/usr/bin/env Rscript
# alpha = J o P_perp is total inversion, and nothing about it is chosen.
#
# Section 2.1 asserts three things in one line and computes only the third. The wedge reflection's
# point map J fixes the bifurcation sphere, the transverse parity P_perp is -Id, and the
# composition is the free antipodal involution alpha. The uniqueness of -Id among free involutive
# isometries of a round sphere is checked in perp_is_unique.R. What has never been computed is the
# composition itself: that J o P_perp is total inversion X -> -X on the embedding, hence that
# writing the fold as J o P_perp rather than J is not a modelling choice but the only way to get a
# map with no fixed points. This does it, and prices each of the alternatives.
#
# de Sitter space as the hyperboloid  -X0^2 + X1^2 + X2^2 + X3^2 + X4^2 = L^2  in R^{1,4}. The
# static patch's boost acts in the (X0, X1) plane, so its bifurcation surface is X0 = X1 = 0 with
# X2^2 + X3^2 + X4^2 = L^2, a round two-sphere. On the embedding,
#     J        : (X0, X1) -> (-X0, -X1),   transverse fixed
#     P_perp   : (X2, X3, X4) -> -(X2, X3, X4),   boost plane fixed
#     J o P_perp = -Id.
# Both factors have fixed points on the hyperboloid and the product has none. That is the whole of
# the difference the manuscript calls "the whole of what is proposed", made a number.

L <- 1.4
TOL <- 1e-13
fail <- 0
note <- function(ok, what) {
  if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 }
  invisible(ok)
}
Q <- diag(c(-1, 1, 1, 1, 1))                      # the R^{1,4} form
form <- function(X) as.numeric(t(X) %*% Q %*% X)

# A sample of hyperboloid points: pick X0 and a direction, then solve for the spatial radius.
set.seed(20260927)
NS <- 4000
samp <- t(sapply(1:NS, function(i) {
  x0 <- rnorm(1, 0, 2)
  u <- rnorm(4); u <- u / sqrt(sum(u^2))
  c(x0, u * sqrt(L^2 + x0^2))
}))
cat(sprintf("=== 0. the sample: %d points, worst |form - L^2| = %.2e ===\n",
            NS, max(abs(apply(samp, 1, form) - L^2))))
note(max(abs(apply(samp, 1, form) - L^2)) < 1e-12, "sample lies on the hyperboloid")

Jm  <- diag(c(-1, -1, 1, 1, 1))
Pm  <- diag(c(1, 1, -1, -1, -1))
Rot <- diag(c(1, 1, -1, -1, 1))                   # a pi rotation in the (2,3) plane, not -Id
maps <- list(J = Jm, P_perp = Pm, "J o P_perp" = Jm %*% Pm, "J o rot" = Jm %*% Rot)

cat("\n=== 1. each candidate is an isometry of the hyperboloid, and which are involutions ===\n")
cat("      map              preserves the form   maps hyperboloid to itself   M^2 = Id   det\n")
for (nm in names(maps)) {
  M <- maps[[nm]]
  dq <- max(abs(t(M) %*% Q %*% M - Q))
  dh <- max(abs(apply(samp %*% t(M), 1, form) - L^2))
  di <- max(abs(M %*% M - diag(5)))
  cat(sprintf("   %-16s %18.2e %28.2e %11.1e %6.0f\n", nm, dq, dh, di, det(M)))
  note(dq < TOL && dh < 1e-12 && di < TOL, paste(nm, "is an involutive isometry"))
}

cat("\n=== 2. the composition is total inversion, exactly ===\n")
comp <- Jm %*% Pm
cat(sprintf("   max |J P_perp - (-Id)| over the matrix          = %.1e\n",
            max(abs(comp + diag(5)))))
cat(sprintf("   max |J P_perp X + X| over %d sampled points    = %.1e\n",
            NS, max(abs(samp %*% t(comp) + samp))))
note(max(abs(comp + diag(5))) < TOL, "J P_perp is -Id as a matrix")
note(max(abs(samp %*% t(comp) + samp)) < 1e-12, "and pointwise on the hyperboloid")
cat("   Neither factor alone is:\n")
for (nm in c("J", "P_perp", "J o rot")) {
  cat(sprintf("      max |%s X + X| = %.4f\n", nm, max(abs(samp %*% t(maps[[nm]]) + samp))))
  note(max(abs(samp %*% t(maps[[nm]]) + samp)) > 0.1, paste(nm, "is not -Id"))
}

cat("\n=== 3. fixed points: the reason the composition and not either factor ===\n")
cat("   A fold pairs two sheets, so no point may be its own image. Counting points of the\n")
cat("   hyperboloid fixed by each map, as the dimension of the fixed subspace intersected with\n")
cat("   the hyperboloid, and measuring the closest approach over the sample:\n")
cat("      map              +1 eigenvalues   fixed set on the hyperboloid   min |X' - X| sampled\n")
descr <- c(J = "the bifurcation 2-sphere X0=X1=0",
           P_perp = "the 1-dim hyperbola X2=X3=X4=0",
           "J o P_perp" = "empty",
           "J o rot" = "two points X0=X1=X2=X3=0")
for (nm in names(maps)) {
  M <- maps[[nm]]
  np <- sum(abs(eigen(M, symmetric = TRUE)$values - 1) < 1e-12)
  md <- min(sqrt(rowSums((samp %*% t(M) - samp)^2)))
  cat(sprintf("   %-16s %13d   %-30s %12.4f\n", nm, np, descr[[nm]], md))
}
cat("   The sampled minimum is not the right instrument for the last column, because a random\n")
cat("   sample misses a measure-zero fixed set. The fixed sets are found exactly instead, by\n")
cat("   intersecting each +1 eigenspace with the hyperboloid:\n")
for (nm in names(maps)) {
  M <- maps[[nm]]
  ev <- eigen(M, symmetric = TRUE)
  V <- ev$vectors[, abs(ev$values - 1) < 1e-12, drop = FALSE]
  if (ncol(V) == 0) { cat(sprintf("   %-16s no +1 eigenspace, so no fixed point\n", nm)); next }
  # the restriction of the form to the eigenspace: nonempty on the hyperboloid iff it takes L^2
  Qr <- t(V) %*% Q %*% V
  sg <- eigen(Qr, symmetric = TRUE)$values
  reach <- any(sg > 1e-12)                      # some direction with positive form
  cat(sprintf("   %-16s eigenspace dim %d, form signature (%s), reaches +L^2: %s\n",
              nm, ncol(V), paste(sprintf("%+.0f", sign(sg)), collapse = ","),
              if (reach) "YES, fixed points exist" else "no"))
  note(if (nm == "J o P_perp") !reach else reach, paste("fixed-point verdict for", nm))
}
cat("   So J fixes a two-sphere, P_perp fixes a hyperbola, a rotation in place of the parity\n")
cat("   fixes two points, and only the total inversion is free. The composition is not chosen\n")
cat("   over J for elegance; J is not available.\n")

cat("\n=== 4. the differential, and the orientation on the sphere ===\n")
cat("   The map is linear, so its differential is itself: -Id on the 5 embedding directions and\n")
cat("   -Id on the 4-dimensional tangent space of the hyperboloid. Checked at three points by\n")
cat("   projecting off the normal:\n")
for (i in c(1, 17, 333)) {
  X <- samp[i, ]
  n <- Q %*% X; n <- n / sqrt(abs(as.numeric(t(n) %*% solve(Q) %*% n)))
  P <- diag(5) - (X %*% t(Q %*% X)) / as.numeric(t(X) %*% Q %*% X)   # tangent projector
  dT <- P %*% comp %*% P
  cat(sprintf("      point %4d: max |dTheta + P| on the tangent space = %.2e, trace = %+.6f\n",
              i, max(abs(dT + P)), sum(diag(dT))))
  note(max(abs(dT + P)) < 1e-12, "dTheta is -Id on the tangent space")
  note(abs(sum(diag(dT)) + 4) < 1e-10, "and its trace is -4, so all four directions turn over")
}
cat(sprintf("   det P_perp on the transverse R^3 = %+.0f, so it reverses orientation there and\n",
            det(diag(c(-1, -1, -1)))))
cat("   lies outside SO(3), which is why no rotation absorbs it.\n")
note(det(diag(c(-1, -1, -1))) < 0, "P_perp reverses transverse orientation")

cat("\n=== 5. the (-1)^ell that the parity puts on the cross term ===\n")
cat("   P_perp is the antipodal map of the sphere, so Y_lm(-n) = (-1)^l Y_lm(n). Checked with\n")
cat("   the Legendre form of the addition theorem, P_l(-c) = (-1)^l P_l(c), at 400 random pairs\n")
cat("   of directions per multipole:\n")
Pl <- function(l, x) {                             # Legendre by recurrence
  p0 <- rep(1, length(x)); if (l == 0) return(p0)
  p1 <- x; if (l == 1) return(p1)
  for (k in 2:l) { p2 <- ((2*k-1)*x*p1 - (k-1)*p0)/k; p0 <- p1; p1 <- p2 }
  p1
}
ua <- matrix(rnorm(1200), 400); ua <- ua / sqrt(rowSums(ua^2))
ub <- matrix(rnorm(1200), 400); ub <- ub / sqrt(rowSums(ub^2))
cs <- rowSums(ua * ub)
cat("      l    max |P_l(-c) - (-1)^l P_l(c)|\n")
for (l in 0:6) {
  d <- max(abs(Pl(l, -cs) - (-1)^l * Pl(l, cs)))
  cat(sprintf("   %4d %32.2e\n", l, d))
  note(d < 1e-12, sprintf("parity of P_%d", l))
}
cat("   Odd multipoles change sign and even ones do not, which is the (-1)^l on the cross-sheet\n")
cat("   term and the reason a stress evaluated at one squeeze comes out with opposite signs at\n")
cat("   the pole and at the equator.\n")

cat("\n=== 6. plants ===\n")
cat("   (a) A map that is not an isometry, X2 scaled by 1.1 inside the reflection:\n")
Bad <- Jm; Bad[3, 3] <- 1.1
cat(sprintf("       max |M^T Q M - Q| = %.3f\n", max(abs(t(Bad) %*% Q %*% Bad - Q))))
note(max(abs(t(Bad) %*% Q %*% Bad - Q)) > 0.1, "PLANT (a) fires: not an isometry")
cat("   (b) Claiming J alone is the fold. Its +1 eigenspace meets the hyperboloid, so the\n")
cat("       claim is refuted by section 3 above rather than by taste:\n")
ev <- eigen(Jm, symmetric = TRUE); V <- ev$vectors[, abs(ev$values - 1) < 1e-12, drop = FALSE]
cat(sprintf("       J has a %d-dimensional +1 eigenspace and the form restricted to it is\n",
            ncol(V)))
cat(sprintf("       positive definite, so it carries a whole %d-sphere of fixed points.\n",
            ncol(V) - 1))
note(ncol(V) == 3 && all(eigen(t(V) %*% Q %*% V, symmetric = TRUE)$values > 0),
     "PLANT (b) fires: J has a 2-sphere of fixed points")
cat("   (c) An odd multipole with the parity dropped, to show section 5 can fail:\n")
d <- max(abs(Pl(3, -cs) - Pl(3, cs)))
cat(sprintf("       max |P_3(-c) - P_3(c)| with the sign omitted = %.4f\n", d))
note(d > 0.1, "PLANT (c) fires: dropping (-1)^l breaks the identity")

cat("\n=== 7. what this settles ===\n")
cat("   Not that the fold exists, which is Section 2's hypothesis. What it settles is that once\n")
cat("   the fold is a free involution of the de Sitter cover, its form is forced: total inversion,\n")
cat("   which factorises through the wedge reflection and the transverse parity in exactly one\n")
cat("   way. So the transverse parity is not a cost sitting beside the construction. It is what\n")
cat("   the absence of fixed points costs, and the (-1)^l it puts on the cross term is a\n")
cat("   prediction to be carried rather than an assumption to be priced.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}
