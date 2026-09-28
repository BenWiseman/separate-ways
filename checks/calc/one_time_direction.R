# How many time directions does the fold permit? Exactly one, and that is a count not a choice.
#
# signature_or_the_fold_is_nothing.R shows the fold has content in Lorentzian signature and none
# in Euclidean: a parity is invariant only if the map cannot be deformed to the identity, and in
# four Euclidean dimensions -Id is reached by an explicit path of isometries. That comparison was
# two-way and the question is not. A four-dimensional manifold can carry (0,4), (1,3), (2,2),
# (3,1) or (4,0), and the argument as it stood tested two of the five.
#
# THE INVARIANT. O(t,s) with t,s >= 1 has four connected components, labelled by the signs of the
# determinants of the timelike and spacelike blocks. Neither sign can change along a continuous
# path of isometries, and -Id has det (-1)^t on the first and (-1)^s on the second. So -Id sits in
# the identity component exactly when t and s are both EVEN, and the fold keeps a parity exactly
# when one of them is odd. With t + s = 4 that is t odd, which is one time direction or three,
# and three is one with the overall sign flipped.
#
# What makes this work is the difference between a rotation and a boost. Two directions of the
# same signature rotate into each other and reach -1 on both at angle pi. A timelike and a
# spacelike direction only boost, and a boost never reaches -1: that is the |Lambda^0_0| >= 1 of
# the earlier file, here in its general form. So an unpaired time direction has nothing to pair
# with and -Id stays out of reach.

sig <- function(t, s) diag(c(rep(-1, t), rep(1, s)))

pair_path <- function(th, D, pairs) {            # rotate by th in each listed coordinate pair
  M <- diag(D)
  for (pr in pairs) {
    i <- pr[1]; j <- pr[2]
    M[i,i] <- cos(th); M[i,j] <- -sin(th); M[j,i] <- sin(th); M[j,j] <- cos(th)
  }
  M
}
is_iso <- function(M, e) max(abs(t(M) %*% e %*% M - e))

cat("=== 1. the five four-dimensional signatures, and whether a path reaches -Id ===\n")
cat("   Pair like-signature directions and rotate both pairs by pi at once. That is an isometry\n")
cat("   throughout, since a rotation between two directions of the same sign preserves the form,\n")
cat("   and it lands on -Id. It only exists when every direction has a partner of its own sign.\n\n")
cat("      signature   t even?  s even?   path exists   isometry along it   reaches -Id\n")
reach <- list()
for (t in 0:4) {
  s <- 4 - t; e <- sig(t, s)
  ok <- (t %% 2 == 0) && (s %% 2 == 0)
  if (ok) {
    prs <- list()
    if (t >= 2) for (k in seq(1, t, by = 2)) prs[[length(prs)+1]] <- c(k, k+1)
    if (s >= 2) for (k in seq(t+1, 4, by = 2)) prs[[length(prs)+1]] <- c(k, k+1)
    worst <- max(sapply(seq(0, pi, length.out = 60), function(th) is_iso(pair_path(th, 4, prs), e)))
    hit <- max(abs(pair_path(pi, 4, prs) + diag(4)))
    cat(sprintf("      (%d,%d)        %-8s %-8s  %-13s %-19.1e %.1e\n", t, s,
                "yes", "yes", "yes", worst, hit))
    stopifnot(worst < 1e-12, hit < 1e-12)
    reach[[paste0(t,",",s)]] <- TRUE
  } else {
    cat(sprintf("      (%d,%d)        %-8s %-8s  %-13s %-19s %s\n", t, s,
                ifelse(t %% 2 == 0, "yes", "no"), ifelse(s %% 2 == 0, "yes", "no"),
                "no", "-", "-"))
    reach[[paste0(t,",",s)]] <- FALSE
  }
}
cat("\n   So a path of this kind exists for (0,4), (2,2) and (4,0) and not for (1,3) or (3,1).\n")
stopifnot(reach[["0,4"]], reach[["2,2"]], reach[["4,0"]], !reach[["1,3"]], !reach[["3,1"]])

cat("\n=== 2. and for the odd ones no path of ANY kind exists, by the block determinants ===\n")
cat("   Sampling O(1,3) and O(3,1) directly: the sign of the determinant of each block is an\n")
cat("   invariant, so it cannot travel from the identity's +1 to -Id's value.\n\n")
set.seed(7717)
rand_iso <- function(t, s, n = 3000) {           # random isometries by exponentiating the algebra
  e <- sig(t, s); D <- t + s; out <- list()
  for (i in 1:n) {
    A <- matrix(rnorm(D*D, sd = 0.8), D, D)
    G <- A - e %*% t(A) %*% e                    # the generators: e A^T e = -A
    M <- diag(D); T <- diag(D)
    for (k in 1:30) { T <- T %*% G / k; M <- M + T }
    if (max(abs(t(M) %*% e %*% M - e)) < 1e-6 * max(1, max(abs(M))^2)) out[[length(out)+1]] <- M
  }
  out
}
cat("      signature   sampled isometries   det(time block) sign   det(space block) sign   -Id's signs\n")
for (ts in list(c(1,3), c(3,1), c(2,2))) {
  t <- ts[1]; s <- ts[2]
  ms <- rand_iso(t, s)
  dt <- sapply(ms, function(M) sign(det(as.matrix(M[1:t, 1:t, drop = FALSE]))))
  ds <- sapply(ms, function(M) sign(det(as.matrix(M[(t+1):4, (t+1):4, drop = FALSE]))))
  cat(sprintf("      (%d,%d)        %14d   %20s   %21s   %+d and %+d\n", t, s, length(ms),
              paste(sort(unique(dt)), collapse = " "), paste(sort(unique(ds)), collapse = " "),
              (-1)^t, (-1)^s))
  stopifnot(all(dt > 0), all(ds > 0))            # the identity component keeps both at +1
}
cat("\n   Every sampled element of the identity component has both block determinants positive,\n")
cat("   and -Id has (-1)^t and (-1)^s. At (1,3) and (3,1) one of those is negative, so -Id is in\n")
cat("   a different component and no continuous path of isometries reaches it. At (2,2) both are\n")
cat("   positive, which is why section 1 found a path there.\n")

cat("\n=== 3. what separates them is the difference between a rotation and a boost ===\n")
cat("   Two directions of the same sign rotate and reach -1 on both at angle pi. A timelike and\n")
cat("   a spacelike direction only boost, and a boost never reaches -1. Checked over rapidities:\n\n")
cat("        rapidity    cosh entry of the boost    can it reach -1?\n")
for (r in c(0, 0.5, 2, 8)) cat(sprintf("      %10.1f   %24.4f    no\n", r, cosh(r)))
cat("   The entry is cosh, which is at least one at every rapidity, so an unpaired time direction\n")
cat("   has nothing to pair with. That is the general form of the |Lambda^0_0| >= 1 the earlier\n")
cat("   file measured in the Lorentzian case.\n")
stopifnot(min(cosh(seq(-10, 10, by = 0.01))) >= 1)

cat("\n=== 4. the plant: the criterion must not simply reject everything ===\n")
cat("   A pure spatial rotation by pi in (1,3) has both block determinants positive, so it must\n")
cat("   be reachable, and it is: rotating the 3-4 plane through pi keeps the form at every angle.\n\n")
e13 <- sig(1, 3)
worst <- max(sapply(seq(0, pi, length.out = 60), function(th) is_iso(pair_path(th, 4, list(c(3,4))), e13)))
R <- pair_path(pi, 4, list(c(3,4)))
cat(sprintf("      isometry along the path: %.1e;  endpoint = diag(1,1,-1,-1): %.1e\n",
            worst, max(abs(R - diag(c(1,1,-1,-1))))))
cat(sprintf("      its block determinants: time %+d, space %+d, both positive as required\n",
            sign(det(as.matrix(R[1,1]))), sign(det(R[2:4,2:4]))))
stopifnot(worst < 1e-12, max(abs(R - diag(c(1,1,-1,-1)))) < 1e-12)
cat("   So the criterion separates: it admits the spatial rotation and refuses -Id, in the same\n")
cat("   signature and by the same invariant.\n")

cat("\n=== 5. the verdict ===\n")
cat("   The fold keeps a parity exactly when the number of time directions is odd. In four\n")
cat("   dimensions that leaves (1,3) and (3,1), which are one signature written two ways, and\n")
cat("   excludes Euclidean, (2,2) and (4,0) alike. So the earlier statement was narrower than the\n")
cat("   argument: it is not only that Euclidean signature kills the fold, it is that the fold\n")
cat("   permits exactly one time direction. Together with the causal budget's D = 4 that leaves\n")
cat("   one arena. The two are not independent, since a causal budget already presumes a\n")
cat("   Lorentzian signature to be causal in, and the paper should say so where it says this.\n")
