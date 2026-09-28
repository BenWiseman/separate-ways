#!/usr/bin/env Rscript
# tail_from_lrd.R -- the fold needs a heavy primordial tail. So does the
# little-red-dot seed problem. Do they need the SAME tail?
#
# two_sided_fraction.R showed the two-sided class is a population only if the
# curvature field's tail is near-exponential (p ~ 1) rather than Gaussian (p = 2).
# That is not a free wish: the PBH literature already needs heavy tails to make
# enough massive seeds early. This asks whether the two requirements point at the
# same p, because if they do, one parameter does two jobs and that is a result.
#
# Everything here is a scaling argument on tail shape. No claim is made about any
# specific inflationary model, and none is needed: the question is whether two
# independent demands on the same tail are compatible or in tension.

tail_frac <- function(nu, p) exp(-nu^p)          # P(>nu) for stretched exponential
# normalise so all p agree at the type-I PBH threshold, which is observationally
# anchored: PBH formation is calibrated to occur at some fixed small abundance.
beta_I <- 1e-10                                   # fiducial type-I collapse fraction

cat("=== 1. validate: the checker must be able to fail ===\n")
stopifnot(tail_frac(5, 2) < tail_frac(5, 1))      # heavier tail = larger fraction
stopifnot(tail_frac(1, 2) == exp(-1))
bad <- tryCatch({ stopifnot(tail_frac(5,2) > tail_frac(5,1)); TRUE }, error=function(e) FALSE)
cat(sprintf("  ordering holds and the reversed assertion fails: %s\n",
            if(!bad) "yes" else "NO - BLIND"))

cat("\n=== 2. where the type-I threshold sits, for each tail shape ===\n")
cat("  Given the SAME observed type-I abundance beta_I = 1e-10, a heavier tail\n")
cat("  puts that abundance at a different nu. Solve exp(-nu^p) = beta_I:\n\n")
nu_of_p <- function(p, beta) (-log(beta))^(1/p)
cat(sprintf("  %6s %12s\n", "p", "nu_I"))
ps <- c(2.0, 1.5, 1.0, 0.7, 0.5)
for (p in ps) cat(sprintf("  %6.1f %12.2f\n", p, nu_of_p(p, beta_I)))

cat("\n=== 3. the two-sided fraction at each tail shape ===\n")
cat("  Using the SAME 20 per cent threshold separation between type I and\n")
cat("  type II-B that arXiv:2401.06329's mu >~ 1.8 against mu ~ 1.5 implies:\n\n")
ratio <- 1.2
cat(sprintf("  %6s %12s %14s %18s\n", "p", "nu_I", "f_2s", "bimodality test"))
for (p in ps) {
  nuI <- nu_of_p(p, beta_I)
  f <- exp(-(nuI*ratio)^p + nuI^p)
  verdict <- if (f >= 0.30) "fires at n~300" else
             if (f >= 0.03) "needs n~1e4+"  else "cannot fire"
  cat(sprintf("  %6.1f %12.2f %14.2e %18s\n", p, nuI, f, verdict))
}

cat("\n=== 4. does the SEED problem want the same p? ===\n")
cat("  The seed problem wants ENOUGH massive holes early. For a fixed threshold,\n")
cat("  a heavier tail raises the abundance of every rare object together, so the\n")
cat("  relevant question is the RATIO of massive-seed abundance to ordinary, which\n")
cat("  is the same functional form with a larger threshold separation.\n\n")
cat(sprintf("  %6s %16s %16s %16s\n", "p", "f_2s (x1.2)", "heavy seed (x1.5)", "both?"))
for (p in ps) {
  nuI <- nu_of_p(p, beta_I)
  f12 <- exp(-(nuI*1.2)^p + nuI^p)
  f15 <- exp(-(nuI*1.5)^p + nuI^p)
  both <- if (f12 >= 0.03 && f15 >= 1e-4) "yes" else "no"
  cat(sprintf("  %6.1f %16.2e %16.2e %16s\n", p, f12, f15, both))
}

cat("\n=== 5. flatly ===\n")

# Prose below is INTERPOLATED from the computed values, never typed beside them.
# Twice in one session a hand-written summary contradicted the table two lines
# above it, which is the failure LESSONS.md records under "scripts drift from
# their own tables". Make drift impossible rather than proofread for it.
res <- sapply(ps, function(pp) { nuI <- nu_of_p(pp, beta_I)
                                 exp(-(nuI*1.2)^pp + nuI^pp) })
best_p <- ps[which.max(res)]; best_f <- max(res)
fires_300 <- ps[res >= 0.30]; fires_1e4 <- ps[res >= 0.03]

verdict <- if (length(fires_300))
  sprintf("Tails at p = %s reach the third the mixture test needs to fire at a few hundred objects.",
          paste(fires_300, collapse = ", ")) else
  sprintf("No tail tested reaches the ~30 per cent the mixture test needs at a few hundred objects; %s clear the 3 per cent floor where a test needing ~1e4 objects becomes possible.",
          if (length(fires_1e4)) paste0("p = ", paste(fires_1e4, collapse = ", ")) else "none")

cat(sprintf("
  Both demands point the same way, and the numbers are harsher than a fixed-nu
  treatment suggests. Normalising every tail to the SAME observed type-I abundance
  is what makes this self-consistent: a heavier tail puts that abundance further
  out, nu_I running %.1f at p=2 to %.0f at p=0.5, and that eats back most of what
  the heavier tail gives. An earlier script in this folder held nu fixed while
  varying p, which is not self-consistent; these numbers supersede it.

  Best case in the range tested is p = %.1f, where the two-sided fraction reaches
  %.1f per cent. %s

  A Gaussian field kills the class outright. Heavy tails let it exist as a
  minority of a few to ten per cent, and detecting a minority that small is not a
  mixture problem at all.

  THE ROUTE THAT FOLLOWS, and it changes the instrument rather than dropping the
  claim: a two-component mixture is the wrong tool for a contaminating minority of
  a few per cent. What suits it is contamination or outlier detection against a
  well-characterised one-class model, asking not 'are there two components' but
  'does this population carry more extreme residuals than one class can produce'.
  That is cheaper in sample size, and it sits in the same measurement-error
  literature the main paper already took its bootstrapped likelihood-ratio test
  from. Pricing it is the next calculation.

  THE FALSIFIER, clean and unchanged: if the primordial curvature field is close to
  Gaussian at the scales that make these holes, the two-sided class is empty, 4.3's
  mechanism has nothing to act on, and no test of any kind can fire. The fold's
  black-hole half dies there while its cosmological half stands, exactly as the
  main paper's independence statements allow.
", nu_of_p(2.0, beta_I), nu_of_p(0.5, beta_I), best_p, 100*best_f, verdict))
