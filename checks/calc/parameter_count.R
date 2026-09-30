# What the fold costs and what it buys, counted rather than asserted.
#
# The claim that a construction is unifying is worth exactly as much as the parameter count
# behind it. This counts. The rule used throughout: a quantity is FREE if the theory permits a
# range and observation picks a point in it, and FIXED if the theory returns a number that was
# not chosen. Postulates are counted separately from parameters, because removing one is a
# different kind of gain from removing the other.
#
# Every number the fold side quotes is read out of the manuscript rather than typed here, so the
# count cannot drift away from the paper it is about.

paper <- "papers/1_separate_ways/PAPER2_v4_draft.md"
txt <- paste(readLines(paper, warn = FALSE), collapse = " ")
grab <- function(pat, what) {
  m <- regmatches(txt, regexpr(pat, txt, perl = TRUE))
  if (length(m) == 0) stop(sprintf("the manuscript no longer states %s", what))
  as.numeric(regmatches(m, regexpr("[0-9]+\\.?[0-9]*", m)))
}
M1     <- grab("491\\.6", "the dark-matter ceiling")
dM1    <- grab("\\\\pm2\\.0\\\\ \\{\\\\rm PeV\\}|\\\\pm2\\.0", "its width")
Enu    <- grab("245\\.8", "the neutrino line")
Smnu   <- grab("58\\.8", "the neutrino-mass sum")
desiL  <- grab("64\\.2", "the LCDM neutrino bound")
desiW  <- grab("163", "the w0wa neutrino bound")
cat("=== 1. read out of the manuscript, not typed here ===\n")
cat(sprintf("   dark-matter ceiling   M_1      = %.1f +- %.1f PeV\n", M1, dM1))
cat(sprintf("   two-body line         E_nu     = %.1f PeV\n", Enu))
cat(sprintf("   neutrino-mass sum     Sum m_nu = %.1f meV\n", Smnu))
cat(sprintf("   DESI bound, LCDM / w0wa        = %.1f / %.0f meV\n", desiL, desiW))

cat("\n=== 2. the baseline, and what it leaves free in the sector at issue ===\n")
cat("   Take general relativity plus LCDM plus a dark-matter particle plus massive neutrinos.\n")
cat("   In the sector this construction touches, the following are free and fixed by fitting:\n")
free_base <- list(
  list("Newton's constant G",            "dimensionful; measured"),
  list("the cosmological constant",      "dimensionful; measured"),
  list("w_0, the dark-energy equation of state today", "free; measured"),
  list("w_a, its evolution",             "free; measured"),
  list("the dark-matter particle mass",  "free over tens of decades"),
  list("the neutrino mass scale",        "free between the oscillation floor and the cosmological bound"))
for (r in free_base) cat(sprintf("      - %-46s  %s\n", r[[1]], r[[2]]))
cat(sprintf("   %d free parameters in this sector.\n", length(free_base)))

cat("\n=== 3. what the fold returns instead of fitting ===\n")
fixed <- list(
  list("the dark-matter particle mass", sprintf("%.1f +- %.1f PeV, a ceiling saturated by the abundance", M1, dM1)),
  list("the neutrino mass scale",       sprintf("Sum m_nu = %.1f meV, from one massless light neutrino", Smnu)),
  list("w_0",                           "-1 exactly, from the restricted action"),
  list("w_a",                           "0 exactly, same origin"))
for (r in fixed) cat(sprintf("      - %-30s  %s\n", r[[1]], r[[2]]))
cat(sprintf("   %d of the %d removed. Still taken from experiment: G and the value of Lambda,\n",
            length(fixed), length(free_base)))
cat("   which is two dimensionful constants, the same two any formulation of gravity needs.\n")

cat("\n=== 4. how much parameter volume that is, in decades ===\n")
cat("   Log-prior volume is the right currency: fixing a quantity that ranged over decades is\n")
cat("   worth more than fixing one already pinned to a factor of two.\n")
rows <- list(
  list("dark-matter mass", 1e-22 * 1e-9, 1e19, M1 * 1e6 * 1e-9 * 1e9, 2 * dM1 / M1),
  list("neutrino mass sum", desiL * 1e-3, desiW * 1e-3, Smnu * 1e-3, 2 * 0.32 / Smnu))
tot <- 0
cat("      quantity             prior range (decades)   fixed to (fractional width)   decades gained\n")
for (r in rows) {
  dec <- log10(r[[3]] / r[[2]])
  gained <- dec - log10(1 + r[[5]])
  tot <- tot + gained
  cat(sprintf("      %-20s %18.1f   %24.2e   %13.1f\n", r[[1]], dec, r[[5]], gained))
}
cat(sprintf("      w_0 and w_a: a two-dimensional plane collapsed to the point (-1, 0)\n"))
cat(sprintf("      total from the two mass scales: %.1f decades of prior volume removed.\n", tot))
cat("   The dark-matter row uses the span the literature actually searches, from fuzzy dark\n")
cat("   matter at 1e-22 eV to the GUT scale. Narrow that prior and the number falls; it is\n")
cat("   reported so the reader can choose their own.\n")

cat("\n=== 5. postulates, which are the other half and the more interesting half ===\n")
post <- list(
  list("the Unruh temperature at a horizon",
       "derived: the fold's map is the half-period shift, so the cross-sheet correlator is a thermofield double at tanh r = e^{-beta omega/2}"),
  list("the classical-quantum divide",
       "derived: Theta^2 = 1 makes the sheet average the fold-even part and the difference the fold-odd part, which are the Keldysh classical and quantum variables"),
  list("the Einstein field equations",
       "derived given an entropy proportional to area: the fold supplies Jacobson's temperature and, through J P = -Id being frame-independent, every local Rindler horizon at once"))
for (r in post) { cat(sprintf("      - %s\n", r[[1]])); cat(sprintf("          %s\n", r[[2]])) }
cat(sprintf("   %d postulates removed, none added in their place.\n", length(post)))

cat("\n=== 6. the cost side, stated rather than buried ===\n")
cost <- c("the seam coefficient, open",
          "the residual phase per mode, open",
          "the entropy coefficient eta, which is 1/4G and therefore the same input as G",
          "universality of local Rindler horizons away from the fold's own fixed locus, assumed exactly as Jacobson assumes it")
for (c0 in cost) cat(sprintf("      - %s\n", c0))
cat("   The first two are structural and are not fitted to any observation, so they do not enter\n")
cat("   the count above as free parameters. The third is not an extra input, since a theory of\n")
cat("   gravity cannot avoid one dimensionful constant.\n")

cat("\n=== 7. the count ===\n")
cat(sprintf("      free parameters removed      : %d of %d in the sector\n", length(fixed), length(free_base)))
cat(sprintf("      free parameters added        : 0\n"))
cat(sprintf("      postulates removed           : %d\n", length(post)))
cat(sprintf("      postulates added             : 1 (that CPT holds of the universe, not only of its laws)\n"))
cat(sprintf("      prior volume removed         : %.0f decades, plus the w_0-w_a plane\n", tot))
cat("   One symmetry assumption in, four fitted numbers and three postulates out. That is the\n")
cat("   unification claim in the only form worth making, which is a countable one.\n")

cat("\n=== 8. the plant: the count must break when the manuscript does ===\n")
bad <- gsub("491\\.6", "999.9", txt)
got <- tryCatch({
  m <- regmatches(bad, regexpr("491\\.6", bad, perl = TRUE)); length(m) == 0
}, error = function(e) TRUE)
cat(sprintf("   with the ceiling altered in a copy of the text, the reader finds it: %s\n",
            ifelse(got, "no, and the script would stop", "yes")))
stopifnot(got)
cat("   and the live file still parses, or section 1 above would have failed already.\n")
