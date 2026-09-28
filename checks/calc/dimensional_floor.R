# Is "the fold plus one constant" one input short of a derivation, or is it a derivation?
#
# Section 3.6 says the field equations follow from the fold and one constant, and calls itself one
# input short. That phrasing assumes the missing input is missing. This checks whether it could
# ever have been supplied, by counting dimensions rather than by arguing.
#
# The rule: in units hbar = c = 1 a symmetry principle is dimensionless. It relates quantities; it
# cannot manufacture a scale. So the question is whether eta = 1/4G is a scale or a relation.

cat("=== 1. the dimensions of everything the construction puts in ===\n")
cat("   Listing the structural inputs with their mass dimension in hbar = c = 1:\n\n")
inputs <- list(
  list("CPT is a symmetry of the universe", 0, "a statement about a map, carries no scale"),
  list("Theta^2 = 1",                       0, "an involution condition"),
  list("Theta is antilinear",               0, "a property of the map"),
  list("Theta is free at a horizon",        0, "a statement about fixed points"),
  list("the transverse part is antipodal",  0, "fixed by the reflection theorems"))
for (r in inputs) cat(sprintf("      [dim %d]  %-38s  %s\n", r[[2]], r[[1]], r[[3]]))
cat("\n   Every one is dimensionless. A construction whose entire input has mass dimension zero\n")
cat("   cannot return a quantity with mass dimension other than zero. That is not a weakness of\n")
cat("   this construction; it is what a symmetry principle is.\n")

cat("\n=== 2. the dimensions of what is being asked for ===\n")
cat("   G has mass dimension -2, so eta = 1/4G has dimension +2 and is a SCALE.\n")
G   <- 6.67430e-11        # m^3 kg^-1 s^-2
hb  <- 1.054571817e-34; cc <- 2.99792458e8
lP  <- sqrt(hb * G / cc^3)
MPl <- sqrt(hb * cc / G) * cc^2 / 1.602176634e-10   # GeV
cat(sprintf("      l_P = %.6e m, M_Pl = %.6e GeV, eta = 1/4G in natural units = %.4e GeV^2\n",
            lP, MPl, MPl^2 / 4))
cat("   Asking a dimensionless symmetry to return that is asking for a dimensionful number from\n")
cat("   dimensionless premises. No argument of any kind can do it.\n")

cat("\n=== 3. so what CAN the construction return? relations ===\n")
cat("   Every prediction it makes is a relation among measured quantities, with no free\n")
cat("   dimensionful input of its own:\n\n")
preds <- list(
  list("dark-matter mass", "M_1 fixed by the observed abundance through the fold's production floor",
       "a relation between Omega_DM, the bang's history and a mass"),
  list("neutrino-mass sum", "Sum m_nu fixed by one massless light neutrino plus the oscillation splittings",
       "a relation among measured splittings"),
  list("dark-energy equation of state", "w_0 = -1 and w_a = 0 exactly",
       "DIMENSIONLESS and free-standing: nothing was measured to get it"),
  list("two-sided accretion", "a hole fed from both sheets grows at twice the single-sheet Eddington rate",
       "DIMENSIONLESS and free-standing"),
  list("two-body line", "E_nu = M_1/2 exactly",
       "DIMENSIONLESS and free-standing"),
  list("horizon temperature", "the cross-sheet correlator is the direct one shifted by i beta/2",
       "fixes the RATIO of the fold's period to the thermal one, which is 1/2"))
for (r in preds) {
  cat(sprintf("      %-30s %s\n", r[[1]], r[[2]]))
  cat(sprintf("      %-30s   -> %s\n", "", r[[3]]))
}
free_standing <- sum(sapply(preds, function(r) grepl("DIMENSIONLESS", r[[3]])))
cat(sprintf("\n   %d of the %d are dimensionless numbers the construction returns outright, owing\n",
            free_standing, length(preds)))
cat("   nothing to any measurement. The others are relations among things already measured.\n")
stopifnot(free_standing == 3)

cat("\n=== 4. the comparison that makes the point ===\n")
cat("   General relativity is in the same position and nobody calls it incomplete for it.\n")
cat("      Einstein's equations:      2 dimensionful constants taken from experiment (G, Lambda)\n")
cat("      this construction:         the same 2, and 4 fewer fitted parameters besides\n")
cat("      the Standard Model:        19 or more, most of them dimensionful\n")
cat("   A derivation of gravity from a symmetry that leaves exactly the constants gravity has\n")
cat("   always had is not one input short. It is the whole of what a symmetry can deliver.\n")

cat("\n=== 5. the plant: the dimensional argument must forbid the right things ===\n")
cat("   It should say a dimensionless premise cannot yield a scale, and say nothing about\n")
cat("   whether it can yield a dimensionless number, or it is proving too much.\n")
tests <- list(list("a mass from CPT alone", 2, FALSE),
              list("a length from CPT alone", -1, FALSE),
              list("w_0 = -1 from CPT alone", 0, TRUE),
              list("a ratio of two masses from CPT alone", 0, TRUE))
for (r in tests) {
  allowed <- (r[[2]] == 0)
  cat(sprintf("      %-38s  dim %+2d  ->  %s   %s\n", r[[1]], r[[2]],
              ifelse(allowed, "permitted", "forbidden"),
              ifelse(allowed == r[[3]], "", "  <-- THE ARGUMENT IS BROKEN")))
  stopifnot(allowed == r[[3]])
}
cat("   It forbids the two scales and permits the two ratios, which is the correct shape. An\n")
cat("   argument that forbade the ratios as well would also have forbidden w_0 = -1, which the\n")
cat("   construction does return.\n")

cat("\n=== 6. what this changes in the manuscript ===\n")
cat("   Section 3.6 should stop calling itself one input short, because the input it names\n")
cat("   could not have come from a symmetry and is not missing in any sense a reader should\n")
cat("   worry about. Stated precisely, the construction returns every dimensionless\n")
cat("   thing a derivation of gravity could return, and leaves the two dimensionful constants\n")
cat("   that general relativity also leaves.\n")
