# How much of general relativity the construction produces, itemised, with the assumptions named.
#
# A claim to derive a theory can be checked only line by line. This lists what general relativity
# consists of, marks each line derived or assumed, and says by what. No line is marked derived
# unless a script in this release produces it.
#
# 2026-10-06. Each derived line now also says whether it is new. Jacobson's 1995 paper
# (gr-qc/9504004) already took the horizon temperature from the Unruh effect, the area law from
# vacuum entanglement with a fixed cutoff (Bombelli, Koul, Lee and Sorkin 1986), and local Rindler
# horizons through every point from local flatness, and it wrote Lambda as the integration constant
# of the Bianchi step and G as 1/(4 hbar eta). By Bisognano-Wichmann the wedge conjugation is
# CPT times a rotation by pi in the transverse plane, so J composed with P_perp is CPT itself. The
# manuscripts had said those inputs were assumed in 1995 and supplied here. They were supplied in
# 1995, by the same quantum field theory. So every derived line carries a group: "Jacobson 1995"
# (in his argument or immediate from the inputs he cited), "standard" (textbook general relativity
# or quantum field theory in curved spacetime), or "fold" (needs the fold's global structure).

item <- function(what, status, by) list(what = what, status = status, by = by)

L <- list(
  item("the field equations' FORM, R_ab - R g_ab/2 + Lambda g_ab proportional to T_ab",
       "derived", "Clausius on the fold's horizons; Raychaudhuri integrated, the null-vector algebra checked, the Bianchi step written out"),
  item("a temperature at every local Rindler horizon",
       "derived", "the fold's map is the half-period shift, so the cross-sheet correlator is a thermofield double"),
  item("those horizons, at every boost and orientation",
       "derived", "the wedge reflection and the transverse antipode are frame-dependent and their composition -Id is not"),
  item("an entropy proportional to the horizon area",
       "derived", "the entanglement entropy of the state the fold forces, with flat-space transverse mode counting; Bekenstein-Hawking never invoked"),
  item("the SAME entropy coefficient at every horizon, which is what makes the result Einstein's",
       "derived", "the surface gravity cancels between the occupation and the WKB mode density because the wedge metric in boost coordinates has none and Bisognano-Wichmann fixes the modular temperature at 2 pi there; S/A holds to 1e-9 while beta runs over four decades, and the cancellation survives three different transverse weightings. One condition: the cutoff must be a single proper length, since one tracking kappa moves S/A by a factor of fifty"),
  item("the transverse involution P_perp at a black hole",
       "derived", "the only free involutive isometry of a round bifurcation sphere is -Id, and of Kerr's axisymmetric surface is the antipodal map; by the uniqueness theorems those exhaust the stationary vacuum case"),
  item("that Newton's constant does not run with epoch or location",
       "derived", "eta is fixed by the modular temperature and a transverse mode count, neither of which is cosmological, so the construction has no dial that could make G vary"),
  item("the SIGN of Newton's constant, so that gravity attracts",
       "derived", "8 pi G = 2 pi / eta and eta is an entropy density, which is positive"),
  item("Lambda appearing at all, as an integration constant",
       "derived", "the contracted Bianchi identity plus matter conservation leave exactly one constant of integration, and both of those are themselves derived"),
  item("null-focusing surfaces for Raychaudhuri to act on",
       "derived", "contact is conjugacy; the fold's own contact surfaces are caustics, at every charge and dimension"),
  item("four large spacetime dimensions and no others",
       "derived", "the turning available for contact closes only at D = 4; contact needs a turning of pi because the antipode is an involution of a sphere, and the interior supplies pi/(D-3) per leg, which at D = 5 reaches pi only in a limit that is no point of the spacetime and from D = 6 falls short outright. What is derived is that four is the only dimension in which the sheets touch anywhere, so a higher-dimensional fold is allowed and has nothing to add"),
  item("the classical-quantum divide, and its amplitude",
       "derived", "the fold's parity, with weight 4 tanh^2(beta omega/4)"),
  item("black-hole thermodynamics",
       "derived", "the same half-period shift, at a bifurcate Killing horizon"),
  item("the opposite time orientation of the two sheets, which the cosmological accounts assume",
       "derived", "the modular flow of a wedge is the boost, whose generator gives dX0/ds = X1, positive at every point of the right static patch and negative at every point of its antipode; the two limits hold, since the flow preserves the state and the relation doing the work is the commutant's Delta_{A}^{-1} rather than any property of J"),
  item("a hot bang, with no cold component at the fold's fixed slice",
       "derived", "Theta invariance needs a odd in conformal time, which admits a fluid only where -3w is an odd integer and excludes dust; the radiation bang was previously assumed"),
  item("the VALUE of Newton's constant",
       "assumed", "dimensionally impossible for a symmetry to supply: every structural input has mass dimension zero"),
  item("the VALUE of Lambda",
       "assumed", "an integration constant of the derivation, and both places it could have come from are now closed: the contact budget is short at every finite radius in the late universe, and at the bang an a-independent image stress would be of order M_1^4, some 81.4 orders above the observed value"),
  item("Lorentzian signature",
       "assumed", "the fold is built on a spacetime that already has one, and it is the only signature in which the fold has content: -Id sits in the identity component of O(t,s) exactly when t and s are both even, so the fold keeps a parity exactly when the number of time directions is odd, which over the five four-dimensional signatures leaves (1,3) and (3,1) and excludes Euclidean, (2,2) and (4,0) alike"),
  item("one metric, carrying both the fold's horizons and the matter action",
       "assumed", "the metric postulate; the Clausius step needs it already, since T_ab means the variation of the matter action with respect to that metric and has no other meaning. Not a convenience either: with a second metric, conservation holds against the wrong connection and the mismatch is 4e5 times the numerical floor, while a constant rescaling, which changes no connection, costs 7e-15"),
  item("matter following the metric's geodesics, the weak equivalence principle",
       "derived", "the eikonal characteristics of a field on that metric are its geodesics, checked against the Christoffel system to 2.4e-12, and the mass cancels out of the path"),
  item("the contracted Bianchi identity, which is what leaves exactly one integration constant",
       "derived", "a theorem about the Levi-Civita connection of any metric, measured at 2.2e-5 relative on a generic curved metric and falling as h^2; a connection with torsion breaks it at order unity"),
  item("matter conservation, grad^a T_ab = 0",
       "derived", "Noether's second theorem applied to a covariant matter action: the divergence of the stress tensor is the matter field equation contracted with the field, checked off shell in a generic curved metric"),
  item("local Rindler horizons away from the fold's own fixed locus",
       "derived", "the geodesic symmetry at any point has the fold's differential -Id and is an isometry through second order, which is the order the balance is computed at; its failure is third order in grad R and its weight against the balance vanishes linearly in the patch size"))

# The header above says no line is marked derived unless a script in this release produces it.
# Until 2026-09-27 nothing checked that, which for the paper's central claim is the wrong way
# round. Each line now names the file, the files are asserted to exist, and gate 7 runs them all.
SRC <- c(
  "calc/jacobson_rederived.R",            # 1  the field equations' form
  "calc/half_period_match.R",             # 2  a temperature at every local Rindler horizon
  "calc/rindler_fold_identity.R",         # 3  those horizons, at every boost and orientation
  "calc/area_law_from_fold.R",            # 4  an entropy proportional to the area
  "calc/eta_is_universal.R",              # 5  the same coefficient at every horizon
  "calc/perp_is_unique.R",                # 6  the transverse involution at a black hole
  "calc/eta_is_universal.R",              # 7  G does not run: eta carries no epoch and no location
  "calc/area_law_from_fold.R",            # 8  the sign of G, from the entropy density's positivity
  "calc/lambda_as_boundary.R",            # 9  Lambda as an integration constant
  "calc/contact_conjugacy_general.R",     # 10 null-focusing surfaces for Raychaudhuri
  "calc/contact_dimension.R",             # 11 four dimensions and no others
  "calc/quantum_classical_amplitude.R",   # 12 the classical-quantum divide and its amplitude
  "calc/horizon_squeeze.py",              # 13 black-hole thermodynamics
  "calc/opposite_orientation.R",          # 14 the two sheets' opposite time orientation
  "calc/bang_parity_rule.R",              # 15 a hot bang
  "-",                                    # 16 the VALUE of G: an input, and no script can supply it
  "calc/bang_image_scale.R",              # 17 the VALUE of Lambda: both sources closed
  "calc/signature_or_the_fold_is_nothing.R",  # 18
  "calc/one_metric_is_forced.R",          # 19
  "calc/equivalence_from_one_metric.R",   # 20 the weak equivalence principle
  "calc/bianchi_is_a_theorem.R",          # 21 the contracted Bianchi identity
  "calc/conservation_from_diffeo.R",      # 22 matter conservation
  "calc/universality_local_fold.R")       # 23 local Rindler horizons away from the fixed locus
stopifnot(length(SRC) == length(L))
for (i in seq_along(L)) L[[i]]$src <- SRC[i]

# Where each derived line stood before this paper. "input" marks the four assumed lines.
GROUP <- c(
  "Jacobson 1995",  # 1  his eq. (6): Clausius, Raychaudhuri, null-vector algebra, Bianchi
  "Jacobson 1995",  # 2  the Unruh temperature; Bisognano-Wichmann is what makes it exact
  "Jacobson 1995",  # 3  local Poincare symmetry gives the local Rindler horizon in every null direction
  "Jacobson 1995",  # 4  entanglement entropy across the horizon, with a cutoff, citing BKLS 1986
  "Jacobson 1995",  # 5  one fundamental cutoff length, so area in units of l_c^2 everywhere
  "fold",           # 6  the transverse involution at a black hole is the fold's choice of map
  "Jacobson 1995",  # 7  eta a constant fixed by the cutoff, so G = 1/(4 hbar eta) is constant
  "Jacobson 1995",  # 8  G = 1/(4 hbar eta) with eta an entropy density
  "Jacobson 1995",  # 9  "f = -R/2 + Lambda for some constant Lambda"
  "fold",           # 10 the fold's own contact surfaces
  "fold",           # 11 contact closes only at D = 4
  "fold",           # 12 the fold's parity as the Keldysh split, and its weight
  "standard",       # 13 Hawking 1975, Israel 1976, Sewell 1982 (from PCT), Kay-Wald 1991
  "fold",           # 14 the sheets' opposite time orientation
  "fold",           # 15 Theta invariance at the bang excludes dust
  "input", "input", "input", "input",
  "standard",       # 20 eikonal limit: a field on a metric follows its geodesics
  "standard",       # 21 the contracted Bianchi identity is a theorem
  "standard",       # 22 Noether's second theorem; Jacobson assumed it outright
  "Jacobson 1995")  # 23 local flatness gives local Rindler horizons at every point
stopifnot(length(GROUP) == length(L))
for (i in seq_along(L)) L[[i]]$group <- GROUP[i]
GROUPS <- c("Jacobson 1995", "standard", "fold")
stopifnot(all(sapply(L, function(x) (x$status == "assumed") == (x$group == "input"))),
          all(sapply(L, function(x) x$group %in% c(GROUPS, "input"))))

cat("=== the ledger ===\n\n")
der <- Filter(function(x) x$status == "derived", L)
asu <- Filter(function(x) x$status == "assumed", L)
cat(sprintf("DERIVED (%d)\n", length(der)))
for (x in der) { cat(sprintf("   %s\n", x$what)); cat(sprintf("        by: %s\n", x$by))
                 cat(sprintf("        in: %s\n", x$src)); cat(sprintf("        before: %s\n", x$group)) }
cat(sprintf("\nASSUMED (%d)\n", length(asu)))
for (x in asu) { cat(sprintf("   %s\n", x$what)); cat(sprintf("        why: %s\n", x$why <- x$by)) }

cat("\n=== every derived line names a file, and every file is there ===\n")
miss <- character(0)
for (x in L) {
  if (x$src == "-") { stopifnot(x$status == "assumed"); next }
  if (!file.exists(paste0("checks/", x$src))) miss <- c(miss, x$src)
}
if (length(miss)) { cat("   MISSING:", paste(unique(miss), collapse=", "), "\n"); stop("a ledger line names a file that is not there") }
cat(sprintf("   %d derived lines, %d distinct files, all present\n",
            length(der), length(unique(sapply(der, function(x) x$src)))))
cat("   the plant: a line naming a file that is not in the release must stop this script.\n")
cat(sprintf("      asking for calc/no_such_script.R -> %s\n",
            ifelse(file.exists("checks/calc/no_such_script.R"), "would pass", "STOPS, as it must")))
stopifnot(!file.exists("checks/calc/no_such_script.R"))

cat(sprintf("\n=== the count: %d derived, %d assumed ===\n", length(der), length(asu)))
nstruct <- 2
cat(sprintf("   Of the %d assumed, %d are not assumptions about gravity at all but about what it\n",
            length(asu), nstruct))
cat("   means to write a metric theory: a Lorentzian signature and one metric carrying both the\n")
cat("   fold's horizons and a covariant matter action. The remaining two are the dimensionful\n")
cat("   constants, which no symmetry can produce, for the reason the irreducibility argument\n")
cat("   gives: every structural input has mass dimension zero.\n")
cat("   Two entries that used to sit here have moved. The equivalence principle went, because\n")
cat("   its weak form is what the eikonal limit of a field on the shared metric does. And\n")
cat("   diffeomorphism invariance went, because its two uses have different status: the\n")
cat("   contracted Bianchi identity is a theorem about any metric, and the covariance of the\n")
cat("   matter action is part of the metric postulate rather than beside it. Nothing\n")
cat("   substantive about gravity is left on the assumed side.\n")

cat("\n=== what is new, line by line ===\n")
ng <- sapply(GROUPS, function(g) sum(sapply(der, function(x) x$group == g)))
for (g in GROUPS) {
  cat(sprintf("   %-14s %2d\n", g, ng[[g]]))
  for (x in der) if (x$group == g) cat(sprintf("      %s\n", x$what))
}
stopifnot(sum(ng) == length(der), all(ng > 0))
cat("   Jacobson's 1995 argument already had the temperature, the area law, the horizons and the\n")
cat("   integration constant, from the quantum field theory he cited. What the fold does with\n")
cat("   those lines is trace them to one map: by Bisognano-Wichmann the wedge conjugation\n")
cat("   composed with the transverse antipode is CPT, so the same involution that defines the\n")
cat("   folded cosmology is the one behind Jacobson's two inputs. The fold's own lines are the\n")
cat("   ones that need a second sheet: the black-hole involution, the contact surfaces, four\n")
cat("   dimensions, the classical-quantum split, the sheets' time orientation and the hot bang.\n")
cat("   The non-local image stress is the fold's too, and sits in the field equations, not here.\n")

cat("\n=== the reading ===\n")
# TWO STALE LINES FIXED 2026-09-28. The count said three of the four assumed lines are what a
# metric theory is, against the two the list actually carries and the two the section above
# states; and the sentence naming what is untouched still said the equivalence principle is
# assumed, after its weak form moved onto the derived list. Both contradicted this file's own
# output, and both are counted rather than typed now so the first cannot recur.
n_units <- sum(grepl("^the VALUE of", sapply(asu, function(a) a$what)))
n_kin   <- length(asu) - n_units
cat(sprintf("   %d lines of general relativity come out; %d go in, of which %d are what a metric\n",
            length(der), length(asu), n_kin))
cat(sprintf("   theory is and %d are units. Of the %d, %d were in Jacobson's 1995 argument and %d are\n",
            n_units, length(der), ng[["Jacobson 1995"]], ng[["standard"]]))
cat(sprintf("   standard results; %d need the fold. What it does not touch: special relativity is\n",
            ng[["fold"]]))
cat("   input rather than output, and nothing here bears on the gauge structure of matter.\n")
stopifnot(n_kin == 2, n_units == 2)

cat("\n=== the plant ===\n")
cat("   Every line marked derived must name a script in this release, or the ledger is a wish\n")
cat("   list. Checking that each 'by' clause is non-empty and specific:\n")
thin <- sum(sapply(der, function(x) nchar(x$by) < 40))
cat(sprintf("      derived lines with a thin justification: %d\n", thin))
cat("   And the assumed list must be non-empty, or the ledger is not a ledger:\n")
cat(sprintf("      assumed lines: %d\n", length(asu)))
stopifnot(thin == 0, length(asu) >= 3, length(der) >= 10)
cat("   And a derived line filed under a group the paper does not name must stop the count:\n")
bad <- GROUP; bad[1] <- "new here"
cat(sprintf("      line 1 relabelled 'new here' -> %s\n",
            ifelse(all(bad %in% c(GROUPS, "input")), "would pass", "STOPS, as it must")))
stopifnot(!all(bad %in% c(GROUPS, "input")))
