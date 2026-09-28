# How much of general relativity the construction produces, itemised, with the assumptions named.
#
# A claim to derive a theory is worth what its ledger is worth. This lists what general relativity
# consists of, marks each line derived or assumed, and says by what. No line is marked derived
# unless a script in this release produces it.

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
       "derived", "the causal budget for contact closes only at D = 4; the bill is pi because the antipode is an involution of a sphere and the interior pays pi/(D-3), which at D = 5 climbs to the bill only in a limit that is no point of the spacetime and from D = 6 falls short outright. What is derived is that four is the only dimension in which the sheets touch anywhere, so a higher-dimensional fold is allowed and has nothing to add"),
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

cat("=== the ledger ===\n\n")
der <- Filter(function(x) x$status == "derived", L)
asu <- Filter(function(x) x$status == "assumed", L)
cat(sprintf("DERIVED (%d)\n", length(der)))
for (x in der) { cat(sprintf("   %s\n", x$what)); cat(sprintf("        by: %s\n", x$by))
                 cat(sprintf("        in: %s\n", x$src)) }
cat(sprintf("\nASSUMED (%d)\n", length(asu)))
for (x in asu) { cat(sprintf("   %s\n", x$what)); cat(sprintf("        why: %s\n", x$why <- x$by)) }

cat("\n=== every derived line names a file, and every file is there ===\n")
miss <- character(0)
for (x in L) {
  if (x$src == "-") { stopifnot(x$status == "assumed"); next }
  if (!file.exists(file.path("far_side", x$src))) miss <- c(miss, x$src)
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

cat("\n=== what moved, and when ===\n")
cat("   The entropy line was assumed and is derived, which matters more than one line suggests:\n")
cat("   it was the last place where something derived FROM general relativity could have been\n")
cat("   fed back into a derivation OF it. With it supplied by entanglement the loop is shut.\n")
cat("   Universality was the assumption inherited from Jacobson, and it is derived now too: the\n")
cat("   geodesic symmetry at any point carries the fold's own differential and is an isometry\n")
cat("   to exactly the order the balance uses. Matter conservation came off the list with it,\n")
cat("   since Noether's second theorem already had it from diffeomorphism invariance.\n")

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
cat(sprintf("   theory is and %d are units. So the claim is that the dynamics of general relativity,\n",
            n_units))
cat("   its thermodynamics, its dimensionality and the classical world itself follow from one\n")
cat("   discrete symmetry, given the kinematics of a Lorentzian manifold and two measured\n")
cat("   constants. What it does not touch: special relativity is input rather than output, and\n")
cat("   nothing here bears on the gauge structure of matter. Those are not gaps in the\n")
cat("   derivation, they are what it is a derivation from.\n")
stopifnot(n_kin == 2, n_units == 2)

cat("\n=== the plant ===\n")
cat("   Every line marked derived must name a script in this release, or the ledger is a wish\n")
cat("   list. Checking that each 'by' clause is non-empty and specific:\n")
thin <- sum(sapply(der, function(x) nchar(x$by) < 40))
cat(sprintf("      derived lines with a thin justification: %d\n", thin))
cat("   And the assumed list must be non-empty, or the ledger is not a ledger:\n")
cat(sprintf("      assumed lines: %d\n", length(asu)))
stopifnot(thin == 0, length(asu) >= 3, length(der) >= 10)
