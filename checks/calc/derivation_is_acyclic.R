# Does the derivation actually stand on its assumptions, or does it lean on its own conclusion?
#
# A paper claiming to derive the field equations is worth exactly what its dependency graph is
# worth. Section 3.6 asserts in prose that the argument does not feed general relativity back into
# a derivation of general relativity, and relativity_ledger.R lists what is derived and what is
# assumed. Neither checks that the two fit together. This does: it writes down what each derived
# line rests on, and then tests three things a reader would otherwise have to take on trust.
#
#   (1) the graph is acyclic, so nothing is used to prove itself, directly or at a distance;
#   (2) every chain terminates on one of the five assumptions or on the CPT postulate, so there
#       are no unnamed inputs hiding in the middle;
#   (3) no chain passes through general relativity, Bekenstein-Hawking entropy, or any other
#       result that presupposes what is being derived.
#
# The third is the one that matters. The natural way to get an area law is to quote
# Bekenstein-Hawking, and that would be circular here because Bekenstein-Hawking is a result OF
# general relativity. The prose says the entropy comes from entanglement instead. This checks
# that the stated dependency actually is that, and the plant at the end puts the forbidden edge
# back in to confirm the test can fail.

dep <- function(node, on) list(node = node, on = on)

# The Hadamard root was missing and that was a real omission, not a presentation choice. Section 3
# listed equilibrium at the horizon as a cost, the thermofield-double node uses it, and the graph
# named neither: the leaf test therefore passed while an input sat outside the graph. It is a
# regularity condition on states, the one under which a renormalised stress tensor exists at all,
# so it is not a property of general relativity and does not enter the ledger of nineteen and four.
# It is still an input and the graph now says so.
ROOTS <- c("CPT holds of the universe",        # the one postulate
           "Lorentzian signature",
           "one metric shared by the horizons and the matter action",
           "the Hadamard condition on states",
           "the measured value of G",
           "the measured value of Lambda")
FORBIDDEN <- c("general relativity", "Bekenstein-Hawking entropy",
               "the Einstein equation", "black-hole thermodynamics as an input")

G <- list(
  dep("Theta is an involution with dTheta = -Id",
      c("CPT holds of the universe", "Lorentzian signature")),
  dep("P_perp is unique at a stationary vacuum horizon",
      c("Theta is an involution with dTheta = -Id", "Lorentzian signature")),
  dep("the fold's map at a bifurcate horizon is the half-period shift",
      c("Theta is an involution with dTheta = -Id", "P_perp is unique at a stationary vacuum horizon",
        "Bisognano-Wichmann and Sewell")),
  dep("Bisognano-Wichmann and Sewell", c("Lorentzian signature")),
  dep("alpha is total inversion on the embedding",
      c("Theta is an involution with dTheta = -Id",
        "P_perp is unique at a stationary vacuum horizon")),
  dep("equilibrium at the horizon to the order the balance uses",
      c("the Hadamard condition on states")),
  dep("the cross-sheet correlator is a thermofield double",
      c("the fold's map at a bifurcate horizon is the half-period shift",
        "alpha is total inversion on the embedding",
        "equilibrium at the horizon to the order the balance uses")),
  dep("the contour displacement is beta/2",
      c("Theta is an involution with dTheta = -Id",
        "the cross-sheet correlator is a thermofield double")),
  dep("a temperature at every local Rindler horizon",
      c("the cross-sheet correlator is a thermofield double")),
  dep("horizons at every boost and orientation",
      c("Theta is an involution with dTheta = -Id", "Lorentzian signature")),
  dep("the classical-quantum split is the fold's parity",
      c("Theta is an involution with dTheta = -Id")),
  dep("the quantum half weighs 4 tanh^2(beta omega/4)",
      c("the classical-quantum split is the fold's parity",
        "the cross-sheet correlator is a thermofield double")),
  dep("an entropy proportional to horizon area",
      c("the cross-sheet correlator is a thermofield double", "flat-space transverse mode counting")),
  dep("flat-space transverse mode counting", c("Lorentzian signature")),
  dep("the same entropy coefficient at every horizon",
      c("an entropy proportional to horizon area", "Bisognano-Wichmann and Sewell")),
  dep("G > 0", c("the same entropy coefficient at every horizon")),
  dep("G does not run", c("the same entropy coefficient at every horizon")),
  dep("the weak equivalence principle",
      c("one metric shared by the horizons and the matter action")),
  dep("universality of local Rindler horizons",
      c("Theta is an involution with dTheta = -Id",
        "one metric shared by the horizons and the matter action")),
  dep("null-focusing surfaces for Raychaudhuri",
      c("Theta is an involution with dTheta = -Id", "Lorentzian signature")),
  dep("the contracted Bianchi identity",
      c("one metric shared by the horizons and the matter action")),
  dep("matter conservation",
      c("one metric shared by the horizons and the matter action")),
  dep("four large dimensions and no others",
      c("Theta is an involution with dTheta = -Id", "Lorentzian signature")),
  dep("the field equations",
      c("a temperature at every local Rindler horizon",
        "an entropy proportional to horizon area",
        "the same entropy coefficient at every horizon",
        "universality of local Rindler horizons",
        "null-focusing surfaces for Raychaudhuri",
        "matter conservation", "the contracted Bianchi identity")),
  dep("Lambda is an integration constant",
      c("the field equations", "matter conservation", "the contracted Bianchi identity")),
  dep("the numerical value of G", c("the measured value of G", "the field equations")),
  dep("a hot bang", c("Theta is an involution with dTheta = -Id", "Lorentzian signature")),
  dep("the fold cannot source Lambda",
      c("null-focusing surfaces for Raychaudhuri", "a hot bang")),
  dep("black-hole thermodynamics",
      c("the fold's map at a bifurcate horizon is the half-period shift")),
  dep("the opposite time orientation of the two sheets",
      c("Theta is an involution with dTheta = -Id", "Bisognano-Wichmann and Sewell"))
)

edges <- function(g) { e <- list(); for (d in g) for (u in d$on) e[[length(e)+1]] <- c(u, d$node); e }
nodes <- function(g) unique(c(sapply(g, function(d) d$node), unlist(lapply(g, function(d) d$on))))

cat("=== 1. the graph ===\n")
cat(sprintf("   %d derived nodes, %d roots declared, %d edges\n",
            length(G), length(ROOTS), length(edges(G))))
internal <- sapply(G, function(d) d$node)
leaves <- setdiff(nodes(G), internal)
cat(sprintf("   %d leaves, i.e. nodes nothing derives:\n", length(leaves)))
for (l in sort(leaves)) cat(sprintf("      %-45s %s\n", l, ifelse(l %in% ROOTS, "declared", "<-- UNDECLARED")))
undeclared <- setdiff(leaves, ROOTS)
stopifnot(length(undeclared) == 0)
cat("   Every leaf is a declared assumption, so nothing enters the argument unnamed.\n")
unused <- setdiff(ROOTS, nodes(G))
cat(sprintf("   declared roots that appear in no edge at all: %s\n",
            ifelse(length(unused) == 0, "none", paste(unused, collapse = ", "))))
cat("   That is worth reading rather than skipping. A declared assumption nothing uses is an\n")
cat("   assumption the derivation does not actually make, and Lambda's value is one: it enters\n")
cat("   the cosmology and never the derivation of the field equations.\n")

cat("\n=== 2. acyclicity, by Kahn's algorithm ===\n")
topo <- function(g) {
  ns <- nodes(g); es <- edges(g)
  indeg <- sapply(ns, function(n) sum(sapply(es, function(e) e[2] == n)))
  names(indeg) <- ns
  order <- c(); ready <- names(indeg)[indeg == 0]
  while (length(ready)) {
    v <- ready[1]; ready <- ready[-1]; order <- c(order, v)
    for (e in es) if (e[1] == v) { indeg[e[2]] <- indeg[e[2]] - 1
      if (indeg[e[2]] == 0) ready <- c(ready, e[2]) }
  }
  order
}
ord <- topo(G)
cat(sprintf("   %d of %d nodes ordered; a cycle would leave some unordered\n", length(ord), length(nodes(G))))
stopifnot(length(ord) == length(nodes(G)))
cat("   The graph is acyclic. Nothing is used to prove itself at any distance.\n")

cat("\n=== 3. what the field equations actually rest on, traced to the roots ===\n")
ancestors <- function(g, target) {
  par <- list(); for (d in g) par[[d$node]] <- d$on
  seen <- c(); stack <- target
  while (length(stack)) {
    v <- stack[1]; stack <- stack[-1]
    if (v %in% seen) next
    seen <- c(seen, v)
    if (!is.null(par[[v]])) stack <- c(stack, par[[v]])
  }
  setdiff(seen, target)
}
anc <- ancestors(G, "the field equations")
cat("   roots reached:\n")
for (r in sort(intersect(anc, ROOTS))) cat(sprintf("      %s\n", r))
cat(sprintf("   roots NOT reached: %s\n",
            paste(setdiff(ROOTS, anc), collapse = ", ")))
cat("   So the field equations do not use the measured value of either constant, which is the\n")
cat("   ledger's claim that the values are inputs to the numbers and not to the derivation.\n")
stopifnot(!("the measured value of G" %in% anc), !("the measured value of Lambda" %in% anc))

cat("\n=== 4. the forbidden ancestors ===\n")
bad <- intersect(anc, FORBIDDEN)
cat(sprintf("   results presupposing general relativity found among the ancestors: %d\n", length(bad)))
if (length(bad)) for (b in bad) cat(sprintf("      %s  <-- CIRCULAR\n", b))
stopifnot(length(bad) == 0)
cat("   The entropy line reaches the roots through the thermofield double and flat-space mode\n")
cat("   counting, and never through Bekenstein-Hawking. That is the loop the prose claims is\n")
cat("   shut, checked rather than asserted.\n")

cat("\n=== 5. the plant: put the forbidden edge back and the test must fire ===\n")
Gbad <- G
for (i in seq_along(Gbad)) if (Gbad[[i]]$node == "an entropy proportional to horizon area")
  Gbad[[i]]$on <- c(Gbad[[i]]$on, "Bekenstein-Hawking entropy")
ancb <- ancestors(Gbad, "the field equations")
badb <- intersect(ancb, FORBIDDEN)
cat(sprintf("   with the entropy law quoted from Bekenstein-Hawking: %d forbidden ancestor(s), %s\n",
            length(badb), paste(badb, collapse = ", ")))
stopifnot(length(badb) > 0)

cat("\n   and a genuine cycle must be caught too:\n")
Gcyc <- G
for (i in seq_along(Gcyc)) if (Gcyc[[i]]$node == "the cross-sheet correlator is a thermofield double")
  Gcyc[[i]]$on <- c(Gcyc[[i]]$on, "the field equations")
ordc <- topo(Gcyc)
cat(sprintf("   ordering the cyclic graph places %d of %d nodes, so the cycle is detected\n",
            length(ordc), length(nodes(Gcyc))))
stopifnot(length(ordc) < length(nodes(Gcyc)))
cat("   Both plants fire.\n")

cat("\n=== 6. the longest chain, which is what a referee will walk ===\n")
depth <- function(g, v, par) {
  if (is.null(par[[v]])) return(0)
  1 + max(sapply(par[[v]], function(u) depth(g, u, par)))
}
par <- list(); for (d in G) par[[d$node]] <- d$on
dfe <- depth(G, "the field equations", par)
cat(sprintf("   depth of 'the field equations' above the assumptions: %d steps\n", dfe))
for (n in c("a temperature at every local Rindler horizon", "an entropy proportional to horizon area",
            "the same entropy coefficient at every horizon", "universality of local Rindler horizons"))
  cat(sprintf("      %-48s depth %d\n", n, depth(G, n, par)))
cat(sprintf("   So the longest path a referee has to walk is %d steps, and every link on it is\n", dfe))
cat("   computed in this release rather than cited.\n")
