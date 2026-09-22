# ==========================================================================================
# 2.4 states the tree-level branching h nu : Z nu : W l = 1 : 1 : 2 without deriving it, and
# an editorial pass objected that this is the one quantitative input making the line sharp.
# It follows from the Goldstone equivalence theorem plus counting the Higgs multiplet, and
# the counting is worth printing because it is two lines.
#
# The complex Higgs doublet H = (G+, (v + h + i G0)/sqrt2) has FOUR real components. For
# M_N >> m_W the heavy Majorana N decays through its Yukawa to each with equal strength, and
# the longitudinal gauge bosons ARE the eaten Goldstones:
#
#     N -> nu h     (1 channel, neutral)
#     N -> nu G0    (1 channel, neutral)      G0  is Z_L
#     N -> l- G+    (1 channel, charged)      G+- is W+-_L
#     N -> l+ G-    (1 channel, charged)
#
# so h : Z : W = 1 : 1 : 2, and the neutral (hard-neutrino) fraction is 2/4.
chan <- c(h=1, Z=1, W=2)
cat("  channels:", paste(names(chan), chan, sep="="), "\n")
cat(sprintf("  ratio h : Z : W = %d : %d : %d\n", chan['h'], chan['Z'], chan['W']))
cat(sprintf("  neutral fraction (hard neutrinos per decay) = %d/%d = %.3f\n",
            chan['h']+chan['Z'], sum(chan), (chan['h']+chan['Z'])/sum(chan)))
cat("  which is the 0.5 per decay 3.2 uses, so the two statements are one statement.\n")

# --- how good is the equivalence theorem here? corrections are O((m_V/M_N)^2) -------------
M1 <- 491.6e6      # GeV
cat("\n  The theorem is exact only as m_V/M_N -> 0. Size of the correction here:\n\n")
cat("      boson      mass (GeV)     (m_V/M_N)^2\n")
for (nm in c("W","Z","h")) {
  mV <- switch(nm, W=80.377, Z=91.1876, h=125.25)
  cat(sprintf("   %8s %14.3f %18.3e\n", nm, mV, (mV/M1)^2))
}
cat("\n  So the ratio is 1:1:2 to fourteen decimal places at this mass, and the objection\n")
cat("  that it 'depends on conventions for the Z-nu-nu vector and axial couplings' does not\n")
cat("  arise: in the equivalence regime the decay is to the GOLDSTONES, whose couplings are\n")
cat("  the Yukawa itself, not to transverse gauge bosons whose couplings differ.\n")
cat("\n  WHAT IT DOES ASSUME, and 2.4 already says so: one Yukawa column, a single heavy\n")
cat("  Majorana N, and the decay treated at tree level. Summing over charged-lepton flavour\n")
cat("  does not change the RATIO because every flavour enters both the charged and neutral\n")
cat("  channels through the same column; it changes only which flavour the charged channel\n")
cat("  produces, which is the flavour freedom 3.2 already declines to use as a discriminant.\n")
