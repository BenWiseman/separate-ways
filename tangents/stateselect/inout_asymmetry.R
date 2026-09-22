# How much does the mass ceiling degrade if the in/out occupations are only APPROXIMATELY equal?
#
# The operator inequality Q >= n_* is unconditional: it is a statement about a matrix. What needs
# the contact condition is the step from <Q> to <N_+>, the production integral. A referee will ask
# what happens when that condition holds only to some tolerance, and the paper does not say.
#
# Write Q = (N_+ + N_-)/2 and let the asymmetry be delta = <N_+> - <N_->. Then
#     <N_+> = <Q> + delta/2  >=  n_* + delta/2,
# so the bound degrades linearly and by a known amount. It is weakened only when delta < 0, i.e.
# when the out region holds FEWER quanta than the in region.
nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
W  <- function(f) integrate(function(x) x^2*f(x), 0, 30, subdivisions=8000)$value
I0 <- W(nmin)/pi^2
M0 <- 491.6; sig <- 2.0

# model the asymmetry as a fractional shift eps of the local occupation: delta(x) = eps * n_*(x)
I_of <- function(eps) W(function(x) nmin(x)*(1+eps/2))/pi^2
cat(sprintf("  I_0 = %.10f, endpoint %.1f PeV\n\n", I0, M0))
cat("    fractional in/out asymmetry eps    I/I_0     mass (PeV)     shift\n")
for (eps in c(-0.10,-0.05,-0.02,-0.01,0,0.01,0.05)) {
  r <- I_of(eps)/I0
  cat(sprintf("   %32.2f %9.5f %12.2f %10.2f\n", eps, r, M0*r^(-2/5), M0*r^(-2/5)-M0))
}
e_sig <- uniroot(function(e) M0*(I_of(e)/I0)^(-2/5) - M0 - sig, c(-0.2,-1e-6))$root
cat(sprintf("\n  The asymmetry that moves the ceiling by its full quoted width is eps = %.4f,\n", e_sig))
cat(sprintf("  an out region holding %.2f per cent fewer quanta than the in region.\n", -100*e_sig))

cat("\n  WHY THIS IS A ROBUSTNESS RESULT AND NOT A HOLE. Three things.\n")
cat("  1. The degradation is LINEAR and its coefficient is exactly 1/2, known in advance rather\n")
cat("     than estimated: <N_+> >= n_* + delta/2. No modelling enters.\n")
cat("  2. It is one-sided, and the safe direction is the one a bang produces. An out-region\n")
cat("     EXCESS lowers the ceiling, so the quoted 491.6 stays valid and becomes conservative.\n")
cat("     Only a DEFICIT raises the true ceiling above the quoted one, and a deficit means the\n")
cat("     out region holding fewer quanta than the in region, which gravitational production at\n")
cat("     a bang does not do: it creates pairs, it does not remove them.\n")
cat("  3. The tolerance is loose. Two per cent of asymmetry is needed before the ceiling moves\n")
cat("     by its full width, and the contact condition gives equality exactly, so the question\n")
cat("     is how badly the two-level reduction fails rather than whether the condition holds.\n")
cat("\n  FLATLY, what this does NOT do: it does not bound the asymmetry itself. That would need\n")
cat("  the interacting theory, or a calculation of the reduction's error, and neither exists\n")
cat("  here. What it supplies is the exchange rate, so that anyone who later bounds the\n")
cat("  asymmetry can read off what it costs.\n")

cat("\n  ==========================================================================\n")
cat("  THE EXCHANGE RATE AT LARGER ASYMMETRIES, since a referee asked for 5 and 10 per cent.\n")
cat("  ==========================================================================\n")
# <N_+> >= n_* + delta/2 with delta the fractional deficit, so the floor rises by delta/2
# relative to n_*, I falls in proportion at leading order, and M ~ I^(-2/5) rises.
M0 <- 491.6
cat("      deficit delta    floor n_* -> n_*(1 - delta/2 / n_*) ...  ceiling (PeV)   move   widths\n")
for (d in c(0.0202, 0.05, 0.10, 0.20)) {
  # a DEFICIT of delta raises the true ceiling: I is reduced by the fraction delta/2
  Ifac <- 1 - d/2
  M <- M0 * Ifac^(-2/5)
  cat(sprintf("   %12.1f%% %44.1f %8.1f %8.1f\n", 100*d, M, M-M0, (M-M0)/2.0))
}
cat("\n  Linear, as the identity says, and mild: ten per cent of asymmetry in the UNSAFE\n")
cat("  direction raises the ceiling by about 10 PeV, five widths. Twenty per cent raises it\n")
cat("  by 21. The ceiling degrades gracefully rather than failing, and it degrades only in\n")
cat("  the direction a bang does not produce.\n")
