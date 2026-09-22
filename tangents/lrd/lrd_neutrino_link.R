# Is the JWST / neutrino-line link real at OBSERVED little-red-dot densities, or only formal?
#
# The paper says a JWST result and a neutrino-telescope result are not independent tests of this
# framework, because a compact-object fraction f of the dark budget lowers the mass as
# M_1 = M_1(0) * (1-f)^(2/5), and the two-body line follows at half of it. True, but never checked
# for reach. Compute the f that today's LRD seed population actually implies, and the f needed to
# move the line by one quoted error bar.
M0 <- 491.6; Enu0 <- 245.8; sigE <- 1.0          # PeV
OmDM_h2 <- 0.1200; h <- 0.674
rho_crit <- 2.775e11 * h^2                       # Msun / Mpc^3
rho_DM   <- OmDM_h2/h^2 * rho_crit               # Msun / Mpc^3
cat(sprintf("  dark-matter density: %.3e Msun/Mpc^3\n\n", rho_DM))

f_of <- function(n, Ms) n*Ms/rho_DM              # comoving number density n per Mpc^3
line <- function(f) Enu0*(1-f)^(2/5)

cat("   seed mass   number density      f (fraction of budget)   neutrino line (PeV)   shift\n")
cases <- list(c(1e5,1e-4), c(1e5,1e-3), c(1e6,1e-4), c(1e7,1e-4), c(1e5,1e-1))
for (cc in cases) {
  f <- f_of(cc[2], cc[1])
  cat(sprintf("  %9.0e %17.0e %24.3e %20.4f %9.2e\n", cc[1], cc[2], f, line(f), Enu0-line(f)))
}

cat("\n  What f would move the line by its own quoted width of 1.0 PeV?\n")
f_need <- uniroot(function(f) Enu0 - line(f) - sigE, c(1e-6, 0.5))$root
cat(sprintf("     f = %.5f, i.e. %.2f per cent of the dark budget\n", f_need, 100*f_need))
base <- f_of(1e-4, 1e5)
cat(sprintf("     the fiducial LRD seed budget is f = %.3e\n", base))
cat(sprintf("     so the aggregate seed mass would have to be %.1e times larger\n", f_need/base))
cat(sprintf("     e.g. %.0e seeds of 1e5 Msun per Mpc^3, against the ~1e-4 observed\n",
    f_need*rho_DM/1e5))

cat("\n  FLATLY: the link is formal, not observational. At the seed budget the little red dots\n")
cat("  actually imply, the neutrino line moves by 2.95e-8 PeV, seven and a half orders below\n")
cat("  its own quoted width of 1.0 PeV. A JWST measurement cannot move the\n")
cat("  neutrino prediction and a neutrino measurement cannot constrain the seeds. The paper's\n")
cat("  'not independent tests' is true as algebra and empty as a test.\n")
cat("\n  NEXT ROUTE, and it inverts the question. The interesting number is not the shift but the\n")
cat("  CEILING: the closed budget caps the aggregate primordial seed mass at f = 1, which is\n")
cat(sprintf("     n*M_s <= %.2e Msun/Mpc^3.\n", rho_DM))
cat("  Observed LRD seeds sit nine orders below that, so the fold is in no tension with them and\n")
cat("  would only be strained if the seed population were found to be nine orders denser. That\n")
cat("  is a statement about what CANNOT happen, which is the shape a falsifier takes, and it is\n")
cat("  worth stating because most dark-sector models have no such ceiling at all.\n")
