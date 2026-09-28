# Is the entropy coefficient the same at every horizon?
#
# Jacobson's second input is "an entropy proportional to the horizon's area with a UNIVERSAL
# coefficient eta". Section 3.6 supplies the proportionality, from the entanglement of the
# thermofield double the fold forces, and universality_local_fold.R supplies the horizons. What
# neither supplies is the universality of the coefficient, and that is the half that decides
# whether the field equation is Einstein's. If eta depended on the surface gravity then
# 2 pi T_kk = eta R_kk would have a different constant at every point and the resulting theory
# would not be general relativity at all.
#
# WHY IT COULD HAVE GONE EITHER WAY. The entropy is sum over modes of s(n) with n the Bose factor
# at beta omega, and beta = 2 pi / kappa differs from horizon to horizon. So kappa is visibly
# present in the occupation. The claim to test is that it cancels against the mode density.
#
# WHERE THE CANCELLATION COMES FROM. Write the wedge in boost coordinates,
#     ds^2 = -rho^2 deta^2 + drho^2 + dx_perp^2,
# with rho the proper distance to the horizon and eta the dimensionless boost parameter. No kappa
# appears anywhere in that metric: kappa enters only when eta is identified with kappa times the
# proper time of one particular accelerated observer. Bisognano-Wichmann then makes the modular
# temperature 2 pi in eta, for every wedge. So with Omega the boost frequency conjugate to eta,
# beta omega = 2 pi Omega exactly, whatever kappa is. Everything below checks that the entropy
# per unit area really does come out free of kappa, and what it is instead a function of.

s_of_n <- function(n) ifelse(n <= 0, 0, (1 + n)*log(1 + n) - n*log(n))
n_bose <- function(x) 1/(expm1(x))                       # x = beta omega

cat("=== 1. the occupation in boost variables carries no kappa ===\n")
cat("   omega is the proper frequency of an observer with surface gravity kappa, Omega = omega/kappa\n")
cat("   the boost frequency, and beta = 2 pi / kappa, so beta omega = 2 pi Omega identically.\n\n")
cat("      kappa     Omega     omega = kappa Omega    beta omega    n\n")
for (kap in c(0.025, 1, 7.3, 250)) {
  Om <- 0.4; om <- kap*Om; bo <- (2*pi/kap)*om
  cat(sprintf("   %9.2f   %6.2f   %18.4f   %11.6f   %.8f\n", kap, Om, om, bo, n_bose(bo)))
}
cat("   The occupation is the same number at every surface gravity, which is the cancellation\n")
cat("   stated rather than the one computed. Sections 2 and 3 compute it.\n")

cat("\n=== 2. the mode density, by WKB, in the OBSERVER'S proper variables ===\n")
cat("   The cancellation is only a result if kappa is carried through rather than left out, so\n")
cat("   everything below is written in the proper frequency omega of an observer with surface\n")
cat("   gravity kappa, in the wedge metric ds^2 = -kappa^2 rho^2 dtau^2 + drho^2 + dx_perp^2.\n")
cat("   The radial turning point is at kappa rho k = omega, and the WKB count below omega is\n")
cat("      N(<omega) = (1/pi) int_eps^{omega/(kappa k)} sqrt(omega^2/(kappa^2 rho^2) - k^2) drho,\n")
cat("   which carries kappa in three separate places.\n\n")
Ncount <- function(om, k, eps, kap) {
  top <- om/(kap*k)
  if (top <= eps) return(0)
  f <- function(r) sqrt(pmax(om^2/(kap^2*r^2) - k^2, 0))
  integrate(f, eps, top, rel.tol = 1e-10, subdivisions = 2000)$value / pi
}
dN <- function(om, k, eps, kap, h = 1e-6) {
  (Ncount(om + h*kap, k, eps, kap) - Ncount(max(om - h*kap, 0), k, eps, kap)) / (2*h*kap)
}
cat("      kappa    omega     k      N(<omega)     and N at the same Omega = omega/kappa\n")
for (kap in c(0.5, 1, 4)) {
  om <- 1.2*kap                                   # hold Omega = 1.2 fixed while kappa moves
  cat(sprintf("   %8.2f  %7.3f  %5.2f   %11.6f   %s\n", kap, om, 0.6,
              Ncount(om, 0.6, 1e-4, kap), "Omega = 1.20"))
}
cat("   The count depends on omega and kappa only through Omega = omega/kappa, which is the\n")
cat("   first half of the cancellation and is measured here, not assumed.\n")

cat("\n=== 3. the entropy per unit area, at four surface gravities, kappa carried throughout ===\n")
cat("   S/A = int d^2k/(2 pi)^2 int domega (dN/domega) s(n(beta omega)), beta = 2 pi / kappa.\n")
cat("   The occupation carries kappa one way and the density the other. If they fail to cancel\n")
cat("   the coefficient is not universal and the field equation is not Einstein's.\n\n")
SA <- function(eps, kap, kfac = 1) {
  kmax <- 1/eps
  beta <- 2*pi/kap
  fk <- function(k) {
    fO <- function(om) sapply(om, function(o) dN(o, k, eps, kap) * s_of_n(n_bose(beta*o*kfac)))
    integrate(fO, 1e-6*kap, 6*kap, rel.tol = 1e-8, subdivisions = 800)$value
  }
  ks <- seq(0.02*kmax, 0.98*kmax, length.out = 40)
  inner <- sapply(ks, function(k) 2*pi*k*fk(k))
  sum((inner[-1] + inner[-length(inner)])/2 * diff(ks)) / (2*pi)^2
}
eps <- 0.02
base <- SA(eps, 1)
cat("      kappa        beta          S/A            relative to kappa = 1\n")
vals <- c()
for (kap in c(0.025, 1, 7.3, 250)) {
  v <- SA(eps, kap); vals <- c(vals, v)
  cat(sprintf("   %10.2f   %11.4f   %12.6f   %20.2e\n", kap, 2*pi/kap, v, abs(v/base - 1)))
}
cat("\n   Kappa runs from 0.025 to 250, four decades exactly, and S/A does not move. The\n")
cat("   coefficient is universal,\n")
cat("   which is the half of Jacobson's assumption that decides whether the field equation is\n")
cat("   Einstein's. The cancellation is between the Jacobian domega = kappa dOmega and the\n")
cat("   1/kappa the density carries, and neither the wedge metric in boost coordinates nor the\n")
cat("   modular flow knows kappa at all.\n")
stopifnot(max(abs(vals/base - 1)) < 1e-6)

cat("\n=== 4. what it IS a function of, which is the other half of the story ===\n")
cat("   S/A is not a pure number: it runs with the cutoff, and that is where the one dimensionful\n")
cat("   constant of the construction lives.\n\n")
cat("        eps        S/A          S/A x eps^2\n")
es <- c(0.04, 0.03, 0.02, 0.015)
pr <- c()
for (e in es) { v <- SA(e, 1); pr <- c(pr, v*e^2)
  cat(sprintf("   %9.3f   %12.5f   %14.6f\n", e, v, v*e^2)) }
ex <- unname(coef(lm(log(sapply(es, function(e) SA(e,1))) ~ log(es)))[2])
cat(sprintf("\n   fitted exponent %.4f against -2, so S/A goes as 1/eps^2\n", ex))
cat("   An area law with a coefficient of mass dimension two, diverging with the cutoff. That is\n")
cat("   the same object the dimensional argument says no symmetry could ever supply: eta = 1/4G\n")
cat("   carries dimension two and every structural input of the construction carries zero. The\n")
cat("   two statements section 3.6 makes separately are one statement. What the fold fixes is\n")
cat("   that the coefficient is the same everywhere; what it cannot fix is its value.\n")
stopifnot(abs(ex + 2) < 0.1)

cat("\n=== 5. the plant: break the universal modular temperature and kappa must reappear ===\n")
cat("   Bisognano-Wichmann is what makes beta omega = 2 pi Omega at every wedge. Suppose instead\n")
cat("   the modular temperature ran with the surface gravity, beta omega = 2 pi Omega kappa/kappa_0.\n")
cat("   Then the occupation depends on kappa and the entropy must follow.\n\n")
cat("      kappa       S/A with a running modular temperature     relative to kappa = 1\n")
b0 <- SA(eps, 1, kfac = 1)
spread <- 0
for (kap in c(0.3, 1, 3)) {
  v <- SA(eps, 1, kfac = kap)
  spread <- max(spread, abs(v/b0 - 1))
  cat(sprintf("   %9.2f   %30.6f   %18.4f\n", kap, v, abs(v/b0 - 1)))
}
cat(sprintf("\n   spread %.3f, so the check is reading the temperature and the agreement in section 3\n", spread))
cat("   is a cancellation rather than an insensitivity of the quadrature.\n")
stopifnot(spread > 0.1)

cat("\n=== 5b. what the universality actually rests on, found by trying to break it ===\n")
cat("   The cutoff above is a fixed PROPER distance eps from the horizon, and that choice was\n")
cat("   made without comment. It is doing work. A cutoff is a property of the theory, a length\n")
cat("   like the Planck length, and a proper distance is the same length at every horizon; but\n")
cat("   nothing in the arithmetic forces that reading, so it is worth seeing what a horizon-\n")
cat("   dependent cutoff would do before a reader asks.\n\n")
cat("      prescription                        kappa = 0.1      kappa = 1      kappa = 7.3     spread\n")
row <- function(lab, f) {
  v <- sapply(c(0.1, 1, 7.3), function(k) SA(f(k), k))
  cat(sprintf("      %-34s %11.5f %14.5f %14.5f  %9.2e\n", lab, v[1], v[2], v[3],
              max(abs(v/v[2] - 1))))
  max(abs(v/v[2] - 1))
}
s_fixed <- row("eps fixed, a proper length",        function(k) 0.02)
s_kappa <- row("eps proportional to 1/kappa",        function(k) 0.02/k)
s_sqrt  <- row("eps proportional to 1/sqrt(kappa)",  function(k) 0.02/sqrt(k))
cat("\n   So the cancellation is not unconditional. It holds when the cutoff is a fixed proper\n")
cat("   length and fails when the cutoff is allowed to track the surface gravity, which is what\n")
cat("   one would expect from S/A going as 1/eps^2. Stating it the other way round is the\n")
cat("   useful form: eta is the same at every horizon provided the ultraviolet cutoff is one\n")
cat("   length rather than one per horizon. That is a mild condition and it is not nothing, and\n")
cat("   the paper should carry it rather than leave a referee to find it.\n")
stopifnot(s_fixed < 1e-6, s_kappa > 1, s_sqrt > 0.1)

cat("\n=== 5c. and it is not the WKB prescription either ===\n")
cat("   The mode count above came from one turning-point formula. If universality were an\n")
cat("   artefact of that formula it should break when the count is changed. Weighting the\n")
cat("   density by an arbitrary smooth function of the transverse momentum, which is what a\n")
cat("   different regularisation of the transverse integral looks like, must leave the kappa\n")
cat("   independence alone while moving the value.\n\n")
SAw <- function(eps, kap, w) {
  kmax <- 1/eps; beta <- 2*pi/kap
  fk <- function(k) {
    fO <- function(om) sapply(om, function(o) dN(o, k, eps, kap) * s_of_n(n_bose(beta*o)))
    w(k) * integrate(fO, 1e-6*kap, 6*kap, rel.tol = 1e-8, subdivisions = 800)$value
  }
  ks <- seq(0.02*kmax, 0.98*kmax, length.out = 40)
  inner <- sapply(ks, function(k) 2*pi*k*fk(k))
  sum((inner[-1] + inner[-length(inner)])/2 * diff(ks)) / (2*pi)^2
}
cat("      weight w(k)                kappa = 0.1      kappa = 1      kappa = 7.3     spread\n")
for (nm in list(list("1 (as above)", function(k) 1),
                list("exp(-k eps)",  function(k) exp(-k*0.02)),
                list("1/(1 + k^2 eps^2)", function(k) 1/(1 + (k*0.02)^2)))) {
  v <- sapply(c(0.1, 1, 7.3), function(kk) SAw(0.02, kk, nm[[2]]))
  cat(sprintf("      %-24s %11.5f %14.5f %14.5f  %9.2e\n", nm[[1]], v[1], v[2], v[3],
              max(abs(v/v[2] - 1))))
  stopifnot(max(abs(v/v[2] - 1)) < 1e-6)
}
cat("\n   The value moves with the weight and the kappa independence does not, so the\n")
cat("   cancellation belongs to the boost structure and not to one way of counting modes.\n")

cat("\n=== 6. the verdict ===\n")
cat("   Jacobson's second input has two halves and the fold supplies both. Proportionality to\n")
cat("   area comes from the transverse mode count, which section 3.6 already had. Universality\n")
cat("   of the coefficient comes from the wedge's modular temperature being 2 pi in boost time\n")
cat("   at every horizon, which is Bisognano-Wichmann and is the same theorem the fold's own map\n")
cat("   is read from. The surface gravity cancels between the occupation and the mode density\n")
cat("   because neither the wedge metric in boost coordinates nor the modular flow knows it.\n")
cat("   What survives is one divergent number, which is the dimensionful constant.\n")
cat("   One condition attaches and section 5b found it: the ultraviolet cutoff has to be a\n")
cat("   single proper length rather than one per horizon, because S/A goes as 1/eps^2 and a\n")
cat("   cutoff tracking the surface gravity moves it by a factor of fifty. That is what a\n")
cat("   cutoff is, a length belonging to the theory, so the condition is mild; it is stated\n")
cat("   because it is load-bearing and not because it is doubtful.\n")
