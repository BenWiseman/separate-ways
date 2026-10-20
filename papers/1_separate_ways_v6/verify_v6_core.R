## verify_v6_core.R
## Comparison rule, register of checks and planted errors, inputs, and small numerical tools for
## VERIFY_v6.Rmd. Base R only.

REG <- new.env()
REG$checks <- list()
REG$plants <- list()

## ---- comparison rule ----------------------------------------------------------------------

## A number as the paper prints it, and half a unit in its last printed digit.
printed <- function(s) {
  s <- trimws(s)
  mant <- sub("[eE].*$", "", s)
  ex <- if (grepl("[eE]", s)) as.integer(sub("^.*[eE]", "", s)) else 0L
  dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0L
  list(value = as.numeric(s), half_ulp = 0.5 * 10^(ex - dec))
}

## The one comparison every check goes through. A missing or infinite value fails.
judge <- function(computed, value, tol) {
  isTRUE(is.finite(computed) && is.finite(value) && is.finite(tol) && tol >= 0 &&
         abs(computed - value) <= tol)
}

## Record a check. `paper` is either the string the paper prints (tolerance: half a unit in its
## last digit, unless `tol` or `rel` is given) or an exact number with `label` for display.
chk <- function(id, sec, what, computed, paper, tol = NULL, rel = NULL, label = NULL,
                quote = NA_character_, where = "") {
  if (is.character(paper)) {
    p <- printed(paper); value <- p$value; shown <- paper
    tol_note <- "half a unit in the last digit"
    if (is.null(tol) && is.null(rel)) tol <- p$half_ulp
  } else {
    value <- as.numeric(paper)
    shown <- if (is.null(label)) formatC(value, digits = 10, format = "g") else label
    tol_note <- ""
  }
  if (!is.null(rel)) { tol <- rel * abs(value); tol_note <- sprintf("relative %g", rel) }
  if (is.null(tol)) stop("check ", id, " has no tolerance")
  if (!nzchar(tol_note)) tol_note <- "absolute"
  ok <- judge(computed, value, tol)
  REG$checks[[id]] <- list(id = id, sec = sec, what = what, computed = computed, value = value,
                           shown = shown, tol = tol, tol_note = tol_note, ok = ok,
                           quote = quote, where = where)
  invisible(ok)
}

## A planted error: a wrong input pushed through the same calculation, compared with the same
## paper value and tolerance as check `check_id`. It must fail.
plant <- function(id, check_id, what, planted) {
  c0 <- REG$checks[[check_id]]
  if (is.null(c0)) stop("plant ", id, " refers to unknown check ", check_id)
  ok <- judge(planted, c0$value, c0$tol)
  REG$plants[[id]] <- list(id = id, check_id = check_id, what = what, planted = planted,
                           verdict = if (ok) "PASS" else "FAIL", caught = !ok)
  invisible(!ok)
}

## The generic planted error for every check: the computed value moved by at least 1 per cent of
## the paper's value and by at least 2.5 times the tolerance, so it must fail.
generic_plant <- function(c0) {
  shift <- max(0.01 * abs(c0$value), 2.5 * c0$tol, 1e-300)
  planted <- c0$computed + shift
  list(planted = planted, caught = !judge(planted, c0$value, c0$tol))
}

## ---- display ------------------------------------------------------------------------------

fmt <- function(x, d = 7) {
  if (is.character(x)) return(x)
  if (!is.finite(x)) return(as.character(x))
  if (x != 0 && (abs(x) < 1e-3 || abs(x) >= 1e5)) formatC(x, digits = d - 1, format = "e")
  else formatC(x, digits = d, format = "g", flag = "#") |> sub(pattern = "\\.$", replacement = "")
}

## x as LaTeX scientific notation with d significant digits, for inline text
tex_sci <- function(x, d = 3) {
  e <- floor(log10(abs(x))); m <- signif(x / 10^e, d)
  if (abs(m) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("$%s\\times10^{%d}$", formatC(m, digits = d - 1, format = "f"), e)
}

html_escape <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

## An HTML table passed through pandoc untouched, so that quoted LaTeX stays as the paper writes it.
html_table <- function(df, size = "90%", escape = FALSE, ...) {
  tab <- knitr::kable(df, format = "html", escape = escape, row.names = FALSE, ...,
                      table.attr = sprintf('class="table table-condensed" style="font-size:%s"', size))
  knitr::asis_output(paste0("\n```{=html}\n", paste(tab, collapse = "\n"), "\n```\n"))
}

## Paper text, for the column that says whether each quoted fragment occurs in it.
find_paper <- function() {
  cand <- c("SEPARATE_WAYS_v6.md", "SEPARATE_WAYS_v6_src.md")
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(hit[1])
  md <- list.files(".", pattern = "\\.md$")
  for (f in md) {
    first <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) "")
    if (any(grepl("Separate Ways and the Upside Down", first))) return(f)
  }
  NA_character_
}
PAPER_FILE <- find_paper()
PAPER_TEXT <- if (is.na(PAPER_FILE)) NA_character_ else
  gsub("\\s+", " ", paste(readLines(PAPER_FILE, warn = FALSE), collapse = " "))
in_paper <- function(q) {
  if (is.na(q) || !nzchar(q)) return("")
  if (is.na(PAPER_TEXT)) return("no text")
  if (grepl(gsub("\\s+", " ", q), PAPER_TEXT, fixed = TRUE)) "yes" else "NO"
}

## Table of the checks of one section.
show_checks <- function(sec) {
  rows <- Filter(function(r) r$sec == sec, REG$checks)
  if (!length(rows)) return(invisible(NULL))
  df <- data.frame(
    ID = vapply(rows, `[[`, "", "id"),
    Check = vapply(rows, function(r) html_escape(r$what), ""),
    Computed = vapply(rows, function(r) fmt(r$computed), ""),
    Paper = vapply(rows, `[[`, "", "shown"),
    Tolerance = vapply(rows, function(r) paste0(fmt(r$tol, 2), " (", r$tol_note, ")"), ""),
    Verdict = vapply(rows, function(r) if (r$ok) "PASS" else "<b>FAIL</b>", ""),
    `Paper's words` = vapply(rows, function(r)
      if (is.na(r$quote)) "" else paste0("<code>", html_escape(r$quote), "</code>",
                                         if (nzchar(r$where)) paste0(" (", r$where, ")") else ""), ""),
    `In text` = vapply(rows, function(r) in_paper(r$quote), ""),
    check.names = FALSE, stringsAsFactors = FALSE)
  html_table(df)
}

## ---- inputs ------------------------------------------------------------------------------
## Every physical constant used below, with the place it comes from. The author's scripts are in
## checks/calc (working tree) or checks/calc (public release).

CONST <- list(
  M_P      = list(2.435e18,  "GeV",   "reduced Planck mass; paper Section 4 and v6_bath.py (MP)"),
  M_PL     = list(1.22089e19, "GeV",  "Planck mass from CODATA G; v6_abundance.py (M_PL)"),
  hbarc    = list(1.97327e-14, "GeV cm", "v6_abundance.py (GEV_PER_INVCM3)"),
  rho_crit = list(1.05368e-5, "GeV cm^-3", "critical density over h^2 (PDG); v6_abundance.py (RHO_CRIT_H2)"),
  s0_cm3   = list(2891.2,     "cm^-3", "entropy density today for T0 = 2.7255 K; paper Section 3.2 and v6_abundance.py (S0_CM3)"),
  s0_GeV3  = list(2.2215e-38, "GeV^3", "the same in GeV^3 as the paper prints it; v6_bath.py (RHO_OVER_S)"),
  rho_DM   = list(9.7e-48,    "GeV^4", "Boyle, Finn and Turok's rho_DM,0; paper Section 3.2, v6_abundance.py, v6_bath.py"),
  gstar    = list(106.75,     "",      "relativistic degrees of freedom; paper Section 3.2, v6_abundance.py (GSTAR)"),
  Omega_c  = list(0.1200,     "",      "Planck 2018 Omega_c h^2; paper Section 3.2 and v6_abundance.py (OMEGA_C_H2)"),
  Omega_err= list(0.0012,     "",      "its error; v6_abundance.py (OMEGA_ERR)"),
  s0_BFT   = list(2.3e-38,    "GeV^3", "Boyle, Finn and Turok's entropy density; paper Section 3.2, v6_abundance.py"),
  T0_K     = list(2.7255,     "K",     "CMB temperature (Fixsen); paper Section 3.2"),
  kB       = list(8.617333262e-5, "eV/K", "Boltzmann constant (CODATA 2018); used only to rebuild s0 from T0, not in the scripts"),
  f2_BTV   = list(7e-5,       "",      "Boyle, Turok and Vaibhav's bound on the Weyl coupling; paper Section 4, v6_grav_rate.py")
)
K <- function(name) CONST[[name]][[1]]

const_table <- function() {
  df <- data.frame(Name = names(CONST),
                   Value = vapply(CONST, function(z) fmt(z[[1]], 8), ""),
                   Unit = vapply(CONST, `[[`, "", 2),
                   Source = vapply(CONST, `[[`, "", 3), stringsAsFactors = FALSE)
  df$Source <- html_escape(df$Source)
  html_table(df)
}

## ---- numerical tools ------------------------------------------------------------------------

## Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
gauss_legendre <- function(n) {
  k <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- J[cbind(k + 1, k)] <- k / sqrt(4 * k^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
}
## the same rule on [a, b], split into equal panels
gl_nodes <- function(a, b, n, panels = 1) {
  g <- gauss_legendre(n); e <- seq(a, b, length.out = panels + 1)
  list(x = unlist(lapply(seq_len(panels), function(i) (e[i + 1] - e[i]) / 2 * g$x + (e[i + 1] + e[i]) / 2)),
       w = unlist(lapply(seq_len(panels), function(i) (e[i + 1] - e[i]) / 2 * g$w)))
}

## Riemann zeta for real s > 1: direct sum plus the Euler-Maclaurin tail.
zeta_num <- function(s, N = 30) {
  n <- seq_len(N - 1)
  sum(n^-s) + N^(1 - s) / (s - 1) + 0.5 * N^-s + s * N^(-s - 1) / 12 -
    s * (s + 1) * (s + 2) * N^(-s - 3) / 720 +
    s * (s + 1) * (s + 2) * (s + 3) * (s + 4) * N^(-s - 5) / 30240
}

## Best rational approximation p/q of x with q <= qmax (continued fractions).
as_fraction <- function(x, qmax = 1e10, tol = 1e-12) {
  sgn <- sign(x); x <- abs(x)
  h0 <- 0; h1 <- 1; k0 <- 1; k1 <- 0; y <- x
  repeat {
    a <- floor(y)
    h2 <- a * h1 + h0; k2 <- a * k1 + k0
    if (k2 > qmax) break
    h0 <- h1; h1 <- h2; k0 <- k1; k1 <- k2
    if (abs(x - h1 / k1) <= tol * x || y == a) break
    y <- 1 / (y - a)
  }
  paste0(if (sgn < 0) "-" else "", format(h1, scientific = FALSE),
         if (k1 == 1) "" else paste0("/", format(k1, scientific = FALSE)))
}

## mu = (4 pi G / 3)^(-1/2) = sqrt(3 / 4 pi) M_Pl, Boyle, Finn and Turok's mu-hat
MU <- sqrt(3 / (4 * pi)) * K("M_PL")
