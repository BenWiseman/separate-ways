# helpers.R -- shared reporting utility for the Paper 2 reproducibility companion.
#
# This file contains NO physics: it only formats a side-by-side comparison of an
# "expected" number (quoted verbatim from a named source file) against a "reproduced"
# number (computed in the calling script from the same formula) and prints a
# machine-parseable RESULT line that run_all.sh greps to build the pass/fail table.
#
# Modes:
#   "rel" (default) -- PASS if relative error <= tol (tol given as a fraction, e.g. 0.005 = 0.5%)
#   "oom"            -- PASS if reproduced is within a multiplicative factor `tol` of expected;
#                        use ONLY for numbers the source itself labels ORDER OF MAGNITUDE.
#   "abs"            -- PASS if |reproduced - expected| <= tol (absolute); use for quantities
#                        that can be exactly zero (e.g. components of a phase).

fmt_num <- function(x, digits = 7) {
  formatC(x, format = "g", digits = digits)
}

report <- function(label, expected, reproduced, tol = 0.005,
                    mode = c("rel", "oom", "abs"), note = "") {
  mode <- match.arg(mode)
  relpct <- 100 * abs(reproduced - expected) / abs(expected)  # always computed, never hidden
  if (mode == "rel") {
    pass <- (relpct / 100) <= tol
    tolstr <- sprintf("rel<=%.3g%%", tol * 100)
    errstr <- sprintf("%.4g%%", relpct)
  } else if (mode == "oom") {
    ratio <- reproduced / expected
    factor <- max(ratio, 1 / ratio)
    pass <- factor <= tol
    tolstr <- sprintf("factor<=%.2g", tol)
    errstr <- sprintf("factor %.3g (%.4g%% rel)", factor, relpct)
  } else { # abs
    abserr <- abs(reproduced - expected)
    pass <- abserr <= tol
    tolstr <- sprintf("abs<=%.3g", tol)
    errstr <- sprintf("%.3g (%.4g%% rel)", abserr, relpct)
  }
  status <- if (pass) "PASS" else "FAIL"
  # FLAT, UN-HIDEABLE DISCLOSURE: the default rule is "tolerance 0.5% unless
  # the source gives more digits". Whenever this call uses anything OTHER than a
  # plain 0.5%-relative test (i.e. a loosened tolerance, used only where documented
  # in the calling script -- coarser source precision, or an explicit ORDER OF
  # MAGNITUDE label), this line ALSO prints whether the number clears the strict
  # 0.5% default, so a loosened tolerance can never quietly hide a number that
  # would fail it.
  strict_flag <- ""
  is_plain_default <- mode == "rel" && isTRUE(all.equal(tol, 0.005))
  if (!is_plain_default) {
    strict_pass <- relpct <= 0.5
    strict_flag <- sprintf("  [[vs. plain 0.5%% default: %s, %.4g%% relative]]",
                            if (strict_pass) "PASS" else "FAIL", relpct)
  }
  cat(sprintf("  %-58s expected=%-16s reproduced=%-16s err=%-12s tol=%-12s [%s]%s\n",
              label, fmt_num(expected), fmt_num(reproduced), errstr, tolstr, status, strict_flag))
  if (nzchar(note)) cat(sprintf("      note: %s\n", note))
  # RESULT lines are pipe-delimited for run_all.sh to parse; defensively strip any
  # literal "|" out of free-text fields (label/note can never affect the numbers)
  # so a stray character never desyncs the columns.
  safe_label <- gsub("|", "/", label, fixed = TRUE)
  strict_col <- if (is_plain_default) status else if (relpct <= 0.5) "PASS" else "FAIL"
  cat(sprintf("RESULT|%s|%s|%s|%s|%s|%s|%s\n",
              safe_label, fmt_num(expected), fmt_num(reproduced), errstr, tolstr, status, strict_col))
  invisible(pass)
}

section <- function(title) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}
