#!/usr/bin/env Rscript

# Render every figure both papers use.
#
# This file was stale: it listed the EARLIER paper's figures and none of the nine
# PAPER2_v3.md actually places, so a reader following the repository's own instructions
# regenerated the wrong set and the current figures silently kept whatever was on disk.
# Check 22 in separate_ways/consistency_pass.py now fails if a figure the manuscript
# includes has no generator listed here, so this cannot drift again.
#
# The two groups resolve paths differently. The earlier scripts locate their output from
# their own --file= path and do not care where they are run from; the current ones write to
# "pub/paper2/figs" relative to the repository root. Setting the working directory to the
# root satisfies both.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
script_dir  <- dirname(normalizePath(script_file))
repo_root   <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

# PAPER2_v3.md, in the order the manuscript places them
current <- c(
  "graphical_abstract_v3.R",
  "fig_fold.R",        # Figure 1
  "fig_halfhalf.R",    # Figure 2
  "fig_crossing.R",    # Figure 3
  "fig_floor.R",       # Figure 4
  "fig_band.R",        # Figure 5
  "fig_line.R",        # Figure 6
  "fig_budget.R",      # Figure 8
  "fig_power.R",       # Figure 9
  "fig_twosided.R"     # Figure 10
)

# Figure 7 (fig5_data) is NOT here on purpose. It is built by the de-AI'd Python generator
# (pub/paper2/compute_figs.py and calc/generate_fig5_v3.py); fig5_data.R in this directory is
# a parameter block, not a renderer. Running an R "fig5" here would overwrite a good figure
# with a different build, which is the failure SUBMISSION_CHECKLIST.md already warns about.

# The earlier paper (PAPER2_v2.md / SUPPLEMENT_v3.md)
earlier <- c(
  "fig1_dictionary.R",
  "fig2_janus.R",
  "fig3_curvature.R",
  "fig4_tilt.R"
)

old_wd <- getwd()
setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

missing <- character(0)
for (figure in c(current, earlier)) {
  path <- file.path(script_dir, figure)
  if (!file.exists(path)) {            # fig4_tilt.R was untracked while its three siblings
    missing <- c(missing, figure)      # were committed; skip rather than die on someone
    next                               # else's oversight, but say so loudly.
  }
  message("Rendering ", figure)
  status <- system2("Rscript", c(path, args))
  if (status != 0L) stop("Render failed: ", figure)
}
if (length(missing))
  warning("generator(s) not on disk, nothing rendered for them: ",
          paste(missing, collapse = ", "), call. = FALSE)
message("Rendered ", length(current), " current + ", length(earlier), " earlier figures.")
