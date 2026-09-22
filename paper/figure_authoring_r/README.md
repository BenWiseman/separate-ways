# Paper 2 figure authoring in R

This folder gives each V3 figure a small R entry point. The scripts edit the
current SVG masters deterministically, then render matching PDF and PNG files
with local Chrome. They do not change the release assets unless `--publish` is
present.

## Quick use

From the project root:

```bash
Rscript pub/paper2/figure_authoring_r/graphical_abstract_v3.R
```

The preview appears in:

```text
tmp/figure_authoring_r/graphical_abstract_v3/
```

Open the R file, change the named `text_edit()` values, and rerun it. The most
useful controls are:

- `x`, `y`: baseline position in the SVG coordinate system
- `size`: font size
- `value`: replacement text
- `anchor`: `start`, `middle`, or `end`
- `subscript_size` and `subscript_shift`: equation alignment
- `max_width`: force a label into a fixed width by adjusting glyph spacing;
  use this only after moving, wrapping, or shortening the label has failed

To inspect every text node in a figure:

```bash
Rscript pub/paper2/figure_authoring_r/graphical_abstract_v3.R --list-text
```

To render all six figures into the preview tree:

```bash
Rscript pub/paper2/figure_authoring_r/render_all.R
```

To write one approved revision into `pub/paper2/`:

```bash
Rscript pub/paper2/figure_authoring_r/graphical_abstract_v3.R --publish
```

The publish path first copies the current SVG, PDF and PNG into a timestamped
folder under `tmp/figure_authoring_r/backups/`. Publishing a figure changes the
release package, so rebuild the manuscript PDFs, the arXiv archive and
`SHA256SUMS.txt` afterwards.

## Files

- `graphical_abstract_v3.R`: hero image, equations, branch and horizon labels
- `fig1_dictionary.R`: crease dictionary
- `fig2_janus.R`: regular and irregular Bel-Robinson curves
- `fig3_curvature.R`: discrete curvature constraints
- `fig4_tilt.R`: primordial-tilt comparison
- `fig5_data.R`: neutrino mass, geometry-only profiles and exact-Lambda panel
- `_helpers.R`: SVG selection, text adjustment, backup and rendering
- `render_all.R`: one command for all previews

The original numerical generators remain authoritative for data and curve
values. These R files alter presentation only. They preserve the existing
marks, labels, values, ordering and caveats unless you deliberately change a
text value.
