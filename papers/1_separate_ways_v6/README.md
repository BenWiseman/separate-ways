# Separate Ways and the Upside Down, version 6

Version 6 of paper 1, under the subtitle "corrected graviton-exchange rates and the dark matter of a CPT-symmetric universe". Version 5.6 stays in `papers/1_separate_ways/`, unchanged. Version 6 is a large revision
and sits beside it rather than replacing it.

Version 6 keeps only the dark matter. It calculates the rate at which graviton exchange in the hot
plasma after the bang makes the dark-matter particle $N_1$, $1.20\times10^{-2}\,T^8/M_P^4$, finds that
a line of published rates counts the Standard Model's 45 Weyl fermions as 45 Dirac fermions, gives
corrected coefficients for scalar, fermion and vector dark matter, and sets each published rate
beside the corrected one, equation by equation. Without inflation, every CPT-invariant state makes at
least as many particles as Boyle, Finn and Turok's, so their mass, 491.6 PeV, is a ceiling, and in
Einstein gravity a plasma that ever reached $8.0\times10^{13}$ GeV makes all the dark matter.

The case with inflation after the bang is a separate paper in this directory, "Upside Down Inflation": there the choice of state drops out and the reheating history sets the mass.

## Files

- `Separate_Ways_and_the_Upside_Down_v6.pdf`: the paper, 18 pages, built from `SEPARATE_WAYS_v6.md`.
- `SEPARATE_WAYS_v6_src.md`: the source as written, citing by key, as in `[@Garny2016]`.
- `SEPARATE_WAYS_v6.md`: the same with the citations numbered and the reference list written, made
  from the `_src` file by `v6_refs.py`. Edit the `_src` file and rebuild this one.
- `fig_v6_occupations.pdf`, `fig_v6_rates.pdf` and their `.png` copies: the two figures, drawn by
  `v6_figures.py`. The markdown points to the `.png` files and the PDF uses the `.pdf` ones.
  `fig_v6_occupations_vac.npz` holds the occupation numbers of the massless vacuum in Fig. 2, so that it can be redrawn without solving the modes again. Delete it to solve them again.
- `graphical_abstract_NPB_v6_toc.pdf` and `.png`: the graphical abstract for Nuclear Physics B, sized to read
  at the 200 pixels high that contents pages show; `graphical_abstract_NPB_v6.pdf` and `.png` are a fuller version.
  `v6_graphical_abstract.py` writes `graphical_abstract_NPB_v6.html` from the printed output of the
  scripts in `checks/calc/`, so it adds no number of its own, and renders it with headless Chrome.
- `v6_numbers_check.py`: the check described below.
- `VERIFY_v6.Rmd` and `VERIFY_v6.html`: an R Markdown report that recomputes the rates, the abundance
  and the ceiling in R, independently of the Python, and sets each beside the paper's value. Render
  it with `Rscript -e 'rmarkdown::render("VERIFY_v6.Rmd")'` from this directory.
- `Upside_Down_Inflation_v1.pdf`: the companion paper on the case with inflation, from
  `INFLATON_v1_src.md` (references by `inflaton_refs.py`, figures by `fig_inflaton.py`, Table 1 by
  `mass_table.py`). `inflaton_check.py` checks its numbers the same way, and Table 1 cell by cell;
  three of them quote an unpublished test that is not in this release, and the check names them as
  skipped. `fig_infl_states.npz` holds the occupations behind its first figure.
- `v6_style.py`: the shared figure style (Computer Modern type, the purple-green palette).

The scripts behind the numbers are the `v6_*.py` files in `checks/calc/`.

To rebuild the PDF from this directory:

    python3 v6_refs.py
    sed 's/\.png)/.pdf)/g' SEPARATE_WAYS_v6.md > build.md
    pandoc build.md -o Separate_Ways_and_the_Upside_Down_v6.pdf --pdf-engine=xelatex \
        -V geometry:margin=1in -V fontsize=10pt

The released PDF also used a LaTeX preamble that is not in this repository, so a rebuild differs
from it slightly in layout and type.

## Checking the numbers

From the repository root, with Python 3.12, numpy, scipy, mpmath and matplotlib:

    python3 papers/1_separate_ways_v6/v6_numbers_check.py              # every claim
    python3 papers/1_separate_ways_v6/v6_numbers_check.py --selftest   # planted errors must fail
    python3 papers/1_separate_ways_v6/v6_refs.py --selftest            # so must planted citations
    python3 papers/1_separate_ways_v6/inflaton_check.py                # the inflaton paper, and its Table 1

The check takes each number the paper quotes, runs the script that makes it, and fails if the
number is not in that script's output.

Nine of the scripts take from minutes to hours, so the check reads their saved output from
`checks/calc/logs/` instead of running them. It refuses a log if the script, or any `v6_` module
the script imports, changed after the log was made. A git checkout gives every file the time of the
checkout, so file times cannot show what came first. `checks/calc/logs/SOURCES.sha256` holds the
hash of each script as released, and a script that still has that hash counts as unchanged. A log
written after the last edit to its script also counts. To rerun one, for example `v6_hinf.py`:

    ( python3 checks/calc/v6_hinf.py; echo "exit $?" ) > checks/calc/logs/v6_hinf.log 2>&1

`v6_mass_smooth.py` runs four ways, and the log name carries the flags: `--w13` writes
`v6_mass_smooth_w13.log`, `--ec` writes `v6_mass_smooth_ec.log`, and both write
`v6_mass_smooth_w13_ec.log`.
