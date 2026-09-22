# Separate Ways and the Upside Down

**A CPT fold, a ceiling on the dark matter, and the neutrino line that would break it**

B. H. Wiseman · Linnet Labs, Sydney, Australia

Repository: <https://github.com/BenWiseman/separate-ways> *(placeholder — update once pushed)*

## What this is

If CPT — the symmetry that swaps matter for antimatter and reverses time, one of the
best-tested symmetries in physics — holds of the universe as a whole rather than only of
the laws inside it, the Big Bang has a far side: a mirror universe running the other way.
Treating the two sides as related by an antilinear fold, rather than assuming the single
minimum-energy state the original construction (Boyle, Finn and Turok) prescribes, derives
a hard ceiling on the dark-matter particle's mass — **491.6 ± 2.0 PeV** — from an operator
inequality over every state the fold permits, not a fit. The same argument fixes a
two-body neutrino line at half that mass (**245.8 ± 1.0 PeV**, sharp enough that one
well-measured event above it refutes the model), pins the neutrino-mass sum at its
kinematic floor (**58.78 ± 0.32 meV**), and gives JWST's overweight early black holes a
growth channel for free.

Every quantitative claim in the paper is reproduced by a script in this repository. There
are 178 of them, in base R with no package loaded — a reader needs only an R installation.
A cross-checking pass (`tools/consistency_pass.py`) runs 32 checks across the manuscript,
verifying that every number repeated between sections agrees, every derived number obeys
the law that derives it, and every cross-reference lands on a section that exists.

## Layout

| Path | What it is |
|---|---|
| `paper/PAPER2_v3.md` | The manuscript, Markdown source (internal filename; the manuscript itself carries no printed version number) |
| `paper/Separate_Ways_and_the_Upside_Down_v4.2.pdf` | The manuscript, built — this is the file deposited to Zenodo, version 4.2 |
| `paper/SUPPLEMENT_v3.md` | Supplement, sections S1–S13 |
| `paper/Separate_Ways_Supplement_v4.2.pdf` | Supplement, built |
| `paper/fig_*.png`, `paper/fig_*.pdf`, `paper/graphical_abstract_v3.*` | Figure masters |
| `paper/figure_authoring_r/` | R scripts that render the figures (17 scripts) |
| `paper/cross_checks/` | R scripts verifying claims used across both manuscripts (22 scripts) |
| `tangents/` | The calculations behind individual claims, one script per question asked of the paper (139 scripts, organized by topic: `landau/`, `seam/`, `stateselect/`, `stats/`, `desi/`, `lrd/`, and others) |
| `far-side-of-the-horizon/` | Companion paper — **in preparation**, see its own README |
| `data/` | Acquisition instructions, pinned upstream commit and SHA-256 hashes for the external Pantheon+ inputs (not redistributed here) |
| `tools/consistency_pass.py` | The 32-check cross-verification pass, runnable against the manuscript directly |

## Reproducing a number

Each script in `tangents/`, `paper/cross_checks/` and `paper/figure_authoring_r/` names in
its header the section and the number it produces. Run it with base R:

```bash
Rscript tangents/landau/lz_bang.R
```

No script loads a package. A handful source a sibling `helpers.R` or `_helpers.R` in their
own directory — nothing reaches outside the directory it sits in.

To re-run the full cross-check pass against the manuscript:

```bash
python3 tools/consistency_pass.py
```

## The companion

`far-side-of-the-horizon/` holds the paper's companion, on what the same fold does (and
mostly does not) at a black hole horizon. It is cited in the main paper, marked "in
preparation" throughout, and **is not yet finished** — see its own README for exactly what
that means and, more importantly, what in the main paper does and does not depend on it
(short answer: nothing in the main paper's headline results does; see §1.1 and §2 of
`paper/PAPER2_v3.md` for the explicit statement).

## License

Code (`.R`, `.py`, `.sh`, `.lua`): MIT. Manuscript, supplement, figures and documentation:
CC BY 4.0. See `LICENSE`. External data (Pantheon+) is not redistributed here; see
`data/README.md`.

## Citation

Machine-readable metadata is in `CITATION.cff`.
