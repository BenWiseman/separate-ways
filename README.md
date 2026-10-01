# Separate Ways and the Upside Down

**The classical world and the field equations from a CPT fold, and the neutrino line that
would break it**

B. H. Wiseman · Linnet Labs, Sydney, Australia

Repository: <https://github.com/BenWiseman/separate-ways>

This repository holds the papers in this line of work, their figures, and every script
behind every number in them. Start at `papers/1_separate_ways/`.

## What this is

CPT swaps matter for antimatter and reverses time, and it is one of the best-tested
symmetries in physics. Take it to hold of the universe as a whole rather than only of the
laws inside it. The Big Bang then has a far side, a mirror sheet running the other way,
and the two sheets are related by an antilinear fold.

The main result is that Einstein's field equations follow from that fold. Jacobson showed
in 1995 that the field equations can be derived from the Clausius relation applied to
local Rindler horizons, given three things he had to assume: a temperature at every such
horizon, an entropy proportional to that horizon's area with the same coefficient
everywhere, and a horizon of that kind at every point in every null direction. The fold
supplies all three. The entropy law comes from the entanglement of the state the fold
forces rather than from Bekenstein-Hawking, so nothing is borrowed back from the theory
being derived. The metric is a fixed background at every step, and no gravitational
degree of freedom is quantised anywhere. Nineteen properties of general relativity come
out against four inputs, among them the sign of Newton's constant, four large spacetime
dimensions, and the cosmological constant as an integration constant rather than a
coupling.

The same fold fixes matter, which is what makes the construction breakable. Treating the
two sides as related by an antilinear fold, rather than assuming the single
minimum-energy state the original Boyle, Finn and Turok construction prescribes, gives a
hard ceiling on the dark-matter particle's mass, **491.6 ± 2.0 PeV**, from an operator
inequality over every state the fold permits rather than from a fit. It fixes a two-body
neutrino line at half that mass, **245.8 ± 1.0 PeV**, sharp enough that one well-measured
event above it refutes the model. It floors the neutrino-mass sum at **58.8 meV** and
commits to constant dark energy, which shuts the $w_0w_a$ escape that relaxes the present
cosmological bound to 163 meV. And it gives JWST's overweight early black holes a growth
channel for free.

Four observations would end it: a securely measured variation in Newton's constant at any
level, a robust detection of evolving dark energy, a cosmological neutrino-mass sum below
58.8 meV, and an identified decay above the two-body endpoint.

## The papers

Three. Paper 1 is the foundation and paper 2 depends on it; paper 3 is a Letter carrying
one result of paper 2 on its own. Read them in that order.

| | |
|---|---|
| **1. Separate Ways and the Upside Down** | The fold, the field equations, the matter sector and the tests. 50 pages. `papers/1_separate_ways/` · [doi:10.5281/zenodo.22888119](https://doi.org/10.5281/zenodo.22888119) |
| **2. Somewhere Over the Horizon** | Where the two halves of a folded universe can touch, why only in four dimensions, and where the classical world comes from. Carries the interior results the first paper cites. 88 pages. `papers/2_over_the_horizon/` · [doi:10.5281/zenodo.23030633](https://doi.org/10.5281/zenodo.23030633) |
| **3. Road to Nowhere** | A Letter. Where a point can reach its own image on the far side, and why the answer picks out four dimensions. Paper 2 states the result and keeps the extensions; this derives it. 5 pages. `papers/3_road_to_nowhere/` · [doi:10.5281/zenodo.23056601](https://doi.org/10.5281/zenodo.23056601) |

Each directory holds one paper and everything that belongs to it: the PDF, the Markdown
source it was built from, the figures that source references, and the arXiv metadata
actually deposited. Nothing else is in there, and nothing of a paper's is anywhere else.

## Every number is checked

Every quantitative claim in the manuscripts is produced by a script in this repository,
and a checking pass refuses to pass if a quoted number and the script that produces it
come apart. There are <!--N:calc-->175<!--/N--> such scripts under `checks/calc/`.

<!--N:calcR-->149<!--/N--> of them are R, and not one loads a package, so base R runs them
all. The remaining <!--N:calcpy-->26<!--/N--> are Python and do use numerical libraries:

```bash
python3 -m pip install -r checks/requirements.txt
```

The pass itself runs <!--N:gates-->28<!--/N--> gates over both manuscripts:

```bash
bash checks/check_all.sh
```

It reproduces <!--N:claims-->157<!--/N--> claims against the scripts named for them,
checking every number each one names, and it also measures prose. It
catches an abstract over the arXiv character cap, a cross-reference to a section that does
not exist, a citation that does not say what the sentence says it says, a figure label
printed over by its own figure's ink, a term of art used before it is glossed, and a
sentence that appears twice. A final gate makes every checker fail on purpose, because a
check that reports success while matching nothing is worse than no check.

## Layout

| Path | What it is |
|---|---|
| `papers/1_separate_ways/` | Paper 1: PDF, Markdown source, its six figures, arXiv metadata |
| `papers/2_over_the_horizon/` | Paper 2: PDF, Markdown source, its seventeen figures, arXiv metadata |
| `papers/3_road_to_nowhere/` | Paper 3: PDF, Markdown source, its figure |
| `checks/` | The checking pass and the prose measures it runs |
| `checks/calc/` | One script per claim, named for the claim it supports |
| `checks/CLAIMS.tsv` | Which script backs which passage, and which numbers it must reproduce |
| `tangents/` | Calculations behind individual questions asked of the paper, organised by topic (<!--N:tangents-->139<!--/N--> scripts) |
| `archive/superseded_figures/` | Figures from earlier releases that no current manuscript references |
| `archive/release_4_2_scripts/` | The 4.2 release's cross-checks and figure authoring, kept so that version stays reproducible |
| `data/` | Acquisition instructions, pinned upstream commit and hashes for the external Pantheon+ inputs, which are not redistributed here |
| `archive/` | Superseded manuscripts and PDFs from earlier versions |
| `SHA256SUMS.txt` | Hashes of every paper and script file in this release |

Figures sit with the paper that uses them, and a figure no manuscript references is in
`archive/superseded_figures/` rather than beside the live ones. Twelve were, until
2026-09-30.

## Reproducing a number

Each script names in its header the passage and the number it produces. Run it directly:

```bash
Rscript checks/calc/relativity_ledger.R
```

No script loads a package. A handful source a sibling helpers file in their own directory,
and nothing reaches outside the directory it sits in.

## License

Code (`.R`, `.py`, `.sh`): MIT. Manuscripts, figures and documentation: CC BY 4.0. See
`LICENSE`. External data (Pantheon+) is not redistributed here; see `data/README.md`.

## Citation

Machine-readable metadata is in `CITATION.cff`.

## A note on how this was written

Contemporary research tools were used throughout, large language models among them. No
such tool is an author. Every number in the papers is produced by a script in this
repository and checked against it, which is the point of the apparatus above: the claims
stand on the code rather than on anyone's word, including the tools'.
