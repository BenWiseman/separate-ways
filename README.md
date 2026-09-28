# Separate Ways and the Upside Down

**The classical world and the field equations from a CPT fold, and the neutrino line that
would break it**

B. H. Wiseman · Linnet Labs, Sydney, Australia

Repository: <https://github.com/BenWiseman/separate-ways>

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

## The three documents

| Document | What it is |
|---|---|
| `paper/Separate_Ways_and_the_Upside_Down_v4.5.pdf` | The main paper. The fold, the field equations, the matter sector, and the tests. |
| `paper/Somewhere_Over_the_Horizon_v1.pdf` | The companion. Where the two halves can touch, why only in four dimensions, and where the classical world comes from. It carries the interior results the main paper cites. |
| `paper/Letter_Einstein_from_CPT_v1.pdf` | A Letter-length version of the field-equation result alone, for a venue with a hard length limit. Everything in it is in the main paper. |

Markdown sources sit beside each PDF. The arXiv metadata actually deposited is in
`paper/ARXIV_METADATA.txt` and `paper/ARXIV_METADATA_COMPANION.txt`.

## Every number is checked

Every quantitative claim in the manuscripts is produced by a script in this repository,
and a checking pass refuses to pass if a quoted number and the script that produces it
come apart. There are <!--N:calc-->169<!--/N--> such scripts under `checks/calc/`.

<!--N:calcR-->147<!--/N--> of them are R, and not one loads a package, so base R runs them
all. The remaining <!--N:calcpy-->22<!--/N--> are Python and do use numerical libraries:

```bash
python3 -m pip install -r checks/requirements.txt
```

The pass itself runs <!--N:gates-->26<!--/N--> gates over both manuscripts and the Letter:

```bash
bash checks/check_all.sh
```

It reproduces <!--N:claims-->152<!--/N--> claims against the scripts named for them,
checking every number each one names, and it also measures prose. It
catches an abstract over the arXiv character cap, a cross-reference to a section that does
not exist, a citation that does not say what the sentence says it says, a figure label
printed over by its own figure's ink, a term of art used before it is glossed, and a
sentence that appears twice. A final gate makes every checker fail on purpose, because a
check that reports success while matching nothing is worse than no check.

## Layout

| Path | What it is |
|---|---|
| `paper/` | Both manuscripts and the Letter, as PDF and as Markdown source, with their arXiv metadata |
| `paper/fig_*.pdf`, `paper/fig_*.png` | Figure masters. Re-running the pass regenerates them in place |
| `checks/` | The checking pass and the prose measures it runs |
| `checks/calc/` | One script per claim, named for the claim it supports |
| `checks/CLAIMS.tsv` | Which script backs which passage, and which numbers it must reproduce |
| `tangents/` | Calculations behind individual questions asked of the paper, organised by topic (<!--N:tangents-->139<!--/N--> scripts) |
| `paper/cross_checks/`, `paper/figure_authoring_r/` | The earlier release's scripts, kept so version 4.2 stays reproducible |
| `data/` | Acquisition instructions, pinned upstream commit and hashes for the external Pantheon+ inputs, which are not redistributed here |
| `archive/` | Superseded manuscripts and PDFs from earlier versions |
| `SHA256SUMS.txt` | Hashes of every paper and script file in this release |

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
