# What a CPT-symmetric universe adds to the Einstein equation of state (letter)

Submitted to Nature Physics on 29 September 2026 (NPHYS-2026-09-04663) as "Einstein's equations
from CPT, with nothing quantised"; declined by the editor for insufficient contribution. v3 corrects
how every earlier version described Jacobson's 1995 derivation, and is retitled to match.

## Files

- `LETTER_CPT_Einstein_equation_of_state_v3.md`: current source.
- `Letter_CPT_Einstein_equation_of_state_v3.pdf`: built from that source with the vector ledger
  figure (pandoc, xelatex, 1 inch margins, US letter, 8 pages). This is the version for Zenodo.
- `fig_ledger.pdf`, `fig_ledger.png`: the ledger figure as redrawn on 6 October 2026, with the
  nineteen outputs grouped by where each stood before. The PDF uses the vector `.pdf`; the markdown
  points to the `.png`.

v2 and v2.1, source and PDF, are in `archive/`.

To rebuild the PDF from this directory:

    sed 's/!\[\](fig_ledger.png)/![](fig_ledger.pdf)/' LETTER_CPT_Einstein_equation_of_state_v3.md > build.md
    pandoc build.md -o Letter_CPT_Einstein_equation_of_state_v3.pdf --pdf-engine=xelatex \
        -V geometry:margin=1in -V papersize=letter

## v3 against v2.1 (6 October 2026)

v2.1 and every version before it said Jacobson's derivation assumed a horizon temperature, an area
law and horizons through every point, that the fold supplies all three, and that the fold's entropy
is not taken from Bekenstein-Hawking. That misdescribes his paper (gr-qc/9504004, checked at
source). He took the temperature from the Unruh effect, the area law from the entanglement of vacuum
fluctuations with a fixed cutoff (citing Bombelli, Koul, Lee and Sorkin 1986), and the horizons
from local flatness. He wrote $\Lambda$ as the integration constant of the Bianchi step, gave
$G=1/4\hbar\eta$, and argued that an equation of state need not be quantised. Neither he nor the fold
uses Bekenstein-Hawking. By Bisognano-Wichmann the wedge conjugation composed with the transverse
antipode is CPT itself, so the fold's map at a flat wedge is the operator behind his inputs.

- Title: was "Einstein's equations from CPT, with nothing quantised".
- Abstract and opening rewritten on that basis: one map underlies Jacobson's quantum inputs and
  defines the folded cosmology; the letter re-derives his argument from it and says what the fold
  adds.
- Each re-derivation paragraph now says where Jacobson had the same input. "What the derivation
  needs" is now "Jacobson's argument, run on the fold"; "Why the entropy law is not circular" is now
  "An entropy proportional to area".
- The nineteen outputs are grouped: nine in Jacobson's argument, four standard results, six that
  need the folded universe. Figure 1 and Appendix A carry the groups. The gate that checks
  Appendix A against `relativity_ledger.R` now also checks each line's group.
- Row 17 no longer calls the weight $4\tanh^2(\beta\omega/4)$ the fold's. The weight follows from
  the standard symmetric thermal contour; what is the fold's is that the two parts are its even and
  odd sectors.
- The closing claim that the field equations need not wait on a quantum theory of gravity is
  credited to Jacobson.
- Twelve references added, each checked on Crossref or arXiv: Unruh 1976, Bombelli, Koul, Lee and
  Sorkin 1986, Srednicki 1993, Jacobson 2016, Padmanabhan 2010, Verlinde 2011, Faulkner, Guica,
  Hartman, Myers and Van Raamsdonk 2014, Sakharov 1967, Susskind and Uglum 1994, Hawking 1975,
  Sewell 1982, Kay and Wald 1991.
- Appendix A's tables carry weighted column widths, so the first column no longer takes a third of
  the page.

No number changed. Every number still matches the manuscript of ref. 5:
`checks/letter_numbers_check.py` checks 152 numeric expressions and finds none missing. That gate
no longer reads citation brackets as numbers, and its appendix plant was re-anchored after v2.1
deleted the sentence it used to plant into.
