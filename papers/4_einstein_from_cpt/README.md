# Einstein's equations from CPT, with nothing quantised (letter)

Submitted to Nature Physics on 29 September 2026 (NPHYS-2026-09-04663); declined by the editor for
insufficient contribution. Moved here for a Zenodo deposit and an alternative venue.

## Files

- `LETTER_Einstein_from_CPT_v2.1.md`: current source.
- `Letter_Einstein_from_CPT_v2.1.pdf`: built from that source with the vector ledger figure
  (pandoc, xelatex, 1 inch margins, US letter, 8 pages). This is the version for Zenodo.
- `Letter_Einstein_from_CPT_v1_as_sent_to_NaturePhysics_20260929.pdf`: the PDF exactly as submitted,
  for the record.
- `fig_ledger.pdf`, `fig_ledger.png`: the ledger figure as redrawn on 4 October 2026 (no diagonal
  connector; centre box "CPT of the universe, fused with algebraic QFT"). The PDF uses the vector
  `.pdf`; the markdown points to the `.png`.

v2, source and PDF, is in `archive/`. Its source is `nobel_shot/claude/pub/paper2/LETTER_PRL_v1.md`
as last edited on 1 October 2026.

To rebuild the PDF from this directory:

    sed 's/!\[\](fig_ledger.png)/![](fig_ledger.pdf)/' LETTER_Einstein_from_CPT_v2.1.md > build.md
    pandoc build.md -o Letter_Einstein_from_CPT_v2.1.pdf --pdf-engine=xelatex -V geometry:margin=1in

## v2 against the version sent

- Figure 1: the redrawn ledger.
- "The first property is worth a moment" now reads "needs a moment".
- "the contact budget in the late universe" now reads "the turning available in the late universe".
- The AI disclosure now points to Appendix C for where the scripts live, instead of to the release of
  ref. 5.

Nothing else in the text changed.

## v2.1 against v2 (6 October 2026)

Prose only. No number, equation, reference or claim changed, and every number still matches the
manuscript of ref. 5.

- The last financial words are out. "The result would be worth nothing" now reads "the derivation
  would assume what it set out to prove", and Appendix A row 11 now reads "Contact needs a turning
  of $\pi$ and the interior supplies $\pi/(D-3)$ per leg, which closes only at $D=4$".
- Headings and labels no longer open on "The": "A temperature at every horizon", "Horizons at every
  boost and orientation", "Why the entropy law is not circular", "Inputs and outputs", "Nineteen
  that come out", "Four that go in", "Two objections to the 1995 argument".
- The sentence calling the entropy step load-bearing, and saying the danger is stated before the
  answer, is gone.
- Terms a reader met before their definition are now defined where they are used: $\Phi_c$ and
  $\Phi_q$, $W(0)$ and $W_{\rm cross}$, $M_1$ and the production integral $I$, the function $f$,
  and $J$ and $P_\perp$ in the horizons paragraph. "Exact stabilisation", never defined in the
  letter, is gone, and the abstract no longer uses "the fold" before saying what it is.
- Fewer contrastive "X rather than Y" and "X and not Y" constructions (23 to 15), and no "is not X.
  It is Y."
- "differing by $3\times10^{15}$" now reads "larger than its tenth by a factor of
  $3\times10^{15}$", which is what `checks/calc/jacobson_rederived.R` prints.
- Appendix A's closing line ("two are units") is gone; the caption and the main text already say two
  of the inputs are measured constants.
- Long sentences split at their joins, and the source reflowed to 100 columns.

Checked with the benlm prose tools, tidyhumanize, a paragraph-by-paragraph pass with M3 and DeepSeek
each reading a paragraph after everything before it, and the letter gates in
`nobel_shot/claude/far_side`: numbers against the manuscript, length, editorial voice, and the
ledger table against `relativity_ledger.R`.
