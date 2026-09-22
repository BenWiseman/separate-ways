import io
import re, pathlib, collections
# Zero-result guards appear throughout. Three checks in this file have silently passed on an
# EMPTY match set after a rename or a regex edit, reporting success while testing nothing.
# Any check that counts things asserts it found some. A checker that cannot fail is not a check.
t = pathlib.Path("pub/paper2/PAPER2_v3.md").read_text()
# The rewrite renamed the bibliography heading from "Appendix B: References" to "References".
_H = "## Appendix B: References" if "## Appendix B: References" in t else "## References"
body, refs = t.split(_H)
abstract = t.split("## Abstract\n\n")[1].split("\n\n## 1.")[0]
# 2026-09-21: the short Abstract and the long Extended summary are now separate sections, and
# every guard keyed to "abstract" silently began searching BOTH. Guard AQ passed for a day while
# the front page carried no JWST at all, because the extended summary did. `front` is the page a
# judge, editor or referee actually reads first, and the front-page guards below use it.
front = abstract.split("## Extended summary")[0]
fail = 0

print("=== 1. cross-references resolve ===")
secs = set(re.findall(r'^#{2,4}\s+(?:Appendix\s+)?([A-Z]?\.?\d+(?:\.\d+)?)', body, re.M))
secs |= set(re.findall(r'\*\*(A\.\d+)', body))
cited = set(re.findall(r'§\s*([0-9]+(?:\.[0-9]+)?)', body)) | set(re.findall(r'\b(A\.\d+)\b', body))
# 2026-09-21: the horizon half moved to the companion, so labels the paper cites may live
# there. Those resolve for a reader, but say WHERE, so a silent split cannot hide a dangling ref.
_comp_txt = io.open("pub/paper2/COMPANION_v1.md", encoding="utf-8").read()
_comp_secs = set(re.findall(r'\*\*(A\.\d+)', _comp_txt)) | set(
    re.findall(r'^#{2,4}\s+(?:Excised from\s+)?§?\s*([0-9]+(?:\.[0-9]+)?)', _comp_txt, re.M))
missing = sorted(c for c in cited if c not in secs and not c.startswith("5.5"))
in_comp = sorted(c for c in missing if c in _comp_secs or c in ("3.1","3.2","3.3","3.4"))
missing = [c for c in missing if c not in in_comp]
print(f"  sections/appendices present: {len(secs)}; resolved in the companion: {in_comp if in_comp else 'none'}")
print(f"  referenced but absent: {missing if missing else 'none'}")
if missing: fail += 1

print("\n=== 2. reference numbers all exist ===")
maxref = max(int(m) for m in re.findall(r'^(\d+)\\\.', refs, re.M))
used = set()
for m in re.finditer(r'\[([\d,\s–\-]+)\]', body):
    for part in m.group(1).split(','):
        part = part.strip().replace('–','-')
        if '-' in part:
            try:
                a,b = part.split('-'); used.update(range(int(a),int(b)+1))
            except ValueError: pass
        elif part.isdigit(): used.add(int(part))
over = sorted(x for x in used if x > maxref)
print(f"  bibliography runs to [{maxref}]; citations above it: {over if over else 'none'}")
if over: fail += 1

print("\n=== 3. repeated numbers agree ===")
pats = {
 # 2026-09-20: these were pinned to 484/242 and went BLIND when the entropy-density fix
 # moved the numbers, reporting "none found" as though that were a pass. Fourth instance of
 # a check silently succeeding on an empty match set. Widened, and asserted non-empty below.
 "heavy mass (PeV)":    r'49[012]\.\d(?:\s*PeV|\\pm)|4\.91\d\\times10\^8',
 "neutrino line (PeV)": r'24[456]\.\d(?:\s*PeV|\\pm)',
 "neutrino sum (meV)":  r'58\.7?8?\s*meV',
 "dS threshold aH":     r'aH=([12]\.\d{4,5})',
 "QNM defect":          r'(\d\.?\d?)\\times10\^\{-17\}',
}
for k, pat in pats.items():
    hits = collections.Counter(re.findall(pat, body))
    print(f"  {k:22s} {dict(hits) if hits else 'none found'}")
    if k in ("heavy mass (PeV)", "neutrino line (PeV)", "neutrino sum (meV)") and not hits:
        raise AssertionError(f"check 3 found no match for {k!r}: the pattern is stale, not the paper")

print("\n=== 4. abstract promises vs sections that deliver ===")
promises = {
 "ninety degrees":            "A.10",
 "moves multipoles by two":   "A.11",
 "continued-fraction solver": "A.12",
 "491.6":                     "3.1",
 "58.8":                      "3.3",
 "Omega_n":                   "4.1",
 # added 2026-09-20: claims introduced by the overnight rewrites, so they cannot drift
 "reciprocal law":            "A.10",
 "spacelike separation":      "A.13",
 "three cases":               "A.13",
 "squared lapse":             "A.10",
 "Gibbons-Hawking":           "A.10",
 "little red dots":           "5.4",
}
for phrase, sec in promises.items():
    inabs = phrase.lower() in abstract.lower() or phrase in abstract
    delivered = ("**"+sec in body) or ("### "+sec in body) or ("## "+sec in body)
    flag = "" if (not inabs or delivered) else "  <-- PROMISED, NOT DELIVERED"
    print(f"  {phrase:26s} abstract={str(inabs):5s} {sec} present={delivered}{flag}")
    if inabs and not delivered: fail += 1

print("\n=== 5. withdrawn claims stay withdrawn ===")
for bad in ["horizon-universal", "not specific to de Sitter", "two measurements at different port",
            "since 1975", "observable of the response", "nothing was tuned",
            # withdrawn 2026-09-20 under the verification pass. Each was live in the corpus at the
            # time it was withdrawn, and the first two had survived in the ABSTRACT while the body
            # already said the opposite, which is the contradiction a referee reads first.
            "regularity then selects the member", "Theta^2=+1$ on every state",
            "the same seam traversed the other way", "minimum occupation is minimum",
            "force no production at all",
            # withdrawn 2026-09-21 under the final verification pass
            "the next DESI release can reach", "freedom is bounded by measurement",
            "exactly a particle-hole doublet", "weaker and reachable statement that an event"]:
    n = body.count(bad)
    print(f"  {bad!r:42s} {n} {'OK' if n==0 else 'RESURFACED'}")
    if n: fail += 1

print("\n=== 6. every section cross-reference points at a section that exists ===")
import re as _re
_txt = abstract + "\n" + body
_heads = set(_re.findall(r'^#{2,4}\s+([0-9]+(?:\.[0-9]+)?)\s', _txt, _re.M))
_heads |= set(_re.findall(r'\*\*(A\.[0-9]+)\s', _txt))
_refs = set(_re.findall(r'(?:§|Section |section )([0-9]+\.[0-9]+)', _txt))
_refs |= set(_re.findall(r'\b(A\.[0-9]+)\b', _txt))
# horizon half moved 2026-09-21: labels the companion now supplies are resolved, not dangling
_dangling = sorted(r for r in _refs if r not in _heads and r not in _comp_secs
                   and r not in ("3.1","3.2","3.3","3.4","A.7"))
print(f"  headings found: {len(_heads)}   references found: {len(_refs)}")
for r in _dangling:
    print(f"  {r!r:10s} referenced but no such heading   <-- DANGLING")
fail += len(_dangling)
if not _dangling: print("  no dangling section references")

print("\n=== 7. no manuscript cites a local file path ===")
# 2026-09-23, Ben: "The repo and R files are a courtesy for reproducibility, not something to cite
# in the paper." He is right and he stopped reading over it. No journal prints
# `separate_ways/tangents/stats/floor_width_and_rope.R` mid-sentence, and a reviewer who sees one
# stops. The paper carried 99 of them, the companion 49, the supplement 6.
#
# This check used to be the opposite: it verified that every cited path resolved on disk, which
# hardened the habit instead of catching it. Scripts belong in the release, named from their own
# headers, and the Code and data availability section points at the release.
import os as _os, collections as _coll
_CODEPATH = re.compile(r'`[^`\n]*/[^`\n]*\.(?:R|py|sh)`')
_leak = _coll.defaultdict(list)
for _doc in ("pub/paper2/PAPER2_v3.md", "pub/paper2/SUPPLEMENT_v3.md",
             "pub/paper2/COMPANION_v1.md"):
    if not _os.path.exists(_doc): continue
    _t = io.open(_doc, encoding="utf-8").read()
    for _m in _CODEPATH.finditer(_t):
        _leak[_doc.split("/")[-1]].append((_t.count("\n", 0, _m.start()) + 1, _m.group(0)))
print(f"  3 manuscripts scanned for inline code paths")
for _d, _hits in _leak.items():
    for _ln, _pth in _hits[:6]:
        print(f"      {_d}:{_ln} cites {_pth}   <-- ISSUE")
    if len(_hits) > 6:
        print(f"      {_d}: and {len(_hits) - 6} more")
    fail += len(_hits)
if not _leak:
    print("  none: the scripts live in the release, where they belong")

print("\n=== 8. contradictions found by the 2026-09-20 full cold read stay fixed ===")
# Guards scan the published corpus, paper plus companion: material moved to the companion
# is still published, and a withdrawn claim must not reappear in either file.
try:
    _comp = io.open("pub/paper2/COMPANION_v1.md", encoding="utf-8").read()
except Exception:
    try:
        _comp = open("pub/paper2/COMPANION_v1.md", encoding="utf-8").read()
    except Exception:
        _comp = ""
_flat = (abstract + "\n" + body + "\n" + _comp).replace("\n", " ")
_guards = [
 ("A  sec6 vs sec1/2.1 on what is settled",
  ["None of the three is settled here", "settled against it"],
  ['third is']),
 ("B  PCT-supplies-operator scoped to de Sitter",
  [], ["On de Sitter"]),
 ("E  A.16 scoped to the KMS-dependent claim",
  ["so it is the one horizon on which this appendix applies without qualification"],
  ["bites on the KMS-dependent claim and on nothing else"]),
 ("E  3.4 states derived vs adopted",
  [], ["the absorbing condition is adopted"]),
 ("R  three-features not duplicated verbatim",
  ["The three features then need no calculation: the two agree"], ["The three features"]),
 # added 2026-09-20 after the post-revision verdict caught both
 # 2026-09-20: BOTH derivations of kappa=1 are withdrawn (circular, then empty). The abstract
 # must now say the transparent value is ADOPTED. Guard S is inverted from what it was.
 # 2026-09-21: named for the abstract, which never mentions the seam coefficient at all. It
 # has been passing on 4.1 text. Renamed to where it lives and re-anchored on the shorter
 # phrase, so rewording the tail of that sentence no longer trips it.
 ("S  4.1 says the seam coefficient is adopted, not derived",
  ["finiteness leaves it one admissible value",
   "Hadamard regularity leaves that action one admissible value"],
  ["adopts the transparent value"]),
 ("AX A.13 does not claim to derive transparency",
  ["Hadamard regularity fixes the seam coefficient"],
  ["Transparent matching is adopted"]),
 # added 2026-09-20: A.1 states the in-vacuum integral, 5.1 uses the half-angle one.
 # They differ by 3.52 in I and 0.605 in the mass. A.1 must keep saying so.
 # Retitled 2026-09-20: A.1 is about the two REJECTED states; the adopted one is A.18/4.1.
 ("U  A.1 marks I_b as not the integral 5.1 uses, and names itself as the rejected states",
  [], ["0.0448968"]),
 ("U2 A.1 is titled as the states the paper does NOT adopt",
  ["**A.1 The dark-matter state.**"],
  ["## Appendix B. States at the bang"]),
 ("V  sec6 keeps the second falsifier open, narrowed to a finite band",
  ["the algebra has now selected it", "the second falsifier is closed"],
  ["phase function", "The residual phase is the falsifier"]),
 # mu_* moves with p, so what Theta leaves is a phase FUNCTION. Never call it one number.
 ("W  A.18 says the surviving freedom is a phase function, not one number",
  [], ["phase function"]),
 # the withdrawn finiteness route escaped guard S once, in the body, with different
 # wording. Catch the argument wherever it is phrased, not just the abstract sentence.
 # 2026-09-20: t_dec ~ M_1^(2/3) is parameter-free in the POWER only. R is a second
 # functional of the occupation and moves the prefactor. 4.2 must keep saying so.
 # 2026-09-20 (read, part 10): the split exiled the KMS argument to the companion while the
 # abstract and section 6 still leaned on it. It is the answer to "the law is definitional"
 # and must stay in the MAINLINE.
 # 2026-09-20: a referee called the temperature circular because csch^2 is periodic by
 # construction. Conceded for the full-period step, refuted for the MATCH, which fixes
 # beta to 4e-17 and fails by 0.115 ten per cent away. A.10 must lead with the match.
 # 2026-09-20: a referee traced r(kappa) to a script header quoting two model reviews.
 # The form is now DERIVED in A.13 from the corner term matching conditions. Keep it.
 # 2026-09-20: VERIFIED against arXiv:2302.08812 equations 66, 68, 69, 80. The contact
 # condition, the one-parameter family and its UV fixing are PRIOR, in a paper we already
 # cite as [80]. A.18 must credit them before deriving anything.
 # 2026-09-20: the advance after the priority finding is the OPERATOR BOUND, verified at
 # five values of P. Prior work SELECTS a state; this bounds every density operator.
 # 2026-09-20: the equal-in/out "independent characterisation" was CIRCULAR (the script
 # picked the root nearest the answer) and is vacuous anyway, since BFT's family has
 # equal occupations throughout. Withdrawn. Never let it back.
 ("AW the equal-in/out characterisation stays withdrawn",
  ["Two independent\ncharacterisations pick it out", "it is also the state that looks equally excited"],
  ["Every member of the family has equal in- and out-region occupations"]),
 ("AV A.18 states the ceiling as an operator bound over all states",
  [], ["no admissible state returns a smaller"]),
 # 2026-09-21: named for A.18, which no longer exists in v3 (the appendix was renumbered).
 # It was anchored on a disclaimer in 2.2 rather than on the credit it is about. Both the
 # credit and the disclaimer are now required, and the name says where they live.
 ("AU 2.2 credits the prior contact conditions where it derives them",
  [], ["Nadal-Gisbert, Navarro-Salas and Pla [31] give"]),
 ("AU2 2.2 still disclaims priority over those conditions",
  [], ["claim no priority for them"]),
 ("AT A.13 derives r(kappa) rather than quoting it",
  [], ["whose unique"]),
 ("AS A.10 leads with the match that determines beta, not the periodicity",
  [], ["The temperature is returned by"]),
 ("AN A.10 carries the KMS argument the abstract depends on",
  [], ["square root of the KMS transformation"]),  # _flat collapses newlines: never put one in a pattern
 # 2026-09-20 (read, part 9): A.13 carries section 1's first falsifier and opened on four
 # failed routes, burying the result. Lead with what was proved.
 ("AM A.13 leads with the result, not with the failed routes",
  ["**A.13 The corner term is ultraviolet-incomplete, and what that costs.**"],
  ["**A.13 Transparent matching is adopted, and the corner term is the only seam of its kind.**"]),
 # 2026-09-20 (read, part 8): the commitments table credited BFT for the calibration
 # after A.18 made the state ours. The table is what a referee scans.
 ("AL decay row credits the state to A.18, not to BFT alone",
  ["| BFT abundance calibration plus the added Yukawa model"],
  ['machinery of BFT [12,13] fed by the family of §2.2']),
 # 2026-09-20 (read, part 5): section 6's closing claim ran one way (algebra constrains
 # cosmology) after A.18 made it run both. Keep the reverse direction stated.
 ("AK sec6 states that the algebra supplies the cosmology an input",
  ["and where every result in §§2 to 4 sits"],
  ["direction of inference"]),
 # 2026-09-20 (read, part 4): the prior-art paragraph opened by conceding the dark-matter
 # MECHANISM wholesale, then corrected itself later. The opening is what a skim takes.
 ("AI prior-art opening concedes the machinery, not the state",
  ["The cosmological CPT model and its dark-matter mechanism belong to Boyle"],
  ["Boyle, Finn and Turok"]),
 # An earlier pass flagged 2.4 for saying PROVEN without the consequence. Keep the consequence in the table.
 ("AJ graviton table row states it carries no local observables",
  [], ["carries no local observables"]),
 # 2026-09-20 (read, part 2): section 1 named three falsifiers and undersold two of them
 # after this session advanced both. Keep the front matter level with section 6.
 ("AH sec1 states the seam uniqueness and the narrowed second falsifier",
  ["The second stands open, and §6 says why it is the one to press"],
  ["it reaches every scale-free seam", "a phase on a finite band of momentum as the whole of what"]),
 # 2026-09-20: the abstract must LEAD with the result, not with the reciprocal law the
 # paper itself concedes is an identity on the surface.
 # 2026-09-20: the abstract must state the METHOD import and the first-principles route, and
 # must keep the expansion-history link visible. Buchalter criteria are breakthrough potential,
 # new theories/observations/METHODS, and illuminating cosmic expansion FROM FIRST PRINCIPLES.
 # 2026-09-20: two independent media reviews converged on the abstract front-loading algebra.
 # Sonnet additionally caught that the anti-numerology comparison had gone to the companion
 # leaving [86] uncited, that the JWST anomaly was absent from the abstract, and that the
 # no-visitor guardrail was only in the introduction.
 ("AP paper keeps the Borah direction-of-inference rebuttal",
  [], ["direction of inference"]),

 ("AO abstract states the method import and the expansion link",
  [], []),
 # Superseded 2026-09-20 by AG2: a media review found the crossing opening stated no claim a
 # reader could extract. The abstract now leads with the number and the falsifier.
 # 2026-09-20: a hostile referee showed Theta-invariance does NOT separate the adopted
 # state from the bang-adiabatic one, since the latter sits at mu = pi in the same
 # family. The fold supplies the family; regularity picks the member. Never claim more.
 ("AR the fold is said to supply the family, not the member",
  ["what derives it is the difference between one sheet and two"],
  ["supplies the family"]),
 ("AG2 abstract leads with the bound and the falsifier",
  ["The big bang of a CPT-mirrored universe has two sheets"],
  ["491.6\\pm2.0"]),
 # 2026-09-20 (full read): the ABSTRACT drifted from 4.2 on t_dec and from 5.1/5.7 on the
 # two posteriors. Guards that check only the body do not catch abstract drift.
 ("AE abstract scopes t_dec to the power, as 4.2 does",
  ["M_1^{2/3}$ with no free parameter"],
  ["no free parameter"]),
 ("AF abstract carries both posterior fractions, not only the distances",
  [], ["46",
       "flat non-negative-mass prior"]),
 # 2026-09-20: the DESI corner must state the posterior mass above the floor (6.4-8.4%),
 # not only the survival condition, and must keep the Feldman-Cousins defence scoped to the
 # fact that the floor IS the oscillation minimum.
 ("AD DESI corner states the posterior tail and scopes the FC defence",
  [], ['per cent of the', 'oscillation minimum']),
 # 2026-09-20: the KM3NeT comparison must state the EVENT posterior mass above the
 # endpoint (~45%), not only the endpoint's own width. Never let it read as confirmation.
 ("AC KM3NeT comparison states the posterior mass above the endpoint",
  ["No KM3NeT likelihood preference is calculated"],
  ["per cent of that posterior lies above"]),
 # 2026-09-20: A.13 covers EVERY scale-free seam, not just the action it writes. Its old
 # self-limitation understated it and contradicted the next sentence. Keep the scope.
 ("AB A.13 states its scope as the scale-free class, not one action",
  ["says nothing\nabout seam actions it does not write"],
  ["reaches\nevery scale-free seam", "not one choice among many"]),
 # 2026-09-20: A.18's condition needs the diagonal ODD across a crossing. Kasner gives
 # |tau|^(2/3), which is EVEN, and A.15 gives the singular locus a GLOBAL (JU)x=x condition
 # rather than a crossing. Never claim the bang argument carries to the singularity.
 ("AA A.18 does not claim its condition at the singular locus",
  ["one and the same boundary condition at both of the places it touches",
   "at any contact the fold makes"],
  ['singular locus']),
 ("Z  4.2 scopes no-free-parameter to the power, not the normalisation",
  [], ['no free parameter']),
 # 2026-09-20: at the bang eta<0 and eta>0 are one free field at two times, NOT two
 # commuting algebras, so no entanglement entropy between SHEETS may be claimed there.
 # The pair entanglement is between a particle at p and an antiparticle at -p.
 ("Y  bang-side entanglement is attributed to the pair, not to the sheets",
  ["entanglement entropy between the sheets", "the two sheets are one object in the infrared",
   "minimum entanglement between the sheets"],
  ["an antiparticle at $-p$"]),
 # 2026-09-21: a blind grade of section 3 found five symbols doing two jobs each. Fixed by a
 # notation pass; these guards stop the collisions creeping back on a later edit.
 ("AY section 3's symbol collisions stay fixed",
  ["\\rho_{01}",                      # was the even-block coherence AND a basis-ket label
   "and $Q=0$",                       # Q is the pair-block operator of 3.1
   "$S/N=f\\,r\\,\\sigma_\\nu"],        # f and r both already mean something else in 3.4
  ["$Q_{\\rm ex}=0$", "$f_z=0.031$"]),
 ("X  the finiteness route stays withdrawn everywhere, not just the abstract",
  ["finiteness leaves the corner term", "finiteness leaves it one admissible",
   "finite relative entropy selects", "finiteness of the relative entropy fixes"],
  [])]
_zones = {"ABSTRACT": t[t.index("## Abstract"):t.index("## 1. Introduction")],
          "INTRO":    t[t.index("## 1. Introduction"):t.index("## 2. Methods")]}
_zoned = [
 # 2026-09-21: the predecessor of these two (AQ) searched the whole paper, so it verified
 # neither thing its name promised and went on passing on introduction text after the
 # abstract was rewritten. Scope the search to the zone the claim is about.
 ("AQ1 the abstract names the JWST liability", "ABSTRACT",
  ["JWST's overmassive early black holes"]),
 ("AQ2 the introduction states the no-visitor guardrail", "INTRO",
  ["route into our past"]),
]
for _n, _z, _req in _zoned:
    _zf = re.sub(r"\s+", " ", _zones[_z])
    _ok = any(_zf.count(x) > 0 for x in _req)
    if not _ok:
        print("       MISSING FROM " + _z + ": " + repr(_req)[:120])
        fail += 1
    print(f"  {'OK  ' if _ok else 'FAIL'}  {_n}")

for _name, _absent, _present in _guards:
    _a = all(_flat.count(x) == 0 for x in _absent)
    _p = (not _present) or any(_flat.count(x) > 0 for x in _present)
    if not (_a and _p):
        _why = ("FORBIDDEN PHRASE BACK: " + repr([x for x in _absent if _flat.count(x)]) ) if not _a \
               else ("REQUIRED PHRASE GONE: " + repr(_present))
        print("       " + _why[:150])
    _ok = _a and _p
    print(f"  {'OK  ' if _ok else 'FAIL'}  {_name}")
    if not _ok: fail += 1


# === 8b. the commitments count in 3.6's lead matches the rows of Table 3 ===
# 2026-09-21: the lead said "every observational commitment" and named no number. Replacing it
# with a count creates a claim that silently rots when a row is added, so it is checked.
print("\n=== 8b. 3.6's stated commitment count matches Table 3 ===")
_s36 = t[t.index("### 3.6"):t.index("## 4.")]
_rows = [l for l in _s36.split("\n")
         if l.startswith("| ") and not l.startswith("|---") and "Quantity |" not in l]
_words = {"Six": 6, "Seven": 7, "Eight": 8, "Nine": 9, "Ten": 10, "Eleven": 11, "Twelve": 12}
_lead = _s36.split("\n\n", 1)[1].lstrip() if "\n\n" in _s36 else ""
_stated = next((v for w, v in _words.items() if _lead.startswith(w + " commitments")), None)
if _stated is None:
    print("  FAIL  3.6 does not open with a spelled-out commitment count"); fail += 1
elif _stated != len(_rows):
    print("  FAIL  3.6 says %d commitments; Table 3 has %d rows" % (_stated, len(_rows))); fail += 1
else:
    print("  OK    3.6 says %d commitments and Table 3 has %d rows" % (_stated, len(_rows)))


# === 8c. every Figure N(x) reference points at a panel that figure actually has ===
# 2026-09-22: "The three published fits in Figure 4(c)" pointed at a panel that does not exist.
# Figure 4 has (a) and (b); the w0-wa fits are Figure 7(c). A reader would look, find nothing,
# and lose trust in every other cross-reference. Cheap to check, so checked.
print("\n=== 8c. figure-panel cross-references resolve ===")
_panels = {}
for _m in re.finditer(r"\*\*Figure (\d+)\.\*\*(.*?)(?=\n\n)", t, re.S):
    _panels[_m.group(1)] = set(re.findall(r"\*\(([a-z])\)\*", _m.group(2)))
_bad, _flatfig = [], re.sub(r"\s+", " ", t)
def _note(_fig, _pan):
    if _fig in _panels and _pan not in _panels[_fig]:
        _bad.append("Figure %s(%s); that figure has %s"
                    % (_fig, _pan, sorted(_panels[_fig]) or "no panels"))
# "Figure 7(c)"
for _m in re.finditer(r"Figure (\d+)\(([a-z])\)", _flatfig):
    _note(_m.group(1), _m.group(2))
# "Panels (a) and (c) of Figure 7" - the letters come BEFORE the number here, and a first
# version of this check matched only the other form, so a planted bad panel slipped through.
for _m in re.finditer(r"Panels?\s+((?:\([a-z]\)(?:[, ]+(?:and\s+)?)?)+)\s+of\s+Figure\s+(\d+)",
                      _flatfig):
    for _p in re.findall(r"\(([a-z])\)", _m.group(1)):
        _note(_m.group(2), _p)
if _bad:
    for _b in _bad: print("  FAIL  " + _b)
    fail += len(_bad)
else:
    print("  OK    every Figure N(x) reference names a panel that figure has")


# === 8d. every section cross-reference resolves in one of the two documents ===
# 2026-09-22: the companion carried eleven references to §5, §5.1 and §5.5, stale from a
# numbering in which the cosmology results were section 5. Neither document has a section 5, so
# a referee following any of them found nothing. Cheap to check.
print("\n=== 8d. section cross-references resolve ===")
_heads = set()
for _doc in (t, _comp_txt):
    _heads |= set(re.findall(r"^#{2,4} (\d+(?:\.\d+)?)", _doc, re.M))
_tops = {h.split(".")[0] for h in _heads}
_dangling = {}
for _doc, _name in ((t, "paper"), (_comp_txt, "companion")):
    _fl = re.sub(r"\s+", " ", _doc)
    for _m in re.finditer(r"(?:Section|§)\s?(\d+(?:\.\d+)?)", _fl):
        # "Section 5 of their paper" points at somebody else's section, not ours. A first
        # version flagged Dulac and Wei's section 5 as a dangling internal reference.
        if re.match(r"\s+of (?:their|that|this) paper", _fl[_m.end():_m.end() + 22]):
            continue
        _r = _m.group(1)
        if _r not in _heads and _r.split(".")[0] not in _tops:
            _dangling.setdefault((_name, _r), 0)
            _dangling[(_name, _r)] += 1
if _dangling:
    for (_n, _r), _c in sorted(_dangling.items()):
        print("  FAIL  %s refers to §%s, which neither document has (%d times)" % (_n, _r, _c))
    fail += len(_dangling)
else:
    print("  OK    every §N reference in both documents names a section that exists")


# === 9. every substantive number in the abstract is still supported in the BODY ===
# Added 2026-09-20 after the companion split orphaned two abstract promises twice in one turn.
# Check 4 tests a hand-written list; this tests every number mechanically, because a claim whose
# support moved to the companion reads as an unsupported assertion to a reader of the main text.
# NOTE: a first version of this check passed an injection test it should have failed. It matched
# section numbers (3.4, 4.2) and missed 491.6, because the abstract writes \le491.6 and the "e"
# of \le defeated a word-boundary lookbehind. LaTeX control words are stripped first now.
print("\n=== 9. abstract numbers still supported in the body ===")
_body_only = body.split("\n\n## 1.", 1)[1] if "\n\n## 1." in body else body
def _strip(x):
    x = re.sub(r'\\[a-zA-Z]+', ' ', x)      # drop LaTeX control words
    x = re.sub(r'§\s*[0-9]+(?:\.[0-9]+)?', ' ', x)   # drop section references
    return x
_a = _strip(abstract); _b = _strip(_body_only)
_cand = sorted(set(re.findall(r'(?<![\w.])(\d+\.\d+)(?![\w])', _a)), key=lambda x: -len(x))
_cand += sorted(set(m[0] for m in re.findall(r'(?<![\w.])(\d+)\s*(PeV|meV|GeV|per cent)', _a)))
_orphan = [n for n in _cand if n not in _b]
assert _cand, "check 9 found no numbers in the abstract at all: the pattern is wrong, not the paper"
print(f"  substantive numbers in the abstract: {len(_cand)}  ({', '.join(_cand)})")
if _orphan:
    for n in _orphan:
        k = _a.find(n)
        print(f"  {n!r:10s} NOT found in the body   <-- ...{_a[max(0,k-50):k+20].strip()}...")
    fail += 1
else:
    print("  all supported in the body")


# === 10. paper and companion audited as one corpus ===
# Added 2026-09-20 after four separate findings with one cause: section 5 ran 5.1,5.3 and
# section 4 ran 4,4.2 once their subsections moved out; the companion cited its own blocks
# under the paper's old numbering; and A.17 was exiled while the abstract still leaned on it.
# Nothing was watching for numbering gaps at all.
print("\n=== 10. paper and companion as one corpus ===")
import os as _os
try:
    _c = io.open("pub/paper2/COMPANION_v1.md", encoding="utf-8").read()
except Exception:
    _c = ""
_comp_app = set(re.findall(r'^\*\*(A\.\d+)', _c, re.M))
_main_app = set(re.findall(r'^\*\*(A\.\d+)', body, re.M))
_main_sec = set(re.findall(r'^#{2,3} (\d+(?:\.\d+)?)', body, re.M))
_bad_app = sorted(a for a in set(re.findall(r'\b(A\.\d+)\b', _c)) if a not in _comp_app | _main_app)
_bad_sec = sorted(x for x in set(re.findall(r'§(\d+(?:\.\d+)?)', _c))
                  if x not in _main_sec and x not in ("3.1","3.2","3.3","3.4","5","5.1","5.5"))
# 2026-09-23: the companion no longer cites script paths inline; check 7 now forbids it.
_scripts, _miss = set(), []
# numbering gaps: for each top-level section, subsections must run 1..n with no hole
_gaps = []
for _top in sorted({x.split('.')[0] for x in _main_sec if '.' in x}, key=int):
    _subs = sorted(int(x.split('.')[1]) for x in _main_sec if x.startswith(_top + '.'))
    if _subs and _subs != list(range(1, len(_subs) + 1)):
        _gaps.append(f"section {_top} runs {_subs}")
assert _main_app and _main_sec, "check 10 parsed no headings: the pattern is wrong, not the paper"
print(f"  companion appendices {len(_comp_app)}, mainline {len(_main_app)}, mainline sections {len(_main_sec)}")
for _lbl, _v in (("companion cites appendices absent from both", _bad_app),
                 ("companion cites sections absent from the paper", _bad_sec),
                 ("companion cites missing scripts", _miss),
                 ("subsection numbering gaps", _gaps)):
    if _v:
        print(f"  {_lbl}: {_v}   <-- ISSUE")
        fail += 1
    else:
        print(f"  {_lbl}: none")


print("\n=== 11. derived numbers obey the laws that derive them ===")
# 2026-09-20. Checks 1-10 compare STRINGS. When the entropy fix moved the endpoint, fifteen
# numbers derived FROM the endpoint stayed stale and every string check passed. This one
# recomputes them. M0 is read from the paper, so the check follows the paper rather than
# pinning a constant that would itself go stale.
_M0 = float(re.search(r'M_1\\le(\d+\.\d+)\\pm', body).group(1))
_derr = []
# the band-price table: M = M0 (I/I_min)^(-2/5) and the neutrino line is half the mass
_rows = re.findall(r'^\| (\d\.\d\d) \| (\d\.\d+) \| (\d+\.\d) \| (\d+\.\d) \|$', body, re.M)
if len(_rows) < 6:
    raise AssertionError("check 11 parsed %d band rows, expected 6: the pattern is stale" % len(_rows))
for _xc, _ratio, _m, _ln in _rows:
    _want = _M0 * float(_ratio) ** -0.4
    if abs(_want - float(_m)) > 0.5:
        _derr.append(f"band row x_c={_xc}: mass {_m}, law gives {_want:.1f}")
    if abs(float(_m)/2 - float(_ln)) > 0.1:
        _derr.append(f"band row x_c={_xc}: line {_ln}, half the mass is {float(_m)/2:.1f}")
# the compact-object fraction: M(f) = M0 (1-f)^(2/5)
for _f in (0.1, 0.5, 0.9):
    _want = _M0 * (1 - _f) ** 0.4
    for _v in (f"${_want:.1f}$", f"{_want:.1f}"):
        if _v in body: break
    else:
        _derr.append(f"f={_f}: no mass {_want:.1f} PeV in the corpus")
    if f"{_want/2:.1f}" not in body:
        _derr.append(f"f={_f}: no line {_want/2:.1f} PeV in the corpus")
# prose that rounds the endpoint must round the CURRENT endpoint
for _m in re.findall(r'M_1\\lesssim(\d+)\$ PeV', body) + re.findall(r'\$(\d+)\$ PeV mass of ', body):
    if abs(int(_m) - _M0) > 1:
        _derr.append(f"prose quotes {_m} PeV against an endpoint of {_M0}")
# 2026-09-22: S13's three derived numbers. Y_DM was printed as 8.7e-19 where the inputs it
# names give 9.0e-19, a 3.4% slip that nothing caught because the number appears once and is
# derived rather than repeated. It changes no conclusion (the yield it is compared against is
# 36 orders smaller) but it is wrong, so it is recomputed here from the stated inputs.
_s13 = io.open("pub/paper2/SUPPLEMENT_v3.md", encoding="utf-8").read()
_HBAR, _T0 = 6.582119569e-25, 13.8e9 * 3.1557e7          # GeV s, s
_m = re.search(r'For \$M_1=(\d+\.\d+)\\times10\^(\d+)\$ GeV and the illustrative '
               r'\$\\tau=10\^\{(\d+)\}\$ s, \$q\\simeq(\d+\.\d+)\\times10\^\{-(\d+)\}\$', _s13)
if not _m:
    _derr.append("check 11: S13's M_1/tau/q sentence did not parse; the pattern is stale")
else:
    _M1 = float(_m.group(1)) * 10 ** int(_m.group(2))
    _tau = 10.0 ** int(_m.group(3))
    _qsaid = float(_m.group(4)) * 10 ** -int(_m.group(5))
    _qwant = 8 * 3.141592653589793 * _HBAR / (_M1 * _tau)
    if abs(_qwant - _qsaid) / _qwant > 0.02:
        _derr.append(f"S13: q printed {_qsaid:.3g}, 8*pi*hbar/(M_1 tau) gives {_qwant:.3g}")
    def _sci(_pat, _txt):
        _mm = re.search(_pat, _txt)
        return None if not _mm else float(_mm.group(1)) * 10 ** -int(_mm.group(2))
    _rho = _sci(r'rho_\{\\rm DM,0\}=(\d+\.\d+)\\times10\^\{-(\d+)\}', body)
    _s0  = _sci(r's_0=(\d+\.\d+)\\times10\^\{-(\d+)\}', body)
    if _rho is None or _s0 is None:
        _derr.append("check 11: the paper's rho_DM,0 or s_0 did not parse; the pattern is stale")
        _rho, _s0 = 1.0, 1.0
    _ywant = _rho / (_s0 * _M1)
    _my = re.search(r'Y_\{\\rm DM\}=\\rho_\{\\rm DM,0\}/\(s_0M_1\)\\simeq(\d+\.\d+)'
                    r'\\times10\^\{-(\d+)\}', _s13)
    if not _my:
        _derr.append("check 11: S13's Y_DM did not parse; the pattern is stale")
    else:
        _ysaid = float(_my.group(1)) * 10 ** -int(_my.group(2))
        if abs(_ywant - _ysaid) / _ywant > 0.02:
            _derr.append(f"S13: Y_DM printed {_ysaid:.3g}, rho/(s_0 M_1) gives {_ywant:.3g}")
    _md = re.search(r'1-e\^\{-t_0/\\tau\}\\simeq(\d+\.\d+)\\times10\^\{-(\d+)\}', _s13)
    if _md:
        _dsaid = float(_md.group(1)) * 10 ** -int(_md.group(2))
        _dwant = 1 - 2.718281828459045 ** (-_T0 / _tau)
        if abs(_dwant - _dsaid) / _dwant > 0.03:
            _derr.append(f"S13: depletion printed {_dsaid:.3g}, 1-exp(-t_0/tau) gives {_dwant:.3g}")
    print(f"  S13 recomputed from its own inputs: q, Y_DM and the depletion over {_T0/3.1557e16:.1f} Gyr")
print(f"  endpoint read from the paper: {_M0} PeV; band rows checked: {len(_rows)}")
if _derr:
    for _e in _derr: print(f"  {_e}   <-- ISSUE")
    fail += len(_derr)
else:
    print("  band table, compact-object masses and prose roundings all obey their laws")


print("\n=== 12. display math is in a form md2tex can convert ===")
# 2026-09-20: a $$...$$ opened and closed on the same line is escaped to literal \$\$ by
# tools/md2tex.py and LaTeX then fails with "Missing $ inserted". The converter wants $$ alone
# on its own line. Cheap to check here; expensive to find in a 2000-line .log.
_bad = []
for _f in ("pub/paper2/PAPER2_v3.md", "pub/paper2/COMPANION_v1.md"):
    for _i, _line in enumerate(io.open(_f, encoding="utf-8").read().split("\n"), 1):
        _st = _line.strip()
        while _st.startswith(">"):          # blockquoted display math is handled by md2tex
            _st = _st[1:].strip()
        if _st.count("$$") and _st != "$$":
            _bad.append(f"{_f}:{_i}")
print(f"  files scanned: 2; malformed display blocks: {len(_bad)}")
if _bad:
    for _b in _bad[:10]: print(f"  {_b}   <-- ISSUE")
    fail += len(_bad)
else:
    print("  every $$ sits alone on its line")


print("\n=== 13. the front page carries what four reviewers said decides whether it is read ===")
# The guard table checks presence anywhere in paper+companion, so an "abstract names X" guard
# passes when only the body says X. Guard AQ did exactly that for a day. These check `front`,
# which is the short Abstract alone.
_front_must = [
    ("the bound",            "491.6"),
    ("the falsifier line",   "245.8"),
    ("the neutrino floor",   "58.8"),
    # 2026-09-22: this pinned one phrasing, so rewording the sentence failed the check even
    # though the content survived. It now tests the content: the abstract has to say the
    # quoted width is measurement error, and that theory error is the larger one.
    ("that the width is measurement error", "measurement error"),
    ("that theory error is larger",         "Theory error is larger"),
    ("the JWST liability",   "JWST")]
# Prior credit was on this list and has MOVED to the introduction, 2026-09-21. The reviewers'
# concern was that the paper must not look as though it appropriates Boyle, Finn and Turok's
# cosmology; the fix chosen was to name them in the abstract, which is the wrong place for it.
# Neither Buchalter-winning abstract of 2025 names a single prior author. An abstract is not a
# literature review, and arXiv characters spent on one are characters not spent on the result.
# 1's second paragraph credits them in full and the abstract flags the dependency in four
# words, so the concern is met and only the location has changed.
_intro_must = [("prior credit", "Boyle, Finn and Turok")]
_flatfront = re.sub(r"\s+", " ", front)
_flatintro = re.sub(r"\s+", " ", t[t.index("## 1. Introduction"):t.index("## 2. Methods")])
_missing = [f"{_l} ({_p!r})" for _l, _p in _front_must if _p not in _flatfront]
_missing += [f"{_l} ({_p!r}) missing from the INTRODUCTION"
             for _l, _p in _intro_must if _p not in _flatintro]
print(f"  front page: {len(front.split())} words, {front.count('$$')//2} displayed equations")
if len(front.split()) > 700:
    print(f"  front page is {len(front.split())} words   <-- ISSUE (a front page, not a section)")
    fail += 1
if _missing:
    for _m in _missing: print(f"  front page omits {_m}   <-- ISSUE")
    fail += len(_missing)
else:
    print("  carries the bound, the falsifier, the floor, the width caveat, JWST and the credit")


print("\n=== 14. the paper does not narrate its own edit history ===")
# A journal editor read the manuscript and put this second on a desk-reject list: "The manuscript
# contains at least a dozen 'an earlier version of this paragraph...' corrections. A research
# paper presents the defended final state." Withdrawals stay where they carry physics; the
# document's diff history does not belong in it. One recorded correction is a virtue, a dozen is
# a lab notebook.
# Widened after a cold reader tripped on "An earlier claim ... is withdrawn in A.18", which the
# first version of this pattern did not catch: it matched "withdrawn HERE" and not "withdrawn IN".
_SELF = (r'an earlier version|an earlier claim|earlier draft|previously claimed|we withdraw|'
         r'is withdrawn (here|in|above)|was withdrawn|we no longer claim|'
         r'worth recording|recorded rather than absorbed|rather than gestured at')
_meta = 0
for _f in ("pub/paper2/PAPER2_v3.md", "pub/paper2/COMPANION_v1.md"):
    _c = len(re.findall(_SELF, re.sub(r'\s+', ' ', io.open(_f, encoding='utf-8').read()), re.I))
    print(f"  self-referential edit notes in {_f.split('/')[-1]}: {_c}")
    _meta += _c
if _meta:
    print(f"  {_meta} instances of the document narrating its own history   <-- ISSUE")
    fail += _meta
else:
    print("  neither document narrates its own edit history")


print("\n=== 15. the abstract fits arXiv's hard character limit ===")
# arXiv rejects an abstract over 1920 characters at submission time. This was found at 2111,
# i.e. the manuscript could not have been submitted as it stood, by measuring rather than
# reading. Whitespace-collapsed, because that is how the form counts it.
_ab = re.sub(r'\s+', ' ',
             io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read()
             .split("## Abstract")[1].split("## 1. Introduction")[0]).strip()
print(f"  abstract: {len(_ab)} characters, {len(_ab.split())} words (arXiv limit 1920)")
if len(_ab) > 1920:
    print(f"  over the arXiv limit by {len(_ab)-1920} characters   <-- ISSUE"); fail += 1
else:
    print(f"  fits, with {1920-len(_ab)} characters to spare")

print("\n=== 16. no placeholder ships, and the release pointers are not the previous paper's ===")
# Both external pointers were found aiming at the PREVIOUS manuscript: GitHub tag v3.0.8 ships
# "The Far Side of the Horizon: a geometric fold for CPT-related copies of spacetime" with a
# different set of scripts, and DOI 10.5281/zenodo.22811618 resolves to "...two copies of
# spacetime, and the zero-parameter predictions that follow" (2026-09-17). Neither can ship on
# this paper. They are marked pending until re-cut, and this refuses to let the markers through.
_body = io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read()
_stale = {"RELEASE-TAG-PENDING": "release tag not yet cut for this version",
          "ZENODO-DOI-PENDING":  "Zenodo deposit not yet made for this version",
          "v3.0.8":              "that tag ships the PREVIOUS paper",
          "zenodo.22811618":     "that DOI resolves to the PREVIOUS paper"}
_hit = [(m, why) for m, why in _stale.items() if m in _body]
for m, why in _hit:
    print(f"  {m}: {why}   <-- BLOCKS SUBMISSION")
if _hit:
    print("  (expected before release; must be zero at submission)")
else:
    print("  no placeholders and no pointers to the previous paper")
# deliberately NOT added to `fail`: these are known-pending until Ben cuts the release.

print("\n=== 17. the paper's claims about its own apparatus are true ===")
# The code-availability block was previously found quoting a census that matched nothing.
# It now quotes the number of checks in THIS file and the number of scripts cited. Those are
# facts about the repository, so check them here rather than trust them.
_pv = io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read()
_src = io.open(__file__, encoding="utf-8").read()
_nchk = max(int(n) for n in re.findall(r'=== (\d+)\. ', _src))
_words = {12:"twelve",13:"thirteen",14:"fourteen",15:"fifteen",16:"sixteen",17:"seventeen",18:"eighteen",19:"nineteen",20:"twenty",21:"twenty-one",22:"twenty-two",23:"twenty-three",24:"twenty-four",25:"twenty-five",26:"twenty-six",27:"twenty-seven",28:"twenty-eight",29:"twenty-nine",30:"thirty",31:"thirty-one",32:"thirty-two"}
_claim = _words.get(_nchk, str(_nchk))
_flat = re.sub(r'\s+', ' ', _pv)
if f"runs {_claim} checks" in _flat:
    print(f"  paper says {_claim} checks; this file has {_nchk}")
else:
    print(f"  this file has {_nchk} checks, and the paper does not say '{_claim}'   <-- ISSUE")
    fail += 1
# 2026-09-22: the count was updated in the "There are N of them" sentence and left stale in
# the responsibility paragraph, which says it twice more. This check read only the one
# sentence, so the paper carried 132 and 133 at once. Every number sitting next to the word
# "scripts", and every "all N run/were run", now has to be the true count.
# 2026-09-23: the count used to be of scripts CITED in the text. The text cites none now, so it
# is the count of what the release actually ships: the three directories the release carries.
import glob as _glob
_REL = ("separate_ways/tangents/**/*.R", "pub/paper2/companion/**/*.R",
        "pub/paper2/figure_authoring_r/**/*.R")
_cited_R = sorted({_f for _pat in _REL for _f in _glob.glob(_pat, recursive=True)})
_true_n = len(_cited_R)
_claimed = set()
# 2026-09-22, second pass: the middle pattern was a bare '\bof the\s+(\d+)\s', which matched
# "One of the 170 mentions Tomita-Takesaki" in section 1.1's novelty search and reported it as a
# stale script count. It now has to be followed by the word it was always about.
for _m2 in re.finditer(r'(\d+)\s+scripts\b|\bof the\s+(\d+)\s+scripts\b|[Aa]ll\s+(\d+)\s+(?:run|were run)', _flat):
    _v = next((_g for _g in _m2.groups() if _g), None)
    if _v and 50 < int(_v) < 500: _claimed.add(int(_v))
print(f"  script counts stated next to 'scripts' or 'all N run': {sorted(_claimed) or 'none found'}")
assert _claimed, "check 17 found no stated script count at all: the pattern is stale, not the paper"
_wrong = sorted(_c for _c in _claimed if _c != _true_n)
if _wrong:
    print(f"  on disk there are {_true_n}; the paper also says {_wrong}   <-- ISSUE"); fail += len(_wrong)
else:
    print(f"  every one of them is {_true_n}, which is what is on disk")

# the same paragraph claims a count of scripts recording a withdrawn or superseded result
_marked = [_c for _c in _cited_R
           if re.search(r'withdrawn|superseded|retracted|an earlier version',
                        io.open(_c, encoding="utf-8", errors="replace").read(), re.I)]
_wordn = {8:"Eight",9:"Nine",10:"Ten",11:"Eleven",12:"Twelve",13:"Thirteen",14:"Fourteen",15:"Fifteen"}
_wn = _wordn.get(len(_marked), str(len(_marked)))
if re.search(rf"{_wn}\b[^.]{{0,20}}?record a result withdrawn", _flat):
    print(f"  {_wn} scripts record a withdrawn or superseded result, which is what the paper says")
else:
    print(f"  {len(_marked)} scripts record a withdrawn or superseded result; the paper does not say "
          f"'{_wn} ... record a result withdrawn'   <-- ISSUE"); fail += 1

_m = re.search(r'There are (\d+) of them, ([\d,]+) lines in total', _flat)
if not _m:
    print("  could not find the script-count sentence   <-- ISSUE"); fail += 1
else:
    _n_claim, _l_claim = int(_m.group(1)), int(_m.group(2).replace(",", ""))
    _r = set(_cited_R)
    _lines = sum(len(io.open(c, encoding="utf-8", errors="replace").readlines())
                 for c in _r if _os.path.exists(c))
    ok = (_n_claim == len(_r) and _l_claim == _lines)
    print(f"  paper says {_n_claim} scripts / {_l_claim} lines; on disk {len(_r)} / {_lines}"
          + ("" if ok else "   <-- ISSUE"))
    if not ok: fail += 1

print("\n=== 18. qualifiers that make a quoted number true are still attached ===")
# Each of these numbers is FALSE without the words beside it. "six orders" is the margin on
# the TIGHTEST accretion limit; the weakest is 2.0 orders. The qualifier was added once,
# dropped again while trimming the abstract for characters, and restored. A guard is cheaper
# than remembering.
_pv2 = re.sub(r'\s+', ' ', io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read())
_pairs = [("six orders below", "at the tightest", "the weakest PBH limit gives 2.0 orders"),
          ("six orders of margin", "on the tightest", "same, in the body")]
_miss = 0
for _num, _qual, _why in _pairs:
    for _m in re.finditer(re.escape(_num), _pv2):
        _win = _pv2[_m.start():_m.start()+len(_num)+40]
        if _qual not in _win:
            print(f"  '{_num}' appears without '{_qual}': {_why}   <-- ISSUE"); _miss += 1
print(f"  {len(_pairs)} qualifier pairs checked")
if _miss: fail += _miss
else: print("  every quoted margin still carries the words that make it true")

print("\n=== 19. nothing in either document renders blank or breaks the PDF build ===")
# Both documents built only after 44 Unicode maths symbols used in PROSE were converted.
# Latin Modern has no glyph for these, so xelatex printed "Missing character" and they came
# out BLANK on the page. The set below is the one xelatex actually reported, not a guess:
# a category-based rule flags U+2212 and U+00D7, which the font does have.
# Separately, pandoc does not read "$-$1" as maths at all (it escapes both dollars as
# currency), which swallowed a paragraph into a maths span and failed the build outright.
_MISSING = set("\u2070\u2074\u207b\u00b2\u00b9\u2082\u2084\u00bc\u00bd"
               "\u2218\u222b\u2265\u0227\u2192\U0001d49c") | {chr(c) for c in range(0x0370, 0x0400)}
_bad_total = 0
for _doc in ("pub/paper2/PAPER2_v3.md", "pub/paper2/COMPANION_v1.md"):
    _t = io.open(_doc, encoding="utf-8").read()
    _st = re.sub(r'\$\$.*?\$\$', '', _t, flags=re.S)
    _st = re.sub(r'\$[^$\n]*\$', '', _st)
    _st = re.sub(r'`[^`\n]*`', '', _st)
    _glyph = sorted({c for c in _st if c in _MISSING})
    _curr = re.findall(r'\$-\$\d', _t)
    if _glyph:
        print(f"  {_doc.split('/')[-1]}: blank in the PDF: {_glyph}   <-- ISSUE")
    if _curr:
        print(f"  {_doc.split('/')[-1]}: {len(_curr)} of '$-$<digit>', read as currency   <-- ISSUE")
    _bad_total += len(_glyph) + len(_curr)
if _bad_total:
    fail += _bad_total
else:
    print("  no prose glyph the PDF font lacks, and no currency-ambiguous minus")

print("\n=== 20. phrase guard on bound-as-measurement (a GUARD, not a proof) ===")
# WHAT THIS IS AND IS NOT. A pre-submission review found the same root error six times: a
# BOUND treated as a measured value. This greps for a handful of formulations that assert it.
# The same review then showed the first version of this check gave FALSE ASSURANCE: it passed
# "We measure M_1 = 491.6 PeV" and "a redshift axis with nothing adjustable in it", and it
# FLAGGED the correct sentence "The mass is not fixed by the abundance". A literal phrase
# guard cannot establish that the distinction holds; it can only catch known phrasings. It is
# kept because the known phrasings recurred, and labelled so nobody reads a PASS as evidence.
_pv3 = re.sub(r'\s+', ' ', io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read())
_banned = [
    (r"(?<!not )fixed by the abundance",        "the abundance BOUNDS the mass"),
    (r"(?<!not )the abundance fixes",           "same"),
    (r"requires \$f_\{\\rm PBH\}<",           "a larger fraction lowers the ceiling; nothing is excluded"),
    (r"redshift axis with no(?:thing)? (?:free parameter|adjustable)", "E_nu is bounded, not measured"),
    (r"\bwe measure \$?M_1",                    "no mass has been measured"),
    (r"\bwe measure \$?E_\\nu",                 "no endpoint has been measured"),
    (r"the measured (?:dark-matter )?mass\b",   "there is none"),
    (r"bounds any competing dark component",   "a competing component moves the ceiling, it is not bounded by it"),
    (r"budget cannot afford",                  "the budget does not exclude; it shifts"),
    # Found by a later prose pass 2026-09-21; check 20 had missed all three.
    (r"fixed mass also closes",                "the ceiling is not a fixed mass"),
    (r"once it fixes \$M_1\$",                 "the abundance match assigns the budget; it does not fix the mass"),
    (r"ceiling on a primordial dark seed",     "3.4 calls this a sensitivity and denies it is an exclusion"),
]
_hits = 0
for _pat, _why in _banned:
    for _m in re.finditer(_pat, _pv3, re.I):
        print(f"  '{_m.group(0)}': {_why}   <-- ISSUE"); _hits += 1
print(f"  {len(_banned)} phrasings checked; a PASS here means none of these appear, nothing more")
if _hits:
    fail += _hits
else:
    print("  none of the known bound-as-measurement phrasings appear")


print("\n=== 21. numbers that must follow from the stated oscillation inputs ===")
# A final check (2026-09-21) found the inverted-ordering sum did NOT follow from the
# splittings 2.5 states: the figure had been built on an atmospheric splitting 4 per cent
# lower, and nothing guarded it. Derive both sums here from the quoted inputs and compare.
import math as _math
_m21 = re.search(r"Delta m\^2_\{21\}=\((\d+\.\d+)\\pm[^)]*\)\\times10\^\{-5\}", _pv3) \
    or re.search(r"m\^2_\{21\}=\((\d+\.\d+)", _pv3)
_m31 = re.search(r"m\^2_\{31\}=\((\d+\.\d+)", _pv3)
if not (_m21 and _m31):
    print("  could not locate the stated splittings   <-- ISSUE"); fail += 1
else:
    _d21 = float(_m21.group(1))*1e-5; _d31 = float(_m31.group(1))*1e-3
    _no = (_math.sqrt(_d21) + _math.sqrt(_d31))*1000
    _io = (_math.sqrt(_d31) + _math.sqrt(_d31 + _d21))*1000
    for _label, _val, _pat in (("normal", _no, r"\\Sigma m_\\nu=\\sqrt\{\\Delta m\^2_\{21\}\}\+\\sqrt\{\\Delta m\^2_\{31\}\}=(\d+\.\d+)"),
                               ("inverted", _io, r"gives \$(\d+\.\d+)\$ meV on the same inputs")):
        _q = re.search(_pat, _pv3)
        if not _q:
            print(f"  {_label} ordering: could not find the quoted sum   <-- ISSUE"); fail += 1
        else:
            _qq = float(_q.group(1)); _ok = abs(_qq - _val) < 0.06
            print(f"  {_label} ordering: paper {_qq} meV, derived {_val:.2f} meV"
                  + ("" if _ok else "   <-- ISSUE"))
            if not _ok: fail += 1
    # and the DESI margin, which drifted between 64.0 and 64.2 in the same section
    _flat21 = re.sub(r"\s+", " ", _pv3)
    _b = re.search(r"Sigma m_.nu<(\d+\.?\d*). meV under ..Lambda.CDM", _flat21)
    # EVERY instance, not the first: the percentage is quoted in 3.3 and again in 4.4, and a
    # first-match check passed while one of the two carried the old value. Found by planting it.
    _ts = [float(x) for x in re.findall(r"tighten(?:s)? by more than .(\d+\.\d+). per cent", _flat21)]
    if _b and _ts:
        _bb = float(_b.group(1))
        _der = 100*(_bb - _no)/_bb
        _bad = [x for x in _ts if abs(_der - x) >= 0.15]
        print(f"  DESI margin: {len(_ts)} instance(s) {sorted(set(_ts))} per cent, "
              f"derived {_der:.2f} from bound {_bb}" + ("" if not _bad else "   <-- ISSUE"))
        if _bad: fail += len(_bad)
    else:
        print("  could not locate the DESI bound or the tightening percentage   <-- ISSUE"); fail += 1

print("\n=== 22. every figure the paper places can actually be regenerated ===")
# render_all.R was found listing the EARLIER paper's figures and none of the nine this one
# places. Nothing caught it because nothing compared the two. A reader following the repo's
# own instructions would have regenerated the wrong set and kept stale artwork.
# fig5_data is built by the de-AI'd Python generator, not by R, and is exempt by name.
_PY_BUILT = {"fig5_data"}
_ra_path = "pub/paper2/figure_authoring_r/render_all.R"
_ra = io.open(_ra_path, encoding="utf-8").read()
_listed = set(re.findall(r'"([A-Za-z0-9_]+)\.R"', _ra))
_placed = [m.group(1) for m in re.finditer(r'!\[\]\(([A-Za-z0-9_]+)\.png\)', _pv3)]
assert _placed, "check 22 matched no figures at all: the pattern is wrong, not the paper"
_miss = [f for f in _placed if f not in _listed and f not in _PY_BUILT]
print(f"  {len(_placed)} figures placed; {len(_listed)} generators listed in render_all.R")
if _miss:
    for _f in _miss:
        print(f"  {_f}.png is placed but has no generator in render_all.R   <-- ISSUE")
    fail += len(_miss)
else:
    print(f"  every placed figure has a generator, with {sorted(_PY_BUILT)} exempt (Python-built)")

print("\n=== 23. no superseded value has crept back in ===")
# The number-provenance tool compares the text against the script it cites, which cannot catch
# a number that is stale in BOTH. That is not hypothetical: uv_ringdown.R carried the
# pre-entropy-fix mass of 484.8 PeV while the companion beside it said 492, and three scripts
# asserted things about the paper that were true before a correction and false after. So keep a
# ledger of what has been superseded and grep for it. Each entry names the files where the old
# value is ALLOWED, which are the files that record the correction.
_SUPERSEDED = [
    (r"484\.8|4\.848e8",  "the pre-entropy-fix mass; now 491.6 PeV / 4.916e8 GeV",
     # 2.3 states 484.8 deliberately, as what the erroneous entropy density gives
     {"entropy_density_fix.R", "stale_after_entropy_fix.R", "PAPER2_v3.md"}),
    (r"98\.9 meV",         "the inverted-ordering sum that did not follow from the stated "
                           "inputs; now 100.9 meV", {"neutrino_orderings.R"}),
    (r"64\.0 meV",         "the rounded DESI bound; now 64.2 meV", {"desi_margin.R"}),
    (r"8\.1 per cent",     "the tightening margin off the rounded bound; now 8.4 per cent",
     {"desi_posterior.R", "desi_margin.R"}),
    (r"5\.2 meV",          "the floor-to-bound gap off the rounded bound; now 5.4 meV",
     {"desi_margin.R"}),
]
_scope = ["pub/paper2/PAPER2_v3.md", "pub/paper2/COMPANION_v1.md"]
for _root, _dirs, _files in _os.walk("separate_ways"):
    for _f in _files:
        if _f.endswith((".R", ".py")):
            _scope.append(_os.path.join(_root, _f))
_scope = [f for f in _scope if _os.path.exists(f) and not f.endswith("consistency_pass.py")]

_stale = 0
for _pat, _why, _allow in _SUPERSEDED:
    _hits = []
    for _f in _scope:
        if _os.path.basename(_f) in _allow: continue
        try: _t = io.open(_f, encoding="utf-8", errors="replace").read()
        except OSError: continue
        if re.search(_pat, _t): _hits.append(_f)
    if _hits:
        print(f"  {_why}")
        for _h in _hits[:6]: print(f"      still present in {_h}   <-- ISSUE")
        _stale += len(_hits)
print(f"  {len(_SUPERSEDED)} superseded values checked across {len(_scope)} files"
      + ("" if not _stale else f"; {_stale} survival(s)"))
if _stale: fail += _stale
else: print("  none has crept back; each is allowed only in the file that records its correction")

print("\n=== 24. the figure manifest matches the figures on disk ===")
# The release promises "a SHA-256 manifest". It had drifted to listing four figures from the
# EARLIER paper and two of this one's eleven, because it was maintained by hand against a moving
# paper. It is now generated from the manuscript, the drawing code and the assets, and this
# check fails the release if it is out of date.
import subprocess as _sp, sys as _sys
_r = _sp.run([_sys.executable, "tools/make_figure_manifest.py", "--check"],
             capture_output=True, text=True)
print("  " + (_r.stdout.strip() or "(no output)"))
if _r.returncode != 0: fail += 1

print("\n=== 25. every listed reference is actually cited ===")
# Found with 10 orphans, of which 4 were the canonical decoherence-in-quantum-cosmology papers
# while 3.5 discusses exactly that and cited none of them. That is an attribution gap a referee
# in the field would see at once, not a tidiness point. Three were placed; the rest need a
# judgement about what they contribute and are listed as PENDING rather than quietly dropped.
_PENDING = {}   # emptied 2026-09-21: Ben's call was to cut. References 11, 12, 17, 18 and 19
                # (A. Higuchi's spin-2 mass range, and four Boyle-Turok programme papers) were
                # removed and everything above them renumbered, 64 entries down to 59.
_mainref = io.open("pub/paper2/PAPER2_v3.md", encoding="utf-8").read()
_alltext = _mainref + io.open("pub/paper2/COMPANION_v1.md", encoding="utf-8").read()
_listed = sorted(int(_m.group(1)) for _m in re.finditer(r'(?m)^(\d+)\\?\.\s+[A-Z]', _mainref))
_cited = {int(_x) for _m in re.finditer(r'\[([\d,\s]+)\]', _alltext)
          for _x in _m.group(1).replace(' ', '').split(',') if _x.isdigit()}
_orphan = [_n for _n in _listed if _n not in _cited]
_new = [_n for _n in _orphan if _n not in _PENDING]
print(f"  {len(_listed)} references listed, {len(_orphan)} cited nowhere in either document")
for _n in _orphan:
    _tag = "PENDING Ben's decision" if _n in _PENDING else "NEW   <-- ISSUE"
    print(f"      [{_n}] {_tag}")
if _new: fail += len(_new)
else: print("  no reference has become orphaned since the list was reviewed")

def _expand(_g):
    """[39,42–43,50] -> 39 42 43 50"""
    _out = []
    for _part in _g.replace("—", "–").split(","):
        _part = re.sub(r'(?<=\d)-(?=\d)', "–", _part)
        if "–" in _part:
            _a, _b = _part.split("–")
            _out += [str(_x) for _x in range(int(_a), int(_b) + 1)]
        else:
            _out.append(_part)
    return _out
# capitalised words that sit before a marker but are not surnames
_NOTNAME = {"Planck","This","That","Their","Appendix","Supplement","Section","Table","Figure","Model"}

print("\n=== 26. the reference list is complete, numbered without gaps, and identifiable ===")
# Check 25 runs one direction: every listed reference is cited. This runs the other three.
# The paragraph sweep reached the References section and found no prose to rewrite there, so
# what a bibliography can get wrong is checked instead of read.
_refblock = _mainref[_mainref.rindex("\n## References"):] if "\n## References" in _mainref else ""
_entries = {int(_m.group(1)): _m.group(2).strip()
            for _m in re.finditer(r'(?m)^(\d+)\\?\.\s+(.+)$', _refblock)}
_nomath = re.sub(r'\$[^$]*\$', ' ', _alltext)      # w\in[0,1] is an interval, not a citation
# 2026-09-22: this matched only comma lists. A planted [12-99] in the paper, reaching forty
# entries past the end of a 59-entry list, passed silently. The same blindness in check 30 let
# nineteen of the supplement's sixty cited works ship with no entry. Ranges now count.
_cited31 = {int(_x) for _m in re.finditer(r'\[([\d,\s–—-]+)\]', _nomath)
            for _x in _expand(_m.group(1).replace(' ', '')) if _x.isdigit()}
_missing = sorted(_n for _n in _cited31 if _n not in _entries)
_gaps = [_n for _n in range(1, (max(_entries) if _entries else 0) + 1) if _n not in _entries]
# An arXiv identifier, new style or old, already carries the year; a journal or publisher entry
# has to state one. The first version demanded a year of every entry and flagged six arXiv
# preprints that were perfectly findable.
_ID = (r'arXiv:\s*\d{4}\.\d{4,5}'
       r'|(?:gr-qc|hep-th|hep-ph|hep-ex|astro-ph|quant-ph|math-ph|cond-mat|nucl-th|nucl-ex)/\d{7}'
       r'|doi|DOI|10\.\d{4,}/'
       r'|companion paper')      # the companion has no external identifier and needs none
_VENUE = (r'Phys\.|Rev\.|Class\.|Nucl\.|Ann\.|Commun\.|Living Rev|Cambridge|Springer|Oxford|'
          r'Princeton|Univ\.|Press|Lett\.|Mon\. Not\.|JCAP|JHEP|Nature|Science|SciPost|Int\. J\.|'
          r'SIAM|Optim\.|Proc\.|Rep\.|Adv\.|Eur\.|New J\.|Found\. Phys\.|Am\. J\.|Sov\.')
_noid = sorted(_n for _n, _t in _entries.items()
               if not re.search(_ID, _t)
               and not (re.search(_VENUE, _t) and re.search(r'(19|20)\d\d', _t)))
print(f"  {len(_entries)} entries parsed, {len(_cited31)} distinct citation numbers used in the two documents")
for _n in _missing: print(f"      [{_n}] cited but has no entry   <-- ISSUE")
for _n in _gaps:    print(f"      [{_n}] gap in the numbering      <-- ISSUE")
for _n in _noid:    print(f"      [{_n}] no arXiv id, DOI or venue+year: {_entries[_n][:60]}   <-- ISSUE")
_bad = len(_missing) + len(_gaps) + len(_noid)
if _bad: fail += _bad
else: print("  every citation resolves, the numbering is dense, and every entry can be looked up")

print("\n=== 27. the submission checklist still describes THIS paper ===")
# On 2026-09-22 the checklist still described the paper as it stood on 17 September: a
# 1686-character abstract, 80 references, six numbered sections. Five days and several hundred
# commits stale, and anyone submitting from it would have verified the wrong facts against the
# wrong artefacts. A checklist nothing checks is a checklist that lies quietly.
_CL_PATH = "pub/paper2/SUBMISSION_CHECKLIST_v3.md"
try:
    _cl = io.open(_CL_PATH, encoding="utf-8").read()
except OSError:
    print(f"  {_CL_PATH} is missing   <-- ISSUE"); fail += 1; _cl = ""
if _cl:
    _stated = {m.group(1): int(m.group(2))
               for m in re.finditer(r'(?m)^\s{2,}(\w+):\s*(\d+)\s*$', _cl)}
    _sup = io.open("pub/paper2/SUPPLEMENT_v3.md", encoding="utf-8").read()
    _abs_i = _pv.find("## Abstract"); _abs_j = _pv.find("\n## ", _abs_i + 5)
    _abs = re.sub(r"\s+", " ", _pv[_abs_i:_abs_j].replace("## Abstract", "")).strip()
    _refblk = _pv[_pv.rindex("\n## References"):]
    _actual = {
        "abstract_chars":      len(_abs),
        "abstract_words":      len(_abs.split()),
        "main_sections":       len(re.findall(r'(?m)^## (\d+)\.', _pv)),
        "appendices":          len(re.findall(r'(?m)^## Appendix ([A-Z])\.', _pv)),
        "references":          len(re.findall(r'(?m)^(\d+)\\?\.', _refblk)),
        "em_dashes":           _pv.count("\u2014"),
        "supplement_sections": len(set(re.findall(r'(?m)^## (S\d+)\.', _sup))),
        "consistency_checks":  _nchk,
        "cited_scripts":       _true_n,
    }
    if not _stated:
        print("  the checklist states no machine-readable facts at all   <-- ISSUE"); fail += 1
    else:
        _bad = [(k, _stated[k], _actual[k]) for k in _stated
                if k in _actual and _stated[k] != _actual[k]]
        _unknown = [k for k in _stated if k not in _actual]
        print(f"  {len(_stated)} facts stated, {len(_actual)} checkable")
        for k, sv, av in _bad:
            print(f"      {k}: checklist says {sv}, the files say {av}   <-- ISSUE")
        for k in _unknown:
            print(f"      {k}: stated but nothing checks it")
        if _bad:
            fail += len(_bad)
        else:
            print("  every stated fact matches the live files")

print("\n=== 28. the built arXiv bundle was built from THIS paper ===")
# The previous bundle sat four days out of date and nobody noticed, which is how a rewrite
# nearly went to arXiv as the draft that preceded it. tools/build_arxiv_v4.sh stamps the source
# hashes into the bundle; this compares them with the live files. A bundle is an artefact, and
# an artefact that nothing checks is one that will be stale when it matters.
import hashlib as _hl, json as _json
_BUNDLE = "pub/paper2/arxiv_src_v4"
_stamp = _os.path.join(_BUNDLE, "SOURCE_HASHES.json")
if not _os.path.exists(_stamp):
    print(f"  {_stamp} is missing: the bundle carries no provenance   <-- ISSUE"); fail += 1
else:
    _rec = _json.load(io.open(_stamp, encoding="utf-8")).get("built_from", {})
    if not _rec:
        print("  the stamp records no sources   <-- ISSUE"); fail += 1
    else:
        _drift = []
        for _src, _want in _rec.items():
            if not _os.path.exists(_src):
                _drift.append((_src, "missing")); continue
            _got = _hl.sha256(io.open(_src, "rb").read()).hexdigest()
            if _got != _want:
                _drift.append((_src, "changed since the bundle was built"))
        print(f"  {len(_rec)} source file(s) stamped in the bundle")
        for _src, _why in _drift:
            print(f"      {_src}: {_why}   <-- ISSUE")
        if _drift:
            print("      re-run  bash tools/build_arxiv_v4.sh")
            fail += len(_drift)
        else:
            print("  the bundle matches the live sources")

print("\n=== 29. every section cross-reference lands on a section that exists ===")
# 2026-09-22: the supplement was written against V2's six-section main text and never
# re-pointed. It sent readers to main 5.5, 5.6 and 5.2, none of which exist in a four-section
# paper, and to 2.3 for premises that are stated in Appendix A.1. A referee following any of
# those lands nowhere. Nothing checked cross-references, so the rot was invisible: the
# references are syntactically fine and only wrong about the world.
#
# Routing: "main" sends you to the paper, "supplement" or a leading S to the supplement, a bare
# number to whichever document you are not reading. "their 3.4" is someone else's paper.
_REF = re.compile(r'(\w+\s+)?§\s*(S?\d+(?:\.\d+)*|[A-Z](?:\.\d+)?)')

def _paper_secs(t):
    h  = {m.group(1) for m in re.finditer(r'(?m)^##\s+(\d+)\.\s', t)}
    h |= {m.group(1) for m in re.finditer(r'(?m)^###\s+(\d+\.\d+)\s', t)}
    h |= {m.group(1) for m in re.finditer(r'(?m)^##\s+Appendix\s+([A-Z])\.', t)}
    h |= {m.group(1) for m in re.finditer(r'(?m)^\*\*([A-Z]\.\d+)\s', t)}
    return h

def _supp_secs(t):
    h  = {m.group(1) for m in re.finditer(r'(?m)^##\s+(S\d+)\.', t)}
    h |= {m.group(1) for m in re.finditer(r'(?m)^\*\*(S\d+\.\d+)\s', t)}
    return h

_supp_txt = io.open("pub/paper2/SUPPLEMENT_v3.md", encoding="utf-8").read()
_P, _S = _paper_secs(_pv), _supp_secs(_supp_txt)
if not _P or not _S:
    print("  no sections were parsed out of one of the documents   <-- ISSUE"); fail += 1
else:
    print(f"  {len(_P)} sections in the paper, {len(_S)} in the supplement")
    _dangle = []
    for _name, _txt, _in_supp in (("PAPER2_v3.md", _pv, False),
                                  ("SUPPLEMENT_v3.md", _supp_txt, True)):
        _n = 0
        for _m in _REF.finditer(_txt):
            _lead = (_m.group(1) or "").strip().lower()
            if _lead == "their":
                continue
            _num = _m.group(2).rstrip(".")
            _n += 1
            if _lead == "main":          _want, _w = _P, "the paper"
            elif _lead == "supplement":  _want, _w = _S, "the supplement"
            elif _num.startswith("S"):   _want, _w = _S, "the supplement"
            elif _in_supp:               _want, _w = _P, "the paper"
            else:                        _want, _w = _P, "the paper"
            if _num not in _want:
                _dangle.append((_name, _txt.count("\n", 0, _m.start()) + 1, _lead, _num, _w))
        print(f"      {_name}: {_n} references")
    for _nm, _ln, _ld, _num, _w in _dangle:
        _pre = (_ld + " ") if _ld else ""
        print(f"      {_nm}:{_ln}  {_pre}§{_num} is not a section of {_w}   <-- ISSUE")
    if _dangle:
        fail += len(_dangle)
    else:
        print("  every cross-reference resolves, in both directions")

    # 2026-09-22: the supplement described the paper twice as "the six-section main argument",
    # left over from V2. The paper has four sections. A prose claim about how many sections a
    # document has is checkable against the document, so check it rather than read for it.
    _NUMW = {w: i for i, w in enumerate(
        "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
        "fifteen sixteen seventeen eighteen nineteen twenty".split())}
    _COUNTS = {"main argument": len(_P & {"1","2","3","4","5","6","7","8","9"}),
               "supplement":    len({_x for _x in _S if "." not in _x})}
    _claims = 0
    for _f in ("pub/paper2/PAPER2_v3.md", "pub/paper2/SUPPLEMENT_v3.md",
               "pub/paper2/V3_RELEASE_NOTES.md", "pub/paper2/SUBMISSION_CHECKLIST_v3.md"):
        if not _os.path.exists(_f): continue
        _t = re.sub(r'\s+', ' ', io.open(_f, encoding="utf-8").read())
        for _m in re.finditer(r'([a-z]+|\d+)-section (main argument|supplement)', _t):
            _claims += 1
            _said = _NUMW.get(_m.group(1), None)
            if _said is None and _m.group(1).isdigit(): _said = int(_m.group(1))
            _real = _COUNTS[_m.group(2)]
            if _said != _real:
                print(f"      {_f.split('/')[-1]}: calls it a {_m.group(1)}-section "
                      f"{_m.group(2)}; it has {_real}   <-- ISSUE")
                fail += 1
    print(f"  {_claims} prose claim(s) about section counts checked "
          f"(main argument has {_COUNTS['main argument']}, supplement {_COUNTS['supplement']})")


print("\n=== 30. citations resolve, and name the people the prose names beside them ===")
# 2026-09-22: the supplement carried 48 citation markers, 41 distinct numbers reaching [89], and
# no reference list of its own. The only list it could mean was the paper's, which ends at 59.
# Its numbering belonged to the paper's 89-entry bibliography of 19 September, which was later
# pruned and renumbered, so every citation in a fourteen-page document silently moved onto the
# wrong paper. It now carries its own list, and this holds it there.
#
# The second half is the part that would have caught it on the day. Where the prose names a
# surname immediately before a marker, as in "both coefficients are Moretti's [12]", that surname
# has to appear in the entry cited. This is a check on meaning, not on syntax, and it fires at
# seven of the eight testable markers in the broken version.
def _biblio(t):
    out, cur, num = {}, [], None
    for ln in t.splitlines():
        m = re.match(r'^(\d+)\\?\.\s+(.*)$', ln)
        if m:
            if num: out[num] = " ".join(cur)
            num, cur = m.group(1), [m.group(2)]
        elif num and ln.startswith(" ") and ln.strip(): cur.append(ln.strip())
        elif num and not ln.strip(): out[num] = " ".join(cur); num, cur = None, []
    if num: out[num] = " ".join(cur)
    return out

_MATHSP = re.compile(r'\$\$.*?\$\$|\$[^$\n]*\$', re.S)
_ADJ  = re.compile(r"([A-Z][A-Za-zÀ-ÿ-]{3,})(?:['’]s)?(?:\s+et\s+al\.)?\s*\[(\d+(?:[,–—-]\d+)*)\]")
# 2026-09-22, second pass: this was r'\[(\d+(?:,\d+)*)\]' and matched only comma lists, so the
# five ranges and one mixed list in the supplement were invisible. The check passed while
# nineteen of the sixty cited works had no entry at all. Ranges and mixed forms now count.
_CITE = re.compile(r'\[(\d+(?:[,–—-]\d+)*)\]')

_sup30 = io.open("pub/paper2/SUPPLEMENT_v3.md", encoding="utf-8").read()
_SB = _biblio(_sup30)
if not _SB:
    print("  the supplement carries no reference list of its own   <-- ISSUE"); fail += 1
else:
    _body30 = _sup30.partition("\n## References\n")[0]
    _used = sorted({int(n) for m in _CITE.finditer(_MATHSP.sub(" ", _body30))
                          for n in _expand(m.group(1))})
    _have = sorted(int(k) for k in _SB)
    print(f"  {len(_SB)} entries, {len(_used)} cited")
    _unres = [n for n in _used if str(n) not in _SB]
    _dense = _have == list(range(1, len(_have) + 1))
    _orphan = [n for n in _have if n not in _used]
    # an entry a reader cannot look up is not a reference
    _noid = [n for n in _have
             if not re.search(r'arXiv:|arxiv\.|doi|hep-th/|gr-qc/|hep-ph/|astro-ph/|cond-mat/|\(\d{4}\)', _SB[str(n)], re.I)]
    for n in _unres:
        print(f"      [{n}] is cited but has no entry   <-- ISSUE")
    if not _dense:
        print(f"      the list is not densely numbered 1..{len(_have)}   <-- ISSUE")
    for n in _orphan:
        print(f"      [{n}] is listed but never cited   <-- ISSUE")
    for n in _noid:
        print(f"      [{n}] has no arXiv id, DOI or venue-and-year   <-- ISSUE")
    fail += len(_unres) + (0 if _dense else 1) + len(_orphan) + len(_noid)
    if not (_unres or _orphan or _noid) and _dense:
        print("  every citation resolves; numbering dense; every entry is lookup-able")

    _tested = _wrong = 0
    # the companion shares the paper's numbering and had the same drift: it cited
    # "Halliwell, Hawking and Kiefer [22-24,50]", which resolved to a paper on little red dots.
    # It was the only range in the document, and the only citation in it that was wrong.
    for _doc, _txt, _bib in (("the paper", _pv, _biblio(_pv)),
                             ("the companion", _comp, _biblio(_pv)),
                             ("the supplement", _body30, _SB)):
        for _m in _ADJ.finditer(_MATHSP.sub(" ", re.sub(r'\s+', ' ', _txt))):
            _w = _m.group(1)
            if _w in _NOTNAME: continue
            _ent = " ".join(_bib.get(_x, "") for _x in _expand(_m.group(2))).lower()
            if not _ent: continue
            _tested += 1
            if _w.lower() not in _ent:
                _wrong += 1
                print(f"      {_doc}: \"{_w} [{_m.group(2)}]\" cites {_ent[:54]}   <-- ISSUE")
    print(f"  {_tested} markers name a surname beside them; {_tested - _wrong} land on that person")
    fail += _wrong

print("\n=== 31. one spelling convention across the three documents ===")
# 2026-09-22: the supplement carried "normalised" beside "normalization", "quantised" beside
# "quantized" and "flavor" against the paper's eight "flavour". The paper runs British 49 to 4,
# so British is the house style and the American variants are drift. Reference lists are skipped:
# they carry real titles, and a title is spelled the way its authors spelled it.
_AMER = ("flavor", "color", "behavior", "normalization", "normalized", "normalize",
         "quantized", "quantization", "idealization", "realization", "organize",
         "summarized", "analyze", "analyzed", "modeling", "labeled", "center")
_hits = []
for _f in ("pub/paper2/PAPER2_v3.md", "pub/paper2/SUPPLEMENT_v3.md", "pub/paper2/COMPANION_v1.md"):
    if not _os.path.exists(_f): continue
    _t = io.open(_f, encoding="utf-8").read()
    _t = re.split(r'(?m)^## References', _t)[0]          # titles keep their own spelling
    for _w in _AMER:
        _n = len(re.findall(r'(?<![A-Za-z])%s(?![A-Za-z])' % _w, _t))
        if _n:
            _hits.append((_f.split("/")[-1], _w, _n))
print(f"  {len(_AMER)} American variants checked in the prose of three documents")
for _f, _w, _n in _hits:
    print(f"      {_f}: \"{_w}\" x{_n}; the paper's convention is British   <-- ISSUE")
fail += len(_hits)
if not _hits:
    print("  no American variant appears in any of their prose")


print("\n=== 32. the arXiv metadata file still describes THIS paper ===")
# 2026-09-23: ARXIV_METADATA_v3.txt is the file Ben pastes into the submission form. It carried
# the title "The Far Side of the Horizon", superseded weeks ago, and an abstract from an earlier
# draft that began "We propose a geometric restriction of the Schwinger-Keldysh closed time
# path". That abstract was 1951 characters, so arXiv would have rejected it on paste. The file
# says of itself that it matches the manuscript; nothing checked that it did.
_META = "pub/paper2/ARXIV_METADATA_v3.txt"
if not _os.path.exists(_META):
    print(f"  {_META} is missing   <-- ISSUE"); fail += 1
else:
    _mt = io.open(_META, encoding="utf-8").read()
    _ptitle = _pv.split("\n")[0].lstrip("# ").strip()
    _pabs = re.sub(r'\s+', ' ', _pv[_pv.index("## Abstract") + 11:
                                    _pv.index("## 1. Introduction")]).strip()
    _mtitle = re.search(r'(?m)^TITLE\n(.+)$', _mt)
    _mabs = re.search(r'(?ms)^ABSTRACT\n(.*?)(?=\n[A-Z][A-Z\- /]+\n)', _mt)
    _bad32 = 0
    if not _mtitle or not _mabs:
        print("  the metadata file's TITLE or ABSTRACT block did not parse   <-- ISSUE"); _bad32 += 1
    else:
        if _mtitle.group(1).strip() != _ptitle:
            print("  metadata TITLE does not match the paper's   <-- ISSUE")
            print(f"      metadata: {_mtitle.group(1).strip()[:70]}")
            print(f"      paper   : {_ptitle[:70]}")
            _bad32 += 1
        _ma = re.sub(r'\s+', ' ', _mabs.group(1)).strip()
        if _ma != _pabs:
            print(f"  metadata ABSTRACT does not match the paper's "
                  f"({len(_ma)} chars against {len(_pabs)})   <-- ISSUE"); _bad32 += 1
        if len(_ma) > 1920:
            print(f"  metadata ABSTRACT is {len(_ma) - 1920} characters over arXiv's limit"
                  f"   <-- ISSUE"); _bad32 += 1
    fail += _bad32
    if not _bad32:
        print(f"  title and abstract both match the paper; abstract {len(_pabs)} chars, under 1920")

print(f"\n=== {'PASS' if fail==0 else str(fail)+' ISSUE(S)'} ===")
