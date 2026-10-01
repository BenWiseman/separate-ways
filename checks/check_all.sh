#!/bin/bash
# check_all.sh -- every gate the companion has to pass, in one command.
# Run from the repository root:   bash checks/check_all.sh
# Exit status 0 means all of them passed. Each section prints what it checked.
set -u
PAPER="${1:-papers/2_over_the_horizon/COMPANION_v1.md}"
BENLM="${BENLM:-/home/ben/benlm}"
fail=0

# Say what is missing before running anything. checks/reproduce.sh has had this preflight for
# a while and this pass did not, so on 2026-09-30 a python3 without mpmath produced fourteen
# claims reading "number not in the output of ..." and an hour went into looking for a
# numerical discrepancy that was never there. claims_check now names a script that died, and
# this stops the run before it can mislead anyone in the first place.
_miss=$(python3 - <<'PYDEP'
import importlib.util
print(" ".join(m for m in ("numpy", "scipy", "sympy", "mpmath")
               if importlib.util.find_spec(m) is None))
PYDEP
)
if [ -n "$_miss" ]; then
  echo "This interpreter is missing:$_miss" >&2
  echo "  python3 -m pip install -r checks/requirements.txt" >&2
  echo "Claims behind the Python scripts cannot be checked without them; stopping." >&2
  exit 2
fi
# About a third of the failure assignments below are bare shell comparisons that print nothing
# distinctive, so a failing run showed only self-test plants and the gate had to be found by
# sweeping every checker's exit code by hand. That cost an hour on 2026-09-29. hdr compares the
# count against its value at the previous gate and names whichever gate raised it.
#
# It is a COUNT and not a flag, and that is the whole point. It was a flag until 2026-09-29, so
# the comparison below could fire exactly once: the first gate to fail was named and every gate
# failing after it was silent. Gate 5 and gate 9 both failed that afternoon, gate 9 was invisible,
# and the run was reported as one failure. A flag cannot say "and also". The footer now lists
# every gate that raised a failure, and the logic was tested by planting three failures among
# five gates and requiring all three names and neither of the two passing ones.
_prev_fail=0
_gate="(startup)"
_failed=""
hdr () {
  if [ "$fail" != "$_prev_fail" ]; then
    printf '\033[1m   ^^ GATE RAISED THE FAILURE: %s\033[0m\n' "$_gate"
    _failed="$_failed
   $_gate"
    _prev_fail=$fail
  fi
  _gate="$1"
  printf '\n\033[1m== %s ==\033[0m\n' "$1"
}

hdr "1. arXiv abstract cap"
python3 checks/calc/abstract_len.py "$PAPER" || fail=$((fail+1))

hdr "2. internal and external cross-references"
python3 checks/calc/xref_check.py || fail=$((fail+1))
# And across the two papers, in both directions. This existed as a function from 2026-09-23,
# sat below the __main__ block's sys.exit so it could not run, and was called by nothing.
python3 checks/calc/xref_check.py --external || fail=$((fail+1))

# A reference to a document the release does not carry is worse than a dangling internal one,
# because a reader can see at once that it is missing. Three had accumulated: two to a "Supplement"
# and one to "Supplement S11.8", none of which ships. The gate names the words that must not appear
# and carries a plant, since a grep for something absent passes trivially.
python3 - papers/2_over_the_horizon/COMPANION_v1.md papers/1_separate_ways/PAPER2_v4_draft.md <<'PYSUP'
import io, re, sys
# "the supplement" in lower case slipped past a capital-S ban for a day and pointed at a
# document the release does not carry. The angle sense, "an angle with its supplement", is
# legitimate and appears three times, so what is banned is the article, not the word.
BAD = [r"Supplement\s*(?:§\s*)?S?\d", r"\bSupplement\b", r"\bthe\s+supplement\b"]
bad = 0
for path in sys.argv[1:]:
    t = io.open(path, encoding="utf-8").read()
    hits = [m.group(0) for p in BAD for m in re.finditer(p, t)]
    print("   %-28s references to a document not in the release: %d%s"
          % (path.split("/")[-1], len(hits), "  <-- ISSUE" if hits else ""))
    for h in sorted(set(hits))[:4]:
        print("        %r" % h)
    bad |= bool(hits)
plant = "see Supplement S4 for the rest"
fired = any(re.search(p, plant) for p in BAD)
print("   plant: a planted Supplement reference is caught: %s" % ("yes" if fired else "NO"))
if not fired:
    bad = 1
clean = "see the companion paper for the rest"
spared = not any(re.search(p, clean) for p in BAD)
print("   plant: an ordinary companion reference is spared: %s" % ("yes" if spared else "NO"))
if not spared:
    bad = 1
sys.exit(1 if bad else 0)
PYSUP
[ $? -ne 0 ] && fail=$((fail+1))

# And a number the cosmology paper ATTRIBUTES to the companion has to be in the companion. On
# 2026-09-27 it said "the companion retains I = 0.0127596673634 as a quadrature check"; that
# number appears nowhere in the companion. It lives in archive/release_4_2_scripts/cross_checks/dm_clock.R, a
# directory whose name is not the paper's, which is almost certainly how it got there. A referee
# who follows the reference finds nothing, and no gate could see it.
python3 - papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY2B'
import io, re, sys
cos = " ".join(io.open(sys.argv[1], encoding="utf-8").read().split())
com = " ".join(io.open(sys.argv[2], encoding="utf-8").read().split())
NUM = r"\d+\.\d{3,}|\d+\.\d+\\times10\^\{-?\d+\}"
def misses(text, other):
    out = []
    for s in re.split(r"(?<=[.!?]) +", text):
        if not re.search(r"\bcompanion\b", s): continue
        for n in set(re.findall(NUM, s)):
            if n not in other: out.append((n, s))
    return out
bad = misses(cos, com)
for n, s in bad:
    print('   the cosmology paper credits the companion with %s, which is not in it  <-- ISSUE' % n)
    print('      "%s"' % s[:150])
if not bad:
    n = sum(1 for s in re.split(r"(?<=[.!?]) +", cos) if re.search(r"\bcompanion\b", s))
    print("   %d sentences credit the companion, and every number in them is in it" % n)
# the plant, on the real text
planted = cos + " The companion computes 0.1234567 for that."
assert misses(planted, com), "plant: a number credited to the companion but absent from it is not caught"
assert not misses(cos + " The companion computes 3.51656 for that.", com), \
       "plant: a number that IS in the companion is wrongly flagged"
print("   plant: a credited number absent from the companion is caught, one present is spared: yes")
sys.exit(1 if bad else 0)
PY2B
[ $? -ne 0 ] && fail=$((fail+1))

hdr "3. citations point at the papers the prose names"
python3 checks/calc/citation_check.py "$PAPER" || fail=$((fail+1))

hdr "4. every claim in CLAIMS.tsv reproduces, and still lands in the paper"
# Ben, 2026-09-22 and again 2026-09-24: "you can't have script paths in a manuscript at all."
# The provenance those paths carried had to survive the strip rather than go with it, so it moved
# to checks/CLAIMS.tsv and this checks both halves of every row: the passage is still there, and
# the script named for it still prints the numbers it claims.
python3 checks/claims_check.py "$PAPER"
[ $? -ne 0 ] && fail=$((fail+1))

hdr "5. machine-prose tics"
python3 "$BENLM/tools/llm_tics.py" "$PAPER" 2>&1 | tail -3
python3 "$BENLM/tools/tic_count.py" "$PAPER" 2>&1 | head -1
echo "   (the bar is 4.0 per 1000, which is where the cosmology paper sat when it was set)"
# The bar above was printed and not enforced for weeks, and the density drifted to 4.1
# without anything saying so. Enforce it. The rate is recomputed from the raw counts
# because the printed figure is rounded and 4.04 prints as 4.0.
ticrate=$(python3 "$BENLM/tools/tic_count.py" "$PAPER" 2>&1 | head -1 \
          | awk '{printf "%.4f", $2 / $7 * 1000}')
echo "   density $ticrate per 1000, bar 4.0000"
awk -v r="$ticrate" 'BEGIN { exit !(r > 4.0) }' && {
  echo "   ABOVE THE BAR  <-- ISSUE"; fail=$((fail+1)); }

# Register frequency, added 2026-09-30. Ben: "instead of computed say calculated more often,
# you say computed a LOT". Not a ban, a mix: tic_count carries the watchlist and the rate bar.
# Enforced for the companion and the Letter. PAPER2_v4_draft.md sits at 0.68 against a bar of
# 0.60 and is NOT edited to fit: it is the submitted manuscript, frozen at commit b1054c3.
for _rf in "$PAPER" papers/3_road_to_nowhere/LETTER2_CONTACT_v2.md; do
  [ -f "$_rf" ] || continue
  _line=$(python3 "$BENLM/tools/tic_count.py" "$_rf" 2>&1 | grep '^   register' || true)
  [ -n "$_line" ] && echo "   $(basename "$_rf"): ${_line#   register: }"
  case "$_line" in *OVER*) echo "   ABOVE THE REGISTER BAR  <-- ISSUE"; fail=$((fail+1));; esac
done
_line=$(python3 "$BENLM/tools/tic_count.py" papers/1_separate_ways/PAPER2_v4_draft.md 2>&1 \
        | grep '^   register' || true)
echo "   PAPER2_v4_draft.md (submitted, not edited): ${_line#   register: }"

# The bar above was only ever applied to $PAPER, so the cosmology paper's own density was
# ungated and drifted from the 4.0 it sat at when the bar was written to 4.4. Most of its
# constructions are scope statements doing real work ("a posterior upper limit and not a hard
# edge"), so the number to hold is its own and not the companion's; what this stops is any
# further drift. Cutting it toward 4.0 means cutting precision and is not worth it.
cosrate=$(python3 "$BENLM/tools/tic_count.py" papers/1_separate_ways/PAPER2_v4_draft.md 2>&1 | head -1 \
          | awk '{printf "%.4f", $2 / $7 * 1000}')
echo "   cosmology paper density $cosrate per 1000, its own bar 4.4100"
awk -v r="$cosrate" 'BEGIN { exit !(r > 4.41) }' && {
  echo "   THE COSMOLOGY PAPER IS ABOVE ITS BAR  <-- ISSUE"; fail=$((fail+1)); }
# Drama fragments. llm_tics only catches the "Their result? Wrong." form and missed
# "Then the useful question.", "Three limits." and "Two consequences." sitting in the body.
# Baseline is ONE allowed fragment: "Outside a horizon, nothing.", which opens 5.1 under the
# heading "The condition, and the answer". The heading asks and the line answers in three words.
# It is on the banned list and it is the one place the ban would cost something, so it is named
# here rather than silently rewritten or silently ignored. Any second fragment fails.
nfrag=$(python3 "$BENLM/tools/fragment_check.py" "$PAPER" 2>&1 | grep -oE '^  [0-9]+ short' | grep -oE '[0-9]+')
echo "   $nfrag verbless short sentence(s) (1 allowed, 5.1's three-word answer)"
[ "${nfrag:-9}" -gt 1 ] && { python3 "$BENLM/tools/fragment_check.py" "$PAPER" | tail -n +3; fail=$((fail+1)); }
# Repeated sentence openings. Both rulesets said llm_tics measured this; it never did.
# Baseline is ONE allowed run: A.12's "If the two anticommute... If the swap is trivial...
# If the swap IS the time reversal...", which is a three-case enumeration and not cadence.
# Any run beyond that is new and fails.
nrun=$(python3 "$BENLM/tools/opener_runs.py" "$PAPER" 2>&1 | grep -oE 'same word: [0-9]+' | grep -oE '[0-9]+')
echo "   $nrun run(s) of 3+ sentences opening on the same word (1 allowed, A.12's three cases)"
[ "${nrun:-9}" -gt 1 ] && { python3 "$BENLM/tools/opener_runs.py" "$PAPER" | sed -n '3,20p'; fail=$((fail+1)); }
# The same check on the cosmology paper, which it had never been applied to. It carried TEN runs
# on 2026-09-27, among them six consecutive sentences opening on "the" in 4.5. Six were rewritten
# by varying the opening word only. The four that remain are enumerations where the parallel
# opening is the right choice and not cadence: 1's roadmap over sections 2 to 4, 3.6's "four
# things / the first / the other three", 4.4's list of returned quantities, and 4.5's three
# refuting observations. Four is therefore the bar, and a fifth run is new.
# Four caveats in a row. REORIENT names the drift as premature reasonableness and Ben's form of
# it is that a result under four caveats reads as an apology; nothing measured it. Overall caveat
# density says nothing useful, both bodies sitting near one sentence in five with no baseline to
# judge that against, but consecutive runs are the failure mode and neither body has one, so the
# bar is zero. The two that exist are in appendices, where a numbered series of distinctions is
# the content.
python3 checks/caveat_runs.py --selftest || fail=$((fail+1))
python3 checks/caveat_runs.py papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md || fail=$((fail+1))

cosrun=$(python3 "$BENLM/tools/opener_runs.py" papers/1_separate_ways/PAPER2_v4_draft.md 2>&1 \
         | grep -oE 'same word: [0-9]+' | grep -oE '[0-9]+')
echo "   $cosrun run(s) in the cosmology paper (4 allowed, all of them enumerations)"
[ "${cosrun:-9}" -gt 4 ] && { python3 "$BENLM/tools/opener_runs.py" papers/1_separate_ways/PAPER2_v4_draft.md \
                              | sed -n '3,24p'; fail=$((fail+1)); }

hdr "6. paragraph lengths"
# This ran on $PAPER alone, so the cosmology paper, which is the one a reader meets first, had
# nine paragraphs between 257 and 368 words and nothing to say so. Both manuscripts are measured
# now, and the appendix anchor is matched on the prefix because one says "Appendices" and the
# other "Appendix A".
python3 - "$PAPER" papers/1_separate_ways/PAPER2_v4_draft.md <<'PY'
import io,re,sys
BAR=250
def longest(lines):
    ps=[q for q in re.split(r"\n\s*\n","\n".join(lines))
        if q.strip() and not q.strip().startswith(("#","$$","|","-","*","![",">"))]
    return (max(len(q.split()) for q in ps) if ps else 0), ps
def measure(path):
    L=io.open(path,encoding="utf-8").read().split("\n")
    i=next(k for k,x in enumerate(L) if x.startswith("## Appendi"))
    j=next((k for k,x in enumerate(L) if x.startswith("## References")), len(L))
    out=[]
    for nm,ls in (("body",L[35:i]),("appendix",L[i:j])):
        w,_=longest(ls)
        out.append((nm,len(" ".join(ls).split()),w))
    return out
bad=0
for path in sys.argv[1:]:
    print("   %s" % path.split("/")[-1])
    for nm,words,w in measure(path):
        print("     %-9s %5d words, longest paragraph %d%s"%(nm,words,w,
              "   <-- OVER %d"%BAR if w>BAR else ""))
        bad |= w>BAR
# The plant. A gate that has just been widened to a second document has to be shown to fire on it,
# because the widening is the part with no history of catching anything.
synth=["x"]*35+["## Appendi","", " ".join(["word"]*(BAR+1)), "", "## References"]
w,_=longest(synth[35:36+3])
print("   plant: a %d-word paragraph is flagged: %s" % (w, "yes" if w>BAR else "NO"))
if w<=BAR: bad=1
short,_=longest(["", " ".join(["word"]*(BAR-1)), ""])
print("   plant: a %d-word paragraph is spared: %s" % (short, "yes" if short<=BAR else "NO"))
if short>BAR: bad=1
sys.exit(1 if bad else 0)
PY
[ $? -ne 0 ] && fail=$((fail+1))

PY=python3; [ -x .venv/bin/python ] && PY=.venv/bin/python
hdr "7. every script runs (note: the fig_*.R generators rewrite their PDFs)"
n=0; bad=0
for f in checks/calc/*.R checks/calc/*.py checks/*.R; do
  [ -f "$f" ] || continue
  n=$((n+1))
  case "$f" in
    *.py) "$PY" "$f" >/dev/null 2>&1 || { echo "   FAIL $f"; bad=1; } ;;
    *)    Rscript "$f" >/dev/null 2>&1 || { echo "   FAIL $f"; bad=1; } ;;
  esac
done
[ $bad -eq 0 ] && echo "   $n scripts, all pass" || fail=$((fail+1))

hdr "8. the arXiv metadata has not drifted from the manuscript"
python3 - "$PAPER" <<'PY4'
import io, re, sys, os
meta = "papers/2_over_the_horizon/ARXIV_METADATA_COMPANION.txt"
if not os.path.exists(meta):
    print("   %s missing" % meta); sys.exit(1)
c = io.open(sys.argv[1], encoding="utf-8").read()
m = re.search(r"(?m)^##\s+Abstract\s*$", c); rest = c[m.end():]
n = re.search(r"(?m)^#{1,6}\s+\S", rest)
paper_abs = " ".join(rest[:n.start()].split())
t = io.open(meta, encoding="utf-8").read()
mm = re.search(r"(?ms)^ABSTRACT\n(.*?)\n\n#", t)
meta_abs = " ".join(mm.group(1).split()) if mm else ""
title = c.split("\n")[0].lstrip("# ").strip()
ok = (paper_abs == meta_abs) and (title in t)
print("   abstract matches: %s   title matches: %s   (%d chars)"
      % (paper_abs == meta_abs, title in t, len(paper_abs)))
# The COMMENTS line states word, figure and reference counts, and every one of them had gone
# stale: 3 figures where there are 14, 32 references where there are 33. A number a submission
# form states about the paper is a number that drifts like any other.
import re as _re
body_all, _, refs = c.partition("\n## References")
body, _, app = body_all.partition("## Appendices")
strip = lambda x: _re.sub(r"!\[.*?\]\([^)]*\)", " ", x, flags=re.S)
nw   = len(strip(body).split()) + len(strip(app).split())
napp = len(strip(app).split())
nfig = len(_re.findall(r"(?m)^!\[", c))
nref = len(_re.findall(r"(?m)^\d+\\?\.\s+[A-Z]", refs))
cm = " ".join(t.split())
# Word counts move on every prose edit, and a gate that fails on every edit is a gate somebody
# deletes, so the COMMENTS line states them as approximations and this checks the approximation
# is close. Figure and reference counts do not churn and stay exact.
for _v, _what in ((nw, "words"), (napp, "are the appendices")):
    _m = _re.search(r"about (\d+) %s" % _re.escape(_what), cm)
    if not _m:
        print('   COMMENTS does not state %s as "about N"' % _what); ok = False
    elif abs(int(_m.group(1)) - _v) > 500:
        print("   COMMENTS says about %s %s; the manuscript has %d" % (_m.group(1), _what, _v))
        ok = False
for _v, _what in ((nfig, "figures"), (nref, "references")):
    if f"{_v} {_what}" not in cm and f"{_v} {_what.split()[0]}" not in cm:
        print("   COMMENTS does not say %d %s" % (_v, _what)); ok = False
if ok:
    print("   COMMENTS: %d words, %d appendix, %d figures, %d references, all stated" %
          (nw, napp, nfig, nref))
sys.exit(0 if ok else 1)
PY4
[ $? -ne 0 ] && fail=$((fail+1))

hdr "17. the positions the release takes are not contradicted somewhere else"
# Every other gate here checks a number, a reference or a word. None of them can see two paragraphs
# taking opposite sides of the same question, which is what A.11 did about the seam for a day while
# the whole suite passed. checks/POSITIONS.tsv lists the contested points with the phrase that has
# to appear and the phrases that must not, and it carries the standing bans too: no draft-history
# narration, no reference to a Supplement neither release has, no script paths. Both failure modes
# were made to fire on the real manuscripts before this was installed.
python3 - papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY17'
# gate: the positions the release takes are stated once and never contradicted
import io, re, sys
COS, COM = sys.argv[1], sys.argv[2]
flat = lambda p: " ".join(io.open(p, encoding="utf-8").read().split())
docs = {"PAPER2": flat(COS), "COMPANION": flat(COM)}
rows = [l.rstrip("\n").split("\t") for l in io.open("checks/POSITIONS.tsv", encoding="utf-8")
        if l.strip() and not l.startswith("#")]
assert rows, "POSITIONS.tsv is empty: the gate would pass by checking nothing"
ok, n = True, 0
for topic, which, must, never in rows:
    t = docs[which]
    if must != "-":
        n += 1
        if not re.search(must, t):
            print(f"   {topic} ({which}): the position is no longer stated  <-- ISSUE")
            print(f'      wanted: "{must[:70]}"'); ok = False
    if never != "-":
        for bad in never.split("|"):
            n += 1
            hit = re.search(bad, t)
            if hit:
                i = max(0, hit.start() - 55)
                print(f"   {topic} ({which}): says the opposite  <-- ISSUE")
                print(f'      "{t[i:hit.end()+55]}"'); ok = False
if ok:
    print(f"   {len(rows)} positions, {n} assertions, none contradicted")
# the plants: both directions have to be catchable
pl = docs["COMPANION"] + " This section previously reported nothing of the kind."
assert re.search("previously reported", pl), "plant: a forbidden phrase is not caught"
assert not re.search("adopted, with the corner term's uniqueness",
                     docs["COMPANION"].replace("adopted, with the corner term's uniqueness", "")), \
       "plant: a missing required phrase is not caught"
print("   plant: a planted draft-history phrase is caught: yes")
print("   plant: a deleted position is caught: yes")
sys.exit(0 if ok else 1)
PY17
[ $? -ne 0 ] && fail=$((fail+1))

# And three of those positions have to come LAST, not merely appear. Ben, 2026-09-27: the abstract
# and the conclusions should finish on the field equations following without gravity being
# quantised, because that is the stronger claim and it was landing on a DESI neutrino bound
# instead. A phrase check cannot see that; this reads the final sentence and the final paragraph.
python3 - papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY17B'
import io, re, sys
cos, com = (io.open(f, encoding="utf-8").read() for f in sys.argv[1:3])

def abstract(c):
    m = re.search(r"(?m)^##\s+Abstract\s*$", c); rest = c[m.end():]
    n = re.search(r"(?m)^#{1,6}\s+\S", rest)
    return " ".join(rest[:n.start()].split())

def last_para(c, head, nxt):
    body = c[c.index(head):c.index(nxt)]
    return " ".join([b for b in body.split("\n\n") if b.strip()][-1].split())

cases = [("the abstract's last sentence",
          re.split(r"(?<=[.!?]) +", abstract(cos))[-1],
          # "owes" was Ben's word to cut on 2026-09-28: a thing a theory "owes" is Claude-ish
          # and the plain verb is "needs". The position this gate defends is unchanged.
          "a quantum theory of gravity still needs"),
         ("section 5's last paragraph",
          last_para(cos, "## 5. Conclusions", "## Appendix A"),
          "waiting on a quantum theory of gravity"),
         ("the companion's last paragraph",
          last_para(com, "## 9. Conclusions", "## Appendices"),
          "gravity did not have to be quantised")]
ok = True
for what, text, want in cases:
    if want not in text:
        print('   %s does not land on the field equations  <-- ISSUE' % what)
        print('      wanted "%s"' % want)
        print('      ends   "...%s"' % text[-90:]); ok = False
    else:
        print("   %-32s lands on it" % what)
assert "gravity did not have to be quantised" not in "a paragraph about something else entirely", \
       "plant: the landing check would pass on unrelated text"
print("   plant: a last paragraph that does not land is caught: yes")
sys.exit(0 if ok else 1)
PY17B
[ $? -ne 0 ] && fail=$((fail+1))

hdr "9. superseded scripts announce it when run"
# A script declares its OWN status with one of these forms. Merely mentioning a
# withdrawal elsewhere in the file is not a declaration: restricted_exponent_ratio.R
# describes the withdrawn 19 per cent delay and is itself perfectly current.
nsup=0
for f in checks/calc/*.R; do
  head -30 "$f" | grep -qE "^# *#* *(SUPERSEDED|WITHDRAWN)|IS WITHDRAWN" || continue
  nsup=$((nsup+1))
  out=$(timeout 90 Rscript "$f" 2>&1)
  h=$(printf '%s' "$out" | head -6 | grep -ciE "superseded|withdrawn")
  t=$(printf '%s' "$out" | tail -12 | grep -ciE "superseded|withdrawn")
  if [ "$h" -ge 1 ] && [ "$t" -ge 1 ]; then
    echo "   $(basename "$f"): banner at both ends"
  else
    echo "   $(basename "$f"): HEADER SAYS SUPERSEDED BUT THE OUTPUT DOES NOT"; fail=$((fail+1))
  fi
done
# The cosmology paper states that count about itself, and a number a paper states about its own
# repository is a number that goes stale. It said 178 scripts when there were 166, so this is
# gated now rather than trusted.
if grep -q "^$nsup of those files carry a banner\|. $nsup of those files carry a banner" papers/1_separate_ways/PAPER2_v4_draft.md; then
  echo "   the cosmology paper's count of $nsup is current"
else
  echo "   THE COSMOLOGY PAPER DOES NOT SAY $nsup OF THOSE FILES CARRY A BANNER  <-- ISSUE"; fail=$((fail+1))
fi
if grep -q "$((nsup+1)) of those files carry a banner" papers/1_separate_ways/PAPER2_v4_draft.md; then
  echo "   plant: an off-by-one count would be caught: no, it is already there"; fail=$((fail+1))
else
  echo "   plant: an off-by-one count would be caught: yes"
fi

hdr "10. every appendix is reachable from the body"
python3 - "$PAPER" <<'PY3'
import io, re, sys
c = io.open(sys.argv[1], encoding="utf-8").read()
i = c.index("## Appendices"); body, app = c[:i], c[i:]
heads = set(m.group(1) for m in re.finditer(r"(?m)^\*\*(A\.\d+)\s", app))
refs  = set("A." + n for n in re.findall(r"\bA\.(\d+)\b", body))
orph  = sorted(heads - refs, key=lambda x: int(x[2:]))
print("   %d appendices, %d referenced from the body" % (len(heads), len(heads & refs)))
if orph: print("   NEVER REFERENCED: %s" % ", ".join(orph))
sys.exit(1 if orph else 0)
PY3
[ $? -ne 0 ] && fail=$((fail+1))

hdr "11. figure type is above the legibility floor"
python3 - <<'PY2'
import re, io, glob, sys
floor = 9.0; ps = 12.0   # R pdf() default pointsize
bad = 0
for f in sorted(glob.glob("checks/fig_*.R")):
    s = io.open(f, encoding="utf-8").read()
    c = [float(m.group(1)) for m in
         re.finditer(r"(?:text|mtext)\([^)]*cex\s*=\s*([\d.]+)", s, re.S)]
    if not c: continue
    pt = min(c) * ps
    print("   %-34s smallest type %.2f pt%s" % (f.split("/")[-1], pt,
          "   <-- BELOW %.1f" % floor if pt < floor else ""))
    bad |= pt < floor
sys.exit(1 if bad else 0)
PY2
[ $? -ne 0 ] && fail=$((fail+1))

hdr "12. the availability note's counts match the tree"
# These drifted: the note said 33 files under checks/calc/ when there were 42, and never
# mentioned the figure generators at all. Every number a paper states about itself is a number
# that can go stale, including the ones about its own repository.
python3 - "$PAPER" papers/1_separate_ways/PAPER2_v4_draft.md <<'PY13'
import io, re, sys, glob
t = io.open(sys.argv[1], encoding="utf-8").read()
cosmo = io.open(sys.argv[2], encoding="utf-8").read()
calc = len(glob.glob("checks/calc/*.R")) + len(glob.glob("checks/calc/*.py"))
figs = len(glob.glob("checks/fig_*.R"))
# The manuscript cites nothing inline now, so a count of "uncited" files would be all of them.
# What the note claims is curated: ten files support no result in this paper, being the three
# checkers, four superseded scripts, one adjudicator, one that supports the cosmology paper and
# one that prices a withdrawn test. That is an author's accounting and not a derived number, so
# the gate checks the two counts that ARE mechanical and leaves the ten to the prose.
un = None
ok = True
# "one of the fourteen" sat in the note while fifteen generators drew the figures, because the
# phrasing carried no noun this gate could match. The count of generators USED here is now stated
# in a form that can be checked, and is.
used = 0
imgs = set(re.findall(r"!\[[^\]]*\]\(([^)]+)\)", t))
for g in sorted(glob.glob("checks/fig_*.R")):
    gs = io.open(g, encoding="utf-8").read()
    if any(im.split("/")[-1].rsplit(".", 1)[0] in gs for im in imgs): used += 1
for n, what in ((calc, "calculation files"), (figs, "figure generators"),
                (used, "of those generators")):
    if f"the {n} {what}" not in t and f"{n} {what}" not in t:
        print(f"   the note does not say {n} {what}"); ok = False
# The cosmology paper carries the same two counts and had a stale total of 178 for years, so it is
# checked here rather than left to the prose.
for n, what in ((calc, "calculation files"), (figs, "figure generators")):
    if f"{n} {what}" not in cosmo:
        print(f"   the cosmology paper does not say {n} {what}"); ok = False
# and the plant: an off-by-one must not already be present in either file
for n, what in ((calc + 1, "calculation files"), (figs + 1, "figure generators")):
    if f"{n} {what}" in t or f"{n} {what}" in cosmo:
        print(f"   plant: an off-by-one {what} count is ALREADY in a manuscript"); ok = False
if ok:
    print(f"   {calc} calc scripts, {figs} figure generators, {un} uncited: both manuscripts current")
    print("   plant: an off-by-one count in either manuscript would be caught: yes")
sys.exit(0 if ok else 1)
PY13
[ $? -ne 0 ] && fail=$((fail+1))

hdr "22. a computed constant and its copies agree"
# J2end = 47.561945 is computed in contact_vanvleck.R and typed into eight other files,
# and it sits under fork 9's calibration target and under kappa, which sets the shell.
# Nothing was comparing the copies with the computation.
python3 checks/shared_constants.py --selftest || fail=$((fail+1))
python3 checks/shared_constants.py || fail=$((fail+1))

hdr "14. no sentence appears twice"
python3 checks/repeat_check.py "$PAPER" || fail=$((fail+1))

hdr "15. the maths closes and every figure is where it says"
python3 checks/structure_check.py "$PAPER" || fail=$((fail+1))
# And the display maths must be in the form the TeX build can read. Balanced dollars are not
# enough: on 2026-09-27 the companion carried one display written as a blockquote with a lone "$"
# on each side, which structure_check counted as balanced and pandoc escaped into literal dollar
# signs, giving three TeX errors and a dropped glyph the first time the companion was ever built.
python3 - papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY15'
import io, re, sys
LONE = re.compile(r"^(?:>\s*)*\$\s*$")          # a single dollar alone on a line, blockquoted or not
DELIM = re.compile(r"^(?:>\s*)*\$\$\s*$")       # the form pandoc reads, at column zero or in a quote
HYPHEN = re.compile(r"^(?!#).*[A-Za-z]-$")       # a wrapped word, which pandoc joins with a space
ok = True
for path in sys.argv[1:]:
    name, n = path.split("/")[-1], 0
    for i, l in enumerate(io.open(path, encoding="utf-8").read().split("\n"), 1):
        if LONE.match(l):
            print("   %s:%d a lone $ opens a display; pandoc escapes it and the build breaks."
                  " Use $$  <-- ISSUE" % (name, i)); ok = False
        elif DELIM.match(l):
            n += 1
        elif HYPHEN.match(l):
            # pandoc joins wrapped lines with a space, so a line ending in a hyphen prints as
            # "classical- quantum". Three of these shipped in the 2026-09-27 build.
            print("   %s:%d the line ends in a hyphen, which prints as \"word- word\""
                  "  <-- ISSUE" % (name, i)); ok = False
    # And the damage after the fact. A wrapper that breaks at a hyphen, then a later reflow that
    # joins lines with a space, leaves "factor- of-two" baked into the file. A suspended hyphen,
    # "in- and out-region", is correct English and is what follows the dash that tells them apart.
    flat = " ".join(io.open(path, encoding="utf-8").read().split())
    for m in re.finditer(r"\w+- ([a-z]\w+)", flat):
        if m.group(1) not in ("and", "or"):
            print('   %s: "%s" is a compound broken by a space  <-- ISSUE' % (name, m.group(0)))
            ok = False
    if n % 2:
        print("   %s: an odd number of $$ delimiters, %d  <-- ISSUE" % (name, n)); ok = False
    elif ok:
        print("   %-22s %d display blocks, every delimiter a $$" % (name, n // 2))
for good in ("$$", "> $$", "text $x$ text"):
    assert not LONE.match(good), "plant: %r is wrongly caught" % good
for bad in ("$", "> $", ">   $  "):
    assert LONE.match(bad), "plant: %r is not caught" % bad
assert HYPHEN.match("the classical-"), "plant: a wrapped hyphen is not caught"
import re as _re2
assert _re2.match(r"\w+- ([a-z]\w+)", "factor- oftwo"), "plant: a space-broken compound is not caught"
assert _re2.match(r"\w+- ([a-z]\w+)", "in- and").group(1) == "and", "plant: the suspended form is recognised"
for good in ("an em-dash free line", "# A heading-", "$x$ = -"):
    assert not HYPHEN.match(good), "plant: %r is wrongly caught" % good
print("   plant: a lone $ is caught and a $$ delimiter is spared: yes")
sys.exit(0 if ok else 1)
PY15
[ $? -ne 0 ] && fail=$((fail+1))

hdr "16. the cosmology paper's arXiv metadata has not drifted either"
# This gate covers the OTHER paper. Gate 8 has guarded the companion's metadata since
# the day it was written and nothing guarded Separate Ways', which is why on 2026-09-24
# its metadata still carried the v3 title and a v3 abstract differing from the draft's
# from the first character. That is what would have been pasted into arXiv.
python3 checks/sync_metadata_paper2.py --check || fail=$((fail+1))

hdr "18. the nominator brief has not drifted from the release"
# The brief is outreach material and the public export does not carry it, so this gate has
# nothing to measure there. Skipping beats failing: a repository whose own suite reports 41
# failures teaches a reader to ignore the suite.
if [ ! -f pub/paper2/NOMINATOR_BRIEF.md ]; then
  echo "   not in this repository (outreach material); gate does not apply here"
else
# The brief is what a nominator reads, and it drifts exactly the way the arXiv metadata drifted.
# On 2026-09-27 it still carried an entanglement witness at 270.30 PeV that v4.5 does not claim,
# "178 scripts" against 173, "thirty-two checks" against 17, and led with the mass ceiling rather
# than the field equations. Two halves: every physics number in the brief must appear in one of
# the manuscripts, and every count it states about the tree must match the tree.
python3 - pub/paper2/NOMINATOR_BRIEF.md papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY18'
import io, re, sys, os, glob
brief, cos, com = sys.argv[1], sys.argv[2], sys.argv[3]
flat = lambda p: " ".join(io.open(p, encoding="utf-8").read().split())
papers = flat(cos) + "  " + flat(com)
b = flat(brief)

# the counts the brief states about the tree, checked against the tree
ncalc = len(glob.glob("checks/calc/*.R")) + len(glob.glob("checks/calc/*.py"))
nfig  = len(glob.glob("checks/fig_*.R"))
nline = sum(sum(1 for _ in io.open(f, encoding="utf-8", errors="replace"))
            for f in sorted(glob.glob("checks/calc/*.R") + glob.glob("checks/calc/*.py")
                            + glob.glob("checks/fig_*.R")))
ngate = len(re.findall(r"(?m)^hdr \"", io.open("checks/check_all.sh", encoding="utf-8").read()))
rows  = [l for l in io.open("checks/CLAIMS.tsv", encoding="utf-8") if l.strip() and not l.startswith("#")]
nclaim = len(rows)
sys.path.insert(0, os.path.expanduser("~/benlm/tools"))
import number_provenance as _np          # count numbers the way claims_check counts them,
nnum = sum(1 for r in rows               # so the two can never disagree
           for lit in r.rstrip("\n").split("\t")[2].split(";")
           if lit and _np.nums_from_claim(lit))
prows = [l.rstrip("\n").split("\t") for l in io.open("checks/POSITIONS.tsv", encoding="utf-8")
         if l.strip() and not l.startswith("#")]
npos = len(prows)
nass = sum((r[2] != "-") + (len(r[3].split("|")) if r[3] != "-" else 0) for r in prows)
nrefc = len(re.findall(r"(?m)^\d+\\?\.\s+[A-Z]", io.open(cos, encoding="utf-8").read().partition("\n## References")[2]))
nrefm = len(re.findall(r"(?m)^\d+\\?\.\s+[A-Z]", io.open(com, encoding="utf-8").read().partition("\n## References")[2]))

ok = True
tree = [(ncalc, "calculation files"), (nfig, "figure generators"),
        (ngate, "gates"), (nnum, "numbers"), (nclaim, "claims"), (npos, "positions"),
        (nass, "assertions"), (nrefc, "references in the cosmology paper"),
        (nrefm, "in the companion")]
for v, what in tree:
    pat = "%d %s" % (v, what)
    if pat not in b and "%d %s" % (v, what.split()[0]) not in b:
        print('   the brief does not say "%s"  <-- ISSUE' % pat); ok = False

# The line count is stated as a floor, not a figure, because every edit to any script moves it
# and a gate that fails on every edit is a gate somebody deletes. The floor still has to be true,
# and it has to be close enough to be worth stating.
m = re.search(r"more than (\d+) lines", b)
if not m:
    print('   the brief does not state the line count as "more than N lines"  <-- ISSUE'); ok = False
elif not (int(m.group(1)) <= nline < int(m.group(1)) + 1000):
    print("   the brief says more than %s lines; the tree has %d  <-- ISSUE" % (m.group(1), nline))
    ok = False

# every physics number in the brief must be in a manuscript
SKIP = {str(v) for v, _ in tree} | {str(nline), m.group(1) if m else "", "2026", "30", "1995"}
nums = set(re.findall(r"\d+\.\d+|\b\d{2,}\b", b))
stray = sorted(n for n in nums - SKIP if n not in papers)
if stray:
    print("   numbers in the brief that appear in neither manuscript: %s  <-- ISSUE" % ", ".join(stray))
    ok = False
if ok:
    print("   %d tree counts stated correctly, and every physics number is in a manuscript" % len(tree))

# the plants, both halves
assert "270.30" not in papers, "plant: the withdrawn witness figure is back in a manuscript"
assert "270.30" not in nums, "plant: the withdrawn witness figure is back in the brief"
assert ("%d calculation files" % (ncalc + 1)) not in b, "plant: a wrong tree count would pass"
print("   plant: a number the papers dropped is caught: yes")
print("   plant: a tree count off by one is caught: yes")
sys.exit(0 if ok else 1)
PY18
[ $? -ne 0 ] && fail=$((fail+1))
fi

hdr "21. no body sentence runs past the word cap"
# Nothing measured sentence length and the cosmology paper carried one of 113 words in its
# body, a four-item list run together with semicolons and a parenthetical inside the first
# item. Four of the worst were split on 2026-09-27, taking the maximum to 69 against the
# companion's 79. The cap is 80, which stops the next runaway without forcing a rewrite of
# nineteen sentences that are long because they carry qualifiers a referee wants.
python3 checks/sentence_length.py --selftest || fail=$((fail+1))
python3 checks/sentence_length.py papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md || fail=$((fail+1))
# The word cap is coarse: a long sentence broken by semicolons and colons gives the reader
# somewhere to put the first half down. What forces a re-read is a long run of clauses with no
# hard break in it, so that run is measured and capped too.
python3 checks/clause_stretch.py --selftest || fail=$((fail+1))
python3 checks/clause_stretch.py papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md || fail=$((fail+1))

hdr "20. the release describes its own dependencies correctly"
# Both manuscripts told a reader the code needs "base R and base Python with nothing imported
# beyond the standard library". Nineteen of the 22 Python calculations import numpy, scipy, sympy
# or mpmath, and gate 7 had been passing because it quietly prefers .venv/bin/python. A referee who
# tries to run the code is exactly the reader who finds that out. Gate 12 counts files; this one
# reads what the files actually import and what the manuscripts actually promise.
python3 - papers/1_separate_ways/PAPER2_v4_draft.md papers/2_over_the_horizon/COMPANION_v1.md <<'PY20'
import io, re, sys, glob
flat = lambda p: " ".join(io.open(p, encoding="utf-8").read().split())
docs = {p.split("/")[-1]: flat(p) for p in sys.argv[1:]}
THIRD = {"numpy", "scipy", "sympy", "mpmath", "pandas", "matplotlib", "astropy", "torch", "healpy"}

pys = sorted(glob.glob("checks/calc/*.py"))
used, stdlib_only = set(), 0
for f in pys:
    mods = set()
    for l in io.open(f, encoding="utf-8"):
        m = re.match(r"\s*(?:import|from)\s+([A-Za-z_][\w.]*)", l)
        if m:
            top = m.group(1).split(".")[0]
            if top in THIRD: mods.add(top)
    used |= mods
    if not mods: stdlib_only += 1
nR    = len(glob.glob("checks/calc/*.R"))
nfig  = len(glob.glob("checks/fig_*.R"))
ok = True

for name, t in docs.items():
    if "nothing imported beyond the standard library" in t:
        print("   %s still promises the standard library alone  <-- ISSUE" % name); ok = False
    for m in sorted(used):
        if m not in t:
            print("   %s does not name %s, which %d script(s) import  <-- ISSUE"
                  % (name, m, sum(1 for f in pys
                                  if re.search(r"(?m)^\s*(?:import|from)\s+%s\b" % m,
                                               io.open(f, encoding="utf-8").read()))))
            ok = False
cos = docs["PAPER2_v4_draft.md"]
# The split itself, checked in EACH document against the tree rather than in one of them.
for name, t in docs.items():
    if not re.search(r"\b%d\b[^.]{0,90}\b%d\b[^.]{0,60}\bR\b" % (nR, nfig), t) \
       and not re.search(r"\b%d\b[^.]{0,90}\b%d\b[^.]{0,60}\bR\b" % (nfig, nR), t):
        print("   %s does not state the R split as %d calculations and %d generators  <-- ISSUE"
              % (name, nR, nfig)); ok = False
    if not re.search(r"\b%d\b[^.]{0,60}Python" % len(pys), t):
        print("   %s does not say %d calculations are Python  <-- ISSUE" % (name, len(pys))); ok = False
    if not re.search(r"\b(?:%d|three)\b[^.]{0,80}standard\s+library" % stdlib_only, t, re.I):
        print("   %s does not say how many need only the standard library (%d)  <-- ISSUE"
              % (name, stdlib_only)); ok = False
if ok:
    print("   %d Python calculations, %d needing only the standard library; the packages the rest"
          % (len(pys), stdlib_only))
    print("   import (%s) are all named in both manuscripts" % ", ".join(sorted(used)))

# Every venue in pub/VENUES.md treats undisclosed AI use as the disqualifying thing, and the
# companion had no acknowledgements section at all, so it carried no statement. Each manuscript
# is a separate deposit and needs its own.
for name, t in docs.items():
    if "large language models" not in t:
        print("   %s carries no AI-use statement  <-- ISSUE" % name); ok = False
    elif "responsible for every claim" not in t:
        print("   %s discloses the tools but does not take responsibility  <-- ISSUE" % name)
        ok = False
    else:
        print("   %-22s carries an AI-use statement" % name)

assert "nothing imported beyond the standard library" not in cos, \
       "plant: the retired promise would not be caught"
assert used, "plant: the import scan found nothing, so the gate would pass on an empty set"
print("   plant: the retired standard-library promise is caught, and the scan is not empty: yes")
sys.exit(0 if ok else 1)
PY20
[ $? -ne 0 ] && fail=$((fail+1))

hdr "19. every label fits inside its panel"
# Gate 11 measures type SIZE and passed every figure. R clips text at the plot region and
# says nothing, so width went unmeasured until the companion was built and page 20 carried
# an annotation cut off mid-word. Two figures were affected.
Rscript checks/label_fit_check.R --selftest || fail=$((fail+1))
Rscript checks/label_fit_check.R checks/fig_*.R 2>/dev/null || fail=$((fail+1))

hdr "23. every term of art is glossed where the reader first meets it"
# Ben, 2026-09-27: the paper has to read exceptionally well and get through peer review, so a
# term of art must arrive with its plain-English gloss the first time a reader sees it. This
# pins each gloss to the first use, so a later edit that adds an earlier bare use is caught.
python3 checks/first_use_gloss.py || fail=$((fail+1))

hdr "24. Appendix E carries the same lines as the file that computes the ledger"
# The paper's headline structural claim is a count, nineteen out and four in, and a referee
# could see the nineteen only as compressed labels on a figure. Appendix E lists them at full
# wording; this checks membership both ways so the hand-written list cannot drift from the
# file that produces the count.
python3 checks/ledger_appendix_check.py || fail=$((fail+1))

hdr "25. no label is printed over by its own figure's ink"
# Gate 19 measures a label against its PANEL and catches text running off the edge. It says
# nothing about what is already inside. The companion was read as a PDF on 2026-09-28 and
# Figure 5's "double-precision floor" was printed over by the bars at l = 5, 7 and 9: inside
# the panel, correctly sized, unreadable. Sixteen labels across nine generators were in that
# state and gate 19 passed all of them.
Rscript checks/label_ink_check.R checks/fig_*.R 2>/dev/null || fail=$((fail+1))

hdr "26. the PRL Letter fits, and quotes nothing the paper does not"
# The Letter is carved out of section 3.6 for a venue with a hard length limit, and it is a
# second place for a number to live. Two ways that goes wrong: it grows past 3750 words while
# nobody is counting equations and figures against the cap, and a number gets retyped a digit
# short of the manuscript it came from. Writing it, 1.417e-32 s had already become 1.4e-32.
if [ ! -f pub/paper2/LETTER_PRL_v1.md ]; then
  echo "   not in this repository (duplicates paper 1); gate does not apply here"
else
LETTER=pub/paper2/LETTER_PRL_v1.md
python3 checks/letter_length_check.py || fail=$((fail+1))
python3 checks/letter_numbers_check.py || fail=$((fail+1))
lt=$(python3 "$BENLM/tools/llm_tics.py" "$LETTER" 2>&1 | grep -oE 'total hits: [0-9]+' | grep -oE '[0-9]+')
echo "   machine-prose tics: ${lt:-?} (0 allowed)"
[ "${lt:-9}" -gt 0 ] && { python3 "$BENLM/tools/llm_tics.py" "$LETTER" | tail -12; fail=$((fail+1)); }
lr=$(python3 "$BENLM/tools/opener_runs.py" "$LETTER" 2>&1 | grep -oE 'same word: [0-9]+' | grep -oE '[0-9]+')
echo "   runs of 3+ sentences opening on the same word: ${lr:-?} (0 allowed)"
[ "${lr:-9}" -gt 0 ] && { python3 "$BENLM/tools/opener_runs.py" "$LETTER" | sed -n '3,14p'; fail=$((fail+1)); }
python3 checks/repeat_check.py "$LETTER" || fail=$((fail+1))
fi

hdr "27. the paper does not narrate its own rhetoric, and every figure is cited first"
# Ben, 2026-09-29, on finding "One separation is worth stating before that list, because a referee
# will want it." still in 4.1 after a paragraph-by-paragraph pass: that sentence alone gets a paper
# desk-rejected. It was not one slip. Forty-one sentences across the two manuscripts addressed a
# referee, told the reader how to read, or explained why the paper is arranged as it is. Nothing
# measured any of it. The same pass found three of five figures never cited before they appeared,
# two never cited at all.
python3 checks/editorial_voice.py --selftest || fail=$((fail+1))
python3 checks/editorial_voice.py || fail=$((fail+1))
python3 checks/figure_order_check.py --selftest || fail=$((fail+1))
python3 checks/figure_order_check.py || fail=$((fail+1))

hdr "28. the exported repository is not behind this one"
# Ben, 2026-09-29: "The PDF in seperate_ways says it hasn't been updated since yesterday and it
# doesn't appear to have the refinements in it." It was a 49-page build from the previous
# afternoon against a working tree at 47, so a night's work was missing from the repository both
# manuscripts print on their own pages. Nothing watched it, the same gap that left the Zenodo
# record two versions behind. This is the half that can be closed mechanically.
python3 checks/export_freshness.py --selftest || fail=$((fail+1))
python3 checks/export_freshness.py || fail=$((fail+1))

hdr "13. the checkers can still fail"
python3 checks/calc/xref_check.py --validate-external >/dev/null 2>&1 \
  && echo "  xref_check --external: both directions catch their plant" \
  || { echo "  xref_check --external: A PLANTED CASE WAS MISSED"; fail=$((fail+1)); }
for v in checks/calc/abstract_len.py checks/calc/xref_check.py checks/calc/citation_check.py; do
  out=$(python3 "$v" --validate 2>&1 | tail -1)
  echo "   $(basename "$v"): $out"
done
python3 "$BENLM/tools/llm_tics.py" --validate >/dev/null 2>&1 && echo "   llm_tics: every branch fired" || { echo "   llm_tics: A BRANCH IS BLIND"; fail=$((fail+1)); }
python3 "$BENLM/tools/fragment_check.py" --validate >/dev/null 2>&1 && echo "   fragment_check: catches its plants and spares its exceptions" || { echo "   fragment_check: A PLANTED CASE MISBEHAVED"; fail=$((fail+1)); }
python3 "$BENLM/tools/opener_runs.py" --validate >/dev/null 2>&1 && echo "   opener_runs: catches runs, density and abstract subjects; spares the clean control" || { echo "   opener_runs: A PLANTED CASE MISBEHAVED"; fail=$((fail+1)); }
python3 checks/repeat_check.py "$PAPER" --selftest || fail=$((fail+1))
# Gate 4 is the gate that matters most and nothing validated it until now.
python3 checks/claims_check.py "$PAPER" --selftest || fail=$((fail+1))
# The contrastive bar is a threshold rather than a script, so its self-test is the
# comparison itself: it must spare a rate under the bar and catch one over it.
ticok=1
awk 'BEGIN { exit !(3.9999 > 4.0) }' && ticok=0          # must NOT fire
awk 'BEGIN { exit !(4.0001 > 4.0) }' || ticok=0          # must fire
[ "$ticok" = 1 ] && echo "  contrastive bar: spares 3.9999 and catches 4.0001" \
                 || { echo "  contrastive bar: THE THRESHOLD IS INERT"; fail=$((fail+1)); }
# Gate 16 by planting a drift in the real file, checking it is caught, then restoring it.
metasrc=papers/1_separate_ways/ARXIV_METADATA.txt
metatmp=$(mktemp); cp "$metasrc" "$metatmp"
sed -i 's/pages, [0-9]* figures\./pages, 99 figures./' "$metasrc"
python3 checks/sync_metadata_paper2.py --check >/dev/null 2>&1 && metacaught=0 || metacaught=1
cp "$metatmp" "$metasrc"; rm -f "$metatmp"
python3 checks/sync_metadata_paper2.py --check >/dev/null 2>&1 && metareal=1 || metareal=0
# These two are different failures and the old message called both of them INERT. On 2026-09-29
# the abstract was rewritten, the metadata legitimately drifted, and the suite reported that a
# working checker had gone inert. Say which one happened.
if [ "$metacaught" = 1 ] && [ "$metareal" = 1 ]; then
  echo "  paper2 metadata: catches a planted drift and spares the real file"
elif [ "$metacaught" != 1 ]; then
  echo "  paper2 metadata: THE CHECK IS INERT, it did not catch a planted 99-figure drift"; fail=$((fail+1))
else
  echo "  paper2 metadata: the checker works, but ARXIV_METADATA.txt HAS DRIFTED from the manuscript"
  echo "     run: python3 checks/sync_metadata_paper2.py"; fail=$((fail+1))
fi
python3 checks/structure_check.py "$PAPER" --selftest || fail=$((fail+1))
python3 checks/first_use_gloss.py --selftest || fail=$((fail+1))
python3 checks/ledger_appendix_check.py --selftest || fail=$((fail+1))
Rscript checks/label_ink_check.R --selftest || fail=$((fail+1))
# Gate 26 validates itself when it runs: both Letter checkers print their plants above,
# and both exit non-zero if one fails to fire, so nothing more is needed here.

printf '\n'
[ "$fail" != "$_prev_fail" ] && { printf '\033[1m   ^^ GATE RAISED THE FAILURE: %s\033[0m\n' "$_gate"; _failed="$_failed
   $_gate"; }
if [ $fail -eq 0 ]; then echo "ALL GATES PASSED"
# $fail counts FAILURES, not gates: one gate can increment it once per offending file, so the
# first version of this line read "THE 5 GATE(S)" above a list of two. Report both numbers.
else printf 'SOMETHING FAILED. %d failure(s), raised by these gate(s):%s\n' "$fail" "$_failed"
fi
exit $fail
