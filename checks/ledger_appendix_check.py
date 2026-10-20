#!/usr/bin/env python3
"""Appendix E must carry the same lines as the file that computes the ledger.

    python3 checks/ledger_appendix_check.py
    python3 checks/ledger_appendix_check.py --selftest

WHY. Ben, 2026-09-27: the paper has to get through peer review. Its headline structural claim is
a count, nineteen out and four in, and until now a referee could see the nineteen only as
compressed labels on a figure. Appendix E lists them at full wording. A list written by hand can
drift from the file that produces the count, so this checks membership both ways: every line of
relativity_ledger.R has an entry, every entry has a line, and the two counts agree.

Matching is by content-word overlap rather than by exact string, because the appendix is typeset
and the ledger is plain text. That is looser than an equality test and is why the counts are
checked as well: a silently dropped line fails on the count even if some other entry would have
matched it.

2026-10-06. The ledger now files each derived line under where it stood before this paper:
"Jacobson 1995", "standard" or "fold". Both the appendix and the Letter list the lines under the
same three subheadings, and a line must sit under the subheading its ledger group names. That is
what stops the paper drifting back to presenting Jacobson's lines as the fold's. A subheading the
checker cannot find aborts the run rather than leaving its lines unchecked.
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LEDGER = os.path.join(HERE, "calc", "relativity_ledger.R")
PAPER = os.path.join(ROOT, "papers/1_separate_ways/PAPER2_v4_draft.md")

STOP = set("the a an of at to in and or is are be been that which what with for from its it "
           "on by as this those these every no not so does do than then there their them".split())


def words(s):
    s = s.lower()
    s = re.sub(r"\$[^$]*\$", " ", s)           # typeset maths is not a content word
    s = re.sub(r"\\[a-z]+", " ", s)
    s = re.sub(r"[^a-z ]", " ", s)
    return {w for w in s.split() if len(w) > 3 and w not in STOP}


GROUPS = ["Jacobson 1995", "standard", "fold"]
# The subheading each group sits under, in the appendix and in the Letter alike. Anchored to a
# whole line so a mention of the phrase in running prose cannot be taken for a heading.
HEADS = {"Jacobson 1995": r"(?m)^\*Already in Jacobson's 1995 argument\.\*\s*$",
         "standard":      r"(?m)^\*Standard results\.\*\s*$",
         "fold":          r"(?m)^\*Added by the fold\.\*\s*$"}


def ledger_lines():
    src = io.open(LEDGER, encoding="utf-8").read()
    out = {"derived": [], "assumed": []}
    order = []
    for m in re.finditer(r'item\("((?:[^"\\]|\\.)*)",\s*\n?\s*"(derived|assumed)"', src):
        out[m.group(2)].append(m.group(1))
        order.append(m.group(2))
    g = re.search(r"GROUP <- c\((.*?)\)\n", src, re.S)
    if not g:
        raise SystemExit("   no GROUP vector in %s  <-- ISSUE" % LEDGER)
    groups = re.findall(r'"([^"]+)"', re.sub(r"#[^\n]*", "", g.group(1)))
    if len(groups) != len(order):
        raise SystemExit("   GROUP has %d entries for %d ledger lines  <-- ISSUE" % (len(groups), len(order)))
    out["by_group"] = {k: [] for k in GROUPS}
    i_der = 0
    for st, gr in zip(order, groups):
        if st == "derived":
            if gr not in GROUPS:
                raise SystemExit("   derived line filed under %r  <-- ISSUE" % gr)
            out["by_group"][gr].append(out["derived"][i_der]); i_der += 1
    return out


def split_groups(chunk, where):
    """Cut a region into the three groups by their subheadings, in order, or abort."""
    pos = []
    for k in GROUPS:
        m = re.search(HEADS[k], chunk)
        if not m:
            raise SystemExit("   %s: no subheading for the %r group (pattern %s)  <-- ISSUE"
                             % (where, k, HEADS[k]))
        pos.append((m.start(), m.end(), k))
    pos.sort()
    if [p[2] for p in pos] != GROUPS:
        raise SystemExit("   %s: the groups are out of order  <-- ISSUE" % where)
    out = {}
    for i, (a, b, k) in enumerate(pos):
        out[k] = chunk[b:pos[i + 1][0] if i + 1 < len(pos) else len(chunk)]
    return out


def appendix_entries(text=None):
    """Appendix E's derived lines, by group, and its assumed lines."""
    src = text if text is not None else io.open(PAPER, encoding="utf-8").read()
    if "## Appendix E." not in src:
        raise SystemExit("   Appendix E is not in the manuscript  <-- ISSUE")
    body = src[src.index("## Appendix E."):]
    body = body[:body.index("\n## References")] if "\n## References" in body else body
    e1 = body.index("**E.1"); e2 = body.index("**E.2")
    def items(chunk):
        return [re.sub(r"\s+", " ", m.group(1)).strip()
                for m in re.finditer(r"(?m)^\s*\d{1,2}\.\s+(\*\*.*?)(?=\n\s*\n|\Z)",
                                     chunk, flags=re.S)]
    groups = split_groups(body[e1:e2], "Appendix E")
    return {k: items(v) for k, v in groups.items()}, items(body[e2:])


LETTER = os.path.join(ROOT, "papers/4_einstein_from_cpt/LETTER_CPT_Einstein_equation_of_state_v3.md")
L_OUT, L_IN = "**Nineteen that come out.**", "**Four that go in.**"


def letter_entries(text=None):
    """The same nineteen and four, as the Letter's appendix tables carry them, by group.

    A third copy of the count is a third place for it to go stale, and the Letter is the copy a
    referee reads first. Both columns are taken as the entry text, since the mechanism column
    carries most of the content words.
    """
    src = text if text is not None else io.open(LETTER, encoding="utf-8").read()
    for h in (L_OUT, L_IN):
        if h not in src:
            raise SystemExit("   the Letter has no %s heading; if it was renamed, edit L_OUT/L_IN  <-- ISSUE" % h)
    body = src[src.index(L_OUT):]
    cut = body.index(L_IN)
    def rows(chunk):
        out = []
        for line in chunk.splitlines():
            m = re.match(r"^\|\s*\d{1,2}\s*\|(.+)\|(.+)\|\s*$", line)
            if m:
                out.append((m.group(1) + " " + m.group(2)).strip())
        return out
    groups = split_groups(body[:cut], "the Letter")
    return {k: rows(v) for k, v in groups.items()}, rows(body[cut:])


def compare(led, got, label, quiet=False):
    bad = 0
    if len(led) != len(got):
        print("   %s: the ledger has %d lines and the appendix %d  <-- ISSUE"
              % (label, len(led), len(got)))
        bad += 1
    used = set()
    for name in led:
        lw = words(name)
        best, bi = 0.0, -1
        for i, entry in enumerate(got):
            if i in used:
                continue
            ew = words(entry)
            frac = len(lw & ew) / max(1, len(lw))
            if frac > best:
                best, bi = frac, i
        if best < 0.5:
            print('   %s: no entry matches "%s" (best overlap %.2f)  <-- ISSUE'
                  % (label, name[:60], best))
            bad += 1
        else:
            used.add(bi)
            if not quiet:
                print("   %-8s %.2f  %s" % (label, best, name[:68]))
    spare = [g for i, g in enumerate(got) if i not in used]
    for g in spare:
        print("   %s: an entry matches no ledger line: %s  <-- ISSUE" % (label, g[:60]))
        bad += 1
    return bad


def compare_groups(led, got, label, quiet=True):
    bad = 0
    for k in GROUPS:
        bad += compare(led["by_group"][k], got[k], "%s/%s" % (label, k), quiet=quiet)
    return bad


def main():
    led = ledger_lines()
    d, a = appendix_entries()
    ld, la = letter_entries()
    sizes = {k: len(led["by_group"][k]) for k in GROUPS}
    print("   ledger groups: %s" % ", ".join("%s %d" % (k, sizes[k]) for k in GROUPS))
    if sum(sizes.values()) != len(led["derived"]) or not all(sizes.values()):
        print("   the groups do not cover the derived lines  <-- ISSUE"); return 1
    if "--selftest" in sys.argv:
        ok = True
        src = io.open(PAPER, encoding="utf-8").read()
        # (1) drop one entry: the count and the membership must both complain
        cut = re.sub(r"(?m)^10\. \*\*Black-hole thermodynamics\.\*\*.*?\n\n", "", src, flags=re.S)
        if cut == src:
            print("   plant (1) did not apply: the entry it removes has moved"); return 1
        dd, aa = appendix_entries(cut)
        if compare_groups(led, dd, "derived") == 0:
            print("   plant: a dropped entry is NOT caught"); ok = False
        else:
            print("   plant: a dropped entry is caught: yes")
        # (2) an entry that matches nothing in the ledger
        extra = dict(d); extra["fold"] = d["fold"] + ["**A twentieth line about nothing.** Invented for the self-test."]
        if compare_groups(led, extra, "derived") == 0:
            print("   plant: a spurious entry is NOT caught"); ok = False
        else:
            print("   plant: a spurious entry is caught: yes")
        # (3) a Jacobson line moved under the fold's subheading: the overstatement this guards
        moved = {k: list(v) for k, v in d.items()}
        x = moved["Jacobson 1995"].pop(1); moved["fold"].append(x)
        if compare_groups(led, moved, "derived") == 0:
            print("   plant: a Jacobson line filed under the fold is NOT caught"); ok = False
        else:
            print("   plant: a Jacobson line filed under the fold is caught: yes")
        # (4) a missing subheading must abort, not return an empty group
        nohead = re.sub(HEADS["standard"], "", src)
        try:
            appendix_entries(nohead); print("   plant: a missing subheading is NOT caught"); ok = False
        except SystemExit:
            print("   plant: a missing subheading aborts: yes")
        # (5) the Letter's own tables have to be caught when they drift
        lsrc = io.open(LETTER, encoding="utf-8").read()
        lcut = re.sub(r"(?m)^\| 10 \|.*\n", "", lsrc)
        if lcut == lsrc:
            print("   plant (5) did not apply: the Letter has no row 10"); return 1
        lld, lla = letter_entries(lcut)
        if compare_groups(led, lld, "letter-derived") == 0:
            print("   plant: a row dropped from the Letter's table is NOT caught"); ok = False
        else:
            print("   plant: a row dropped from the Letter's table is caught: yes")
        # (6) the real pair must pass
        if compare_groups(led, d, "derived") == 0 and \
           compare(led["assumed"], a, "assumed", quiet=True) == 0 and \
           compare_groups(led, ld, "letter-derived") == 0 and \
           compare(led["assumed"], la, "letter-assumed", quiet=True) == 0:
            print("   the real appendix and the real Letter both pass: yes")
        else:
            print("   the real appendix FAILS"); ok = False
        return 0 if ok else 1
    bad = compare_groups(led, d, "derived")
    bad += compare(led["assumed"], a, "assumed", quiet=True)
    bad += compare_groups(led, ld, "letter-derived")
    bad += compare(led["assumed"], la, "letter-assumed", quiet=True)
    print("   %d derived lines in %d groups and %d assumed lines, all in Appendix E and the Letter "
          "under the right subheadings and nothing else there"
          % (len(led["derived"]), len(GROUPS), len(led["assumed"])) if bad == 0 else
          "   %d mismatch(es) between the ledger and Appendix E or the Letter" % bad)
    return 1 if bad else 0


sys.exit(main())
