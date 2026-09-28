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
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LEDGER = os.path.join(HERE, "calc", "relativity_ledger.R")
PAPER = os.path.join(ROOT, "pub", "paper2", "PAPER2_v4_draft.md")

STOP = set("the a an of at to in and or is are be been that which what with for from its it "
           "on by as this those these every no not so does do than then there their them".split())


def words(s):
    s = s.lower()
    s = re.sub(r"\$[^$]*\$", " ", s)           # typeset maths is not a content word
    s = re.sub(r"\\[a-z]+", " ", s)
    s = re.sub(r"[^a-z ]", " ", s)
    return {w for w in s.split() if len(w) > 3 and w not in STOP}


def ledger_lines():
    src = io.open(LEDGER, encoding="utf-8").read()
    out = {"derived": [], "assumed": []}
    for m in re.finditer(r'item\("((?:[^"\\]|\\.)*)",\s*\n?\s*"(derived|assumed)"', src):
        out[m.group(2)].append(m.group(1))
    return out


def appendix_entries(text=None):
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
    return items(body[e1:e2]), items(body[e2:])


LETTER = os.path.join(ROOT, "pub", "paper2", "LETTER_PRL_v1.md")


def letter_entries(text=None):
    """The same nineteen and four, as the Letter's End Matter tables carry them.

    A third copy of the count is a third place for it to go stale, and the Letter is the
    copy a PRL referee reads first. Both columns are taken as the entry text, since the
    mechanism column carries most of the content words.
    """
    src = text if text is not None else io.open(LETTER, encoding="utf-8").read()
    if "**The nineteen that come out.**" not in src:
        raise SystemExit("   the Letter has no ledger table  <-- ISSUE")
    body = src[src.index("**The nineteen that come out.**"):]
    cut = body.index("**The four that go in.**")
    def rows(chunk):
        out = []
        for line in chunk.splitlines():
            m = re.match(r"^\|\s*\d{1,2}\s*\|(.+)\|(.+)\|\s*$", line)
            if m:
                out.append((m.group(1) + " " + m.group(2)).strip())
        return out
    return rows(body[:cut]), rows(body[cut:])


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


def main():
    led = ledger_lines()
    d, a = appendix_entries()
    if "--selftest" in sys.argv:
        ok = True
        src = io.open(PAPER, encoding="utf-8").read()
        # (1) drop one entry: the count and the membership must both complain
        cut = re.sub(r"(?m)^13\. \*\*Black-hole thermodynamics\.\*\*.*?\n\n", "", src, flags=re.S)
        dd, aa = appendix_entries(cut)
        if compare(led["derived"], dd, "derived", quiet=True) == 0:
            print("   plant: a dropped entry is NOT caught"); ok = False
        else:
            print("   plant: a dropped entry is caught: yes")
        # (2) an entry that matches nothing in the ledger
        extra = dd + ["**A twentieth line about nothing.** Invented for the self-test."]
        if compare(led["derived"], extra, "derived", quiet=True) == 0:
            print("   plant: a spurious entry is NOT caught"); ok = False
        else:
            print("   plant: a spurious entry is caught: yes")
        # (2b) the Letter's own tables have to be caught when they drift
        lsrc = io.open(LETTER, encoding="utf-8").read()
        lcut = re.sub(r"(?m)^\| 13 \|.*\n", "", lsrc)
        lld, lla = letter_entries(lcut)
        if compare(led["derived"], lld, "letter-derived", quiet=True) == 0:
            print("   plant: a row dropped from the Letter's table is NOT caught"); ok = False
        else:
            print("   plant: a row dropped from the Letter's table is caught: yes")
        # (3) the real pair must pass
        rld, rla = letter_entries()
        if compare(led["derived"], d, "derived", quiet=True) == 0 and \
           compare(led["assumed"], a, "assumed", quiet=True) == 0 and \
           compare(led["derived"], rld, "letter-derived", quiet=True) == 0 and \
           compare(led["assumed"], rla, "letter-assumed", quiet=True) == 0:
            print("   the real appendix and the real Letter both pass: yes")
        else:
            print("   the real appendix FAILS"); ok = False
        return 0 if ok else 1
    bad = compare(led["derived"], d, "derived", quiet=True)
    bad += compare(led["assumed"], a, "assumed", quiet=True)
    ld, la = letter_entries()
    bad += compare(led["derived"], ld, "letter-derived", quiet=True)
    bad += compare(led["assumed"], la, "letter-assumed", quiet=True)
    print("   %d derived and %d assumed lines, all of them in Appendix E and nothing else there"
          % (len(led["derived"]), len(led["assumed"])) if bad == 0 else
          "   %d mismatch(es) between the ledger and Appendix E" % bad)
    return 1 if bad else 0


sys.exit(main())
