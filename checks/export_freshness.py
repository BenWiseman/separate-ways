#!/usr/bin/env python3
"""The exported repository must not be behind the working tree.

    python3 checks/export_freshness.py
    python3 checks/export_freshness.py --selftest

WHY. Ben, 2026-09-29: "The PDF in seperate_ways says it hasn't been updated since yesterday and
it doesn't appear to have the refinements in it." He was right. The export was a 49-page build
from the previous afternoon and the working tree was at 47 pages, so a whole night's work,
the introduction and abstract rewrite, the section renames, the editorial-voice pass and every
fix the blind panel produced, was missing from the repository both manuscripts print on their
own pages.

Nothing checked it, for the same reason nothing checked the Zenodo record: the export is a step
someone has to remember, and the suite only watched the things inside this repository. Two
outward-facing records have now gone stale in two days. This gate closes the one that can be
closed mechanically. The Zenodo record cannot be, since only Ben can deposit, so it is named in
the output rather than gated.

The check is a hash comparison on the sources, not on the PDFs: a PDF differs on every build
because of its timestamp, so comparing built artefacts would cry wolf every time.
"""

import hashlib
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DST = os.environ.get("SEPARATE_WAYS", os.path.expanduser("~/separate_ways"))

# source in this repository -> where the export puts it
# The PRL Letter is not on this list because it is deliberately not exported; see the
# DO_NOT_PUBLISH block in tools/export_to_separate_ways.sh.
PAIRS = [
    ("papers/1_separate_ways/PAPER2_v4_draft.md", "paper/PAPER2_v4_draft.md"),
    ("papers/2_over_the_horizon/COMPANION_v1.md", "paper/COMPANION_v1.md"),
    ("checks/CLAIMS.tsv", "checks/CLAIMS.tsv"),
    ("checks/check_all.sh", "checks/check_all.sh"),
]

# Files the export rewrites on the way across, so their bytes legitimately differ: it turns
# checks/calc/ into checks/calc/ and tangents/ into tangents/ inside every script, the
# pass and the claims table. For these the check is that the destination exists and is no older
# than the source. Comparing their bytes reports a difference that is the export working
# correctly, which is a gate that cries wolf and then gets ignored.
REWRITTEN = {"checks/check_all.sh", "checks/CLAIMS.tsv"}


def sha(path):
    return hashlib.sha256(io.open(path, "rb").read()).hexdigest()


def audit():
    """Return a list of (label, problem) for everything out of step."""
    bad = []
    if not os.path.isdir(os.path.join(DST, ".git")):
        return [(DST, "no exported repository here")]
    for src, dst in PAIRS:
        s = os.path.join(ROOT, src)
        d = os.path.join(DST, dst)
        if not os.path.exists(s):
            continue
        if not os.path.exists(d):
            bad.append((dst, "missing from the export"))
            continue
        if dst in REWRITTEN:
            if os.path.getmtime(d) < os.path.getmtime(s):
                bad.append((dst, "older than the source it is rewritten from"))
        elif sha(s) != sha(d):
            bad.append((dst, "differs from the working tree"))
    return bad


def main():
    if "--selftest" in sys.argv:
        # The gate has to fail on purpose. Point it at a directory that is not an export.
        global DST
        keep = DST
        DST = "/nonexistent-export-path"
        caught = len(audit()) == 1
        DST = keep
        print("   plant: a missing export is caught: %s" % ("yes" if caught else "NO"))
        real = audit()
        print("   plant: the real export passes: %s" % ("yes" if not real else "NO"))
        sys.exit(0 if (caught and not real) else 1)

    bad = audit()
    if bad:
        print("   the export at %s is BEHIND the working tree  <-- ISSUE" % DST)
        for label, why in bad:
            print("      %-34s %s" % (label, why))
        print("      run: bash tools/export_to_separate_ways.sh")
        sys.exit(1)
    print("   %s carries the current manuscripts, claims table and pass" % DST)
    print("   (the Zenodo record is not checkable from here; only Ben can deposit)")


if __name__ == "__main__":
    main()
