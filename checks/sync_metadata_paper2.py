#!/usr/bin/env python3
"""Keep the cosmology paper's arXiv metadata from drifting off the manuscript.

The companion has had `checks/sync_metadata.py` and a gate refusing to let its metadata drift
since the day the gate was written.  This paper had neither, and it drifted all the way: on
2026-09-24 `ARXIV_METADATA.txt` still carried the v3 title and a v3 abstract that differs from
the v4 draft's from the first character.  That is a submission blocker rather than an untidiness,
because the abstract in the metadata file is the one that would have been pasted into arXiv.

This script rewrites three things from the manuscript and touches nothing else:

    TITLE       the manuscript's H1
    ABSTRACT    the ## Abstract section, flattened to one paragraph
    COMMENTS    the numbered-figure count, which said five while the draft numbers three
    CHECKLIST   the count of images to embed, which is those three plus the graphical
                abstract; the two fields want different numbers and used to carry one
    the header  the character count and the date stamp on it

Everything Ben fills in by hand, the categories, the email, the licence and the endorsement
notes, is left exactly as it is.  The one line of the checklist this does touch states a count
taken from the manuscript, and a wrong number in a list somebody is going to tick is worse than
a number nobody maintains.  Run it after any edit to the
abstract or the title:

    python3 checks/sync_metadata_paper2.py

With --check it writes nothing and exits non-zero if it would have changed anything, which is what
the gate calls.  A sync tool nobody is obliged to run is how this file drifted in the first place.

It also refuses to write a metadata file arXiv would reject: the abstract field is capped at 1920
characters and the fields are ASCII only, so both are checked here rather than discovered in the
submission form.
"""
import io
import re
import sys

PAPER = "papers/1_separate_ways/PAPER2_v4_draft.md"
META = "papers/1_separate_ways/ARXIV_METADATA.txt"
CAP = 1920
STAMP = "2026-09-24"


def manuscript_fields(path):
    c = io.open(path, encoding="utf-8").read()
    # Two different counts, because the two fields ask different questions. arXiv's
    # Comments field means the figures a reader can cite, which are the numbered ones.
    # The embed checklist means every image that has to make it into the PDF, which is
    # those plus the unnumbered graphical abstract.
    nfig = len(re.findall(r"(?m)^\*\*Figure\s+\d+\.", c))
    nimg = len(re.findall(r"(?m)^!\[", c))
    title = c.splitlines()[0].lstrip("# ").strip()
    m = re.search(r"(?m)^##\s+Abstract\s*$", c)
    if not m:
        sys.exit(f"  {path}: no '## Abstract' heading found")
    rest = c[m.end():]
    end = re.search(r"(?m)^#{1,6}\s+\S", rest)
    abstract = " ".join(rest[:end.start()].split())
    return title, abstract, nfig, nimg


def check(title, abstract):
    problems = []
    if len(abstract) > CAP:
        problems.append(f"abstract is {len(abstract)} characters against arXiv's cap of {CAP}")
    for name, field in (("title", title), ("abstract", abstract)):
        bad = sorted({ch for ch in field if ord(ch) > 127})
        if bad:
            problems.append(f"{name} carries non-ASCII characters {bad}; arXiv accepts ASCII only")
    return problems


def rewrite(meta_text, title, abstract, nfig=None, nimg=None):
    t = meta_text
    t, n_title = re.subn(r"(?ms)(^TITLE\n).*?(\n\n)", lambda g: g.group(1) + title + g.group(2),
                         t, count=1)
    t, n_abs = re.subn(r"(?ms)(^ABSTRACT[^\n]*\n).*?(\n\nREPORT-NO)",
                       lambda g: g.group(1) + abstract + g.group(2), t, count=1)
    n_fig = 1
    if nfig is not None:
        t, n_fig = re.subn(r"(?m)^(\[PDF page count[^\n]*?pages, )\d+( figures\.)",
                           lambda g: g.group(1) + str(nfig) + g.group(2), t, count=1)
    n_chk = 1
    if nimg is not None:
        words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
                 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
        t, n_chk = re.subn(r"(all )(?:one|two|three|four|five|six|seven|eight|nine|ten|\d+)"
                           r"( images embedded)",
                           lambda g: g.group(1) + words.get(nimg, str(nimg)) + g.group(2),
                           t, count=1)
    t, n_hdr = re.subn(r"(?m)^# Current TeX-form paste paragraph: [^\n]*$",
                       f"# Current TeX-form paste paragraph: {len(abstract)} ASCII characters "
                       f"(sync_metadata_paper2.py, {STAMP});", t, count=1)
    return t, (n_title, n_abs, n_fig, n_chk, n_hdr)


if __name__ == "__main__":
    check_only = "--check" in sys.argv
    title, abstract, nfig, nimg = manuscript_fields(PAPER)
    problems = check(title, abstract)
    if problems:
        for p in problems:
            print(f"  REFUSING: {p}")
        sys.exit(1)

    before = io.open(META, encoding="utf-8").read()
    after, counts = rewrite(before, title, abstract, nfig, nimg)
    if 0 in counts:
        sys.exit(f"  {META}: a section did not match and nothing was written "
                 f"(title {counts[0]}, abstract {counts[1]}, figures {counts[2]}, "
                 f"checklist {counts[3]}, header {counts[4]})")
    if check_only:
        if after != before:
            print(f"  {META} HAS DRIFTED from {PAPER}  <-- ISSUE")
            print(f"  run: python3 {__file__.split('/')[-1]} (without --check) to resync")
            sys.exit(1)
        print(f"  {META}: current with the manuscript "
              f"(abstract {len(abstract)} of {CAP} chars, {nfig} figures, {nimg} images)")
        sys.exit(0)

    io.open(META, "w", encoding="utf-8").write(after)
    print(f"  {META}: abstract {len(abstract)} chars of {CAP}, title {len(title)} chars, "
          f"{nfig} numbered figures, {nimg} images"
          + ("" if after != before else "   (already current)"))
    print(f"  title: {title}")
