#!/usr/bin/env python3
"""The paper must not narrate its own rhetoric, and must never address a referee.

    python3 checks/editorial_voice.py [files...]
    python3 checks/editorial_voice.py --selftest

WHY. Ben, 2026-09-29, on finding "One separation is worth stating before that list, because a
referee will want it." still sitting in section 4.1 after a paragraph-by-paragraph pass: a
sentence like that gets a paper desk-rejected on its own. It is the author talking about the
submission process inside the submission. Thirteen sentences of the same family were in the
manuscript when this was written, so it is a pattern and not a slip, and nothing measured it.

Three families, in descending order of how fatal they are.

  1. ADDRESSING THE REFEREE. Any mention of a referee, a reviewer or what they will want. There
     is no acceptable use of this inside a paper. Zero tolerance.

  2. NARRATING THE RHETORIC. "is worth stating because", "the point of section 3.6 is",
     "so it is checked rather than repeated". The paper explaining why it is arranged as it is,
     instead of just being arranged that way. A very small number of these are load-bearing,
     for instance a sentence that says where a limit applies, so this family has an allowance
     rather than a ban, and the allowance is recorded here where it can be argued with.

  3. INSTRUCTING THE READER. "a reader should press", "the reader will want to know". Telling
     someone how to read the paper. Reserve for the one or two places where the paper genuinely
     hands over a test.

The allowance for families 2 and 3 is a number, not a judgement, so it can only go down.
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT = [os.path.join(ROOT, "pub", "paper2", "PAPER2_v4_draft.md"),
           os.path.join(ROOT, "pub", "paper2", "COMPANION_v1.md"),
           os.path.join(ROOT, "pub", "paper2", "LETTER_PRL_v1.md")]

REFEREE = re.compile(r"\b(referee|referees|reviewer|reviewers|desk[- ]reject|peer review)\b", re.I)

# "worth <act-of-writing>" is the tic. "worth" about a physical quantity is not, so the verb
# list is what does the work: stating, saying, noting, making, recording, following, naming,
# setting out, turning round, computing-for-the-reader's-sake. The first version of this listed
# only the four forms already found in the manuscript, reported zero, and left nineteen more in
# the paper and fifty-two in the companion. A checker built from the instances you have already
# seen is a checker that encodes your list instead of the class.
_ACTS = (r"stating|saying|noting|making|remarking|recording|naming|having|keeping|establishing"
         r"|setting out|spelling out|turning round|reading the other way|following|answering"
         r"|separating|checking|believing|computing|a sentence|a line|a number|the paragraph"
         r"|one more line|a paragraph")
RHETORIC = re.compile(
    r"\b(worth (?:%s)\b|worth it\b"
    r"|is cheaper than having|checked rather than repeated|deserves to be checked"
    r"|the point of (?:§|section )?\d|that is the point of it|the point is worth"
    r"|before that list|for completeness|it should be noted|it is worth noting"
    r"|as informative as what it does|bears repeating"
    r"|as noted above|as we shall see|in what follows we|we now turn to)\b" % _ACTS, re.I)

# The paper describes physics, not its own audience. Any mention of a reader is the tic.
READER = re.compile(r"\b(a|the|any|every|some)\s+readers?\b|\breaders\b|\bthe reader\b", re.I)

# Sentences that trip family 2 or 3 and are kept on purpose. Each needs a reason, and the list
# is the allowance: adding to it is a decision someone has to defend in review of this file.
ALLOWED = [
    # none yet; the 2026-09-29 pass removed every one it found
]

# 4. THE PAPER AS AN ACTOR. "the construction returns a region", "what it does say is",
#    "this Letter says nothing". The work is made the subject of a reporting verb, so the
#    sentence is about the document instead of about the physics. Added 2026-09-30 after
#    Ben found "the construction returns a region rather than a possibility ... None of that
#    was put in by hand" in a Letter this gate had just passed 0/0/0. Families 2 and 3 did
#    not cover it: nothing was being justified and no reader was addressed.
_SELF = r"(?:the construction|the argument|the analysis|the calculation|the paper|this letter)"
_REPORT = (r"says|shows|argues|claims|returns|tells|asks|answers|demonstrates|reaches|puts"
           r"|selects|closes|offers|gives|finds|concludes|establishes|proves")
ACTOR = re.compile(
    r"\b%s\s+(?:%s)\b"
    r"|\bwhat it does (?:say|show|claim)\b"
    # "nothing is put in by hand" is a claim about the physics and good writing. The tic is
    # the vague-referent closing version, where "that" stands for the whole paper.
    r"|\bnone of (?:that|this|it) (?:was|is) put in by hand\b"
    r"|\bshould be read as\b" % (_SELF, _REPORT), re.I)

CAP_RHETORIC = 0
CAP_READER = 0
CAP_ACTOR = 0


def sentences(text):
    body = text
    for head in ("## References", "## Acknowledgements", "## Code and data"):
        if head in body:
            body = body.split(head)[0]
    body = " ".join(body.split())
    return [x.strip() for x in re.split(r"(?<=[.!?]) +", body) if x.strip()]


def scan(path):
    text = io.open(path, encoding="utf-8").read()
    hits = {"referee": [], "rhetoric": [], "reader": [], "actor": []}
    for s in sentences(text):
        if any(a in s for a in ALLOWED):
            continue
        if REFEREE.search(s):
            hits["referee"].append(s)
        elif READER.search(s):
            hits["reader"].append(s)
        elif RHETORIC.search(s):
            hits["rhetoric"].append(s)
        if ACTOR.search(s):
            hits["actor"].append(s)
    return hits


def report(paths):
    bad = 0
    for p in paths:
        if not os.path.exists(p):
            continue
        h = scan(p)
        name = os.path.basename(p)
        print("   %-24s referee %d (0 allowed), rhetoric %d (%d), reader %d (%d)"
              % (name, len(h["referee"]), len(h["rhetoric"]), CAP_RHETORIC,
                 len(h["reader"]), CAP_READER))
        for kind, cap in (("referee", 0), ("rhetoric", CAP_RHETORIC), ("reader", CAP_READER),
                          ("actor", CAP_ACTOR)):
            if len(h[kind]) > cap:
                bad += 1
                for s in h[kind]:
                    print("      %-9s %s" % (kind, s[:150]))
    return bad


def _selftest():
    """Plant one of each family and check they all fire; spare a clean control."""
    ok = True
    cases = [
        ("referee", "A referee will want the separation stated before that list."),
        ("rhetoric", "The division of labour is worth stating, because it is easy to get backwards."),
        ("reader", "The residual phase is the falsifier a reader should press."),
        ("actor", "Asked where its halves can touch, the construction returns a region."),
        ("actor", "What it does say is that none of that was put in by hand."),
    ]
    for kind, sent in cases:
        h = {"referee": [], "rhetoric": [], "reader": [], "actor": []}
        if REFEREE.search(sent):
            h["referee"].append(sent)
        elif READER.search(sent):
            h["reader"].append(sent)
        elif RHETORIC.search(sent):
            h["rhetoric"].append(sent)
        if ACTOR.search(sent):
            h["actor"].append(sent)
        caught = len(h[kind]) == 1
        print("   plant: a %s sentence is caught: %s" % (kind, "yes" if caught else "NO"))
        ok = ok and caught
    clean = "The surface gravity cancels between the occupation and the WKB mode density."
    spared = not (REFEREE.search(clean) or READER.search(clean) or RHETORIC.search(clean)
                  or ACTOR.search(clean))
    print("   plant: an ordinary sentence is spared: %s" % ("yes" if spared else "NO"))
    return ok and spared


def main():
    if "--selftest" in sys.argv:
        sys.exit(0 if _selftest() else 1)
    paths = [a for a in sys.argv[1:] if not a.startswith("-")] or DEFAULT
    sys.exit(1 if report(paths) else 0)


if __name__ == "__main__":
    main()
