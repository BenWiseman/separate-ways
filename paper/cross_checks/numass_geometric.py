#!/usr/bin/env python3
"""Report the exploratory-fit source/input inventory; no numerical fit is run.

The fifteen R scripts count arithmetic comparisons. The geometry-only profiles
are separate numerical fits; see WORKING_CALCULATIONS.md for scope and inputs.
"""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[3]
DATA = Path(os.environ.get("PAPER2_PANTHEON_DIR", ROOT / "data/pantheon_plus"))

def main():
    print("Geometry-only profile provenance (existence checks, not a fit)")
    for relative in ["calc/tangents/numass_mirror/numass_mirror.py",
                     "calc/tangents/numass_mirror/numass_mirror_output.txt",
                     "calc/desi_dr2_geometry.py", "WORKING_CALCULATIONS.md",
                     "data/README.md"]:
        print(f"[{'OK' if (ROOT / relative).is_file() else 'MISSING'}] {relative}")
    for name in ["PantheonPlusSH0ES.dat", "PantheonPlusSH0ES_STAT+SYS.cov"]:
        print(f"[{'OK' if (DATA / name).is_file() else 'EXTERNAL INPUT NEEDED'}] {name}")
    print("Download and verify the declared data before rerunning the Python fit.")
    print("No fit result is counted among the 96 R comparisons.")

if __name__ == "__main__":
    main()
