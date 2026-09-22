#!/usr/bin/env python3
"""Verify optional, separately acquired Pantheon+ inputs against frozen hashes."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'calc'))
from pantheon_data import pantheon_directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=None)
    args = parser.parse_args()
    directory = args.directory.expanduser().resolve() if args.directory else pantheon_directory()
    expected = json.loads(Path(__file__).with_name('pantheon_inputs.json').read_text())
    passed = True
    for item in expected['files']:
        path = directory / item['filename']
        if not path.is_file():
            print('MISSING', item['filename'])
            passed = False
            continue
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        ok = digest.hexdigest() == item['sha256'] and path.stat().st_size == item['bytes']
        print('PASS' if ok else 'MISMATCH', item['filename'], digest.hexdigest())
        passed &= ok
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
