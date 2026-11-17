"""Check all report files against the saved reference run."""
import argparse
import hashlib
from pathlib import Path


def hashes(folder):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("actual", type=Path)
    parser.add_argument("expected", type=Path, nargs="?", default=Path("expected/report"))
    args = parser.parse_args()
    actual, expected = hashes(args.actual), hashes(args.expected)
    if actual != expected:
        changed = sorted(k for k in actual.keys() | expected.keys() if actual.get(k) != expected.get(k))
        raise SystemExit("Mismatch: " + ", ".join(changed))
    print(f"PASS: all {len(expected)} report files match (SHA-256)")
