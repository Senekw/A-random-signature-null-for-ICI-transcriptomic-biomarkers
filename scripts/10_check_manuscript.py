#!/usr/bin/env python3
"""Step 10 -- reconcile the regenerated numbers against the manuscript.

Reads ``results/headline_numbers.json`` (step 08) and compares each value
against the number quoted in the submitted manuscript, listed in
``config/manuscript_claims.yaml``. Prints a table of agreements and
disagreements and exits non-zero if any headline claim has moved.

Tolerances are per-claim, not global: a proportion quoted to three decimal
places is held to a tighter bound than a percentage quoted as a round
number, and Monte-Carlo quantities (the null's empirical rates) get a band
that reflects the draw count rather than a spurious exact match.

The point is that "does the repository reproduce the paper" is a question
with a mechanical answer, not a reading exercise.

Usage:  python scripts/10_check_manuscript.py
        python scripts/10_check_manuscript.py --json report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def dig(obj, dotted: str):
    """Fetch a dotted path out of nested dicts; None if any hop is absent."""
    cur = obj
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="", help="also write a JSON report here")
    args = ap.parse_args()

    hn_path = ROOT / "results" / "headline_numbers.json"
    if not hn_path.exists():
        raise SystemExit("results/headline_numbers.json missing -- run "
                         "scripts/08_headline_numbers.py first.")

    claims = yaml.safe_load(
        (ROOT / "config" / "manuscript_claims.yaml").read_text()
    )
    hn = json.loads(hn_path.read_text())

    rows = []
    for c in claims["claims"]:
        got = dig(hn, c["path"])
        pub = c["published"]
        tol = float(c["tolerance"])

        if got is None:
            status = "MISSING"
            delta = None
        else:
            got = float(got)
            delta = got - float(pub)
            status = "ok" if abs(delta) <= tol else "MOVED"

        rows.append({
            "claim": c["claim"], "path": c["path"], "published": pub,
            "regenerated": got, "delta": delta, "tolerance": tol,
            "status": status, "note": c.get("note", ""),
        })

    w = max(len(r["claim"]) for r in rows)
    print(f"{'claim':<{w}}  {'paper':>8}  {'rerun':>8}  {'delta':>8}  status")
    print("-" * (w + 40))
    for r in rows:
        got = "--" if r["regenerated"] is None else f"{r['regenerated']:.3f}"
        dl = "--" if r["delta"] is None else f"{r['delta']:+.3f}"
        print(f"{r['claim']:<{w}}  {float(r['published']):>8.3f}  {got:>8}  "
              f"{dl:>8}  {r['status']}")

    moved = [r for r in rows if r["status"] == "MOVED"]
    missing = [r for r in rows if r["status"] == "MISSING"]
    ok = [r for r in rows if r["status"] == "ok"]

    print(f"\n{len(ok)}/{len(rows)} headline claims reproduce within tolerance")

    if missing:
        print(f"\n{len(missing)} not present in headline_numbers.json "
              f"(upstream step not run?):")
        for r in missing:
            print(f"  {r['claim']}  ({r['path']})")

    if moved:
        print(f"\n{len(moved)} MOVED beyond tolerance -- these change a "
              f"manuscript claim and must be reconciled before submission:")
        for r in moved:
            print(f"  {r['claim']}: paper {float(r['published']):.3f} -> "
                  f"rerun {r['regenerated']:.3f} "
                  f"(delta {r['delta']:+.3f}, tol {r['tolerance']})")
            if r["note"]:
                print(f"    {r['note']}")

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=2))
        print(f"\n-> {args.json}")

    sys.exit(1 if (moved or missing) else 0)


if __name__ == "__main__":
    main()
