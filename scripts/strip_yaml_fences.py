#!/usr/bin/env python3
"""
Strip stray markdown code fences (```yaml / ```) from spec YAML files.

Many specs in evals/generated_specs/ were captured from LLM output that
included surrounding markdown code fences. The fences make the files fail
yaml.safe_load. This script:

  1. Reads each *.yaml file under the target directory.
  2. If the file starts with a ```yaml/```yml fence and/or ends with a ```
     fence, strips them (idempotent: files without fences are untouched).
  3. Verifies the result parses with yaml.safe_load before writing.
  4. Reports counts of fixed / clean / still-broken files.

Run with --dry-run to preview without writing.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent

FENCE_OPEN_PREFIXES = ("```yaml", "```yml", "```")


def strip_fences(text: str) -> tuple[str, bool]:
    """Return (cleaned_text, changed).

    If the file contains a ```yaml / ``` fenced block, keep only the lines
    between the opener and closer (drop any prose before the opener and
    after the closer). Otherwise fall back to a no-op."""
    lines = text.splitlines(keepends=True)

    open_idx = None
    for i, line in enumerate(lines):
        if line.lstrip().startswith("```"):
            open_idx = i
            break
    if open_idx is None:
        return text, False

    close_idx = None
    for j in range(open_idx + 1, len(lines)):
        if lines[j].strip() == "```":
            close_idx = j
            break

    if close_idx is None:
        # Only opener, no closer — drop the opener and any prose before it
        kept = lines[open_idx + 1:]
    else:
        kept = lines[open_idx + 1:close_idx]

    return "".join(kept), True


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir", type=Path,
                   default=ROOT / "evals" / "generated_specs")
    p.add_argument("--dry-run", action="store_true",
                   help="Don't write; only report what would change.")
    args = p.parse_args()

    if not args.dir.exists():
        print(f"directory not found: {args.dir}", file=sys.stderr)
        return 2

    files = sorted(args.dir.glob("*.yaml"))
    n_total = len(files)
    n_clean = n_fixed = n_multi_doc = n_still_broken = 0
    still_broken: list[tuple[Path, str]] = []

    for path in files:
        text = path.read_text()
        cleaned, changed = strip_fences(text)
        candidate = cleaned if changed else text

        # Try as single doc first
        parsed_ok = False
        try:
            yaml.safe_load(candidate)
            parsed_ok = True
        except yaml.YAMLError as exc:
            err = str(exc).splitlines()[0]
            # Multi-document fallback: keep the first non-null document
            if "expected a single document" in err:
                try:
                    docs = list(yaml.safe_load_all(candidate))
                    first = next((d for d in docs if d is not None), None)
                    if first is not None:
                        candidate = yaml.safe_dump(first, sort_keys=False)
                        parsed_ok = True
                        n_multi_doc += 1
                        changed = True
                except yaml.YAMLError as exc2:
                    err = f"multi-doc retry: {str(exc2).splitlines()[0]}"

            if not parsed_ok:
                n_still_broken += 1
                prefix = "after strip: " if changed else ""
                still_broken.append((path, f"{prefix}{err}"))
                continue

        if changed and not args.dry_run:
            path.write_text(candidate)
        if changed:
            n_fixed += 1
        else:
            n_clean += 1

    print(f"Total spec files:     {n_total}")
    print(f"Already clean:        {n_clean}")
    print(f"Fenced & fixed:       {n_fixed}{' (dry-run)' if args.dry_run else ''}")
    print(f"  (incl. multi-doc collapsed to first doc: {n_multi_doc})")
    print(f"Still failing parse:  {n_still_broken}")
    if still_broken:
        print("\nFirst 10 still-broken:")
        for path, err in still_broken[:10]:
            print(f"  {path.name}: {err}")
    return 0 if n_still_broken == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
