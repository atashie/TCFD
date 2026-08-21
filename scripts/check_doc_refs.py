#!/usr/bin/env python3
"""Fail if a LIVE doc or config cites a repo path that does not exist.

Guards against the drift this repo has actually had: a skill instructing a retired
formatter, a README indexing deleted processors, a catalog pointing at a renamed script.
Historical records are exempt on purpose — WORKFLOW-ISSUES.md and the dated decision
journals in docs/ describe the past truthfully and are never edited retroactively.

Checked path literals: `scripts/...`, `config/...`, `docs/...` with a code-ish extension.
Skipped: anything containing a glob/brace/placeholder ({, *, <, >), and `data/`,
`reports/`, `deliveries/` (gitignored runtime — their contents are not repo facts).

Usage:  python scripts/check_doc_refs.py        # exit 0 clean, 1 with findings
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

LIVE_DOCS = [
    "CLAUDE.md",
    "GUARDRAILS.md",
    "OUTPUT-SPEC.md",
    "ASSET-CATALOG.md",
    "DATASET-ATTRIBUTES.md",
    "README.md",
    "scripts/README.md",
    "docs/README.md",
]
LIVE_GLOBS = [
    "docs/reporting/**/*.md",
    ".claude/skills/*/SKILL.md",
    "config/*.yaml",
]
# Append-only history — a path that no longer exists is not an error there.
EXEMPT_HINTS = ("WORKFLOW-ISSUES", "data-options-", "water-stress-status-", "water-stress-formulations")

PATH_RE = re.compile(
    r"(?:scripts|config|docs)/[A-Za-z0-9_\-./]+?\.(?:py|sh|yaml|yml|json|md|csv|txt)\b"
)
PLACEHOLDER = re.compile(r"[{}*<>]")


def live_files() -> list[Path]:
    files = [ROOT / p for p in LIVE_DOCS]
    for pattern in LIVE_GLOBS:
        files.extend(sorted(ROOT.glob(pattern)))
    return [f for f in files if f.is_file() and not any(h in f.name for h in EXEMPT_HINTS)]


def main() -> int:
    failures: list[str] = []
    for doc in live_files():
        rel_doc = doc.relative_to(ROOT)
        for lineno, line in enumerate(doc.read_text(encoding="utf-8").splitlines(), 1):
            for match in PATH_RE.findall(line):
                if PLACEHOLDER.search(match):
                    continue
                if not (ROOT / match).exists():
                    failures.append(f"{rel_doc}:{lineno}: cites nonexistent {match}")
    if failures:
        print(f"FAIL — {len(failures)} dangling reference(s):")
        for f in failures:
            print(f"  {f}")
        return 1
    print(f"OK — every checked path literal in {len(live_files())} live docs/configs resolves.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
