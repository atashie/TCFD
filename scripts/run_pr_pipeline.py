"""Unattended driver for the precipitation ingest: stage 0 -> stage 1, with auto-retry.

WHY THIS EXISTS. The temperature ingest finished its main pass with FOUR of 36 members
failed on transient network errors (`[Errno 60] Operation timed out`, `[Errno 8] nodename
nor servname provided`), and sat done-but-incomplete until a human noticed and re-ran it.
Precipitation is 42 members x ~9 chunks = ~378 downloads over a run measured in days at the
observed 4-6 MB/s, so transient failures are near-certain and a stall that waits for
attention could cost a day of wall clock for no reason.

Both stages are idempotent and checkpoint per chunk, so a re-run costs only the chunks that
actually failed -- which is what makes an automatic retry safe rather than wasteful. Nothing
here re-downloads anything a previous attempt completed.

WHAT IT WILL NOT DO. It never reduces scope. If a member keeps failing after MAX_ATTEMPTS it
is reported and left out, loudly, rather than quietly shipping a thinner ensemble -- a
partial ensemble that looks complete is exactly the failure this repository has been bitten
by before. It also never skips stage 0: stage 1 refuses to run without the per-cell
thresholds, and that refusal is a feature.

Run:
    .venv/bin/python3 scripts/run_pr_pipeline.py            # resume/complete both stages
    .venv/bin/python3 scripts/run_pr_pipeline.py --attach   # stage 0 already running
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PY = ".venv/bin/python3"
STAGE0 = "scripts/pr_baseline_percentiles.py"
STAGE1 = "scripts/download_reduce_prthresh_isimip3b.py"
BASELINE_FILE = Path("data/interim/prthresh/baseline_percentiles.nc")
INTERIM = Path("data/interim/prthresh")
LOG = Path("data/interim/pr_pipeline.log")
MAX_ATTEMPTS = 6
#: Back off between attempts. A transient DNS or timeout failure clears in minutes; hammering
#: the file server immediately is how a rate limit gets earned (see the search skill).
BACKOFF_S = [60, 300, 900, 1800, 3600, 3600]


def log(msg):
    line = f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def running(script: str) -> bool:
    out = subprocess.run(["pgrep", "-f", script], capture_output=True, text=True)
    return bool(out.stdout.strip())


def run(script: str, extra: list[str]) -> int:
    cmd = [PY, "-u", script, "--run", *extra]
    log(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd).returncode


def wait_for(script: str, poll=120):
    """Attach to an already-running stage rather than starting a second copy of it."""
    if not running(script):
        return
    log(f"  {script} is already running -- attaching, will not start a second copy")
    while running(script):
        time.sleep(poll)
    log(f"  {script} exited")


def members_done() -> int:
    return len(list(INTERIM.glob("*_pr.nc")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--attach", action="store_true",
                    help="a stage may already be running; wait for it instead of starting it")
    a = ap.parse_args()

    log("=" * 72)
    log("precipitation pipeline driver starting")
    log("=" * 72)

    # ---- Stage 0 ------------------------------------------------------------
    if a.attach:
        wait_for(STAGE0)
    for attempt in range(MAX_ATTEMPTS):
        if BASELINE_FILE.exists():
            log(f"stage 0 complete: {BASELINE_FILE}")
            break
        if attempt:
            wait = BACKOFF_S[min(attempt - 1, len(BACKOFF_S) - 1)]
            log(f"stage 0 attempt {attempt + 1} after {wait}s backoff")
            time.sleep(wait)
        rc = run(STAGE0, ["--workers", str(a.workers)])
        log(f"stage 0 attempt {attempt + 1} exited rc={rc}")
    else:
        if not BASELINE_FILE.exists():
            log("STAGE 0 FAILED after every attempt -- stopping. Stage 1 cannot run without "
                "per-cell thresholds, and running it anyway would silently produce empty "
                "relative rungs.")
            return 1

    # ---- Stage 1 ------------------------------------------------------------
    if a.attach:
        wait_for(STAGE1)
    target = 42
    for attempt in range(MAX_ATTEMPTS):
        done = members_done()
        if done >= target:
            log(f"stage 1 complete: {done}/{target} members")
            break
        if attempt:
            wait = BACKOFF_S[min(attempt - 1, len(BACKOFF_S) - 1)]
            log(f"stage 1 retry {attempt}: {done}/{target} members present, "
                f"backing off {wait}s (re-runs cost only the chunks that failed)")
            time.sleep(wait)
        rc = run(STAGE1, ["--workers", str(a.workers)])
        log(f"stage 1 attempt {attempt + 1} exited rc={rc}, "
            f"{members_done()}/{target} members present")

    done = members_done()
    if done < target:
        missing = target - done
        log("=" * 72)
        log(f"INCOMPLETE: {done}/{target} members after {MAX_ATTEMPTS} attempts -- "
            f"{missing} still missing.")
        log("NOT proceeding to processing. A partial ensemble that looks complete is worse "
            "than an obvious gap: fix the cause, then re-run this driver (it resumes).")
        log("=" * 72)
        return 1

    log("=" * 72)
    log(f"BOTH STAGES COMPLETE -- {done}/{target} members in {INTERIM}")
    log("Next: scripts/process_prthresh.py (stage 2). NOT run automatically -- the decadal "
        "statistic branch and the slope choice are per-metric MEASUREMENTS that need looking "
        "at, not defaults to be applied while nobody is watching.")
    log("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
