"""Post-publish step: emit the QA/QC report and the map collection for a layer.

Every ingest-and-process run must leave behind reviewable HTML evidence, so
processors call :func:`finalize_layer` immediately after ``publish_processed_layer``
rather than relying on someone remembering to run two more scripts.

Both artifacts land in the version prefix but OUTSIDE ``_COMPLETE.json`` (STORAGE.md):
they are regenerable, so producing or re-producing them never invalidates published
data. That is also why a failure here is reported but NOT raised — the data is
already published and gated at this point, and losing a map collection must not
masquerade as a failed publish. Re-run the two CLIs to recover:

    python scripts/generate_qa_report.py {layer_id} --version {version}
    python scripts/generate_maps.py     {layer_id} --version {version}
"""

import sys
import traceback
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1]
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


def finalize_layer(layer_id, version=None, variable=None, local_only=False,
                   qa=True, maps=True):
    """Generate the QA report and maps for a published layer version.

    Args:
        layer_id: Canonical layer id.
        version: Version to document. Defaults to the layer's current version.
        variable: File stem override; inferred from ``data/`` when omitted.
        local_only: Write to the local cache instead of uploading.
        qa: Generate the QA/QC report.
        maps: Generate the interactive map collection.

    Returns:
        ``{"qa": <report dict or None>, "maps": <location or None>,
        "errors": [(stage, message), ...]}``. Never raises: see module docstring.
    """
    out = {"qa": None, "maps": None, "errors": []}

    if qa:
        try:
            import generate_qa_report
            out["qa"] = generate_qa_report.run(
                layer_id, variable=variable, version=version, local_only=local_only)
        except Exception as e:                                  # noqa: BLE001
            out["errors"].append(("qa", f"{type(e).__name__}: {e}"))
            print(f"WARNING: QA report failed -- {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()

    if maps:
        try:
            import generate_maps
            out["maps"] = generate_maps.run(
                layer_id, variable=variable, version=version, local_only=local_only)
        except Exception as e:                                  # noqa: BLE001
            out["errors"].append(("maps", f"{type(e).__name__}: {e}"))
            print(f"WARNING: map generation failed -- {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()

    print("\n" + "=" * 66, flush=True)
    if out["qa"]:
        s = out["qa"]["summary"]
        print(f"QA/QC   {s['verdict']}  ({s['passed']}/{s['checks_run']} passed, "
              f"{s['failed']} failed, {s['warnings']} warnings)", flush=True)
        print(f"        {out['qa'].get('published_to','')}", flush=True)
    if out["maps"]:
        print(f"MAPS    {out['maps']}", flush=True)
    if out["errors"]:
        print(f"NOTE    data is published and gated; {len(out['errors'])} derived "
              f"artifact(s) failed and can be regenerated:", flush=True)
        for stage, msg in out["errors"]:
            print(f"          {stage}: {msg}", flush=True)
    print("=" * 66, flush=True)

    # A FAILED (or absent) QA verdict must not be a line of log a reader can skim
    # past. This function still does not raise -- the data is already published and
    # gated, so aborting here would misreport a successful publish as a failure --
    # but the layer is NOT fit to use, and saying so quietly is how a bad layer gets
    # treated as a good one.
    verdict = (out["qa"] or {}).get("summary", {}).get("verdict")
    if verdict == "FAIL" or out["qa"] is None:
        why = ("QA reported FAILED invariant checks" if verdict == "FAIL"
               else "the QA report could not be produced, so nothing was verified")
        banner = "!" * 66
        print(f"\n{banner}\n"
              f"!! DO NOT USE THIS LAYER VERSION: {why}.\n"
              f"!! It is published and byte-complete, but NOT verified. Investigate\n"
              f"!! before consuming it, pointing any consumer at it, reporting it as\n"
              f"!! good, or running storage.cleanup_raw (GUARDRAILS S11).\n"
              f"{banner}", flush=True)
        out["fit_for_use"] = False
    else:
        out["fit_for_use"] = True
    return out
