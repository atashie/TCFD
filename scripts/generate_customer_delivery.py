#!/usr/bin/env python
"""Generate a customer-specific climate hazard extract (deterministic CSV star schema).

    # 1. See what WOULD be extracted -- this is the default, and it touches no data.
    python scripts/generate_customer_delivery.py --customer "Acme Timber" \
        --input location-analyses/acme-sites.csv

    # 2. Confirm the resolved asset -> layer mapping with the user, THEN extract.
    python scripts/generate_customer_delivery.py --customer "Acme Timber" \
        --input location-analyses/acme-sites.csv --run

Planning is the default and extraction requires --run ON PURPOSE. The workflow requires
that the specific variables resolved for each asset be shown to the user before a run, and
a flag is the only way to make that a property of the tool rather than of whoever is
driving it. See ASSET-CATALOG.md.

INPUT: one row per location-asset combination.
    required   Location, Lat, Lon, Asset_Type
    optional   Sub_Asset_Unit, Country, State, City, Region, Subregion, Layers

`Layers` (semicolon-separated layer_ids) overrides the asset catalog for that row.
`--layers` overrides it for every row.

OUTPUT: deliveries/{customer-slug}/{YYYYMMDD}/
    locations.csv  assets.csv  layers.csv  values.csv  manifest.json  README.md
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.generate_delivery_dashboard import build_dashboard  # noqa: E402
from scripts.utils.delivery import (  # noqa: E402
    DeliveryError,
    build_plan,
    delivery_dir,
    discover_scenarios,
    load_asset_catalog,
    load_input,
    load_registry,
    record_stage,
    run_delivery,
    scenario_path,
)
from scripts.utils.spatial_extract import as_period_dataset  # noqa: E402


def cmd_list_layers(registry, catalog) -> int:
    """Show what is available to deliver, and what each asset type resolves to."""
    print("\nREGISTERED LAYERS")
    print("-" * 78)
    for layer_id, spec in sorted(registry.layers.items()):
        try:
            scenarios = ", ".join(discover_scenarios(registry, spec))
        except DeliveryError as exc:
            scenarios = f"UNAVAILABLE -- {str(exc).splitlines()[0]}"
        print(f"  {layer_id:<14} {spec.hazard:<22} [{spec.status}]")
        print(f"  {'':<14} {spec.hazard_measure}")
        print(f"  {'':<14} scenarios: {scenarios}")
        print(f"  {'':<14} read slope: {spec.recommended_slope}")
        print()

    print("ASSET CATALOG")
    print("-" * 78)
    for name, entry in sorted(catalog.entries.items()):
        confirmed = entry.get("confirmed_on") or "UNCONFIRMED -- needs user sign-off"
        print(f"  {name:<20} -> {', '.join(entry.get('layers') or [])}")
        print(f"  {'':<20}    aliases: {', '.join(entry.get('aliases') or []) or '-'}")
        print(f"  {'':<20}    {confirmed}")
        print()
    return 0


def cmd_measure_slopes(registry) -> int:
    """Re-measure which slope each layer should be read on.

    Judged on ACTIVE cells only (either slope non-zero) at the final decade, per
    OUTPUT-SPEC.md -- the all-cell view is diluted by ocean and permanently-zero land and
    yields the opposite conclusion.
    """
    print(f"\n{'layer':<16}{'active':>9}{'sen==0':>9}{'agree':>8}  recommend")
    print("-" * 62)
    for layer_id, spec in sorted(registry.layers.items()):
        try:
            scenario = discover_scenarios(registry, spec)[0]
            with xr.open_dataset(scenario_path(registry, spec, scenario)) as ds:
                # An observational layer publishes no slopes and has no decade axis. It is
                # not a layer whose slope choice is "unmeasured" -- it is one where the
                # question does not arise, so say that rather than crashing on a KeyError.
                if "ols_slope" not in ds or "sen_slope" not in ds:
                    print(f"{layer_id:<16}  no slopes -- observational layer, not applicable")
                    continue
                decade = int(ds.decade.values[-1])
                ols = ds["ols_slope"].sel(decade=decade).values
                sen = ds["sen_slope"].sel(decade=decade).values
        except DeliveryError as exc:
            print(f"{layer_id:<16}  {str(exc).splitlines()[0]}")
            continue

        finite = np.isfinite(ols) & np.isfinite(sen)
        active = finite & ((ols != 0) | (sen != 0))
        n = int(active.sum())
        if not n:
            print(f"{layer_id:<16}{0:>9}  no active cells")
            continue
        sen_zero = float((sen[active] == 0).sum()) / n
        agree = float((np.sign(ols[active]) == np.sign(sen[active])).mean())
        recommend = "ols_slope" if sen_zero > 0.10 else "sen_slope"
        flag = "" if recommend == spec.recommended_slope else "  <-- REGISTRY DISAGREES"
        print(f"{layer_id:<16}{n:>9}{sen_zero:>9.3f}{agree:>8.3f}  {recommend}{flag}")
    print("\nMeasured on the final decade of each layer's first scenario.")
    return 0


def print_plan(customer, input_path, locations_df, assets_df, work, registry) -> None:
    print(f"\nDELIVERY PLAN -- {customer}")
    print("=" * 78)
    print(f"input     {input_path}")
    print(f"output    {delivery_dir(customer)}")
    print(f"locations {len(locations_df)}    assets {len(assets_df)}")

    print("\nRESOLVED ASSET -> LAYERS")
    print("-" * 78)
    loc_names = locations_df.set_index("location_id")["name"].to_dict()
    for _, a in assets_df.iterrows():
        sub = f" / {a['sub_asset_unit']}" if a["sub_asset_unit"] else ""
        print(f"  {a['asset_id']}  {loc_names[a['location_id']]}")
        print(f"           {a['asset_type']}{sub}  [{a['layer_source']}]")
        for layer_id in a["layer_ids"]:
            spec = registry.get(layer_id)
            print(f"             - {layer_id:<14} {spec.hazard_measure}")

    print("\nLAYERS TO READ")
    print("-" * 78)
    total_rows = 0
    for item in work:
        spec = registry.get(item["layer_id"])
        scenarios = discover_scenarios(registry, spec)
        with xr.open_dataset(scenario_path(registry, spec, scenarios[0])) as raw:
            # Observational layers carry one observed period rather than a decade axis.
            ds = as_period_dataset(raw)
            n_decades = len(ds.decade)
            units = ds.attrs.get("units", "?")
        rows = item["n_assets"] * len(scenarios) * n_decades
        total_rows += rows
        print(f"  {item['layer_id']:<14} {item['n_assets']:>3} asset(s) x "
              f"{len(scenarios)} scenario(s) x {n_decades} decades = {rows:>5} rows")
        print(f"  {'':<14} units {units}   scenarios {', '.join(scenarios)}")
        print(f"  {'':<14} read {spec.recommended_slope}")
        # The relative-baseline note comes FIRST and is labelled loudest. It is the one
        # property of a layer a reader can get wrong while every number in front of them is
        # correct, and this plan is where the mapping is agreed with the customer -- the
        # skill's Stage 1 explicitly relies on these NOTE lines carrying it. Printing only
        # `delivery_note` would hide it, because the substance lives in its own field so
        # that `generate_delivery_caveats.py` can promote it to must-disclose.
        if spec.relative_baseline:
            note = " ".join(spec.relative_baseline_note.split())
            print(f"  {'':<14} READ AS RELATIVE: {note}")
        if spec.delivery_note:
            note = " ".join(spec.delivery_note.split())
            print(f"  {'':<14} NOTE: {note}")

    print("-" * 78)
    print(f"  values.csv will hold {total_rows} rows")
    print("\nNothing has been extracted. Confirm the mapping above with the customer,")
    print("then re-run with --run.")


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--customer", help="Customer name; becomes the delivery folder slug")
    p.add_argument("--input", type=Path, help="Location-asset CSV/TSV")
    p.add_argument("--run", action="store_true",
                   help="Actually extract. Without this, only the plan is printed.")
    p.add_argument("--no-dashboard", action="store_true",
                   help="Skip the dashboard. Leaves the delivery INCOMPLETE and is "
                        "recorded as such in the manifest; for debugging only.")
    p.add_argument("--reports", action="store_true",
                   help="Also build Stage 4 (caveats) and Stage 3a (compliance report). "
                        "The bespoke report is NOT chained -- it needs facet profiles "
                        "chosen and a narrative written.")
    p.add_argument("--layers", help="Semicolon-separated layer_ids overriding the catalog "
                                    "for EVERY row")
    p.add_argument("--out", type=Path, help="Override the delivery directory")
    p.add_argument("--list-layers", action="store_true",
                   help="Show the layer registry and asset catalog, then exit")
    p.add_argument("--measure-slopes", action="store_true",
                   help="Re-measure which slope each layer should be read on, then exit")
    args = p.parse_args()

    try:
        registry = load_registry()
        catalog = load_asset_catalog()

        if args.list_layers:
            return cmd_list_layers(registry, catalog)
        if args.measure_slopes:
            return cmd_measure_slopes(registry)

        if not args.customer or not args.input:
            p.error("--customer and --input are required (or use --list-layers)")

        override = (
            [s.strip() for s in args.layers.split(";") if s.strip()]
            if args.layers else None
        )
        df = load_input(args.input)
        locations_df, assets_df, work = build_plan(df, catalog, registry, override)

        if not args.run:
            print_plan(args.customer, args.input, locations_df, assets_df, work, registry)
            return 0

        out_dir = args.out or delivery_dir(args.customer)
        print(f"\nExtracting -> {out_dir}")
        manifest = run_delivery(
            customer=args.customer,
            input_path=args.input,
            locations_df=locations_df,
            assets_df=assets_df,
            work=work,
            registry=registry,
            out_dir=out_dir,
        )
        counts = manifest["counts"]
        print(f"\nWrote {counts['value_rows']} value rows "
              f"({counts['locations']} locations, {counts['assets']} assets, "
              f"{counts['layers']} layers)")
        for status, n in sorted(counts["data_status"].items()):
            print(f"  {status}: {n}")

        # STAGE 2 IS CSVs *AND* DASHBOARD, ALWAYS. A delivery shipped as bare CSVs has no
        # way for anyone to look at what was extracted, and "check the reference sites"
        # then depends on somebody remembering a second command. Build it here so the two
        # cannot drift apart, and record the stage in the manifest so later stages
        # (compliance report, bespoke report, caveats) can see what has run.
        # A partially-built delivery directory looks exactly like a finished one from the
        # outside, and nothing stops it being copied to a customer. Mark it on disk -- the
        # same INVALID-DO-NOT-USE.md convention the layer loader already refuses on -- and
        # exit non-zero so a caller in a script does not treat it as success.
        incomplete_marker = out_dir / "DELIVERY-INCOMPLETE.md"
        incomplete_marker.unlink(missing_ok=True)
        failed = None
        if args.no_dashboard:
            print("\n  dashboard SKIPPED (--no-dashboard) -- this delivery is incomplete")
            record_stage(out_dir, "dashboard", "skipped")
            failed = "dashboard skipped via --no-dashboard"
        else:
            try:
                build_dashboard(out_dir)
                record_stage(out_dir, "dashboard", "built")
            except Exception as exc:  # noqa: BLE001 -- must not leave a silent half-delivery
                record_stage(out_dir, "dashboard", "failed", detail=str(exc))
                failed = f"dashboard build failed: {exc}"
                print(f"\n  ERROR: {failed}", file=sys.stderr)

        if failed:
            incomplete_marker.write_text(
                f"# DELIVERY INCOMPLETE -- do not ship\n\n{failed}\n\n"
                f"The CSVs in this folder may be correct, but the delivery is not finished. "
                f"Rebuild the dashboard, then re-run "
                f"`python scripts/test_customer_delivery.py {out_dir}`.\n")
            print(f"\n{out_dir}")
            print(f"  wrote {incomplete_marker.name} -- this delivery must not be shipped")
            return 1

        # STAGE 4 THEN STAGE 3a. Caveats first, always -- the caveat set is an INPUT to the
        # reports (each is required to carry every must-disclose entry) rather than a
        # summary of them. The bespoke report is NOT chained: it needs facet profiles chosen
        # and a narrative written, neither of which a batch run can supply.
        if args.reports:
            print()
            # Flush before handing the terminal to a child: this process's stdout is block
            # buffered when piped, the child's is not, so without this the child's output
            # appears BEFORE the extraction log it actually followed.
            sys.stdout.flush()
            rc = 0
            for script in ("generate_delivery_caveats.py", "generate_compliance_report.py"):
                # Subprocess rather than an import: each stage owns its own argparse and
                # its own exit code, and a failure in one must not leave this process
                # holding half-initialised module state.
                rc = subprocess.call(
                    [sys.executable, str(Path(__file__).parent / script), str(out_dir)]
                )
                if rc != 0:
                    break
            if rc != 0:
                incomplete_marker.write_text(
                    "# DELIVERY INCOMPLETE -- do not ship\n\n"
                    "Caveats or the compliance report failed to build.\n")
                print(f"\n  wrote {incomplete_marker.name} -- this delivery must not be shipped")
                return 1

        print(f"\n{out_dir}")
        print("Next: verify, then OPEN the dashboard --")
        print(f"  python scripts/test_customer_delivery.py {out_dir}")
        print(f"  open {out_dir}/dashboard.html")
        if not args.reports:
            print("\nStage 4 + 3a (caveats, compliance report) -- add --reports, or run:")
            print(f"  python scripts/generate_delivery_caveats.py {out_dir}")
            print(f"  python scripts/generate_compliance_report.py {out_dir}")
        print("\nStage 3b (bespoke report) needs facets chosen in report_config.yaml first:")
        print(f"  python scripts/generate_bespoke_report.py {out_dir} --scaffold")
        return 0

    except DeliveryError as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
