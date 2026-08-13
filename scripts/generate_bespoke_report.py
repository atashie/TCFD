#!/usr/bin/env python3
"""Stage 3b -- the customer-specific report, composed from facet profiles and a written
narrative.

    python scripts/generate_bespoke_report.py <delivery> --scaffold   # first: write the shell
    python scripts/generate_bespoke_report.py <delivery>              # then: render it

TWO AUTHORS, AND THE BOUNDARY BETWEEN THEM IS ENFORCED
-------------------------------------------------------
The generator owns everything derivable from the delivery: figures, tables, counts, the
caveat set. It cannot get those wrong twice in two documents because the compliance report
computes them through the same functions.

A person owns the narrative -- what these numbers mean for this business, this asset class,
this decision. That is judgement and no template produces it. What the generator does is
make it impossible to publish judgement that is not sourced:

  * `narrative.md` is a set of SLOTS. An unfilled slot fails the build.
  * A paragraph in a slot marked `requires_citation: yes` must carry a citation.
  * A citation must RESOLVE -- `[data:...]` to a row that exists in the delivered CSVs,
    `[dossier:...]` to a source recorded in `dossier.yaml`.

That last check is the point of the whole design. The failure mode here is not laziness, it
is fluency: a confident paragraph about a customer's operations reads identically whether it
was researched or invented, and by the time it reaches a disclosure nobody can tell which.
An unresolvable citation is worse than no citation, because it looks like evidence.

PROFILES DO NOT APPEAR IN THE OUTPUT
------------------------------------
Facet profiles are surfaced into the scaffold as HTML comments -- transmission channels, what
this reader is deciding, what not to claim -- and the writer works from them. Rendering them
directly would make every loblolly report identical while looking bespoke, which is the exact
failure the library exists to prevent. See scripts/utils/report_profiles.py.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.delivery import DeliveryError, record_stage  # noqa: E402
from utils import report_figures as F  # noqa: E402
from utils import report_profiles as P  # noqa: E402
from utils.report_common import (  # noqa: E402
    DOSSIER_FILENAME,
    NOT_ASSESSED,
    SEVERITY_SHOULD,
    VULNERABLE,
    Delivery,
    assert_baseline_tier_equality,
    assert_caveats_present,
    band_cell,
    build_stamp,
    check_narrative,
    coverage_summary,
    customer_evidence,
    drop_leading_heading,
    number_citations,
    render_references,
    document,
    esc,
    figure_block,
    hazard_layer_ids,
    raw,
    load_delivery,
    markdown,
    must_disclose,
    parse_narrative,
    portfolio_score_series,
    read_stamp,
    table,
    vulnerability_frame,
)
from utils.viz_common import SCENARIO_TIER, TIER_LABELS, TIER_ORDER, risk_band  # noqa: E402

REPORT_FILENAME = "report_bespoke.html"
NARRATIVE_FILENAME = "narrative.md"

#: (slot, heading, requires_citation, profile section that guides it).
#: Order is the document's order. The spine is the reader's decision, not the standard's
#: structure -- that is the whole difference from the compliance report.
SLOTS = [
    ("executive_summary", "What we found", True, "What this reader needs decided"),
    ("portfolio_context", "Your portfolio in context", True, "Framing and vocabulary"),
    ("hazard_reading", "What these hazards mean for this asset class", True,
     "Transmission channels"),
    ("site_findings", "Site-by-site", True, "Metrics that lead"),
    ("decision", "What this changes for the decision in front of you", True,
     "What this reader needs decided"),
    ("confounders", "What would change this conclusion", True, "Known confounders"),
    ("next_steps", "What we recommend doing next", False, "Questions for the customer"),
]

DATA_SECTIONS = [
    ("figures", "The numbers behind this"),
    ("gaps", "What we could not assess"),
    ("sources", "Sources and provenance"),
]


def callout(title: str, body_html: str, kind: str = "") -> str:
    return f'<div class="callout {kind}"><h4>{esc(title)}</h4>{body_html}</div>'


# ---------------------------------------------------------------------------------------
# Scaffold
# ---------------------------------------------------------------------------------------


def scaffold_narrative(delivery: Delivery, profiles: Dict[str, Optional[P.Profile]]) -> str:
    chosen = ", ".join(f"{f}={p.id}" for f, p in profiles.items() if p)
    out = [
        f"<!-- narrative for {delivery.customer} — {delivery.path.name} -->",
        f"<!-- profiles: {chosen} -->",
        "<!--",
        "  HOW THIS FILE WORKS",
        "    * Fill every slot. An empty slot, or one still containing TODO, fails the build.",
        "    * In slots marked requires_citation: yes, EVERY paragraph needs a citation.",
        "    * Citations must resolve:",
        "        [data:values.csv#ASSET/LAYER/SCENARIO/DECADE]",
        "        [data:climate_score.csv#ASSET/TIER/DECADE]",
        "        [data:caveats.json#CAVEAT-ID]  [data:layers.csv#LAYER]",
        "        [dossier:SOURCE-ID]  -- must exist in dossier.yaml",
        "    * The guidance under each slot comes from the selected facet profiles. It is",
        "      for you, not for the customer: do not paste it. Write from it.",
        "-->",
        "",
    ]
    for slot, heading, cite, guide_section in SLOTS:
        out.append(f"<!-- slot: {slot} | requires_citation: {'yes' if cite else 'no'} -->")
        out.append(f"## {heading}")
        out.append("")
        guidance = P.guidance_block(profiles, guide_section)
        if guidance:
            out.append(f"<!-- GUIDANCE — {guide_section}")
            out.append(guidance)
            out.append("-->")
        never = P.guidance_block(profiles, "Do not claim")
        if never:
            out.append("<!-- DO NOT CLAIM")
            out.append(never)
            out.append("-->")
        out.append("TODO")
        out.append("")
    return "\n".join(out) + "\n"


def scaffold_dossier(delivery: Delivery) -> str:
    return f"""# Researched facts about {delivery.customer}, with provenance.
#
# Every non-data claim in narrative.md must cite a source id from this file, and the build
# refuses if a citation does not resolve. That is the mechanism that keeps a fluent,
# plausible, invented sentence out of a customer document.
#
# SOURCE PRECEDENCE, best first. Prefer the highest tier that answers the question, and say
# in the narrative when you are relying on a lower one:
#   customer-direct     something the customer told us, with the date and the person
#   customer-published  their own annual report, website, filing
#   regulator           a filing, register or official dataset
#   peer-reviewed       published research
#   trade-press         industry reporting -- useful for context, weak for facts
#
# TEMPORAL LAG IS A FIELD, NOT A FEELING. `as_of` is when the fact was true; `retrieved_on`
# is when we read it. A 2019 figure retrieved today is a 2019 figure, and a narrative that
# leans on it should say so.

customer: {delivery.customer!r}
compiled_on: "{datetime.now(timezone.utc).date().isoformat()}"

sources: []
  # - id: example-annual-report-2025
  #   title: "Annual Report 2025"
  #   publisher: ""
  #   url: "https://"
  #   source_type: customer-published
  #   retrieved_on: "{datetime.now(timezone.utc).date().isoformat()}"
  #   as_of: "2025-12-31"
  #   what_it_supports: >-
  #     The specific claim this source backs. One source, one claim -- a source listed
  #     without saying what it establishes cannot be checked by anyone else.
  #   verified_by: >-
  #     How the fact was double-checked, and against what. Leave blank ONLY if it was not.
"""


# ---------------------------------------------------------------------------------------
# Data sections (identical objects to the compliance report's, never recomputed)
# ---------------------------------------------------------------------------------------


def sec_figures(delivery: Delivery, vf: pd.DataFrame) -> str:
    hz = delivery.config["horizons"]
    long_dec = int(hz["long"]["decade"]) if "long" in hz else max(delivery.decades)
    tiers = [t for t in TIER_ORDER if t in set(vf["scenario_tier"])]
    tier = tiers[-1] if tiers else "medium"

    sel = vf[(vf["decade"] == long_dec) & (vf["scenario_tier"] == tier)].sort_values(
        "worst_percentile", ascending=False, na_position="last"
    )
    bars = [(delivery.asset_label(r["asset_id"]), r["worst_percentile"])
            for _, r in sel.iterrows()]

    series, panel, kept, dropped = portfolio_score_series(delivery)
    assert_baseline_tier_equality(series)
    score_fig = (
        F.decade_line(series, title="Portfolio Climate Score by forcing tier",
                      y_label="Climate Score — mean risk percentile")
        if series else
        F.empty_figure("Portfolio Climate Score by forcing tier",
                       "No asset is scored in every tier and decade.")
    )
    score_note = (
        f"Balanced panel of {len(panel)} of {len(delivery.assets)} assets."
        if series else
        "No asset carries a score in every tier and decade, so a cross-tier comparison "
        "would compare different sets of assets and is not drawn."
    )

    detail = table(
        ["Asset", "Worst hazard", "Percentile", "Determination", "Hazards assessed"],
        [
            [
                delivery.asset_label(r["asset_id"]),
                r["worst_hazard"] or "—",
                band_cell(r["worst_percentile"]),
                "Vulnerable" if r["status"] == VULNERABLE
                else ("Not assessed" if r["status"] == NOT_ASSESSED else "Not vulnerable"),
                f"{int(r['n_hazards_assessed'])} of {int(r['n_hazards_expected'])}",
            ]
            for _, r in sel.iterrows()
        ],
        align_right=(2,),
        caption=f"{TIER_LABELS[tier]}, {long_dec}s.",
    )
    absent = delivery.layers_absent_at_tier(tier)
    gap = (
        f"<p class='sub'>{', '.join(esc(l) for l in absent)} publishes no "
        f"{esc(TIER_LABELS[tier].lower())} scenario, so assets carrying it show fewer "
        f"hazards assessed at this tier. Those sites were assessed for that hazard at the "
        f"tiers where it exists.</p>"
        if absent else ""
    )
    return f"""
{figure_block(
    F.ranked_bar(bars, title=f"Sites by worst assessed hazard — {TIER_LABELS[tier].lower()}, {long_dec}s"),
    caption="Percentile ranks each site against global land in the 2020s baseline. It is a "
            "level, not a change.",
    table_html=detail,
)}
{gap}
{F.band_legend()}
{figure_block(score_fig, caption=score_note)}
"""


def sec_gaps(delivery: Delivery, cov: dict) -> str:
    out = [
        "<p>Everything below limits what this report can be used for. It is generated from "
        "the delivered data rather than written, so it cannot be softened by accident.</p>"
    ]
    for c in must_disclose(delivery):
        ev = customer_evidence(c["evidence"])
        out.append(
            callout(
                f"{c['id']} — {c['title']}",
                f"<p>{esc(c['text'])}</p>"
                + (f"<p class='sub'>Evidence in this delivery: {esc(ev)}</p>" if ev else ""),
                kind="must",
            )
        )
    out.append("<h3>Hazards outside this assessment</h3>")
    out.append(
        table(
            ["Hazard family", "Class", "Status and what it means for this portfolio"],
            [[u["name"], u["esrs_class"].title(), u["customer_note"]] for u in cov["uncovered"]],
            caption=f"{len(cov['uncovered'])} of {cov['n_families']} families not assessed.",
        )
    )
    should = [c for c in delivery.caveats if c["severity"] == SEVERITY_SHOULD]
    if should:
        out.append("<h3>Notes on reading individual figures</h3>")
        out.append(
            table(["Identifier", "Note"],
                  [[c["id"], f"{c['title']}. {c['text']}"] for c in should])
        )
    return "\n".join(out)


def sec_sources(delivery: Delivery, profiles: Dict[str, Optional[P.Profile]]) -> str:
    sources = list(delivery.dossier.get("sources") or [])
    for p in profiles.values():
        if p:
            sources.extend(p.sources)
    seen, rows = set(), []
    for s in sources:
        sid = str(s.get("id", ""))
        if not sid or sid in seen:
            continue
        seen.add(sid)
        lag = ""
        if s.get("as_of") and s.get("retrieved_on"):
            lag = f"as of {s['as_of']}, read {s['retrieved_on']}"
        elif s.get("as_of"):
            lag = f"as of {s['as_of']}"
        url = str(s.get("url") or "")
        title = str(s.get("title") or sid)
        rows.append(
            [
                sid,
                f'<a href="{esc(url)}" rel="noopener noreferrer">{esc(title)}</a>'
                if url.startswith("http") else esc(title),
                str(s.get("source_type") or ""),
                lag,
                str(s.get("what_it_supports") or ""),
            ]
        )
    rows = [[r[0], raw(r[1]), r[2], r[3], r[4]] for r in rows]

    body = (
        table(["Reference", "Source", "Type", "Currency", "What it supports"], rows,
              caption="Every non-data claim in this report cites one of these.")
        if rows else
        "<p class='muted'>No external sources were used: every statement in this report "
        "derives from the delivered climate data.</p>"
    )
    climate = table(
        ["Layer", "Measure", "Units", "Ensemble", "Scenarios"],
        [
            [
                r["layer_id"], r.get("long_name", ""), r.get("units", ""),
                f"{len([m for m in str(r.get('impact_models') or '').split(',') if m.strip()])} × "
                f"{len([g for g in str(r.get('gcms') or '').split(',') if g.strip()])}",
                str(r.get("scenarios", "")).replace(";", ", "),
            ]
            for _, r in delivery.layers.iterrows()
        ],
        caption="Climate data underlying every figure (ISIMIP impact-model output, 0.5° grid).",
    )
    return f"<h3>Researched sources</h3>{body}<h3>Climate data</h3>{climate}"


# ---------------------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------------------


def build_report(delivery: Delivery, profiles: Dict[str, Optional[P.Profile]]) -> str:
    narrative_path = delivery.path / NARRATIVE_FILENAME
    if not narrative_path.exists():
        raise DeliveryError(
            f"{narrative_path} does not exist. Write the scaffold first:\n"
            f"  python scripts/generate_bespoke_report.py {delivery.path} --scaffold"
        )
    text = narrative_path.read_text()
    problems = check_narrative(delivery, text)
    if problems:
        raise DeliveryError(
            f"{narrative_path} is not ready ({len(problems)} problem(s)):\n"
            + "\n".join(f"  - {p}" for p in problems)
            + "\n\nEvery slot must be filled, every claim in a cited slot must carry a "
              "citation, and every citation must resolve to real data or a recorded source."
        )

    slots = {s["name"]: s for s in parse_narrative(text)}
    cov = coverage_summary(delivery)
    vf = vulnerability_frame(delivery)

    missing_slots = [s for s, _h, _c, _g in SLOTS if s not in slots]
    if missing_slots:
        raise DeliveryError(
            f"{narrative_path} is missing slot(s): {', '.join(missing_slots)}. Regenerate "
            f"the scaffold, or add:\n"
            + "\n".join(f"  <!-- slot: {s} | requires_citation: yes -->" for s in missing_slots)
        )

    # Number citations across the whole narrative BEFORE rendering, so a repeated citation
    # keeps one number and the reference list runs in document order.
    bodies = {s: drop_leading_heading(slots[s]["body"]) for s, _h, _c, _g in SLOTS}
    bodies, refs = number_citations(bodies)

    body_parts: List[str] = []
    toc: List[tuple] = []
    for slot, heading, _cite, _guide in SLOTS:
        toc.append((slot, heading))
        body_parts.append(
            f'<section id="{slot}"><h2>{esc(heading)}</h2>{markdown(bodies[slot])}</section>'
        )

    data_html = {
        "figures": sec_figures(delivery, vf),
        "gaps": sec_gaps(delivery, cov),
        "sources": render_references(delivery, refs) + sec_sources(delivery, profiles),
    }
    for anchor, heading in DATA_SECTIONS:
        toc.append((anchor, heading))
        body_parts.append(
            f'<section id="{anchor}"><h2>{esc(heading)}</h2>{data_html[anchor]}</section>'
        )

    body = "\n".join(body_parts)
    stamp = build_stamp(body, Path(__file__).read_text())
    prof_line = " · ".join(f"{p.name}" for p in profiles.values() if p)
    html = document(
        title="Physical Climate Risk — Portfolio Review",
        subtitle=prof_line or "Prepared for this portfolio",
        delivery=delivery,
        body=body,
        stamp=stamp,
        toc=toc,
    )
    assert_caveats_present(delivery, html)
    return html


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("delivery", type=Path)
    ap.add_argument("--scaffold", action="store_true",
                    help="write narrative.md and dossier.yaml templates, then stop")
    ap.add_argument("--check-stamp", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    path = args.delivery / REPORT_FILENAME
    if args.check_stamp:
        if not path.exists():
            print(f"{path} does not exist")
            return 1
        print(f"on disk: {read_stamp(path.read_text())}")
        return 0

    try:
        delivery = load_delivery(args.delivery)
        profiles = P.load_selected(delivery.config)

        if args.scaffold:
            n = delivery.path / NARRATIVE_FILENAME
            d = delivery.path / DOSSIER_FILENAME
            wrote = []
            if not n.exists():
                n.write_text(scaffold_narrative(delivery, profiles))
                wrote.append(n.name)
            if not d.exists():
                d.write_text(scaffold_dossier(delivery))
                wrote.append(d.name)
            print(f"  scaffold: {', '.join(wrote) if wrote else 'already present, nothing overwritten'}")
            unconf = P.unconfirmed(profiles)
            if unconf:
                print(f"  NOTE: profiles not signed off: {', '.join(unconf)}")
            print(f"  Fill {n.name}, record sources in {d.name}, then rebuild.")
            return 0

        html = build_report(delivery, profiles)
    except DeliveryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    path.write_text(html)
    stamp = read_stamp(html)
    record_stage(args.delivery, "bespoke_report", "built",
                 f"{REPORT_FILENAME} ({len(html) // 1024} KB), build {stamp}")
    if not args.quiet:
        print(f"  bespoke report: {REPORT_FILENAME}  ({len(html) // 1024} KB)  build {stamp}")
        unconf = P.unconfirmed(profiles)
        if unconf:
            print(f"    NOTE: profiles not signed off: {', '.join(unconf)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
