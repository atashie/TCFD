#!/usr/bin/env python3
"""Stage 3a -- the physical climate risk compliance report, on the IFRS S2 spine.

    python scripts/generate_compliance_report.py deliveries/<customer>/<date>
    python scripts/generate_compliance_report.py <delivery> --check-stamp

WHY IFRS S2 AND NOT TCFD
------------------------
TCFD was disbanded in 2023 and the IFRS Foundation took over monitoring; IFRS S2 carries its
four pillars forward and is what the successor regimes point at. California SB 261 names the
IFRS Sustainability Disclosure Standards as an accepted framework, CDP has aligned its
corporate questionnaire with S2, and ESRS E1 asks many of the same questions in its own
vocabulary. So one document on the S2 spine, with a mapping appendix, serves several filings;
four documents on four spines would drift apart at the first data refresh.

THIS DOCUMENT IS FULLY DETERMINISTIC
------------------------------------
Every sentence here is either fixed text about method or a number computed from the delivery.
There are no narrative slots and no researched claims, which is deliberate: a compliance
artifact should be reproducible byte-for-byte from its inputs so that an auditor can rebuild
it and get the same document. Judgement, context and recommendation belong in the bespoke
report, which is where the citation machinery lives.

WHAT IT CANNOT DO, STATED UP FRONT RATHER THAN DISCOVERED LATER
---------------------------------------------------------------
IFRS S2 asks for four pillars. This report can populate Strategy (in part) and Metrics (in
part). It cannot populate Governance or Risk Management -- those describe how the ENTITY is
run and no climate model contains them -- and it cannot convert hazard exposure into
financial effect, which needs vulnerability functions, asset values and business assumptions
we do not have. Those gaps are printed in the document as explicit customer-owned sections
rather than quietly omitted, because a disclosure with a silently missing pillar looks
complete.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.delivery import DeliveryError, record_stage  # noqa: E402
from utils import report_figures as F  # noqa: E402
from utils.report_common import (  # noqa: E402
    NOT_ASSESSED,
    read_stamp,
    SEVERITY_MUST,
    SEVERITY_SHOULD,
    VULNERABLE,
    Delivery,
    assert_baseline_tier_equality,
    assert_caveats_present,
    band_cell,
    build_stamp,
    coverage_summary,
    customer_evidence,
    tbd_block,
    document,
    esc,
    figure_block,
    hazard_layer_ids,
    load_delivery,
    load_hazard_taxonomy,
    must_disclose,
    portfolio_score_series,
    raw,
    table,
    vulnerability_frame,
    vulnerability_rollup,
    vulnerability_sensitivity,
)
from utils.viz_common import SCENARIO_TIER, TIER_LABELS, TIER_ORDER, risk_band  # noqa: E402

REPORT_FILENAME = "report_compliance.html"

#: Which native scenario codes are treated as the Paris-aligned member for S2 22(b)(iv),
#: which asks whether a scenario "aligned with the latest international agreement on climate
#: change" was among those used.
PARIS_ALIGNED = {"rcp26", "ssp126"}

SECTIONS = [
    ("scope", "1. Scope and basis of this assessment"),
    ("risks", "2. Physical climate risks identified"),
    ("horizons", "3. Time horizons"),
    ("concentration", "4. Where physical risk concentrates"),
    ("scenarios", "5. Climate resilience and scenario analysis"),
    ("metrics", "6. Assets vulnerable to physical climate risk — method not yet agreed"),
    ("trend", "7. Direction of change"),
    ("entity", "8. Sections the entity must complete"),
    ("limitations", "9. Limitations and required disclosures"),
    ("not-assessed", "10. Hazards not assessed"),
    ("mapping", "Appendix A. Framework mapping"),
    ("data", "Appendix B. Per-asset results"),
]


def pill(text: str) -> str:
    return f'<span class="pill">{esc(text)}</span>'


def callout(title: str, body_html: str, kind: str = "") -> str:
    cls = f"callout {kind}".strip()
    return f'<div class="{cls}"><h4>{esc(title)}</h4>{body_html}</div>'


def dataset_identity(source_dataset: str) -> str:
    """The dataset's identity, without the internal commentary that follows it.

    `source_dataset` in a processed file opens with the dataset path and then continues into
    notes written for whoever maintains the layer -- which round superseded which, a sibling
    variable's broken land mask, why one sector was preferred over another. An auditor needs
    the identity; a customer does not need our notes on a layer that is not in their report.
    """
    head = str(source_dataset or "").split(" -- ")[0]
    return head.split(". ")[0].strip().rstrip(".")


def horizon_pairs(delivery: Delivery) -> List[Tuple[str, int, str]]:
    """[(key, decade, label)] for the three configured horizons, ordered."""
    h = delivery.config["horizons"]
    out = []
    for key in ("short", "medium", "long"):
        if key in h:
            out.append((key, int(h[key]["decade"]), str(h[key].get("label", key))))
    return out


# ---------------------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------------------


def sec_scope(delivery: Delivery, cov: dict) -> str:
    n_loc, n_ast = len(delivery.locations), len(delivery.assets)
    hazards = hazard_layer_ids(delivery)
    non_hazard = list(cov["non_hazard_layers"])
    countries = sorted({c for c in delivery.locations["country"].fillna("") if c})

    body = f"""
<p>This report presents a screening assessment of physical climate hazard exposure for
{n_ast} asset{'s' if n_ast != 1 else ''} at {n_loc} location{'s' if n_loc != 1 else ''}
{('in ' + ', '.join(countries)) if countries else ''}. It is structured on the disclosure
requirements of <strong>IFRS S2 Climate-related Disclosures</strong>; Appendix&nbsp;A maps
each section to the corresponding CDP, ESRS E1 and California SB&nbsp;261 requirements.</p>

<p>The underlying data are impact-model outputs from the Inter-Sectoral Impact Model
Intercomparison Project (ISIMIP), processed into decadal statistics on a 0.5° global grid.
{len(hazards)} hazard layer{'s' if len(hazards) != 1 else ''} were extracted
({', '.join(hazards)}){(', together with one asset-condition layer (' + ', '.join(non_hazard) + ') that measures productivity rather than a hazard') if non_hazard else ''}.</p>
"""

    body += callout(
        "What this assessment is, and what it is not",
        """
<p><strong>It is</strong> an evidence base of modelled hazard exposure and its direction of
change, comparable across sites, hazards, scenarios and decades, suitable for screening a
portfolio and for prioritising where a detailed assessment is warranted.</p>
<p><strong>It is not</strong> a financial impact assessment. Converting hazard exposure into
expected loss requires vulnerability functions, asset values and business assumptions that
are outside this dataset. It is not a site-level engineering assessment: values describe an
area of roughly 1° × 1° and cannot resolve local topography, drainage or flood defences. And
it does not cover every physical hazard — see section 10.</p>
""",
        kind="limit",
    )
    return body


def sec_risks(delivery: Delivery, cov: dict) -> str:
    rows = []
    for fam in cov["covered"]:
        for lid in fam["by"]:
            layer = delivery.layer(lid)
            models = str(layer.get("impact_models") or "")
            gcms = str(layer.get("gcms") or "")
            n_m = len([m for m in models.split(",") if m.strip()])
            n_g = len([g for g in gcms.split(",") if g.strip()])
            rows.append(
                [
                    fam["name"],
                    fam["esrs_class"].title(),
                    layer.get("long_name", ""),
                    layer.get("units", ""),
                    f"{n_m} impact model{'s' if n_m != 1 else ''} × {n_g} climate model{'s' if n_g != 1 else ''}",
                    str(layer.get("scenarios", "")).replace(";", ", "),
                ]
            )
    body = f"""
<p>The following physical hazards were assessed. Classification into acute and chronic
follows the standard disclosure taxonomy used by IFRS S2, CDP and ESRS E1.</p>
{table(
    ["Hazard", "Class", "What is measured", "Units", "Ensemble", "Scenarios"],
    rows,
    caption="Physical hazards assessed for this portfolio.",
)}
"""
    if cov["non_hazard_layers"]:
        names = ", ".join(f"<code>{esc(k)}</code>" for k in cov["non_hazard_layers"])
        body += callout(
            "One layer measures asset condition, not a hazard",
            f"""<p>{names} measures the asset's own productivity under each scenario rather
than a hazard acting on it. It is reported because it is decision-relevant, but it is
excluded from every hazard count and from the vulnerability metric in section 6. It is
included in the Climate Score, which is a presentational choice and is flagged as such.</p>""",
            kind="limit",
        )
    return body


def sec_horizons(delivery: Delivery) -> str:
    hz = horizon_pairs(delivery)
    src = str(delivery.config["horizons"].get("source", "")).lower()
    rows = [[label, f"{dec}s", key] for key, dec, label in hz]
    body = f"""
<p>IFRS S2 requires the entity to define short, medium and long term and to explain how those
definitions link to its own planning horizons. The horizons used throughout this report are:</p>
{table(["Horizon", "Decade reported", "Key"], rows,
       caption="Time horizons used in this report.")}
"""
    if src == "default":
        body += callout(
            "These horizons have not been set by the entity",
            """<p>The horizons above are defaults chosen to span the available projection
period. They are <strong>not</strong> the entity's planning horizons and do not yet satisfy
IFRS S2 paragraph 10(d). Replacing them changes which results lead: a 25-year timber rotation
and a 10-year hold period place "long term" in different decades.</p>""",
            kind="must",
        )
    else:
        body += (
            "<p>These horizons were supplied by the entity and reflect its own planning "
            "and capital-allocation cycle.</p>"
        )
    return body


def tier_gap_note(delivery: Delivery, tier: str) -> str:
    """Explain, in the tables where it shows, why an asset reads "1 of 2 hazards assessed".

    Without this the column is actively misleading: a reader sees an incomplete count next
    to a coastal warehouse and concludes the site was not checked for cyclones, when in fact
    the cyclone layer publishes no high-forcing scenario at all and the site reads a
    perfectly ordinary low exposure at the two tiers where it does exist.
    """
    absent = delivery.layers_absent_at_tier(tier)
    if not absent:
        return ""
    names = ", ".join(f"<code>{esc(l)}</code>" for l in absent)
    return (
        f"<p class='sub'><strong>Why some assets show fewer hazards assessed at this "
        f"tier:</strong> {names} publishes no {esc(TIER_LABELS[tier].lower())} scenario, so "
        f"it contributes to no asset here. Those assets <em>were</em> assessed for that "
        f"hazard at the tiers where it exists; the gap is in the scenario coverage of the "
        f"dataset, not in the assessment of the site.</p>"
    )


def sec_concentration(delivery: Delivery, vf: pd.DataFrame) -> str:
    hz = horizon_pairs(delivery)
    long_dec = hz[-1][1] if hz else max(delivery.decades)
    tier = "high" if "high" in set(vf["scenario_tier"]) else (
        sorted(set(vf["scenario_tier"]))[-1] if len(vf) else "medium"
    )
    sel = vf[(vf["decade"] == long_dec) & (vf["scenario_tier"] == tier)]
    sel = sel.sort_values("worst_percentile", ascending=False, na_position="last")

    bar_rows = [
        (delivery.asset_label(r["asset_id"]), r["worst_percentile"])
        for _, r in sel.iterrows()
    ]
    band_counts: Dict[str, int] = {}
    for _, r in sel.iterrows():
        label = risk_band(r["worst_percentile"])
        if label:
            band_counts[label] = band_counts.get(label, 0) + 1

    tbl = table(
        ["Asset", "Worst hazard", "Percentile", "Hazards assessed"],
        [
            [
                delivery.asset_label(r["asset_id"]),
                r["worst_hazard"] or "—",
                band_cell(r["worst_percentile"]),
                f"{int(r['n_hazards_assessed'])} of {int(r['n_hazards_expected'])}",
            ]
            for _, r in sel.iterrows()
        ],
        align_right=(2,),
        caption=f"Worst assessed hazard per asset, {TIER_LABELS[tier].lower()}, {long_dec}s.",
    )

    return f"""
<p>The figure below ranks each asset by its most exposed hazard at the
{TIER_LABELS[tier].lower()} tier in the {long_dec}s. Percentile is a level, not a change: it
ranks the site against the global land distribution of the 2020s baseline, so 80 means the
site is more exposed than 80% of global land was in the 2020s.</p>
{figure_block(
    F.ranked_bar(bar_rows, title=f"Assets by worst assessed hazard — {TIER_LABELS[tier].lower()}, {long_dec}s"),
    caption=("Bar length and fill both encode the same percentile, so the ranking survives "
             "greyscale printing. Assets with no assessed hazard draw no bar and are "
             "labelled; a blank is not a zero."),
    table_html=tbl,
)}
{tier_gap_note(delivery, tier)}
{figure_block(
    F.band_stack(band_counts, title=f"Portfolio distribution by risk band — {TIER_LABELS[tier].lower()}, {long_dec}s"),
    caption="Bands are ordinal and always run Very Low to Very High, never sorted by size.",
    table_html="",
)}
{F.band_legend()}
"""


def sec_scenarios(delivery: Delivery) -> str:
    by_tier: Dict[str, List[str]] = {}
    for scen in sorted(delivery.values["scenario"].unique()):
        by_tier.setdefault(SCENARIO_TIER.get(scen, "?"), []).append(scen)
    rows = [
        [
            TIER_LABELS.get(t, t),
            ", ".join(by_tier[t]),
            "Yes" if any(s in PARIS_ALIGNED for s in by_tier[t]) else "No",
        ]
        for t in TIER_ORDER
        if t in by_tier
    ]
    # Full provenance goes in a table, not into a sentence. `source_dataset` carries long
    # explanatory notes for some layers and inlining them produced an unreadable paragraph.
    prov_rows = [
        [
            lid,
            "ISIMIP2b (CMIP5 forcing, RCP pathways)"
            if str(delivery.layer(lid).get("scenarios", "")).startswith("rcp")
            else "ISIMIP3b (CMIP6 forcing, SSP pathways)",
            dataset_identity(delivery.layer(lid).get("source_dataset", "")),
        ]
        for lid in delivery.layer_ids
    ]
    rounds = sorted({r[1].split(" (")[0] for r in prov_rows})

    series, panel, kept, dropped = portfolio_score_series(delivery)
    assert_baseline_tier_equality(series)

    if series:
        fig = F.decade_line(
            series,
            title="Portfolio Climate Score by forcing tier",
            y_label="Climate Score — mean risk percentile across an asset's layers",
        )
        note = (
            f"Computed over a balanced panel of {len(panel)} of {len(delivery.assets)} assets "
            f"— those scored in every tier and decade shown."
            + (f" The {', '.join(str(d) + 's' for d in dropped)} "
               f"{'is' if len(dropped) == 1 else 'are'} excluded because not every layer "
               f"publishes {'that decade' if len(dropped) == 1 else 'those decades'}, so a "
               f"score there would average a different number of layers."
               if dropped else "")
            + " All three tiers meet at the 2020s baseline by construction: the baseline "
              "panel is shared across scenarios, so any separation there would be a "
              "composition artifact rather than a climate signal."
        )
    else:
        fig = F.empty_figure(
            "Portfolio Climate Score by forcing tier",
            "No asset is scored in every tier and decade.",
        )
        note = (
            "No balanced panel exists for this portfolio: no asset carries a Climate Score "
            "in every forcing tier and decade, usually because one of its hazards publishes "
            "no high-forcing scenario. A cross-tier comparison would compare different sets "
            "of assets and is therefore not drawn."
        )

    return f"""
<p>Resilience is assessed by comparing the portfolio across a range of forcing scenarios
drawn from {' and '.join(rounds)}. Because the hazards come from
{'two ISIMIP generations and no scenario code spans both' if len(rounds) > 1 else 'a single ISIMIP generation'},
scenarios are grouped into three forcing tiers.</p>
{table(["Forcing tier", "Scenario codes used", "Paris-aligned pathway present"], rows,
       caption="Scenarios analysed, grouped into forcing tiers (IFRS S2 22(b)).")}
{table(["Layer", "Model generation", "Source dataset"], prov_rows,
       caption="Provenance of each layer's projections.")}
<p>The analysis includes a low-forcing pathway consistent with the goals of the Paris
Agreement, a medium pathway and a high pathway, and therefore spans a diverse range as IFRS
S2 22(b)(iii)–(iv) requires. The tiers are approximately comparable only: the RCP and SSP
pathways rest on different generations of climate model and different socioeconomic
assumptions.</p>
{figure_block(fig, caption=note)}
"""


def sec_metrics(delivery: Delivery, vf: pd.DataFrame) -> str:
    """IFRS S2 29(c) — DEFERRED.

    An earlier version of this section published a vulnerable-asset count derived from a
    percentile threshold. That was wrong, and wrong in a way that looked authoritative:
    `percentile` is a global-relative EXPOSURE rank, "vulnerable" is a statement about
    susceptibility to HARM, and no step in this pipeline connects the two. The threshold
    behaved accordingly -- on one example portfolio every asset was vulnerable at 60 and none
    at 80, which is a property of the cut-point rather than of the portfolio.

    So the section reports nothing until the method has been agreed. See
    `report_common.TBD_SECTIONS`.
    """
    hz = horizon_pairs(delivery)
    tiers = [t for t in TIER_ORDER if t in set(vf["scenario_tier"])]
    n_haz = len(hazard_layer_ids(delivery))

    coverage_notes = "".join(
        f"<li><strong>{esc(TIER_LABELS[t])}:</strong> "
        f"{n_haz - len(delivery.layers_absent_at_tier(t))} of {n_haz} hazards — "
        + ", ".join(f"<code>{esc(l)}</code>" for l in delivery.layers_absent_at_tier(t))
        + " publishes no scenario at this forcing level.</li>"
        for t in tiers if delivery.layers_absent_at_tier(t)
    )
    values_supplied = bool(delivery.config.get("asset_values", {}).get("supplied"))

    body = f"""
<p>IFRS S2 paragraph 29(c) requires the <em>amount and percentage of assets vulnerable to
climate-related physical risks</em>. <strong>This assessment does not yet report that
metric.</strong></p>

{tbd_block("vulnerability_metric")}

<p>What <em>is</em> measured, and reported elsewhere in this document:</p>
<ul>
  <li><strong>Exposure level</strong> per site and hazard, ranked against the global 2020s
      land distribution — section 4.</li>
  <li><strong>Direction of change</strong> at each site, with the robustness of each
      trend — section 7.</li>
  <li><strong>Which sites could not be assessed at all</strong>, kept separate from sites
      assessed and found less exposed — sections 4 and 9.</li>
</ul>
"""
    if coverage_notes:
        body += (
            "<p>A further reason this determination is not straightforward for this "
            "portfolio: hazard coverage is uneven across forcing tiers, so assets are not "
            "compared on a common basis.</p>"
            f"<ul>{coverage_notes}</ul>"
        )
    if not values_supplied:
        body += callout(
            "Asset values were not supplied",
            """<p>Even once the method is agreed, the requirement asks for a monetary
<em>amount</em> and not only a count, and no climate model can supply it. Provide each
asset's carrying amount, currency, valuation date and valuation basis with the site list.
The valuation basis matters as much as the figure — book, insured, market and replacement
values give four different answers to the same question.</p>""",
            kind="must",
        )
    return body


def sec_trend(delivery: Delivery, vf: pd.DataFrame) -> str:
    """One trend strip PER HAZARD. Units differ between layers, and putting two units on one
    axis is the dual-axis mistake in another costume."""
    hz = horizon_pairs(delivery)
    long_dec = hz[-1][1] if hz else max(delivery.decades)
    out = [
        f"""<p>Exposure level and its rate of change are different questions and are reported
separately. The strips below show the trend at each site over the expanding window to the
{long_dec}s, one figure per hazard because the units differ between them.</p>
<p>There is no p-value under this data contract. Two slope estimators are fitted — ordinary
least squares and Theil–Sen — because they fail in opposite regimes, and their
<strong>disagreement</strong> is the robustness signal. Where they disagree the bar is
hatched and the trend should not be relied on.</p>"""
    ]
    for lid in hazard_layer_ids(delivery):
        layer = delivery.layer(lid)
        slope_col = str(layer.get("recommended_slope", "ols_slope"))
        units = str(layer.get("slope_units", ""))
        v = delivery.values
        sel = v[(v["layer_id"] == lid) & (v["decade"] == long_dec)]
        # Show each hazard at the highest forcing tier IT publishes -- cyclone tops out at
        # medium, the SSP layers reach high. Forcing every layer to a common tier would
        # simply drop cyclone's trend from the report.
        available = {SCENARIO_TIER.get(s) for s in sel["scenario"].unique()}
        top_tier = next((t for t in reversed(TIER_ORDER) if t in available), None)
        scen = sorted(s for s in sel["scenario"].unique() if SCENARIO_TIER.get(s) == top_tier)
        sel = sel[sel["scenario"].isin(scen)]

        rows = []
        for _, r in sel.iterrows():
            agree = r.get("slopes_agree")
            agree = None if pd.isna(agree) else bool(agree)
            rows.append((delivery.asset_label(r["asset_id"]), r.get(slope_col), agree))
        if not rows:
            continue
        out.append(
            figure_block(
                F.trend_strip(
                    rows,
                    title=f"{layer.get('hazard', lid)} — trend to the {long_dec}s",
                    units=units,
                ),
                caption=(
                    f"Scenario {', '.join(scen)} ({esc(TIER_LABELS.get(top_tier, top_tier))}, "
                    f"the highest this layer publishes). Estimator: "
                    f"<code>{esc(slope_col)}</code>, selected for this layer by measurement "
                    f"— {esc(str(layer.get('recommended_slope_rationale', '')))}"
                ),
            )
        )
    return "\n".join(out)


def sec_entity() -> str:
    return f"""
<p>IFRS S2 has four pillars. This report can populate parts of <em>Strategy</em> and
<em>Metrics and targets</em> from modelled hazard data. The following must be completed by
the entity, because they describe how the entity is run and no climate dataset contains
them.</p>
{table(
    ["IFRS S2 requirement", "What is needed", "Status here"],
    [
        ["Governance (paras 5–7)",
         "The body responsible for oversight of climate risk, how it is informed, how often it considers it",
         "Not supplied — entity-owned"],
        ["Risk management (paras 24–26)",
         "The process for identifying, assessing, prioritising and monitoring climate risk, and how it integrates with overall risk management",
         "Not supplied — entity-owned"],
        ["Strategy: financial position, performance and cash flows (paras 15–21)",
         "How the exposures in this report translate into effects on the financial statements, and the assumptions used",
         "Not supplied — requires vulnerability functions and asset values"],
        ["Strategy: response and adaptation (para 14)",
         "Adaptation actions taken or planned, and resourcing",
         "Not supplied — entity-owned"],
        ["Metrics: targets (paras 33–37)",
         "Any climate-related targets and progress against them",
         "Not supplied — entity-owned"],
    ],
    caption="Requirements this report does not address, and why.",
)}
{callout(
    "A missing pillar is invisible unless it is named",
    "<p>This table exists so that a reader cannot mistake a partial assessment for a "
    "complete disclosure. Publishing sections 1–7 alone, without it, would present hazard "
    "screening as an IFRS S2 filing.</p>",
    kind="must",
)}
"""


def sec_limitations(delivery: Delivery) -> str:
    out = [
        "<p>The following are required disclosures for this assessment. Each carries the "
        "identifier used in <code>caveats.json</code>, so a statement in this report can be "
        "traced to the delivery that produced it.</p>"
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
    should = [c for c in delivery.caveats if c["severity"] == SEVERITY_SHOULD]
    if should:
        out.append("<h3>Further notes on reading specific values</h3>")
        out.append(
            table(
                ["Identifier", "Note"],
                [[c["id"], f"{c['title']}. {c['text']}"] for c in should],
                caption="Caveats affecting how an individual figure should be read.",
            )
        )
    return "\n".join(out)


def sec_not_assessed(delivery: Delivery, cov: dict) -> str:
    # `customer_note`, never `blocker`. The blocker field is our engineering note about what
    # would close the gap and reads as internal chatter in a filing.
    rows = [
        [u["name"], u["esrs_class"].title(), u["customer_note"]]
        for u in cov["uncovered"]
    ]
    n_cov, n_all = len(cov["covered"]), cov["n_families"]
    return f"""
{callout(
    f"{n_all - n_cov} of {n_all} physical hazard families were not assessed",
    "<p>Their absence from this report is <strong>not</strong> a finding that they are "
    "immaterial. They were not examined. A reader comparing sites in this report is "
    "comparing them on the hazards in section 2 and on nothing else — in particular, "
    "riverine flood, coastal flood and extreme heat are not covered, and for many built "
    "assets one of those is the dominant physical risk.</p>",
    kind="must",
)}
{table(["Hazard family", "Class", "Status and what it means for this portfolio"], rows,
       caption=f"Physical hazard families not covered by this assessment ({len(rows)} of {n_all}).")}
"""


def sec_mapping(delivery: Delivery, cov: dict) -> str:
    s2 = table(
        ["IFRS S2", "Requirement", "Where"],
        [
            ["10, 10(a)", "Climate risks that could affect prospects; physical or transition", "Section 2"],
            ["10(d)", "Definitions of short, medium and long term", "Section 3"],
            ["13(b)", "Where risks are concentrated in the business model", "Section 4"],
            ["22(a)–(b)", "Climate resilience assessed with scenario analysis; scenarios, sources, horizons", "Section 5"],
            ["29(c)", "Amount and percentage of assets vulnerable to physical risks",
             "Section 6 — NOT REPORTED, method not yet agreed"],
            ["28–29", "Basis of preparation, method and limitations of the metrics", "Sections 1, 9"],
            ["5–7, 24–26, 14–21, 33–37", "Governance, risk management, financial effects, targets", "Section 8 — not supplied"],
        ],
        caption="IFRS S2 requirement to report section.",
    )
    cdp = table(
        ["CDP 3.1.1 field", "Supplied here", "Source"],
        [
            ["Primary climate risk driver (hazard)", "Yes", "Section 2 hazard list"],
            ["Where in the value chain", "Partly — direct operations only", "Own sites; upstream and downstream not assessed"],
            ["Country / area", "Yes", "locations.csv"],
            ["Time horizon", "Yes", "Section 3 horizons"],
            ["Likelihood", "Proxy only", "Percentile is an exposure level, not an event probability"],
            ["Magnitude of impact", "Proxy only", "Percentile band; not a financial magnitude"],
            ["Financial impact figure", "No", "Entity-owned — requires asset values and vulnerability functions"],
            ["Primary financial effect", "No", "Entity-owned"],
            ["Explanation of the risk", "Partly", "Method here; business consequence is entity-owned"],
            ["Cost of response / management method", "No", "Entity-owned"],
        ],
        caption="CDP question 3.1.1 per-risk fields. Ten fields are required for full credit; "
                "this assessment supplies the hazard and exposure half.",
    )
    esrs = table(
        ["ESRS E1-9 datapoint", "Supported", "Note"],
        [
            ["Assets at material physical risk — count and share", "Not yet",
             "Section 6 — requires an agreed definition of 'at material risk'"],
            ["Assets at material physical risk — monetary amount", "Not yet",
             "Needs both the definition above and Asset_Value with the site list"],
            ["Disaggregation by acute and chronic", "Yes", "Class column in sections 2 and 10"],
            ["Location of significant assets by NUTS 3 code", "No",
             "Decimal coordinates are delivered for every site; NUTS classification is not applied"],
            ["Short, medium and long term breakdown", "Yes", "Section 3"],
            ["Before adaptation actions", "Yes", "No adaptation is modelled anywhere in this assessment"],
            ["Anticipated financial effects", "No", "Entity-owned"],
        ],
        caption="ESRS E1-9 coverage.",
    )
    sb261 = table(
        ["California SB 261", "Where"],
        [
            ["Climate-related financial risk, physical — exposure disclosed", "Sections 2, 4, 7"],
            ["Assets vulnerable — quantified", "Section 6 — not yet reported"],
            ["Framework used", "IFRS S2, as an accepted framework under the statute"],
            ["Measures adopted to reduce and adapt to risk", "Section 8 — entity-owned"],
            ["Gaps in disclosure and plans to address them", "Sections 8, 9, 10"],
        ],
        caption="California SB 261 index.",
    )
    note = ""
    if not cov["cdp_labels_verified"]:
        # Deliberately NOT `cov['cdp_labels_note']` -- that is the internal bookkeeping
        # note, and it renders as an alarm in a filing.
        note = callout(
            "The CDP hazard labels above are a mapping, not a transcription",
            "<p>Hazard families have been matched to CDP's question 3.1.1 vocabulary by "
            "name rather than copied from the current questionnaire. Reconcile the two "
            "lists against the live form before submitting, since CDP revises its hazard "
            "options between cycles.</p>",
            kind="limit",
        )
    return f"<h3>IFRS S2</h3>{s2}<h3>CDP</h3>{cdp}{note}<h3>ESRS E1-9</h3>{esrs}<h3>California SB 261</h3>{sb261}"


def sec_data(delivery: Delivery, vf: pd.DataFrame) -> str:
    hz = horizon_pairs(delivery)
    tiers = [t for t in TIER_ORDER if t in set(vf["scenario_tier"])]
    rows = []
    for _, a in delivery.assets.iterrows():
        loc = delivery.location_of(a["asset_id"])
        for key, dec, label in hz:
            for t in tiers:
                m = vf[
                    (vf["asset_id"] == a["asset_id"])
                    & (vf["decade"] == dec)
                    & (vf["scenario_tier"] == t)
                ]
                if m.empty:
                    continue
                r = m.iloc[0]
                rows.append(
                    [
                        loc["name"],
                        a["asset_type"],
                        f"{loc['lat']:.4f}, {loc['lon']:.4f}",
                        label,
                        TIER_LABELS[t],
                        r["worst_hazard"] or "—",
                        band_cell(r["worst_percentile"]),
                        f"{int(r['n_hazards_assessed'])} of {int(r['n_hazards_expected'])}",
                    ]
                )
    return f"""
<p>Every figure in this report is derived from the delivered CSV extract. The table below is
the per-asset summary; <code>values.csv</code> in the same delivery folder carries the full
result for every hazard, scenario and decade, together with confidence intervals, both slope
estimators and per-cell ensemble depth.</p>
{table(
    ["Location", "Asset", "Coordinates", "Horizon", "Tier", "Most exposed hazard",
     "Percentile", "Hazards assessed"],
    rows,
    align_right=(6,),
    caption="Per-asset exposure at the three reporting horizons. Percentile is a global-relative exposure rank, not a vulnerability determination.",
)}
"""


# ---------------------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------------------


def build_report(delivery: Delivery) -> str:
    taxonomy = load_hazard_taxonomy()
    cov = coverage_summary(delivery, taxonomy)
    vf = vulnerability_frame(delivery, taxonomy=taxonomy)

    parts: Dict[str, str] = {
        "scope": sec_scope(delivery, cov),
        "risks": sec_risks(delivery, cov),
        "horizons": sec_horizons(delivery),
        "concentration": sec_concentration(delivery, vf),
        "scenarios": sec_scenarios(delivery),
        "metrics": sec_metrics(delivery, vf),
        "trend": sec_trend(delivery, vf),
        "entity": sec_entity(),
        "limitations": sec_limitations(delivery),
        "not-assessed": sec_not_assessed(delivery, cov),
        "mapping": sec_mapping(delivery, cov),
        "data": sec_data(delivery, vf),
    }
    body = "\n".join(
        f'<section id="{anchor}"><h2>{esc(heading)}</h2>{parts[anchor]}</section>'
        for anchor, heading in SECTIONS
    )
    stamp = build_stamp(body, Path(__file__).read_text())
    html = document(
        title="Physical Climate Risk Assessment",
        subtitle="Prepared on the IFRS S2 disclosure structure, with CDP, ESRS E1 and "
                 "California SB 261 mapping",
        delivery=delivery,
        body=body,
        stamp=stamp,
        toc=[(a, h) for a, h in SECTIONS],
    )
    assert_caveats_present(delivery, html)
    return html


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("delivery", type=Path)
    ap.add_argument("--check-stamp", action="store_true",
                    help="print the build stamp of the file on disk and exit")
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
        html = build_report(delivery)
    except DeliveryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    path.write_text(html)
    stamp = read_stamp(html)
    record_stage(args.delivery, "compliance_report", "built",
                 f"{REPORT_FILENAME} ({len(html) // 1024} KB), build {stamp}")
    if not args.quiet:
        print(f"  compliance report: {REPORT_FILENAME}  ({len(html) // 1024} KB)  build {stamp}")
        print(f"    carries {len(must_disclose(delivery))} must-disclose caveats")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
