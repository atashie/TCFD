"""Shared machinery for customer reports: loading, shell, citations, vulnerability.

Both Stage 3 reports are built from this module so that the compliance document and the
bespoke document cannot disagree about a number. They differ in SPINE -- one is organised by
IFRS S2 paragraph, the other by the reader's decision -- and in nothing else. Anything a
report states about the data comes through here.

THE LOAD-BEARING RULE, AND EXACTLY HOW FAR THE MACHINERY ENFORCES IT
--------------------------------------------------------------------
The rule a writer must follow: a report contains only sentences derived from the delivery,
or sentences carrying a citation that resolves.

WHAT THE CODE ACTUALLY CHECKS is narrower, and the gap matters:

  * every PARAGRAPH in a slot marked `requires_citation` contains at least one citation --
    not every sentence, and not every claim;
  * a `[data:...]` marker names a row that EXISTS -- not that the row contains the number
    the sentence states, and not that the cited column is the one being quoted;
  * a `[dossier:...]` marker names a source id that exists in `dossier.yaml` -- the source
    is never fetched, and `what_it_supports` is never compared against the claim.

So a paragraph can carry one resolving citation and several unsupported assertions and pass.
The checks raise the cost of inventing something and make a fabricated reference impossible;
they do not make an unsupported sentence impossible. **The remaining gap is closed by
reading, not by the build.**

That gap is worth being precise about, because the failure mode is FLUENCY: a well-written
paragraph about a customer's operations reads exactly the same whether it was researched or
invented, and by the time it reaches a disclosure nobody can tell which. A guarantee
overstated here would be the same category of error.

EXPOSURE IS COMPUTED IN ONE PLACE, AND VULNERABILITY IS NOT PUBLISHED AT ALL
----------------------------------------------------------------------------
`vulnerability_frame()` is the only implementation. What the reports USE from it is factual:
the worst exposure percentile per asset, which hazard produced it, how many hazards were
assessed, and the NOT_ASSESSED distinction. It excludes non-hazard layers (`conifer-npp`
measures the stand's productivity, not a hazard acting on it -- see
config/hazard_taxonomy.yaml).

Its `VULNERABLE` / `NOT_VULNERABLE` statuses are NOT rendered anywhere. The method for
turning a global-relative exposure rank into a claim about susceptibility to harm has not
been agreed -- see `TBD_SECTIONS` and docs/reporting/compliance/vulnerability-definition.md
-- and the verifier fails any report that publishes such a determination.
"""

from __future__ import annotations

import hashlib
import html
import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
import yaml

from . import report_figures as figures
from .delivery import (
    DeliveryError,
    PROJECT_ROOT,
    load_report_config,
)
from .viz_common import (
    DIVERGING_HIGH,
    DIVERGING_LOW,
    DIVERGING_MID_DARK,
    DIVERGING_MID_LIGHT,
    FONT_STACK,
    INK,
    ORDINAL_RISK,
    ORDINAL_RISK_DARK,
    SCENARIO_TIER,
    TIER_COLOR_DARK,
    TIER_COLOR_LIGHT,
    TIER_LABELS,
    TIER_ORDER,
    risk_band,
)

HAZARD_TAXONOMY_PATH = PROJECT_ROOT / "config" / "hazard_taxonomy.yaml"

DOSSIER_FILENAME = "dossier.yaml"
CAVEATS_JSON = "caveats.json"
CAVEATS_MD = "caveats.md"

#: Caveats at this severity MUST appear in every report built for the delivery. The verifier
#: checks it; the report builder refuses without it.
SEVERITY_MUST = "must_disclose"
SEVERITY_SHOULD = "should_note"
SEVERITY_INFO = "informational"
SEVERITIES = (SEVERITY_MUST, SEVERITY_SHOULD, SEVERITY_INFO)


# ---------------------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------------------


@dataclass
class Delivery:
    """Everything Stage 2 produced, plus the report configuration and (once built) caveats."""

    path: Path
    manifest: dict
    config: dict
    locations: pd.DataFrame
    assets: pd.DataFrame
    layers: pd.DataFrame
    values: pd.DataFrame
    climate_score: pd.DataFrame
    caveats: List[dict] = field(default_factory=list)
    dossier: dict = field(default_factory=dict)

    # -- convenience ---------------------------------------------------------------------
    @property
    def customer(self) -> str:
        return self.manifest.get("customer", "")

    @property
    def layer_ids(self) -> List[str]:
        return list(self.layers["layer_id"])

    def layer(self, layer_id: str) -> dict:
        rows = self.layers[self.layers["layer_id"] == layer_id]
        if rows.empty:
            raise DeliveryError(f"Layer {layer_id!r} is not in this delivery.")
        return rows.iloc[0].to_dict()

    def location_of(self, asset_id: str) -> dict:
        aid = self.assets[self.assets["asset_id"] == asset_id].iloc[0]
        loc = self.locations[self.locations["location_id"] == aid["location_id"]].iloc[0]
        return loc.to_dict()

    def asset_label(self, asset_id: str) -> str:
        a = self.assets[self.assets["asset_id"] == asset_id].iloc[0]
        loc = self.location_of(asset_id)
        # `or ""` is not enough: an empty optional column reads back as float('nan'),
        # which is TRUTHY, so it survives to str() and prints "nan" in a customer document.
        unit_raw = a.get("sub_asset_unit")
        unit = "" if pd.isna(unit_raw) else str(unit_raw).strip()
        base = f"{loc['name']} — {a['asset_type']}"
        return f"{base} ({unit})" if unit else base

    def tiers_by_layer(self) -> Dict[str, set]:
        """{layer_id: {tiers it publishes a scenario in}}."""
        return {
            lid: {
                SCENARIO_TIER.get(s)
                for s in self.values[self.values["layer_id"] == lid]["scenario"].unique()
            }
            for lid in self.layer_ids
        }

    def layers_absent_at_tier(self, tier: str, hazards_only: bool = True) -> List[str]:
        """Layers with NO scenario in `tier`.

        A hazard missing from a tier is a scenario-availability fact about the dataset, and
        it looks identical in a results table to a site that was not assessed. Naming it is
        what keeps a reader from concluding the site was never checked for that hazard.
        """
        pool = hazard_layer_ids(self) if hazards_only else self.layer_ids
        by_layer = self.tiers_by_layer()
        return [lid for lid in pool if tier not in by_layer.get(lid, set())]

    @property
    def decades(self) -> List[int]:
        return sorted(int(d) for d in self.values["decade"].unique())


def load_delivery(path: Path) -> Delivery:
    """Read a Stage 2 delivery. Refuses one marked incomplete."""
    path = Path(path)
    if not (path / "manifest.json").exists():
        raise DeliveryError(f"{path} has no manifest.json -- not a delivery folder.")
    if (path / "DELIVERY-INCOMPLETE.md").exists():
        raise DeliveryError(
            f"{path} carries DELIVERY-INCOMPLETE.md. Stage 2 did not finish, so a report "
            f"built from it would be a finished-looking document over an unfinished "
            f"extract. Fix and re-run the extract first."
        )
    manifest = json.loads((path / "manifest.json").read_text())

    def _csv(name: str, required: bool = True) -> pd.DataFrame:
        p = path / name
        if not p.exists():
            if required:
                raise DeliveryError(f"{p} is missing.")
            return pd.DataFrame()
        return pd.read_csv(p)

    caveats: List[dict] = []
    cav_path = path / CAVEATS_JSON
    if cav_path.exists():
        caveats = json.loads(cav_path.read_text()).get("caveats", [])

    dossier: dict = {}
    dos_path = path / DOSSIER_FILENAME
    if dos_path.exists():
        dossier = yaml.safe_load(dos_path.read_text()) or {}

    return Delivery(
        path=path,
        manifest=manifest,
        config=load_report_config(path),
        locations=_csv("locations.csv"),
        assets=_csv("assets.csv"),
        layers=_csv("layers.csv"),
        values=_csv("values.csv"),
        climate_score=_csv("climate_score.csv", required=False),
        caveats=caveats,
        dossier=dossier,
    )


def load_hazard_taxonomy() -> dict:
    if not HAZARD_TAXONOMY_PATH.exists():
        raise DeliveryError(f"{HAZARD_TAXONOMY_PATH} is missing.")
    return yaml.safe_load(HAZARD_TAXONOMY_PATH.read_text())


def hazard_layer_ids(delivery: Delivery, taxonomy: Optional[dict] = None) -> List[str]:
    """The delivery's layers that are HAZARDS, in registry order.

    `conifer-npp` is excluded here and that is the whole point of the function: it is
    `higher_is_better` with an already-inverted percentile, so it reads on the same 1-100
    risk axis as a hazard and would silently become a fourth hazard in every count.
    """
    tax = taxonomy or load_hazard_taxonomy()
    non_hazard = set(tax.get("non_hazard_layers") or {})
    return [lid for lid in delivery.layer_ids if lid not in non_hazard]


def coverage_summary(delivery: Delivery, taxonomy: Optional[dict] = None) -> dict:
    """Which hazard families this delivery evidences, and which it does not.

    Scoped to the delivery: a family covered by a registry layer that this customer's assets
    did not pull is NOT covered for this customer, and reporting otherwise would be the
    exact overclaim the taxonomy exists to prevent.
    """
    tax = taxonomy or load_hazard_taxonomy()
    present = set(delivery.layer_ids)
    covered, uncovered = [], []
    for key, fam in tax["families"].items():
        entry = {
            "id": key,
            "name": fam["name"],
            "esrs_class": fam["esrs_class"],
            "cdp_label": fam.get("cdp_label", ""),
            "by": [l for l in fam.get("covered_by") or [] if l in present],
            "registry_has": list(fam.get("covered_by") or []),
            "coverage_note": fam.get("coverage_note", ""),
            "materiality_note": fam.get("materiality_note", ""),
            # Customer-facing prose. `blocker` is the INTERNAL note -- what it would take to
            # close the gap, written for whoever would do the work -- and must never be
            # rendered into a customer document: it carries repository paths, dataset
            # defects and words like UNVERIFIED that mean something to us and something
            # alarming to a reader.
            "customer_note": fam.get("customer_note", "") or fam.get("materiality_note", ""),
            "isimip_candidate": fam.get("isimip_candidate") or [],
            "blocker": fam.get("blocker", ""),
        }
        (covered if entry["by"] else uncovered).append(entry)
    return {
        "covered": covered,
        "uncovered": uncovered,
        "n_families": len(tax["families"]),
        "cdp_labels_verified": bool(tax.get("cdp_labels_verified")),
        "cdp_labels_note": tax.get("cdp_labels_note", ""),
        "non_hazard_layers": {
            k: v for k, v in (tax.get("non_hazard_layers") or {}).items()
            if k in present
        },
    }


# ---------------------------------------------------------------------------------------
# Deferred decisions
# ---------------------------------------------------------------------------------------
#
# A REPORT SECTION IS ALLOWED TO SAY "NOT YET DECIDED", AND SHOULD.
#
# The pressure on this pipeline is to fill every box a framework defines, because a complete
# -looking report is what a customer expects to receive. That pressure is exactly what
# produces a number nobody chose deliberately -- and once such a number is in a filing it is
# indistinguishable from one that was reasoned about.
#
# So: where the method for translating our data into a framework's requirement has NOT been
# thought through and agreed, the section states that plainly and reports nothing. Err
# toward an explicit gap over a fast answer. A gap is a conversation; a wrong number is a
# liability with the customer's name on it.
#
# Move an entry out of here only after the decision has actually been made WITH the user,
# and record the decision where the code implements it.

TBD_SECTIONS: Dict[str, dict] = {
    "vulnerability_metric": {
        "requirement": "IFRS S2 paragraph 29(c) / ESRS E1-9 — the amount and percentage of "
                       "assets vulnerable to climate-related physical risks",
        "why_deferred": (
            "Our data measures EXPOSURE — where a site sits in the global distribution of a "
            "modelled hazard — and 'vulnerable' is a statement about susceptibility to harm. "
            "Converting one into the other is a judgement about each asset type, not a "
            "calculation, and it has not been made."
        ),
        "decisions_outstanding": [
            "What counts as vulnerable for each asset type. A timber stand at the 85th "
            "percentile for drought may be well inside its normal operating range while a "
            "data centre at the 60th is not, so a single portfolio-wide cut-point encodes an "
            "assumption we have not tested.",
            "Whether a global-relative rank is the right basis at all. The alternatives — "
            "change from the site's own baseline, or an absolute physical threshold such as "
            "a burnt-area fraction — answer materially different questions and would rank "
            "the same portfolio differently.",
            "How to combine hazards. Taking the worst treats them as interchangeable; any "
            "weighting is a materiality judgement that belongs to the asset owner.",
            "Whether a determination is meaningful when hazard coverage is uneven — cyclone "
            "publishes no high-forcing scenario, and different asset types carry different "
            "hazard sets, so assets are not being judged on a common basis.",
            "Whether to report at all without asset values, given that the requirement asks "
            "for a monetary amount and not only a count.",
        ],
        "reported_instead": (
            "Exposure level per site and hazard, its rank within this portfolio, and its "
            "direction of change — all of which are measured rather than inferred."
        ),
    },
}


def tbd_block(key: str) -> str:
    """Render a deferred decision as a report section. Publishes no number."""
    spec = TBD_SECTIONS[key]
    items = "".join(f"<li>{esc(d)}</li>" for d in spec["decisions_outstanding"])
    return (
        '<div class="callout must">'
        "<h4>Not yet determined — method to be agreed</h4>"
        f"<p><strong>Requirement:</strong> {esc(spec['requirement'])}</p>"
        f"<p>{esc(spec['why_deferred'])}</p>"
        "<p><strong>Decisions outstanding:</strong></p>"
        f"<ol>{items}</ol>"
        f"<p><strong>Reported instead:</strong> {esc(spec['reported_instead'])}</p>"
        "<p>This section is deliberately blank rather than populated with a provisional "
        "figure. A count of vulnerable assets is a monotone function of a threshold "
        "somebody chooses, and publishing one before that choice has been made and agreed "
        "would put an unowned number into a disclosure.</p>"
        "</div>"
    )


# ---------------------------------------------------------------------------------------
# Exposure (the measured half of IFRS S2 29(c))
# ---------------------------------------------------------------------------------------

NOT_ASSESSED = "NOT_ASSESSED"
VULNERABLE = "VULNERABLE"
NOT_VULNERABLE = "NOT_VULNERABLE"


def vulnerability_frame(
    delivery: Delivery,
    threshold: Optional[float] = None,
    taxonomy: Optional[dict] = None,
) -> pd.DataFrame:
    """One row per (asset, forcing tier, decade): is this asset vulnerable, and on what.

    Method, and every step of it is a choice worth stating in the report:

      - HAZARD LAYERS ONLY. Asset-condition layers are excluded (see `hazard_layer_ids`).
      - Within a hazard, native scenario codes belonging to the same tier are AVERAGED --
        the same two-stage estimator the Climate Score uses, so the two cannot disagree.
      - Across hazards, the MAXIMUM is taken, not the mean. Vulnerability is not an average
        condition: an asset exposed to one severe hazard and three mild ones is vulnerable,
        and a mean would dilute exactly the case that matters.
      - A hazard with no finite percentile contributes nothing and is counted as not
        assessed. An asset with NO assessed hazard is `NOT_ASSESSED`, never
        `NOT_VULNERABLE` -- the difference between "we looked and it is fine" and "we did
        not look" is the difference between a disclosure and a misstatement.
    """
    if threshold is None:
        threshold = float(delivery.config["vulnerability"]["threshold"])
    hazards = set(hazard_layer_ids(delivery, taxonomy))

    v = delivery.values.copy()
    v = v[v["layer_id"].isin(hazards)]
    if v.empty:
        return pd.DataFrame(
            columns=[
                "asset_id", "scenario_tier", "decade", "status", "worst_percentile",
                "worst_hazard", "n_hazards_assessed", "n_hazards_expected",
            ]
        )
    v["scenario_tier"] = v["scenario"].map(SCENARIO_TIER)

    expected = {
        r["asset_id"]: len(
            [l for l in str(r["layer_ids"]).split(";") if l and l in hazards]
        )
        for _, r in delivery.assets.iterrows()
    }

    # Stage 1: mean across native codes inside a tier, per hazard.
    per_hazard = (
        v.groupby(["asset_id", "layer_id", "scenario_tier", "decade"], dropna=True)[
            "percentile"
        ]
        .mean()
        .reset_index()
    )

    rows: List[dict] = []
    for (asset_id, tier, decade), grp in per_hazard.groupby(
        ["asset_id", "scenario_tier", "decade"]
    ):
        finite = grp[grp["percentile"].notna()]
        n_assessed = int(len(finite))
        if not n_assessed:
            rows.append(
                {
                    "asset_id": asset_id, "scenario_tier": tier, "decade": int(decade),
                    "status": NOT_ASSESSED, "worst_percentile": None, "worst_hazard": "",
                    "n_hazards_assessed": 0,
                    "n_hazards_expected": expected.get(asset_id, 0),
                }
            )
            continue
        top = finite.loc[finite["percentile"].idxmax()]
        rows.append(
            {
                "asset_id": asset_id,
                "scenario_tier": tier,
                "decade": int(decade),
                "status": VULNERABLE
                if float(top["percentile"]) >= threshold
                else NOT_VULNERABLE,
                "worst_percentile": float(top["percentile"]),
                "worst_hazard": str(top["layer_id"]),
                "n_hazards_assessed": n_assessed,
                "n_hazards_expected": expected.get(asset_id, 0),
            }
        )
    return pd.DataFrame(rows).sort_values(["asset_id", "scenario_tier", "decade"])


def vulnerability_rollup(
    delivery: Delivery,
    tier: str,
    decade: int,
    threshold: Optional[float] = None,
    frame: Optional[pd.DataFrame] = None,
) -> dict:
    """Counts, percentages and (when values are supplied) the monetary amount, for 29(c)."""
    if threshold is None:
        threshold = float(delivery.config["vulnerability"]["threshold"])
    vf = frame if frame is not None else vulnerability_frame(delivery, threshold)
    sel = vf[(vf["scenario_tier"] == tier) & (vf["decade"] == decade)]

    n_total = int(len(delivery.assets))
    n_assessed = int((sel["status"] != NOT_ASSESSED).sum())
    n_vuln = int((sel["status"] == VULNERABLE).sum())

    values = pd.to_numeric(delivery.assets.get("asset_value"), errors="coerce") \
        if "asset_value" in delivery.assets.columns else pd.Series(dtype=float)
    have_values = bool(values.notna().any())
    amount_vuln = amount_assessed = amount_total = None
    if have_values:
        by_id = dict(zip(delivery.assets["asset_id"], values))
        amount_total = float(values.sum(skipna=True))
        amount_vuln = float(
            sum(by_id.get(a) or 0.0 for a in sel.loc[sel["status"] == VULNERABLE, "asset_id"])
        )
        amount_assessed = float(
            sum(by_id.get(a) or 0.0 for a in sel.loc[sel["status"] != NOT_ASSESSED, "asset_id"])
        )

    return {
        "tier": tier,
        "decade": decade,
        "threshold": threshold,
        "n_assets_total": n_total,
        "n_assets_assessed": n_assessed,
        "n_assets_not_assessed": n_total - n_assessed,
        "n_assets_vulnerable": n_vuln,
        # Denominator is the ASSESSED population, and it is named in the key so a consumer
        # cannot mistake it for a share of the portfolio. An unassessed asset is not
        # evidence of safety and must not dilute the rate.
        "pct_of_assessed_vulnerable": (100.0 * n_vuln / n_assessed) if n_assessed else None,
        "pct_of_portfolio_vulnerable": (100.0 * n_vuln / n_total) if n_total else None,
        "monetary_available": have_values,
        "amount_vulnerable": amount_vuln,
        "amount_assessed": amount_assessed,
        "amount_total": amount_total,
        "currency": delivery.config.get("asset_values", {}).get("currency"),
    }


def vulnerability_sensitivity(
    delivery: Delivery, tier: str, decade: int
) -> List[dict]:
    """The same count at the chosen threshold and at its neighbours. Always reported."""
    cfg = delivery.config["vulnerability"]
    thresholds = sorted({float(cfg["threshold"]), *(float(t) for t in cfg["sensitivity"])})
    out = []
    for t in thresholds:
        r = vulnerability_rollup(delivery, tier, decade, threshold=t)
        r["is_chosen"] = t == float(cfg["threshold"])
        out.append(r)
    return out


# ---------------------------------------------------------------------------------------
# Balanced panels
# ---------------------------------------------------------------------------------------


def balanced_assets(
    delivery: Delivery,
    tiers: Sequence[str],
    decades: Sequence[int],
    location_id: Optional[str] = None,
) -> List[str]:
    """Assets carrying a Climate Score in EVERY (tier, decade) cell being compared.

    Any cross-tier or cross-decade rollup must be restricted to these. Averaging each tier
    over whatever assets it happens to have makes the tiers differ for compositional reasons
    that look exactly like a climate signal -- `cyclone` publishes no high-forcing scenario,
    so cyclone-carrying assets vanish from the high tier and drag its mean somewhere else.

    Measured twice on this codebase, because fixing it at one level of rollup did not fix it
    at the next: the portfolio baseline read 39.9 at the high tier against 42.1 at the low,
    and after that was corrected the Shasta location read 62.3 against 51.7 for the same
    reason one level down. Hence `location_id` -- the restriction has to be applied at EVERY
    level of aggregation, not at the level where it first bit.
    """
    s = delivery.climate_score
    if s.empty:
        return []
    if location_id is not None:
        ids = set(
            delivery.assets.loc[delivery.assets["location_id"] == location_id, "asset_id"]
        )
        s = s[s["asset_id"].isin(ids)]
    need = {(t, int(d)) for t in tiers for d in decades}
    out = []
    for asset_id, grp in s.groupby("asset_id"):
        have = {(t, int(d)) for t, d in zip(grp["scenario_tier"], grp["decade"])}
        if need <= have:
            out.append(str(asset_id))
    return sorted(out)


def portfolio_score_series(
    delivery: Delivery, location_id: Optional[str] = None
) -> Tuple[Dict[str, List[Tuple[int, float]]], List[str], List[int], List[int]]:
    """Mean Climate Score by tier and decade over a balanced panel of COMPLETE rows.

    Returns (series, panel asset ids, decades plotted, decades dropped).

    Two restrictions, both necessary, and dropping either one produces a chart that shows a
    climate trend which is really a bookkeeping change:

      - ASSETS: only those present in every tier-decade cell (see `balanced_assets`).
      - DECADES: only those where every panel asset scored over its FULL hazard set in every
        tier. ISIMIP3b layers publish no 2010s panel, so a 2010s score can rest on one
        hazard where the 2020s rests on three. Measured on the example portfolio that reads
        30.4 → 39.9, a 31% jump entirely manufactured by the second and third hazards
        arriving.

    An empty panel returns empty series rather than an unbalanced fallback: a chart that
    silently swaps in a different comparison is worse than one that says it cannot be drawn.
    """
    s = delivery.climate_score
    if s.empty:
        return {}, [], [], []
    tiers = [t for t in TIER_ORDER if t in set(s["scenario_tier"])]
    all_decades = sorted(int(d) for d in s["decade"].unique())
    panel = balanced_assets(delivery, tiers, all_decades, location_id)
    if not panel:
        return {}, [], [], all_decades

    sel = s[s["asset_id"].isin(panel)]
    complete = sel[sel["n_hazards"] >= sel["n_hazards_expected"]]
    kept = [
        d
        for d in all_decades
        if len(complete[complete["decade"] == d]) == len(panel) * len(tiers)
    ]
    dropped = [d for d in all_decades if d not in kept]

    series: Dict[str, List[Tuple[int, float]]] = {}
    for tier in tiers:
        rows = complete[(complete["scenario_tier"] == tier) & (complete["decade"].isin(kept))]
        series[tier] = [
            (int(d), float(g["climate_score"].mean())) for d, g in rows.groupby("decade")
        ]
    return series, panel, kept, dropped


def assert_baseline_tier_equality(
    series: Dict[str, List[Tuple[int, float]]], baseline_decade: int = 2020, tol: float = 0.05
) -> None:
    """The baseline panel is bit-identical across scenarios, so a baseline Climate Score
    MUST be equal across tiers on a balanced panel. Unequal means a composition artifact --
    the same defect that shipped twice. Refuse rather than publish it."""
    vals = {
        t: v for t, pts in series.items() for d, v in pts if d == baseline_decade
    }
    if len(vals) < 2:
        return
    spread = max(vals.values()) - min(vals.values())
    if spread > tol:
        raise DeliveryError(
            f"Baseline ({baseline_decade}s) Climate Score differs across forcing tiers by "
            f"{spread:.2f} on a balanced panel: "
            + ", ".join(f"{t}={v:.2f}" for t, v in sorted(vals.items()))
            + ".\nThe 2020s panel is bit-identical across scenarios by contract, so this is "
              "a composition artifact, not a climate signal. Do not publish it."
        )


# ---------------------------------------------------------------------------------------
# Citations
# ---------------------------------------------------------------------------------------

CITATION_RE = re.compile(r"\[(data|dossier):([^\]]+)\]")

#: A slot's content is "unfilled" if it is blank or still carries the placeholder. Checked
#: literally so that a half-written narrative cannot render.
TODO_MARKER = "TODO"

SLOT_RE = re.compile(
    r"<!--\s*slot:\s*(?P<name>[a-z0-9_]+)\s*\|\s*requires_citation:\s*(?P<cite>yes|no)\s*-->",
    re.IGNORECASE,
)


#: MULTI-LINE. The scaffold's guidance blocks span many lines, and a line-by-line filter
#: catches only the opening `<!--` and renders the rest -- which would put the facet
#: profiles' internal guidance, including the "Do not claim" list, into the customer's
#: document. That is the precise failure the profile design exists to prevent.
COMMENT_RE = re.compile(r"<!--.*?-->", re.S)


def strip_comments(text: str) -> str:
    return COMMENT_RE.sub("", text)


def parse_narrative(text: str) -> List[dict]:
    """Split a narrative file into its slots, in document order.

    Slot bodies come back with HTML comments REMOVED. Comments are the channel the scaffold
    uses to hand the writer their profile guidance and the citation grammar; they are not
    content. Leaving them in would both render guidance into the report and make the
    citation checker trip over the worked examples in the file header.
    """
    slots: List[dict] = []
    matches = list(SLOT_RE.finditer(text))
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        body = strip_comments(text[m.end():end]).strip()
        slots.append(
            {
                "name": m.group("name"),
                "requires_citation": m.group("cite").lower() == "yes",
                "body": body,
            }
        )
    return slots


def _values_key_exists(delivery: Delivery, ref: str) -> bool:
    parts = ref.split("/")
    if len(parts) != 4:
        return False
    asset_id, layer_id, scenario, decade = parts
    try:
        dec = int(decade)
    except ValueError:
        return False
    v = delivery.values
    return bool(
        (
            (v["asset_id"] == asset_id)
            & (v["layer_id"] == layer_id)
            & (v["scenario"] == scenario)
            & (v["decade"] == dec)
        ).any()
    )


def _score_key_exists(delivery: Delivery, ref: str) -> bool:
    parts = ref.split("/")
    if len(parts) != 3 or delivery.climate_score.empty:
        return False
    asset_id, tier, decade = parts
    try:
        dec = int(decade)
    except ValueError:
        return False
    s = delivery.climate_score
    return bool(
        ((s["asset_id"] == asset_id) & (s["scenario_tier"] == tier) & (s["decade"] == dec)).any()
    )


def _manifest_path_exists(delivery: Delivery, dotted: str) -> bool:
    node = delivery.manifest
    for part in dotted.split("."):
        if isinstance(node, dict) and part in node:
            node = node[part]
        else:
            return False
    return True


def check_citations(delivery: Delivery, text: str) -> List[str]:
    """Every citation in `text` must RESOLVE. Returns the violations.

    A citation that does not resolve is worse than no citation: it looks like evidence.
    """
    problems: List[str] = []
    dossier_ids = {
        str(s.get("id")) for s in (delivery.dossier.get("sources") or []) if s.get("id")
    }
    for kind, ref in CITATION_RE.findall(text):
        ref = ref.strip()
        if kind == "dossier":
            if ref not in dossier_ids:
                problems.append(
                    f"[dossier:{ref}] does not name a source in {DOSSIER_FILENAME}"
                    + (f" (known: {', '.join(sorted(dossier_ids))})" if dossier_ids else
                       f" ({DOSSIER_FILENAME} has no sources, or does not exist)")
                )
            continue
        if "#" not in ref:
            problems.append(f"[data:{ref}] has no '#' -- expected file.csv#key")
            continue
        fname, key = ref.split("#", 1)
        fname = fname.strip()
        key = key.split(":", 1)[0].strip()  # allow file#key:column
        if fname == "values.csv":
            ok = _values_key_exists(delivery, key)
        elif fname == "climate_score.csv":
            ok = _score_key_exists(delivery, key)
        elif fname == "manifest.json":
            ok = _manifest_path_exists(delivery, key)
        elif fname == "caveats.json":
            ok = any(str(c.get("id")) == key for c in delivery.caveats)
        elif fname in {"layers.csv", "assets.csv", "locations.csv"}:
            frame = {"layers.csv": delivery.layers, "assets.csv": delivery.assets,
                     "locations.csv": delivery.locations}[fname]
            idcol = frame.columns[0]
            ok = bool((frame[idcol].astype(str) == key).any())
        else:
            problems.append(f"[data:{ref}] names an unknown file {fname!r}")
            continue
        if not ok:
            problems.append(f"[data:{ref}] does not resolve to a row in {fname}")
    return problems


def check_narrative(delivery: Delivery, text: str) -> List[str]:
    """Slots filled, placeholders gone, required paragraphs cited, citations resolve.

    Everything is checked against the COMMENT-STRIPPED slot bodies, i.e. against exactly the
    text that will be rendered. Checking the raw file instead would flag the citation
    examples in the scaffold's own instructions, and -- far worse -- would let a citation
    that appears only inside a comment satisfy the requirement for a paragraph that ships.
    """
    problems: List[str] = []
    slots = parse_narrative(text)
    if not slots:
        problems.append("narrative.md contains no slot markers -- nothing to render.")
    for slot in slots:
        body = slot["body"]
        if not body:
            problems.append(f"slot '{slot['name']}' is empty")
            continue
        if TODO_MARKER in body:
            problems.append(f"slot '{slot['name']}' still contains a {TODO_MARKER} placeholder")
        if slot["requires_citation"]:
            for para in [p for p in re.split(r"\n\s*\n", body) if p.strip()]:
                stripped = para.strip()
                if stripped.startswith("#"):  # a heading is not a claim
                    continue
                if not CITATION_RE.search(para):
                    head = " ".join(stripped.split())[:70]
                    problems.append(
                        f"slot '{slot['name']}': uncited paragraph -- \"{head}...\""
                    )
        problems.extend(check_citations(delivery, body))
    return problems


# ---------------------------------------------------------------------------------------
# Caveats
# ---------------------------------------------------------------------------------------


def must_disclose(delivery: Delivery) -> List[dict]:
    return [c for c in delivery.caveats if c.get("severity") == SEVERITY_MUST]


#: Files that ship INSIDE the delivery folder. A caveat's `evidence` may cite these in a
#: customer document, because the customer can open them.
DELIVERED_FILES = (
    "manifest.json", "values.csv", "layers.csv", "assets.csv", "locations.csv",
    "climate_score.csv", "report_config.yaml", "caveats.json", "dashboard.html",
    "input_locations.csv",
)


def customer_evidence(evidence: str) -> str:
    """Render a caveat's evidence only where it points at something the customer holds.

    `evidence` is written for traceability and some of it names configuration or modules
    that live in our repository, not in the delivery. A filing that cites a path the reader
    cannot open is not auditable, it is just noise with the shape of rigour -- so those are
    dropped from the document and remain in caveats.json for our own audit trail.

    Filtering is per SEMICOLON-SEPARATED SEGMENT, not per string. Evidence commonly reads
    "manifest.json qa_review; config/layer_registry.yaml qa_reviewed_on" -- a whole-string
    test passes it on the strength of the first half and publishes the second.
    """
    if not evidence:
        return ""
    kept = []
    for seg in evidence.split(";"):
        # Match on a BARE filename token, not a substring. `internal/archive/manifest.json`
        # and `config/values.csv.notes` both contain a delivered filename and neither is a
        # file the customer holds; a substring test published both.
        tokens = re.findall(r"[A-Za-z0-9_.\-/]+", seg)
        if any(t in DELIVERED_FILES for t in tokens):
            kept.append(seg.strip())
    return "; ".join(kept)


def assert_caveats_present(delivery: Delivery, rendered_html: str) -> None:
    """Every `must_disclose` caveat must actually appear in the rendered document."""
    if not delivery.caveats:
        raise DeliveryError(
            f"{delivery.path / CAVEATS_JSON} not found. Stage 4 runs BEFORE the reports -- "
            f"the caveat set is an input to them, not a summary of them. Run:\n"
            f"  python scripts/generate_delivery_caveats.py {delivery.path}"
        )
    missing = [c["id"] for c in must_disclose(delivery) if c["id"] not in rendered_html]
    if missing:
        raise DeliveryError(
            "The report does not carry every must-disclose caveat: "
            + ", ".join(missing)
            + "\nA caveat that is in caveats.json and not in the document is a caveat the "
              "reader never sees."
        )


# ---------------------------------------------------------------------------------------
# Markdown (minimal, deliberately)
# ---------------------------------------------------------------------------------------

_INLINE = (
    (re.compile(r"\*\*(.+?)\*\*"), r"<strong>\1</strong>"),
    (re.compile(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])"), r"<em>\1</em>"),
    (re.compile(r"`([^`]+)`"), r"<code>\1</code>"),
    (re.compile(r"\[([^\]]+)\]\((https?://[^)\s]+)\)"),
     r'<a href="\2" rel="noopener noreferrer">\1</a>'),
)


def _inline(text: str) -> str:
    """Escape first, then apply markup. Never the other way round -- a customer-supplied
    name containing markup must render as text, not as tags."""
    out = html.escape(text, quote=False)
    for pattern, repl in _INLINE:
        out = pattern.sub(repl, out)
    # Numbered citation markers, substituted upstream by `number_citations()`. A hover
    # title is not a citation in a document that will be printed and filed, so each marker
    # is a numbered link into a reference list that survives on paper.
    out = re.sub(
        r"\{\{cite:(\d+)\}\}",
        lambda m: f'<sup class="cite"><a href="#ref-{m.group(1)}">{m.group(1)}</a></sup>',
        out,
    )
    return out


#: A list marker is a bullet or number FOLLOWED BY WHITESPACE. Testing `startswith("*")`
#: instead classifies `**Bold lead-in.**` as a list item and loses the paragraph.
LIST_RE = re.compile(r"^(?:[-*+]\s+|\d+[.)]\s+)")

LEADING_H2_RE = re.compile(r"\A\s*##\s+.*?$\n?", re.M)


def drop_leading_heading(body: str) -> str:
    """Remove a slot's own `## Heading`.

    The scaffold writes one so narrative.md reads as a document on its own, and the report
    generator emits the section heading itself -- without this the heading appears twice.
    """
    return LEADING_H2_RE.sub("", body, count=1).strip()


def number_citations(bodies: Dict[str, str]) -> Tuple[Dict[str, str], List[dict]]:
    """Replace every citation with a numbered marker and return the reference list.

    Numbering runs in document order across all slots, and a repeated citation keeps its
    first number. Returns (rewritten bodies, references).
    """
    order: Dict[str, int] = {}
    refs: List[dict] = []

    def repl(m):
        key = f"{m.group(1)}:{m.group(2).strip()}"
        if key not in order:
            order[key] = len(order) + 1
            refs.append({"n": order[key], "kind": m.group(1), "ref": m.group(2).strip()})
        return f"{{{{cite:{order[key]}}}}}"

    return {k: CITATION_RE.sub(repl, v) for k, v in bodies.items()}, refs


def render_references(delivery: Delivery, refs: List[dict]) -> str:
    """The reference list the numbered markers point at."""
    if not refs:
        return ""
    sources = {str(s.get("id")): s for s in (delivery.dossier.get("sources") or [])}
    rows = []
    for r in refs:
        if r["kind"] == "dossier":
            s = sources.get(r["ref"], {})
            title = str(s.get("title") or r["ref"])
            url = str(s.get("url") or "")
            what = str(s.get("publisher") or "")
            lag = f"as of {s['as_of']}" if s.get("as_of") else ""
            link = (f'<a href="{esc(url)}" rel="noopener noreferrer">{esc(title)}</a>'
                    if url.startswith("http") else esc(title))
            detail = ", ".join(x for x in (what, lag) if x)
            rows.append([raw(f'<span id="ref-{r["n"]}">{r["n"]}</span>'), "Source",
                         raw(f"{link}{(' — ' + esc(detail)) if detail else ''}")])
        else:
            fname, _, key = r["ref"].partition("#")
            pretty = key.replace("/", " · ")
            rows.append([raw(f'<span id="ref-{r["n"]}">{r["n"]}</span>'), "Delivered data",
                         raw(f"<code>{esc(fname)}</code> — {esc(pretty)}")])
    return table(
        ["#", "Kind", "Reference"],
        rows,
        caption="Every numbered claim above resolves to one of these. Delivered-data "
                "references name a row in this delivery's CSV files; source references name "
                "a document recorded in the engagement dossier.",
    )


def markdown(text: str) -> str:
    """Paragraphs, headings, bullet and numbered lists, blockquotes. Nothing else.

    A full markdown engine is not installed and is not needed: the narrative is written
    against this converter, and a converter whose whole grammar fits on one screen cannot
    surprise a document that will be filed with a regulator.
    """
    out: List[str] = []
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            i += 1
            continue
        if stripped.startswith("<!--"):
            i += 1
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            lvl = min(len(m.group(1)) + 2, 6)
            out.append(f"<h{lvl}>{_inline(m.group(2))}</h{lvl}>")
            i += 1
            continue
        if stripped.startswith(">"):
            block = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                block.append(lines[i].strip().lstrip(">").strip())
                i += 1
            out.append(f"<blockquote><p>{_inline(' '.join(block))}</p></blockquote>")
            continue
        if LIST_RE.match(stripped):
            ordered = bool(re.match(r"^\d+[.)]\s+", stripped))
            items: List[str] = []
            while i < len(lines):
                s = lines[i].strip()
                m2 = re.match(r"^(?:[-*+]|\d+[.)])\s+(.*)$", s)
                if m2:
                    items.append(m2.group(1))
                    i += 1
                    continue
                # A CONTINUATION line: non-blank, not a new item, not another block. The
                # earlier version stopped here and let the rest of the item fall through to
                # the paragraph branch, which split one bullet into a <li> plus a stray <p>.
                if s and not s.startswith(("#", ">", "<!--")) and items:
                    items[-1] = f"{items[-1]} {s}"
                    i += 1
                    continue
                break
            tag = "ol" if ordered else "ul"
            out.append(
                f"<{tag}>" + "".join(f"<li>{_inline(it)}</li>" for it in items) + f"</{tag}>"
            )
            continue
        para = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith(("#", ">", "<!--")) or LIST_RE.match(s):
                break
            para.append(s)
            i += 1
        if para:
            out.append(f"<p>{_inline(' '.join(para))}</p>")
        else:
            # Nothing matched and nothing consumed: emit the line rather than skipping it.
            # The previous version advanced past any line it did not recognise, which
            # silently DELETED every paragraph opening with bold text -- `**Finding.** ...`
            # starts with an asterisk, was rejected as a list, rejected as a paragraph, and
            # dropped. Losing a sentence out of a customer document is the worst failure
            # this converter can have, so the fallback is always to render.
            out.append(f"<p>{_inline(stripped)}</p>")
            i += 1
    return "\n".join(out)


# ---------------------------------------------------------------------------------------
# HTML shell
# ---------------------------------------------------------------------------------------


def esc(value) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def assert_narrative_rendered(bodies: Dict[str, str], rendered_html: str) -> None:
    """Every paragraph the author wrote must appear in the document.

    A markdown converter that silently drops a line produces a document that looks finished
    and is missing an argument -- and nothing else in this pipeline would notice, because
    every remaining number is still correct. This shipped once: paragraphs opening with bold
    text were classified as list items, rejected as lists, and discarded.

    The check is on a distinctive prose fragment from each paragraph rather than the whole
    string, since inline markup is transformed on the way through.
    """
    # Whitespace-INSENSITIVE comparison. Replacing a tag with a space is necessary (otherwise
    # `a</p><p>b` becomes `ab`) but it also inserts a space wherever inline markup closes
    # mid-sentence, so `…two sites</strong>, and…` renders as `two sites , and`. Comparing
    # with spaces removed sidesteps every such difference; a 50-character probe cannot
    # collide by accident.
    def squash(text: str) -> str:
        return re.sub(r"\s+", "", text)

    visible = squash(re.sub(r"<[^>]+>", " ", rendered_html))
    missing = []
    for slot, body in bodies.items():
        # LINE level, not paragraph level. Source lines are joined with single spaces on the
        # way out, so each line's own text stays contiguous in the output -- whereas a
        # paragraph-level probe spans list-item boundaries and reports false misses.
        for line in body.split("\n"):
            s = line.strip()
            if not s or s.startswith(("#", ">", "<!--")):
                continue
            # Probe a run of text BETWEEN citation markers. A probe spanning one cannot
            # match: the marker is stripped from the source but renders as a number, so
            # "…derived , which removes…" would be compared against "…derived 12, which
            # removes…" on a paragraph that is present.
            cleaned = [
                " ".join(
                    html.unescape(
                        re.sub(r"[*`_]|^\s*[-+]\s+|^\s*\d+[.)]\s+", "", frag)
                    ).split()
                )
                for frag in re.split(r"\{\{cite:\d+\}\}", s)
            ]
            probe = max(cleaned, key=len) if cleaned else ""
            if len(probe) < 25:
                continue
            if squash(probe[:50]) not in visible:
                missing.append(f"{slot}: \"{probe[:60]}...\"")
    if missing:
        raise DeliveryError(
            f"{len(missing)} narrative paragraph(s) did not reach the rendered report:\n"
            + "\n".join(f"  - {m}" for m in missing[:8])
            + "\nThe markdown converter dropped them. Do not ship a document that is "
              "missing text the author wrote."
        )


def table(
    headers: Sequence[str],
    rows: Sequence[Sequence],
    *,
    caption: str = "",
    align_right: Sequence[int] = (),
    classes: str = "",
) -> str:
    """A data table. Every figure in a report is paired with one of these.

    A figure communicates shape; a filing needs the number. Wrapped in an overflow-x
    container so a wide table scrolls inside itself rather than making the page scroll.
    """
    right = set(align_right)
    num = ' class="num"'
    head = "".join(
        "<th" + (num if i in right else "") + f">{esc(h)}</th>"
        for i, h in enumerate(headers)
    )
    body = "".join(
        "<tr>"
        + "".join(
            "<td"
            + (num if i in right else "")
            + ">"
            + (str(c) if isinstance(c, _Raw) else esc(c))
            + "</td>"
            for i, c in enumerate(row)
        )
        + "</tr>"
        for row in rows
    )
    cap = f"<caption>{esc(caption)}</caption>" if caption else ""
    return (
        f'<div class="table-wrap"><table class="{classes}">{cap}'
        f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"
    )


class _Raw(str):
    """Mark a cell as already-safe HTML (a band swatch, say). Nothing customer-supplied
    should ever be wrapped in this."""


def raw(html_text: str) -> _Raw:
    return _Raw(html_text)


def band_cell(percentile: Optional[float]) -> _Raw:
    """A percentile with its band swatch -- colour plus the label, never colour alone."""
    label = risk_band(percentile)
    if label is None:
        return raw('<span class="muted">not assessed</span>')
    idx = figures.band_index(percentile)
    return raw(
        f'<span class="swatch" style="background:var(--risk-{idx + 1})"></span>'
        f"{figures.fmt(percentile)} <span class='muted'>{esc(label)}</span>"
    )


def figure_block(svg: str, *, caption: str, table_html: str = "") -> str:
    """A figure, its caption, and optionally the table that carries its numbers."""
    return (
        f'<figure class="fig-block">{svg}'
        f"<figcaption>{caption}</figcaption>"
        f"{table_html}</figure>"
    )


def css_tokens() -> str:
    """Light, dark and print token sets. Every colour a figure uses is defined here."""
    L, D = INK["light"], INK["dark"]

    def block(ink, tiers, risk, div_mid):
        rows = [
            f"--ink-surface:{ink['surface']}", f"--ink-plane:{ink['plane']}",
            f"--ink-primary:{ink['primary']}", f"--ink-secondary:{ink['secondary']}",
            f"--ink-muted:{ink['muted']}", f"--ink-grid:{ink['grid']}",
            f"--ink-axis:{ink['axis']}", f"--ink-border:{ink['border']}",
            f"--div-low:{DIVERGING_LOW}", f"--div-mid:{div_mid}",
            f"--div-high:{DIVERGING_HIGH}",
        ]
        rows += [f"--tier-{t}:{tiers[t]}" for t in TIER_ORDER]
        rows += [f"--risk-{i + 1}:{c}" for i, c in enumerate(risk)]
        return ";".join(rows) + ";"

    light = block(L, TIER_COLOR_LIGHT, ORDINAL_RISK, DIVERGING_MID_LIGHT)
    dark = block(D, TIER_COLOR_DARK, ORDINAL_RISK_DARK, DIVERGING_MID_DARK)
    return f"""
:root{{{light}--font:{FONT_STACK};}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{{dark}}}}}
:root[data-theme="dark"]{{{dark}}}
/* Print always uses the light set: a dark-surface figure on paper is a black rectangle,
   and these documents exist to be printed and attached to a filing. */
@media print{{:root{{{light}}}}}
"""


PAGE_CSS = """
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--ink-surface);color:var(--ink-primary);
  font-family:var(--font);font-size:15px;line-height:1.62;
  -webkit-font-smoothing:antialiased}
.page{max-width:62rem;margin:0 auto;padding:2.5rem 1.5rem 5rem}
h1{font-size:1.9rem;line-height:1.2;margin:0 0 .35rem;letter-spacing:-0.01em}
h2{font-size:1.3rem;margin:2.6rem 0 .7rem;padding-top:.9rem;
  border-top:1px solid var(--ink-border);letter-spacing:-0.005em}
h3{font-size:1.05rem;margin:1.7rem 0 .5rem}
h4{font-size:.95rem;margin:1.2rem 0 .4rem;color:var(--ink-secondary)}
p,li{color:var(--ink-primary)}
a{color:var(--div-low)}
.muted{color:var(--ink-muted)}
.lede{font-size:1.05rem;color:var(--ink-secondary)}
.sub{color:var(--ink-muted);font-size:.85rem}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.86em;
  background:var(--ink-plane);padding:.1em .35em;border-radius:3px}
blockquote{margin:1rem 0;padding:.65rem 1rem;border-left:3px solid var(--ink-axis);
  background:var(--ink-plane);border-radius:0 4px 4px 0}
blockquote p{margin:0;color:var(--ink-secondary)}

/* Figures ------------------------------------------------------------------ */
.fig-block{margin:1.5rem 0;padding:0}
.fig{display:block;overflow:visible}
.fig-title{fill:var(--ink-primary)}
.fig-label{fill:var(--ink-secondary)}
.fig-tick,.fig-axis-title{fill:var(--ink-muted)}
.fig-value{fill:var(--ink-primary)}
.fig-value-on-fill{fill:var(--ink-surface)}
.fig-null{fill:var(--ink-muted);font-style:italic}
figcaption{color:var(--ink-muted);font-size:.83rem;margin:.5rem 0 0;line-height:1.5}

/* Tables ------------------------------------------------------------------- */
.table-wrap{overflow-x:auto;margin:1rem 0;-webkit-overflow-scrolling:touch}
table{border-collapse:collapse;width:100%;font-size:.86rem}
caption{text-align:left;color:var(--ink-muted);font-size:.82rem;padding:0 0 .45rem}
th,td{padding:.4rem .6rem;border-bottom:1px solid var(--ink-border);text-align:left;
  vertical-align:top}
th{color:var(--ink-secondary);font-weight:600;white-space:nowrap;
  border-bottom:1px solid var(--ink-axis)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
tbody tr:last-child td{border-bottom:none}
.swatch{display:inline-block;width:.65rem;height:.65rem;border-radius:2px;
  margin-right:.4rem;vertical-align:baseline}

/* Callouts ----------------------------------------------------------------- */
.callout{margin:1.2rem 0;padding:.85rem 1rem;border-radius:6px;
  border:1px solid var(--ink-border);background:var(--ink-plane)}
.callout h4{margin:0 0 .35rem;color:var(--ink-primary)}
.callout p:last-child{margin-bottom:0}
.callout.must{border-left:3px solid var(--risk-4)}
.callout.limit{border-left:3px solid var(--tier-medium)}
.cite{font-size:.72em;vertical-align:.35em;color:var(--ink-muted);
  border-bottom:1px dotted var(--ink-axis);cursor:help;margin-left:.1em}
.stamp{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78rem;
  color:var(--ink-muted)}
.toc{margin:1.5rem 0 0;padding:.9rem 1rem;background:var(--ink-plane);
  border:1px solid var(--ink-border);border-radius:6px}
.toc ol{margin:.3rem 0 0;padding-left:1.3rem}
.toc li{margin:.15rem 0;font-size:.9rem}
.pill{display:inline-block;font-size:.72rem;padding:.1rem .45rem;border-radius:99px;
  border:1px solid var(--ink-axis);color:var(--ink-secondary);margin-left:.4rem;
  vertical-align:.1em;white-space:nowrap}

/* Print -------------------------------------------------------------------- */
@page{margin:16mm}
@media print{
  body{font-size:10.5pt;background:#fff}
  .page{max-width:none;padding:0}
  .no-print{display:none !important}
  h2{page-break-after:avoid;break-after:avoid}
  h3,h4{page-break-after:avoid;break-after:avoid}
  .fig-block,.callout,.table-wrap{page-break-inside:avoid;break-inside:avoid}
  tr{page-break-inside:avoid;break-inside:avoid}
  thead{display:table-header-group}
  a{text-decoration:none;color:inherit}
  a[href^="http"]::after{content:" (" attr(href) ")";font-size:.8em;color:#666}
  .table-wrap{overflow-x:visible}
}
"""


def build_stamp(*parts: str) -> str:
    """8-hex content stamp, shown in the header.

    A cached page is indistinguishable from a fresh one by eye, and that ambiguity has
    already produced two phantom bug reports against the dashboard. Same fix here.
    """
    h = hashlib.sha256()
    for p in parts:
        h.update(p.encode("utf-8"))
    return h.hexdigest()[:8]


#: Anchored on the stamp span, NOT on a bare "build " substring. Caveat prose legitimately
#: contains the word "build", and a loose split once pulled a whole sentence out as the
#: stamp.
STAMP_RE = re.compile(r'class="stamp">build ([0-9a-f]{8})<')


def read_stamp(html_text: str) -> str:
    m = STAMP_RE.search(html_text)
    return m.group(1) if m else "?"


def document(
    *,
    title: str,
    subtitle: str,
    delivery: Delivery,
    body: str,
    stamp: str,
    toc: Sequence[Tuple[str, str]] = (),
) -> str:
    """The standalone HTML shell shared by both reports."""
    toc_html = ""
    if toc:
        items = "".join(f'<li><a href="#{esc(a)}">{esc(t)}</a></li>' for a, t in toc)
        toc_html = f'<nav class="toc no-print"><strong>Contents</strong><ol>{items}</ol></nav>'
    generated = delivery.manifest.get("generated_utc", "")
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    return f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)} — {esc(delivery.customer)}</title>
<style>{css_tokens()}{PAGE_CSS}</style>
</head><body><main class="page">
<header>
  <h1>{esc(title)}</h1>
  <p class="lede">{esc(subtitle)}</p>
  <p class="sub">{esc(delivery.customer)} &nbsp;·&nbsp; prepared {esc(now)}
     &nbsp;·&nbsp; extract generated {esc(generated)}
     &nbsp;·&nbsp; <span class="stamp">build {esc(stamp)}</span></p>
</header>
{toc_html}
{body}
</main></body></html>
"""
