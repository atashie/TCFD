"""Facet profiles — the library a bespoke report is composed from.

THE COMBINATORIAL PROBLEM, AND WHY THIS IS THE SHAPE OF THE ANSWER
------------------------------------------------------------------
A report should be specific to one asset × region × persona × business vertical × company ×
use case. Written out, that is a document per combination, and the combinations multiply
faster than anyone can write them: five asset types, ten regions, four personas, five
verticals and four use cases is four thousand documents.

So the library is not a document per combination. It is a short profile per VALUE of each
facet — one for "loblolly pine timberland", one for "US Southeast", one for "sustainability
team" — and a report selects one profile from each facet. Six facets with ten values each is
sixty files, and it covers a million combinations.

PROFILES GUIDE THE NARRATIVE. THEY ARE NOT PASTED INTO IT.
----------------------------------------------------------
This is the load-bearing rule and it is easy to get wrong. If profile prose were rendered
into the report, every loblolly report would contain the same paragraphs and the output would
be exactly the generic document the facet library exists to avoid — worse, it would be
generic while *looking* tailored.

So a profile contains no sentences meant for a customer. It contains what the writer needs in
order to write: which hazards actually reach this asset and how, what decision this reader is
making, what vocabulary lands and what misfires, which confounders will make them distrust
the numbers, what to ask them, and what must never be claimed. The generator surfaces all of
that into the narrative scaffold as inline comments, and the writer writes.

`company` is the exception in kind: it is researched per engagement rather than drawn from a
library, so a first engagement writes only `dossier.yaml`. A repeat customer earns a company
profile seeded from that dossier.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from .delivery import PROJECT_ROOT, DeliveryError

PROFILES_ROOT = PROJECT_ROOT / "docs" / "reporting" / "profiles"

#: Directory name per facet. `use_case` is stored hyphenated on disk, like the others.
FACET_DIRS = {
    "asset": "asset",
    "region": "region",
    "persona": "persona",
    "vertical": "vertical",
    "use_case": "use-case",
    "company": "company",
}

#: The fixed section vocabulary. Every profile uses a subset of these headings and nothing
#: else -- a fixed schema is what makes composition mechanical rather than editorial, and a
#: typo'd heading is caught at load rather than silently dropping guidance.
SECTIONS = (
    "Transmission channels",
    "What this reader needs decided",
    "Framing and vocabulary",
    "Metrics that lead",
    "Known confounders",
    "Questions for the customer",
    "Do not claim",
)

REQUIRED_FRONT_MATTER = ("facet", "id", "name", "confirmed_on")

FRONT_MATTER_RE = re.compile(r"\A---\s*\n(.*?)\n---\s*\n(.*)\Z", re.S)
HEADING_RE = re.compile(r"^##\s+(.+?)\s*$", re.M)


@dataclass
class Profile:
    facet: str
    id: str
    name: str
    path: Path
    meta: dict
    sections: Dict[str, str]

    @property
    def confirmed_on(self):
        return self.meta.get("confirmed_on")

    @property
    def sources(self) -> List[dict]:
        return list(self.meta.get("sources") or [])

    def section(self, name: str) -> str:
        return self.sections.get(name, "").strip()


def parse_profile(path: Path) -> Profile:
    text = path.read_text()
    m = FRONT_MATTER_RE.match(text)
    if not m:
        raise DeliveryError(
            f"{path} has no YAML front matter. A profile must open with a --- block "
            f"carrying {', '.join(REQUIRED_FRONT_MATTER)}."
        )
    meta = yaml.safe_load(m.group(1)) or {}
    missing = [k for k in REQUIRED_FRONT_MATTER if k not in meta]
    if missing:
        raise DeliveryError(f"{path} front matter is missing: {', '.join(missing)}")

    body = m.group(2)
    heads = list(HEADING_RE.finditer(body))
    sections: Dict[str, str] = {}
    for i, h in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(body)
        sections[h.group(1)] = body[h.end():end].strip()

    # A region profile without a `matches:` block matches nothing explicitly and therefore
    # EVERYTHING implicitly, which silently defeats `assert_region_coverage` -- the check
    # would pass while the report described the wrong part of the world. Refuse at load.
    if meta["facet"] == "region" and not (meta.get("matches") or {}).get("countries"):
        raise DeliveryError(
            f"{path} is a region profile with no `matches.countries`. Without it the "
            f"profile matches every location and the region-coverage check cannot fail.\n"
            f"  Add, in the front matter:\n"
            f"    matches:\n"
            f"      countries: [USA]\n"
            f"      states: [California]   # empty list = the whole country"
        )

    unknown = [h for h in sections if h not in SECTIONS]
    if unknown:
        raise DeliveryError(
            f"{path} uses heading(s) outside the profile schema: {', '.join(unknown)}.\n"
            f"  Allowed: {', '.join(SECTIONS)}\n"
            f"  A heading the loader does not recognise is guidance nobody will ever see."
        )
    if not sections:
        raise DeliveryError(f"{path} contains no sections.")
    return Profile(
        facet=str(meta["facet"]),
        id=str(meta["id"]),
        name=str(meta["name"]),
        path=path,
        meta=meta,
        sections=sections,
    )


def profile_dir(facet: str) -> Path:
    if facet not in FACET_DIRS:
        raise DeliveryError(f"Unknown facet {facet!r}; known: {', '.join(FACET_DIRS)}")
    return PROFILES_ROOT / FACET_DIRS[facet]


def available(facet: str) -> List[str]:
    d = profile_dir(facet)
    return sorted(p.stem for p in d.glob("*.md")) if d.exists() else []


def load_profile(facet: str, profile_id: str) -> Profile:
    path = profile_dir(facet) / f"{profile_id}.md"
    if not path.exists():
        have = available(facet)
        raise DeliveryError(
            f"No {facet} profile {profile_id!r} at {path}.\n"
            + (f"  Available: {', '.join(have)}" if have else
               f"  No {facet} profiles exist yet. Write one -- see "
               f"docs/reporting/README.md for the schema.")
        )
    prof = parse_profile(path)
    if prof.facet != facet:
        raise DeliveryError(
            f"{path} declares facet {prof.facet!r} but sits in the {facet} directory."
        )
    return prof


def load_selected(config: dict) -> Dict[str, List[Profile]]:
    """Load the profiles named in `report_config.yaml` facets.

    EVERY facet returns a LIST, and `region` in particular is expected to hold more than
    one. A portfolio is not obliged to sit in a single place: the worked example spans
    California, the Gulf coast, the mid-Atlantic and Bavaria, and describing it with one
    regional context would attach a narrative about hurricanes to a Bavarian spruce stand.
    Region is the facet the DATA determines -- see `region_coverage()`, which refuses a
    selection that does not cover every location.
    """
    facets = config.get("facets") or {}
    out: Dict[str, List[Profile]] = {}
    missing = []
    for facet in FACET_DIRS:
        chosen = facets.get(facet)
        if not chosen:
            out[facet] = []
            if facet != "company":
                missing.append(facet)
            continue
        ids = chosen if isinstance(chosen, (list, tuple)) else [chosen]
        out[facet] = [load_profile(facet, str(i)) for i in ids if i]
    if missing:
        listing = "\n".join(
            f"    {f:9} available: {', '.join(available(f)) or '(none written yet)'}"
            for f in missing
        )
        raise DeliveryError(
            "A bespoke report has to know who it is for. These facets are unset in "
            "report_config.yaml:\n"
            f"    {', '.join(missing)}\n\n"
            f"{listing}\n\n"
            "  Set them under `facets:` and rebuild. `region` accepts a list. `company` is "
            "optional -- a first engagement records the customer in dossier.yaml instead."
        )
    return out


def all_profiles(profiles: Dict[str, List[Profile]]) -> List[Profile]:
    return [p for group in profiles.values() for p in group]


def unconfirmed(profiles: Dict[str, List[Profile]]) -> List[str]:
    """Profiles no human has signed off. Same protocol as the asset catalog."""
    return sorted(
        f"{p.facet}/{p.id}" for p in all_profiles(profiles) if not p.confirmed_on
    )


def _blank_safe(value) -> str:
    """Lowercased text, with pandas' float-nan and its "nan" string both read as blank."""
    if value is None:
        return ""
    try:
        if value != value:  # NaN
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "n/a", "-"} else text.lower()


def matches_location(profile: Profile, country: str, state: str) -> bool:
    """Does this region profile describe this location?

    Country must match. If the profile names states, the location's state must be among
    them; a location with no state recorded matches on country alone, because refusing it
    would make an offshore or country-level site impossible to cover at all.
    """
    spec = profile.meta.get("matches") or {}
    countries = {str(c).strip().lower() for c in (spec.get("countries") or [])}
    states = {str(s).strip().lower() for s in (spec.get("states") or [])}
    # An absent optional column reads back from CSV as the FLOAT nan, whose str() is the
    # truthy string "nan". Left unhandled it is treated as a state name, matches nothing,
    # and an offshore site is reported as uncovered for every region on earth.
    cty = _blank_safe(country)
    st = _blank_safe(state)
    if countries and cty not in countries:
        return False
    if not states:
        return True
    return not st or st in states


def region_coverage(locations, region_profiles: List[Profile]):
    """(location_id -> matching profile ids, [uncovered rows]).

    The point of this function is that region is the ONE facet nobody should be free to
    choose: it is a property of the delivered coordinates. A report whose regional context
    does not match the sites it describes is confidently wrong about the thing a reader is
    most likely to check first.
    """
    covered, uncovered = {}, []
    for _, loc in locations.iterrows():
        hits = [
            p.id for p in region_profiles
            if matches_location(p, loc.get("country", ""), loc.get("state", ""))
        ]
        if hits:
            covered[loc["location_id"]] = hits
        else:
            uncovered.append(loc)
    return covered, uncovered


def assert_region_coverage(locations, region_profiles: List[Profile]) -> None:
    covered, uncovered = region_coverage(locations, region_profiles)
    if uncovered:
        rows = "\n".join(
            f"    {r['name']} ({r.get('country') or '?'}"
            + (f" / {r['state']}" if str(r.get('state') or '').strip() else "")
            + ")"
            for r in uncovered
        )
        have = ", ".join(p.id for p in region_profiles) or "(none selected)"
        raise DeliveryError(
            f"{len(uncovered)} location(s) are not covered by any selected region profile:\n"
            f"{rows}\n\n"
            f"  Selected: {have}\n"
            f"  Available: {', '.join(available('region')) or '(none written yet)'}\n\n"
            "  A bespoke report describes the regional context of the sites it reports on. "
            "Add a region profile whose `matches:` block covers these locations, or select "
            "an existing one, and list them all under `facets.region`."
        )


def guidance_block(profiles: Dict[str, List[Profile]], section: str) -> str:
    """Every profile's guidance for one section, formatted for an HTML comment.

    This is what the writer reads while writing. It never reaches the rendered report.
    """
    parts = []
    for facet, group in profiles.items():
        for prof in group:
            text = prof.section(section)
            if text:
                parts.append(f"  [{facet}: {prof.name}]\n" + "\n".join(
                    "    " + line for line in text.splitlines()
                ))
    return "\n\n".join(parts)
