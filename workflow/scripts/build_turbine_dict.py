"""
build_turbine_dict — resolve turbine_dict.json against the OEDB via atlite.

Your wind time-series uses `atlite.resource.get_oedb_windturbineconfig(name, id)`
(Open Energy Database on the Open Energy Platform). This script takes each wind
model used in BC, tries a list of candidate OEDB `turbine_type` names, keeps the
first that atlite can actually resolve, and writes:

    turbine_dict.json :  { "<model>": "<oedb_turbine_type>[*<id>]" }

Run it in your environment (atlite reaches OEDB there):
    python -m workflow.scripts.build_turbine_dict

NOTE: OEDI (data.openei.org) is NOT the source — it has no power-curve library.
OEDB/OEP is. If a model resolves to several OEDB ids, the `*id` is appended so
`pypsa_bc.wind.get_config()` selects the right one.
"""

from __future__ import annotations

import json
from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.logger import get_logger

log = get_logger("build_turbine_dict")

# Code-level inventory baseline — guarantees the dict is never smaller than this,
# even if the on-disk turbine_dict.json is missing/empty. Edit here to curate.
DEFAULT_INVENTORY: dict[str, str] = {
    "E82-3.0":     "E-82/3000",
    "E126-4.2":    "E-126/4200",
    "V90-3.0":     "V90/3000",
    "V90-1.8":     "V90/2000",
    "V100-1.8":    "V100/1800",
    "GE 2.75-120": "GE2.75-120",
    "GE 3.2-103":  "GE3.2-103",
    "3.2-MM114":   "MM114/3200",
    "LTW77-1.5":   "LTW77/1500",
}

# Model (as used in the BC inventory) -> ordered candidate OEDB turbine_type names.
# First candidate that atlite resolves wins. Extend as the fleet grows.
CANDIDATES: dict[str, list[str]] = {
    "E82-3.0":     ["E-82/3000", "E82/3000"],
    "E126-4.2":    ["E-126/4200", "E126/4200"],
    "V90-3.0":     ["V90/3000"],
    "V90-1.8":     ["V90/2000", "V90/1800"],
    "V100-1.8":    ["V100/1800", "V100/2000"],
    "GE 2.75-120": ["GE2.75-120", "GE2.75/120"],
    "GE 3.2-103":  ["GE3.2-103", "GE3.2/103", "GE2.75-120"],   # last = nearest-size fallback
    "3.2-MM114":   ["MM114/3200", "3.2M114/3200", "MM114"],
    "LTW77-1.5":   ["LTW77/1500", "LTW77-1500"],
}


def _resolves(name: str) -> bool:
    """True if atlite can build a turbine config for this OEDB name.

    Returns False (not raise) if atlite/OEDB is unavailable, so the caller
    falls back to the inventory rather than crashing.
    """
    try:
        import atlite
        atlite.resource.get_oedb_windturbineconfig(name)
        return True
    except Exception:
        return False


def build(out_json: str, candidates: dict[str, list[str]] = CANDIDATES):
    """Resolve models via OEDB, falling back to the existing inventory JSON.

    The current `out_json` (the hand-maintained inventory) is the baseline;
    any model atlite can resolve against OEDB overrides its inventory value.
    Models neither resolvable nor in the inventory are reported.
    """
    out_json = Path(out_json)
    file_inv: dict[str, str] = json.loads(out_json.read_text()) if out_json.exists() else {}
    # Baseline = code inventory overlaid with whatever is on disk. Never shrinks below this.
    inventory: dict[str, str] = {**DEFAULT_INVENTORY, **file_inv}

    mapping: dict[str, str] = dict(inventory)   # inventory = fallback baseline
    resolved, fell_back, unmatched = [], [], []
    for model, cands in candidates.items():
        hit = next((c for c in cands if _resolves(c)), None)
        if hit:
            mapping[model] = hit
            resolved.append(model)
        elif model in inventory:
            fell_back.append(model)             # keep inventory value
        else:
            unmatched.append(model)

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(mapping, indent=2))
    log.info(f"turbine_dict: {len(resolved)} resolved via OEDB, "
             f"{len(fell_back)} from inventory fallback, {len(unmatched)} unresolved "
             f"-> {out_json}")
    if fell_back:
        log.warning(f"inventory fallback used (OEDB not resolved): {fell_back}")
    if unmatched:
        log.warning(f"no value at all for: {unmatched} — add to inventory or CANDIDATES")
    return mapping, fell_back, unmatched


def main():
    out = AttributesParser().data_cfg["data"]["wind"]["turbine_dict"]
    build(out)


if __name__ == "__main__":
    main()
