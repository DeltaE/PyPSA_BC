"""
bootstrap_hydro_static.py — build the static hydro assets consumed by
`build_hydro_cascade.py`, and optionally upgrade reservoir storage from GRanD.

Outputs
-------
    hydro_topology.csv    facility -> reservoir -> cascade_group, cascade_order
    hydro_reservoirs.csv  reservoir -> live storage, discharge, sill, spill (+ provenance)

Design
------
1. `extract_facilities`  lift per-facility physics/topology from the withdrawn
   CODERS `hydro_cascade.csv`.
2. `assign_reservoirs`   collapse facilities that share a storage body into one
   reservoir (so storage is never double-counted), numbered by cascade order.
3. `build_topology` / `build_reservoirs`  emit the two tables, stamping every
   reservoir row with provenance (`method`, `verified`, `source`).
4. `apply_storage_override`  swap `live_*_m3` for values from a verified source
   (e.g. GRanD), flipping `method="GRanD"`, `verified=True` per reservoir.
   `grand_to_override` builds that override table from a GRanD export plus a
   small `reservoir -> GRAND_ID` crosswalk.

The legacy freeze is the seed; the override is how you retire it, one verified
reservoir at a time.

Author: Md Eliasinul Islam
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

SENTINEL_SPILL = 1e7  # m3/s; matches build_hydro_cascade.SENTINEL_SPILL

# legacy hydro_cascade column -> hydro_reservoirs column
_RESERVOIR_FIELDS = {
    "min_storage": "live_min_m3",
    "max_storage": "live_max_m3",
    "min_water_discharge": "q_min_m3s",
    "max_water_discharge": "q_max_m3s",
    "power_intake_sill_elevation": "sill_elev_m",
    "max_spill": "spill_max_m3s",
}

# columns that define one physical reservoir (shared storage collapses to one row)
_RESERVOIR_KEY = ["cascade_group", "min_storage", "max_storage"]

TOPOLOGY_COLUMNS = ["facility", "reservoir", "cascade_group", "cascade_order"]
RESERVOIR_COLUMNS = [
    "reservoir", "live_min_m3", "live_max_m3", "q_min_m3s", "q_max_m3s",
    "sill_elev_m", "sill_datum", "spill_max_m3s",
    "flag_spill", "flag_minstor", "flag_dup", "method", "verified", "source",
]


# --------------------------------------------------------------------------- #
#  1. Extract facilities from the legacy table
# --------------------------------------------------------------------------- #

def _facility_of(project_name: str) -> str:
    """Turbine project_name -> facility name: 'Bridge1_04' -> 'Bridge1'."""
    return re.sub(r"_\d+$", "", str(project_name)).strip()


def extract_facilities(legacy_cascade_csv) -> pd.DataFrame:
    """One row per facility with its physics, cascade group, and order."""
    df = pd.read_csv(legacy_cascade_csv)
    df = df.loc[:, [c for c in df.columns if not c.startswith("Unnamed")]]
    df["facility"] = df["project_name"].map(_facility_of)

    cols = ["facility", "cascade_group_name", "number", *_RESERVOIR_FIELDS]
    fac = (
        df[cols]
        .drop_duplicates("facility")
        .rename(columns={"cascade_group_name": "cascade_group", "number": "cascade_order"})
        .reset_index(drop=True)
    )
    return fac


# --------------------------------------------------------------------------- #
#  2. Assign reservoirs (shared storage -> one reservoir, named by order)
# --------------------------------------------------------------------------- #

def assign_reservoirs(fac: pd.DataFrame) -> pd.DataFrame:
    """Add a `reservoir` column, numbered by cascade order within each group.

    Facilities sharing identical (group, min_storage, max_storage) draw the same
    body and get the same reservoir id (e.g. Bridge1 & Bridge2 -> Carpenter).
    """
    fac = fac.copy()
    # representative order per distinct reservoir = its earliest (uppermost) facility
    rep = (
        fac.groupby(_RESERVOIR_KEY, sort=False)["cascade_order"].transform("min")
    )
    fac["_rep_order"] = rep

    names = {}
    for cg, grp in fac.groupby("cascade_group", sort=False):
        # rank distinct reservoirs in this group by their representative order
        distinct = (
            grp.drop_duplicates(_RESERVOIR_KEY)
            .sort_values("_rep_order")
            .reset_index()
        )
        for rank, row in enumerate(distinct.itertuples(), start=1):
            key = (row.cascade_group, row.min_storage, row.max_storage)
            names[key] = f"{cg}_R{rank}"

    fac["reservoir"] = [
        names[(r.cascade_group, r.min_storage, r.max_storage)] for r in fac.itertuples()
    ]
    return fac.drop(columns="_rep_order")


# --------------------------------------------------------------------------- #
#  3. Build the two output tables
# --------------------------------------------------------------------------- #

def build_topology(fac: pd.DataFrame) -> pd.DataFrame:
    return (
        fac[TOPOLOGY_COLUMNS]
        .sort_values(["cascade_group", "cascade_order", "facility"])
        .reset_index(drop=True)
    )


def build_reservoirs(fac: pd.DataFrame) -> pd.DataFrame:
    """One row per reservoir, with quality flags and legacy provenance."""
    res = (
        fac.drop_duplicates("reservoir")
        .rename(columns=_RESERVOIR_FIELDS)[["reservoir", *_RESERVOIR_FIELDS.values()]]
        .copy()
    )

    n_facilities = fac.groupby("reservoir").size()
    res["sill_datum"] = "unknown"
    res["source"] = "CODERS hydro_cascade (withdrawn 2025); provenance not per-field"
    res["method"] = "legacy_CODERS_hydro_cascade"
    res["verified"] = False
    res = _recompute_flags(res, n_facilities)

    return res[RESERVOIR_COLUMNS].sort_values("reservoir").reset_index(drop=True)


def _recompute_flags(res: pd.DataFrame, n_facilities: pd.Series) -> pd.DataFrame:
    """(Re)derive the quality flags from current values."""
    res = res.copy()
    res["flag_spill"] = res["spill_max_m3s"].isna() | (res["spill_max_m3s"] >= SENTINEL_SPILL)
    res["flag_minstor"] = res["live_min_m3"].fillna(0).eq(0)      # dead storage imputed as 0
    res["flag_dup"] = res["reservoir"].map(n_facilities).gt(1)    # >1 facility on this reservoir
    return res


# --------------------------------------------------------------------------- #
#  4. GRanD storage swap (verified-source override)
# --------------------------------------------------------------------------- #

def grand_to_override(grand_csv, crosswalk_csv, out_csv=None,
                      cap_col="CAP_MCM", id_col="GRAND_ID") -> pd.DataFrame:
    """Build a storage-override table from GRanD via a reservoir->GRAND_ID crosswalk.

    Args:
        grand_csv: GRanD/GDW export with a capacity column (million m3) and an id.
        crosswalk_csv: two columns — `reservoir`, `grand_id` (hand-maintained).
        cap_col, id_col: column names in the GRanD export.

    Returns a frame keyed by `reservoir` with `live_max_m3` and `source`.
    """
    grand = pd.read_csv(grand_csv)[[id_col, cap_col]].rename(columns={id_col: "grand_id"})
    xwalk = pd.read_csv(crosswalk_csv)  # reservoir, grand_id
    merged = xwalk.merge(grand, on="grand_id", how="left")
    merged["live_max_m3"] = merged[cap_col] * 1e6          # million m3 -> m3
    merged["source"] = "GRanD/GDW v1.0"
    out = merged[["reservoir", "live_max_m3", "source"]]

    missing = merged.loc[merged["live_max_m3"].isna(), "reservoir"].tolist()
    if missing:
        print(f"!! grand_to_override: no GRanD capacity for {missing}")
    if out_csv:
        out.to_csv(out_csv, index=False)
        print(f">> storage override written: {out_csv} ({out['live_max_m3'].notna().sum()} reservoirs)")
    return out


def apply_storage_override(res: pd.DataFrame, override, method="GRanD") -> pd.DataFrame:
    """Overlay verified storage onto the reservoir table, per reservoir.

    Any of `live_min_m3`, `live_max_m3`, `sill_elev_m`, `spill_max_m3s` present in
    `override` replaces the legacy value; matched reservoirs are stamped
    `method`, `verified=True`, and the override's `source`. Unmatched reservoirs
    keep their legacy freeze.
    """
    ov = override if isinstance(override, pd.DataFrame) else pd.read_csv(override)
    ov = ov.set_index("reservoir")
    res = res.set_index("reservoir")

    fields = [c for c in ("live_min_m3", "live_max_m3", "sill_elev_m", "spill_max_m3s") if c in ov]
    matched = res.index.intersection(ov.index)
    for col in fields:
        vals = ov.loc[matched, col]
        res.loc[matched, col] = res.loc[matched, col].where(vals.isna(), vals)

    res.loc[matched, "method"] = method
    res.loc[matched, "verified"] = True
    if "source" in ov:
        res.loc[matched, "source"] = ov.loc[matched, "source"]

    res = res.reset_index()
    # storage/spill values changed, so refresh the flags that depend on them
    # (flag_dup is topology, not storage, so it is left untouched)
    res["flag_spill"] = res["spill_max_m3s"].isna() | (res["spill_max_m3s"] >= SENTINEL_SPILL)
    res["flag_minstor"] = res["live_min_m3"].fillna(0).eq(0)

    unmatched = sorted(set(res["reservoir"]) - set(matched))
    print(f">> storage override: {len(matched)}/{len(res)} reservoirs verified from {method}; "
          f"unmatched kept legacy: {unmatched or 'none'}")
    return res[RESERVOIR_COLUMNS]


# --------------------------------------------------------------------------- #
#  Orchestration
# --------------------------------------------------------------------------- #

def bootstrap(legacy_cascade_csv, out_dir, storage_override=None):
    """Produce hydro_topology.csv and hydro_reservoirs.csv (optionally GRanD-upgraded)."""
    fac = assign_reservoirs(extract_facilities(legacy_cascade_csv))
    topo = build_topology(fac)
    res = build_reservoirs(fac)
    if storage_override is not None:
        res = apply_storage_override(res, storage_override)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    topo.to_csv(out_dir / "hydro_topology.csv", index=False)
    res.to_csv(out_dir / "hydro_reservoirs.csv", index=False)
    print(f">> hydro_topology.csv  : {len(topo)} facilities")
    print(f">> hydro_reservoirs.csv: {len(res)} reservoirs "
          f"({int(res['verified'].sum())} verified, {int(res['flag_dup'].sum())} shared)")
    return topo, res


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Seed/refresh static hydro assets from the legacy cascade table.")
    p.add_argument("--legacy", required=True, help="withdrawn hydro_cascade.csv")
    p.add_argument("--out-dir", default="data/static")
    p.add_argument("--storage-override", default=None,
                   help="verified-storage CSV (keyed by reservoir); e.g. GRanD output")
    p.add_argument("--grand", nargs=2, metavar=("GRAND_CSV", "CROSSWALK_CSV"),
                   help="build the override from a GRanD export + reservoir->GRAND_ID crosswalk")
    a = p.parse_args()

    override = a.storage_override
    if a.grand:
        override = grand_to_override(a.grand[0], a.grand[1])
    bootstrap(a.legacy, a.out_dir, storage_override=override)
