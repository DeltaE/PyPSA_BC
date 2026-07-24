"""
to_legacy_cascade.py — emit the legacy 34-column `hydro_cascade.csv` that
`create_hydro_assets.main()` consumes, from the curated static inventory plus a
live per-unit `existing_hydro` pull.

Why this exists
---------------
`create_hydro_assets` was written against the withdrawn CODERS `hydro_cascade`
schema (one row per turbine, 33 data columns). The rebuild pipeline stores the
physics as a cleaner inventory (`hydro_topology.csv` + `hydro_reservoirs.csv`).
This adapter closes the loop: it broadcasts the reservoir physics/topology onto
each live turbine row and writes the exact legacy layout — so `create_hydro_assets`
runs unchanged while the storage numbers come from the curated (GRanD-upgradable)
inventory rather than the frozen file.

Join
----
`existing_hydro.generation_facility_name` == `hydro_topology.facility`
(verified 26/26 facilities, 76 turbine rows — the legacy count).

Author: Md Eliasinul Islam
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# exact legacy column order (33 data columns; written with a leading index like the original)
LEGACY_COLUMNS = [
    "project_name", "gen_node_code", "owner", "province", "location",
    "latitude", "longitude", "copper_balancing_area", "operating_region",
    "connecting_node_name", "connecting_node_code", "cascade_group_name",
    "number", "code", "connecting_voltage_in_kv", "start_year", "end_year",
    "gen_type", "gen_type_copper", "install_capacity_in_mw",
    "capacity_adjustment_in_%", "effective_capacity_in_mw",
    "capacity_factor_in_%", "annual_avg_energy_unit_gwh/y",
    "min_storage", "max_storage", "min_water_discharge", "max_water_discharge",
    "power_intake_sill_elevation", "max_spill", "num_of_units",
    "other_references", "notes",
]

# legacy column  <-  existing_hydro (per-unit) column (non-capacity fields)
_FROM_EXISTING = {
    "project_name": "generation_unit_name",
    "gen_node_code": "generation_unit_code",
    "owner": "owner",
    "province": "province",
    "location": "location",
    "latitude": "latitude",
    "longitude": "longitude",
    "copper_balancing_area": "copper_balancing_area",
    "operating_region": "operating_region",
    "connecting_node_name": "network_node_name",
    "connecting_node_code": "connecting_node_code",
    "connecting_voltage_in_kv": "network_node_voltage",
    "start_year": "start_year",
    "end_year": "possible_renewal_year",       # closest live proxy for legacy end_year
    "gen_type_copper": "gen_type_copper",
    "capacity_factor_in_%": "capacity_factor",
    "annual_avg_energy_unit_gwh/y": "unit_average_annual_energy",
    "num_of_units": "total_facility_generation_units",
}

# legacy column  <-  hydro_reservoirs column (broadcast from reservoir to its turbines)
_FROM_RESERVOIR = {
    "min_storage": "live_min_m3",
    "max_storage": "live_max_m3",
    "min_water_discharge": "q_min_m3s",
    "max_water_discharge": "q_max_m3s",
    "power_intake_sill_elevation": "sill_elev_m",
    "max_spill": "spill_max_m3s",
}


def to_legacy_cascade(existing_hydro_csv, static_dir, out_csv=None,
                      capacity_convention: str = "legacy") -> pd.DataFrame:
    """Build the legacy 34-column hydro_cascade table.

    Args:
        existing_hydro_csv: live CODERS existing_hydro (per generation unit).
        static_dir: folder holding hydro_topology.csv and hydro_reservoirs.csv.
        out_csv: optional output path (written with a leading index, as legacy).
        capacity_convention:
            "legacy" (default) — reproduce the old file's convention:
                install = effective = nameplate, capacity_adjustment = 1.0
                (no derate applied). Drop-in; model capacity unchanged.
            "live" — use the real CODERS derate:
                install = nameplate, adjustment = live, effective = derated.
                More correct, but lowers cascade capacity (~30% on derated units).
    """
    assert capacity_convention in ("legacy", "live"), capacity_convention
    static_dir = Path(static_dir)
    e = pd.read_csv(existing_hydro_csv)
    topo = pd.read_csv(static_dir / "hydro_topology.csv")
    res = pd.read_csv(static_dir / "hydro_reservoirs.csv")

    # keep only turbines belonging to cascade facilities
    units = e[e["generation_facility_name"].isin(topo["facility"])].copy()
    dropped_fac = sorted(set(topo["facility"]) - set(e["generation_facility_name"]))
    if dropped_fac:
        print(f"!! {len(dropped_fac)} cascade facilities absent from existing_hydro: {dropped_fac}")

    # attach cascade group/order and the reservoir key, then reservoir physics
    units = units.merge(
        topo, left_on="generation_facility_name", right_on="facility", how="left"
    ).merge(res, on="reservoir", how="left")

    out = pd.DataFrame(index=units.index)
    for legacy_col, src in _FROM_EXISTING.items():
        out[legacy_col] = units[src].values
    for legacy_col, src in _FROM_RESERVOIR.items():
        out[legacy_col] = units[src].values

    # capacity columns — convention-dependent (see docstring)
    if capacity_convention == "legacy":
        out["install_capacity_in_mw"] = units["unit_installed_capacity"].values
        out["capacity_adjustment_in_%"] = 1.0
        out["effective_capacity_in_mw"] = units["unit_installed_capacity"].values
    else:  # "live"
        out["install_capacity_in_mw"] = units["unit_installed_capacity"].values
        out["capacity_adjustment_in_%"] = units["capacity_adjustment"].values
        out["effective_capacity_in_mw"] = units["unit_effective_capacity"].values

    # derived / constant columns
    out["cascade_group_name"] = units["cascade_group"].values
    out["number"] = units["cascade_order"].values
    out["code"] = 0                                   # reservoir-type flag; curate 3 impute plants
    out["gen_type"] = units["generator_fuel_type"].str.capitalize().values  # 'hydro_daily'->'Hydro_daily'
    out["other_references"] = ""
    out["notes"] = units["notes"].fillna("").values

    result = out.reindex(columns=LEGACY_COLUMNS)

    assert list(result.columns) == LEGACY_COLUMNS, "Legacy column schema mismatch."
    unmatched_res = units["reservoir"].isna().sum()
    if unmatched_res:
        print(f"!! {unmatched_res} turbine rows did not resolve a reservoir (check topology).")
    print(f">> legacy hydro_cascade built: {len(result)} turbine rows, {len(result.columns)} columns")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(out_csv, index=True, index_label="")  # leading index like the original
        print(f">> written: {out_csv}")
    return result


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Emit the legacy hydro_cascade.csv for create_hydro_assets.")
    p.add_argument("--existing", required=True, help="live existing_hydro.csv (per unit)")
    p.add_argument("--static-dir", default="data/static")
    p.add_argument("--out", required=True, help="output hydro_cascade.csv")
    a = p.parse_args()
    to_legacy_cascade(a.existing, a.static_dir, a.out)
