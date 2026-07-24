"""
build_hydro_cascade.py — regenerate the deprecated CODERS `hydro_cascade` table.

Context
-------
CODERS deprecated the `hydro_cascade` table that ``create_hydro_assets.main()``
consumes. The generation attributes it held (capacity, energy, years) are still
available from the live `existing_hydro`/`generators` pull; the *cascade-specific*
fields are not:

    cascade_group_name, number, code,
    min_storage, max_storage, min_water_discharge, max_water_discharge,
    power_intake_sill_elevation, max_spill, num_of_units

These are physical constants of the dams and reservoirs — they do not change
release to release. This pipeline therefore **freezes** them from the last-known
cascade file (kept as `cascade_reference.csv`) and **re-attaches fresh
generation attributes** from `existing_hydro`, matched at facility level. It
fabricates nothing: unmatched facilities keep their reference values and are
reported for manual attention.

Join key
--------
``fac_key`` = alphabetic prefix of the code's middle token
(``BC_BR0101_GEN`` -> ``BR``, ``BC_LAJ00_GEN`` -> ``LAJ``,
``BC_ABN01_GEN`` -> ``ABN``). Verified to match 76/76 turbine rows and 25/25
facilities between the legacy cascade file and `existing_hydro` (2026-07-23).

Column policy
-------------
FROZEN  (from cascade_reference): all cascade physics + topology + per-turbine
        identity and per-turbine capacities/energy.
REFRESH (from existing_hydro, facility-level broadcast): owner, latitude,
        longitude, connecting_voltage_in_kv, start_year, gen_type,
        gen_type_copper, capacity_factor. Per-turbine capacities are refreshed
        only if ``refresh_capacity=True`` (proportional rescale that keeps the
        reference's turbine split while matching the new facility total).

Upgrade path
------------
When reservoir storage is sourced properly (GRanD/GDW), replace the frozen
`min_storage`/`max_storage` columns via ``grand_reservoir_adapter.build_table_b``
and drop them from cascade_reference — no other change needed.

Author: Md Eliasinul Islam
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

# Exact target schema (order matters for downstream create_hydro_assets).
TARGET_COLUMNS = [
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

# Fields frozen from the reference (physics + topology + turbine identity).
FROZEN_COLUMNS = [
    "project_name", "gen_node_code", "province", "location",
    "copper_balancing_area", "operating_region",
    "connecting_node_name", "connecting_node_code",
    "cascade_group_name", "number", "code", "end_year",
    "install_capacity_in_mw", "capacity_adjustment_in_%",
    "effective_capacity_in_mw", "annual_avg_energy_unit_gwh/y",
    "min_storage", "max_storage", "min_water_discharge", "max_water_discharge",
    "power_intake_sill_elevation", "max_spill", "num_of_units",
    "other_references", "notes",
]

# Capacity/energy columns eligible for proportional refresh.
CAPACITY_COLUMNS = [
    "install_capacity_in_mw", "effective_capacity_in_mw",
    "annual_avg_energy_unit_gwh/y",
]


def fac_key(code: str) -> str | None:
    """Alphabetic prefix of the middle token: BC_BR0101_GEN -> 'BR'."""
    parts = str(code).split("_")
    if len(parts) < 2:
        return None
    m = re.match(r"[A-Za-z]+", parts[1])
    return m.group(0) if m else None


def extract_reference(legacy_cascade_csv, out_csv=None) -> pd.DataFrame:
    """One-time: pull the frozen physics/topology backbone from the legacy file.

    Keeps the per-turbine identity and every non-CODERS-derivable field, plus a
    ``fac_key`` for joining. Run once, then maintain ``cascade_reference.csv`` by
    hand as dams/reservoirs change (rarely).
    """
    ref = pd.read_csv(legacy_cascade_csv)
    ref = ref.loc[:, [c for c in ref.columns if not c.startswith("Unnamed")]]
    ref["fac_key"] = ref["gen_node_code"].map(fac_key)
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        ref.to_csv(out_csv, index=False)
        print(f">> cascade reference written: {out_csv} ({len(ref)} turbine rows)")
    return ref


def _facility_generation(existing_hydro_csv) -> pd.DataFrame:
    """Aggregate live existing_hydro to one row per facility (generation attrs)."""
    e = pd.read_csv(existing_hydro_csv)
    e["fac_key"] = e["generation_facility_code"].map(fac_key)
    return e.groupby("fac_key").agg(
        owner=("owner", "first"),
        latitude=("latitude", "mean"),
        longitude=("longitude", "mean"),
        connecting_voltage_in_kv=("network_node_voltage", "max"),
        start_year=("start_year", "min"),
        gen_type=("generator_fuel_type", "first"),
        gen_type_copper=("gen_type_copper", "first"),
        capacity_factor=("capacity_factor", "mean"),
        fac_install=("unit_installed_capacity", "sum"),
        fac_effective=("unit_effective_capacity", "sum"),
        fac_energy=("unit_average_annual_energy", "sum"),
    )


def build(
    reference_csv,
    existing_hydro_csv,
    out_csv=None,
    refresh_capacity: bool = False,
    log=None,
) -> pd.DataFrame:
    """Assemble a schema-conforming hydro_cascade table.

    Args:
        reference_csv: frozen cascade backbone (from ``extract_reference``).
        existing_hydro_csv: live CODERS existing_hydro pull.
        out_csv: optional output path.
        refresh_capacity: proportionally rescale per-turbine capacity/energy to
            the new facility totals (preserves the reference split).
        log: optional callable(parameter, value, **kw) e.g. ``log_assumption``.
    """
    ref = pd.read_csv(reference_csv) if not isinstance(reference_csv, pd.DataFrame) else reference_csv.copy()
    if "fac_key" not in ref:
        ref["fac_key"] = ref["gen_node_code"].map(fac_key)
    gen = _facility_generation(existing_hydro_csv)

    matched = ref["fac_key"].isin(gen.index)
    missing = sorted(ref.loc[~matched, "fac_key"].dropna().unique())
    print(f">> facility match: {matched.sum()}/{len(ref)} turbine rows "
          f"({ref['fac_key'].isin(gen.index).groupby(ref['fac_key']).first().sum()} facilities). "
          f"Unmatched: {missing or 'none'}")

    out = ref.copy()
    g = out["fac_key"].map(gen.to_dict("index"))  # dict per row, NaN if unmatched

    def refreshed(col, gen_field):
        return [row.get(gen_field) if isinstance(row, dict) else out.at[i, col] if col in out else np.nan
                for i, row in zip(out.index, g)]

    # --- refresh facility-level scalars (keep reference where unmatched) ---
    for col, gfield in [
        ("owner", "owner"), ("latitude", "latitude"), ("longitude", "longitude"),
        ("connecting_voltage_in_kv", "connecting_voltage_in_kv"),
        ("start_year", "start_year"), ("gen_type", "gen_type"),
        ("gen_type_copper", "gen_type_copper"),
    ]:
        vals = refreshed(col, gfield)
        out[col] = [v if v is not None and not (isinstance(v, float) and np.isnan(v)) else out.at[i, col]
                    for i, v in zip(out.index, vals)]

    # capacity_factor -> capacity_factor_in_%
    cf = refreshed("capacity_factor_in_%", "capacity_factor")
    out["capacity_factor_in_%"] = [v if isinstance(v, (int, float)) and not (isinstance(v, float) and np.isnan(v))
                                   else out.at[i, "capacity_factor_in_%"] for i, v in zip(out.index, cf)]

    # --- optional proportional capacity/energy refresh ---
    if refresh_capacity:
        for turb_col, fac_field in [("install_capacity_in_mw", "fac_install"),
                                    ("effective_capacity_in_mw", "fac_effective"),
                                    ("annual_avg_energy_unit_gwh/y", "fac_energy")]:
            ref_tot = out.groupby("fac_key")[turb_col].transform("sum")
            new_tot = out["fac_key"].map(gen[fac_field])
            scale = (new_tot / ref_tot).where(ref_tot > 0, 1.0).fillna(1.0)
            out[turb_col] = out[turb_col] * scale
        if log:
            log(parameter="hydro_cascade capacities", value="refreshed (proportional rescale)",
                unit="MW/GWh", rationale="Per-turbine split kept from reference; totals set to live existing_hydro",
                source="existing_hydro")

    # --- assemble exact schema ---
    out = out.drop(columns=["fac_key"], errors="ignore")
    result = out.reindex(columns=TARGET_COLUMNS)

    # --- validate ---
    assert list(result.columns) == TARGET_COLUMNS, "Column schema mismatch vs target."
    assert len(result) == len(ref), "Row count changed during build."

    if log:
        log(parameter="hydro_cascade frozen physics",
            value="storage/discharge/spill/sill/cascade-topology carried from last-known file",
            unit="mixed",
            rationale="CODERS deprecated hydro_cascade; physical reservoir constants frozen, generation attrs refreshed",
            source="cascade_reference.csv")
        if missing:
            log(parameter="hydro_cascade unmatched facilities", value=", ".join(missing),
                unit="-", rationale="No facility match in existing_hydro; reference values retained",
                source="build_hydro_cascade")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(out_csv, index=True, index_label="")  # reproduce leading index col
        print(f">> hydro_cascade written: {out_csv} ({result.shape[0]}×{result.shape[1]})")
    return result


def main(legacy_cascade_csv, existing_hydro_csv, out_dir, refresh_capacity=False, log=None):
    """End-to-end: extract reference (if absent) then build hydro_cascade."""
    out_dir = Path(out_dir)
    ref_path = out_dir / "cascade_reference.csv"
    if not ref_path.exists():
        extract_reference(legacy_cascade_csv, ref_path)
    return build(ref_path, existing_hydro_csv, out_dir / "hydro_cascade.csv",
                 refresh_capacity=refresh_capacity, log=log)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Rebuild the deprecated CODERS hydro_cascade table.")
    p.add_argument("--legacy", required=True, help="last-known hydro_cascade.csv")
    p.add_argument("--existing", required=True, help="live existing_hydro.csv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--refresh-capacity", action="store_true")
    a = p.parse_args()
    main(a.legacy, a.existing, a.out_dir, refresh_capacity=a.refresh_capacity)
