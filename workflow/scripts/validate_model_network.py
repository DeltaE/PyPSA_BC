"""Run structural checks on a prepared, unsolved PyPSA network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pypsa

BUS_COLUMNS = {
    "generators": ("bus",),
    "loads": ("bus",),
    "stores": ("bus",),
    "lines": ("bus0", "bus1"),
    "links": ("bus0", "bus1"),
    "transformers": ("bus0", "bus1"),
}


def _validate_dynamic_inputs(network: pypsa.Network) -> list[str]:
    """Return hard errors for non-finite operational time-series inputs."""
    errors: list[str] = []
    checks = {
        "generators_t.p_max_pu": network.generators_t.p_max_pu,
        "generators_t.p_min_pu": network.generators_t.p_min_pu,
        "generators_t.p_set": network.generators_t.p_set,
        "loads_t.p_set": network.loads_t.p_set,
        "stores_t.e_min_pu": network.stores_t.e_min_pu,
        "stores_t.e_max_pu": network.stores_t.e_max_pu,
        "links_t.p_max_pu": network.links_t.p_max_pu,
        "links_t.p_min_pu": network.links_t.p_min_pu,
    }
    for label, frame in checks.items():
        if frame.empty:
            continue
        numeric = frame.astype(float)
        finite = np.isfinite(numeric.to_numpy(dtype=float))
        if finite.all():
            continue
        bad_columns = numeric.columns[(~finite).any(axis=0)].astype(str).tolist()
        errors.append(
            f"{label} contains non-finite values in {len(bad_columns)} component(s): "
            + ", ".join(bad_columns[:10])
        )

    if not network.generators.empty:
        p_set = network.get_switchable_as_dense("Generator", "p_set")
        p_min_pu = network.get_switchable_as_dense("Generator", "p_min_pu")
        p_max_pu = network.get_switchable_as_dense("Generator", "p_max_pu")
        p_nom = network.generators["p_nom"].astype(float)
        lower = p_min_pu.multiply(p_nom, axis=1)
        upper = p_max_pu.multiply(p_nom, axis=1)
        fixed = p_set.notna()
        outside = fixed & ((p_set < lower - 1e-9) | (p_set > upper + 1e-9))
        if outside.any().any():
            bad_columns = outside.columns[outside.any(axis=0)].astype(str).tolist()
            errors.append(
                "generators_t.p_set lies outside dispatch bounds in "
                f"{len(bad_columns)} component(s): " + ", ".join(bad_columns[:10])
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    network = pypsa.Network(args.network)
    errors: list[str] = []
    warnings: list[str] = []
    buses = set(network.buses.index.astype(str))

    if not buses:
        errors.append("Network has no buses.")
    if len(network.snapshots) == 0:
        errors.append("Network has no snapshots.")
    if network.loads.empty:
        errors.append("Network has no loads.")
    errors.extend(_validate_dynamic_inputs(network))

    for table_name, columns in BUS_COLUMNS.items():
        table = getattr(network, table_name)
        for column in columns:
            if column not in table.columns:
                continue
            referenced = set(table[column].dropna().astype(str)) - {""}
            missing = sorted(referenced - buses)
            if missing:
                errors.append(
                    f"{table_name}.{column} references {len(missing)} missing buses: "
                    + ", ".join(missing[:10])
                )

    if {"x", "y"}.issubset(network.buses.columns):
        coordinates = network.buses[["x", "y"]].apply(
            lambda column: np.isfinite(column.astype(float))
        )
        invalid_coordinates = int((~coordinates.all(axis=1)).sum())
        if invalid_coordinates:
            warnings.append(f"{invalid_coordinates} buses have invalid coordinates.")

    summary = {
        "status": "PASS" if not errors else "FAIL",
        "network": str(args.network),
        "checks": {
            "buses": len(network.buses),
            "lines": len(network.lines),
            "generators": len(network.generators),
            "loads": len(network.loads),
            "stores": len(network.stores),
            "links": len(network.links),
            "snapshots": len(network.snapshots),
        },
        "errors": errors,
        "warnings": warnings,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if errors:
        raise RuntimeError("Prepared network failed structural validation: " + "; ".join(errors))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
