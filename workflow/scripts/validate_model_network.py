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
