"""Optimize a previously prepared and validated PyPSA-BC network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pypsa


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    from workflow.scripts.build_model import _choose_solver_name
    from pypsa_bc.studies.cascade.constraints import add_cascade_constraints

    network = pypsa.Network(args.network)
    solver = _choose_solver_name()
    status, termination = network.optimize(
        solver_name=solver,
        extra_functionality=add_cascade_constraints,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    network.export_to_netcdf(args.output)
    result = {
        "status": str(status),
        "termination_condition": str(termination),
        "solver": solver,
        "input_network": str(args.network),
        "output_network": str(args.output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if str(status).lower() not in {"ok", "warning"}:
        raise RuntimeError(f"Optimization failed: status={status}, termination={termination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
