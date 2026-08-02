"""Optimize a previously prepared and validated PyPSA-BC network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pypsa


def _highs_iis(network: pypsa.Network, limit: int = 200) -> dict:
    """Extract a compact HiGHS IIS using Linopy's row-to-label ordering."""
    solver_model = network.model.solver_model
    if solver_model is None or network.model.solver_name != "highs":
        return {"available": False, "reason": "HiGHS solver model unavailable"}
    try:
        status, iis = solver_model.getIis()
        row_indices = [int(value) for value in iis.row_index_]
        labels = network.model.constraints.label_index.clabels
        selected_labels = [int(labels[index]) for index in row_indices if int(labels[index]) >= 0]
        positions = network.model.constraints.get_label_position(selected_labels[:limit])
        constraints = []
        for label, position in zip(selected_labels[:limit], positions, strict=True):
            name, coordinates = position
            constraints.append(
                {
                    "label": label,
                    "name": name,
                    "coordinates": {
                        str(key): str(value) for key, value in (coordinates or {}).items()
                    },
                }
            )
        return {
            "available": bool(iis.valid_ and row_indices),
            "highs_status": str(status),
            "row_count": len(row_indices),
            "column_count": len(iis.col_index_),
            "constraints_truncated": len(selected_labels) > limit,
            "constraints": constraints,
        }
    except (AttributeError, IndexError, RuntimeError, TypeError, ValueError) as exc:
        return {"available": False, "reason": f"IIS extraction failed: {exc}"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--time-limit-seconds", type=float, default=900.0)
    parser.add_argument(
        "--highs-method",
        choices=("choose", "simplex", "ipm"),
        default="choose",
    )
    args = parser.parse_args()

    from workflow.scripts.build_model import _choose_solver_name

    from pypsa_bc.studies.cascade.constraints import add_cascade_constraints

    network = pypsa.Network(args.network)
    solver = _choose_solver_name()
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    running = {
        "status": "running",
        "termination_condition": None,
        "solver": solver,
        "highs_method": args.highs_method,
        "time_limit_seconds": args.time_limit_seconds,
        "input_network": str(args.network),
        "output_network": None,
    }
    args.summary.write_text(json.dumps(running, indent=2) + "\n", encoding="utf-8")
    solver_options = {"time_limit": args.time_limit_seconds}
    if solver == "highs" and args.highs_method != "choose":
        solver_options["solver"] = args.highs_method
    status, termination = network.optimize(
        solver_name=solver,
        extra_functionality=add_cascade_constraints,
        **solver_options,
    )
    success = str(status).lower() == "ok" and str(termination).lower() == "optimal"
    result = {
        "status": str(status),
        "termination_condition": str(termination),
        "solver": solver,
        "highs_method": args.highs_method,
        "time_limit_seconds": args.time_limit_seconds,
        "input_network": str(args.network),
        "output_network": str(args.output) if success else None,
    }
    if not success and str(termination).lower() in {"infeasible", "infeasible_or_unbounded"}:
        result["infeasibility_diagnostic"] = _highs_iis(network)
    if success:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        network.export_to_netcdf(args.output)
    args.summary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not success:
        raise RuntimeError(f"Optimization failed: status={status}, termination={termination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
