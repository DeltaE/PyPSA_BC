"""Solve PyPSA-BC chronologically with explicit reservoir-state hand-offs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from workflow.scripts.build_model import _choose_solver_name

from pypsa_bc.studies.cascade.constraints import add_cascade_constraints
from pypsa_bc.studies.rolling_horizon import (
    build_rolling_windows,
    configure_window_network,
    initial_store_state,
    store_state_at,
)


def _write_summary(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _maximum_nodal_residual(network: pypsa.Network) -> float:
    if network.buses_t.p.empty:
        return float("nan")
    residual = network.buses_t.p.copy()
    # PyPSA's optimized Bus.p is the net injection into passive branches.
    # Subtract each passive branch's terminal flow to obtain the KCL residual.
    for component in ("lines", "transformers"):
        static = getattr(network, component)
        dynamic = getattr(network, f"{component}_t")
        for port in (0, 1):
            frame = getattr(dynamic, f"p{port}")
            if frame.empty:
                continue
            grouped = frame.T.groupby(static[f"bus{port}"]).sum().T
            common = grouped.columns.intersection(residual.columns)
            residual.loc[:, common] -= grouped.loc[:, common]
    values = pd.to_numeric(residual.stack(), errors="coerce").to_numpy(dtype=float)
    return float(np.nanmax(np.abs(values))) if values.size else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--window-hours", type=int, default=168)
    parser.add_argument("--lookahead-hours", type=int, default=24)
    parser.add_argument("--initial-storage-fraction", type=float, default=0.5)
    parser.add_argument("--time-limit-seconds", type=float, default=180.0)
    parser.add_argument("--max-windows", type=int)
    parser.add_argument("--enforce-annual-closure", action="store_true")
    parser.add_argument(
        "--closure-window-hours",
        type=int,
        default=672,
        help="Minimum joint horizon for the final annual-closure solve.",
    )
    args = parser.parse_args()

    source = pypsa.Network(args.network)
    windows = build_rolling_windows(
        source.snapshots,
        window_hours=args.window_hours,
        lookahead_hours=args.lookahead_hours,
        final_window_hours=(args.closure_window_hours if args.enforce_annual_closure else None),
    )
    initial_state, originally_cyclic = initial_store_state(
        source,
        cyclic_fraction=args.initial_storage_fraction,
    )
    selected = windows[: args.max_windows] if args.max_windows else windows
    is_complete_annual_run = len(selected) == len(windows)
    solver = _choose_solver_name()
    payload = {
        "status": "running",
        "input_network": str(args.network),
        "solver": solver,
        "window_hours": args.window_hours,
        "lookahead_hours": args.lookahead_hours,
        "initial_storage_fraction": args.initial_storage_fraction,
        "initial_state_evidence_class": "scenario",
        "enforce_annual_closure": args.enforce_annual_closure,
        "closure_window_hours": args.closure_window_hours,
        "windows_planned": len(windows),
        "windows_selected": len(selected),
        "windows_optimal": 0,
        "claim_boundary": "No annual result is valid until every window and the closure audit pass.",
    }
    _write_summary(args.summary, payload)

    args.output_root.mkdir(parents=True, exist_ok=True)
    carried = initial_state.copy()
    rows: list[dict[str, object]] = []
    water_scale = float(source.meta.get("water_volume_unit_m3", 1.0))
    physical_water_stores = originally_cyclic.intersection(source.stores.index[
        source.stores["bus"].map(source.buses["carrier"]).eq("Water")
    ])

    for position, window in enumerate(selected):
        network = pypsa.Network(args.network)
        closure_target = None
        if (
            args.enforce_annual_closure
            and is_complete_annual_run
            and position == len(selected) - 1
        ):
            closure_target = initial_state.loc[originally_cyclic]
        configure_window_network(
            network,
            window,
            carried,
            annual_closure_target=closure_target,
        )
        status, termination = network.optimize(
            solver_name=solver,
            extra_functionality=add_cascade_constraints,
            time_limit=args.time_limit_seconds,
        )
        optimal = str(status).lower() == "ok" and str(termination).lower() == "optimal"
        row = {
            "window": window.number,
            "solve_start": window.solve_start,
            "solve_end": window.solve_end,
            "commit_end": window.commit_end,
            "solve_hours": len(window.solve_snapshots),
            "commit_hours": len(window.commit_snapshots),
            "status": str(status),
            "termination_condition": str(termination),
            "maximum_nodal_residual_model_units": (
                _maximum_nodal_residual(network) if optimal else np.nan
            ),
        }
        rows.append(row)
        pd.DataFrame(rows).to_csv(args.output_root / "window_audit.csv", index=False)
        if not optimal:
            payload.update(
                {
                    "status": "failed",
                    "failed_window": window.number,
                    "termination_condition": str(termination),
                    "windows_optimal": position,
                }
            )
            _write_summary(args.summary, payload)
            raise RuntimeError(
                f"Rolling window {window.number} failed: status={status}, termination={termination}"
            )

        end_state = store_state_at(network, window.commit_end)
        handoff_gap = float((carried.reindex(network.stores.index) - network.stores.e_initial).abs().max())
        rows[-1]["maximum_state_handoff_gap_model_units"] = handoff_gap
        reservoir_end_state = end_state.reindex(physical_water_stores)
        rows[-1]["committed_end_reservoir_storage_model_units"] = float(
            reservoir_end_state.sum()
        )
        rows[-1]["committed_end_reservoir_storage_m3"] = float(
            reservoir_end_state.sum() * water_scale
        )
        carried = end_state
        window_path = args.output_root / f"window_{window.number:03d}.nc"
        network.export_to_netcdf(window_path)
        payload["windows_optimal"] = position + 1
        _write_summary(args.summary, payload)
        pd.DataFrame(rows).to_csv(args.output_root / "window_audit.csv", index=False)

    closure = carried.loc[originally_cyclic] - initial_state.loc[originally_cyclic]
    max_closure_model = float(closure.abs().max()) if len(closure) else 0.0
    closure_pass = (not args.enforce_annual_closure) or max_closure_model <= 1e-6
    payload.update(
        {
            "status": "annual_complete"
            if is_complete_annual_run and closure_pass
            else "pilot_complete",
            "maximum_annual_closure_gap_model_units": (
                max_closure_model if is_complete_annual_run else None
            ),
            "maximum_annual_closure_gap_m3": (
                max_closure_model * water_scale if is_complete_annual_run else None
            ),
            "maximum_cyclic_state_change_during_pilot_model_units": (
                max_closure_model if not is_complete_annual_run else None
            ),
            "annual_closure_pass": closure_pass if is_complete_annual_run else None,
            "annual_result_valid": bool(is_complete_annual_run and closure_pass),
        }
    )
    _write_summary(args.summary, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
