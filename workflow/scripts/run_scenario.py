"""Run one implemented PyPSA-BC workflow scenario from an immutable case file."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", type=Path, required=True)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--year", type=int, default=2021)
    parser.add_argument("--capacity-choice", default="investment")
    parser.add_argument("--include-vre-investments", action="store_true")
    parser.add_argument("--external-boundary-policy", default="observed_interchange")
    parser.add_argument(
        "--generation-calibration-policy", default="statcan_monthly_nonhydro"
    )
    parser.add_argument("--ror-dispatch-policy", default="fixed_source_profile")
    parser.add_argument(
        "--hydro-dispatch-cost-policy", default="source_generic_vom"
    )
    parser.add_argument(
        "--uniform-hydro-cost-cad-per-mwh", type=float, default=1.97
    )
    parser.add_argument(
        "--observed-generation",
        type=Path,
        default=Path("data/validation/statcan/bc_generation_by_type_2021.csv"),
    )
    parser.add_argument(
        "--named-station-capacity-policy", default="bc_hydro_fiscal_2021"
    )
    parser.add_argument("--named-station-energy-policy", default="none")
    parser.add_argument("--named-station-energy-tolerance", type=float, default=0.25)
    parser.add_argument(
        "--named-station-evidence",
        type=Path,
        default=Path(
            "data/validation/bc_hydro/annual_reports/bc_hydro_fiscal_2021_supply.csv"
        ),
    )
    parser.add_argument(
        "--named-station-mapping",
        type=Path,
        default=Path("studies/chapter2/inputs/bc_hydro_station_validation_map.csv"),
    )
    args = parser.parse_args()

    case = yaml.safe_load(args.scenario.read_text(encoding="utf-8"))
    if not case.get("runnable", False):
        raise RuntimeError(
            f"Scenario {case.get('case_id')} is declared but not implemented. "
            "The workflow will not substitute another reservoir representation."
        )
    representation = case["reservoir_representation"]
    policy = case["policy"]
    storage = case["storage"]

    # Import after validation because build_model imports the full scientific stack.
    from workflow.scripts import build_model

    build_model.main(
        year=args.year,
        capacity_choice=args.capacity_choice,
        include_vre_investments=args.include_vre_investments,
        reservoir_representation=representation,
        water_policy=policy["constraint_mode"],
        storage_policy=storage["constraint_mode"],
        external_boundary_policy=args.external_boundary_policy,
        generation_calibration_policy=args.generation_calibration_policy,
        ror_dispatch_policy=args.ror_dispatch_policy,
        hydro_dispatch_cost_policy=args.hydro_dispatch_cost_policy,
        uniform_hydro_cost_cad_per_mwh=args.uniform_hydro_cost_cad_per_mwh,
        observed_generation_path=args.observed_generation,
        named_station_capacity_policy=args.named_station_capacity_policy,
        named_station_energy_policy=args.named_station_energy_policy,
        named_station_energy_tolerance=args.named_station_energy_tolerance,
        named_station_evidence_path=args.named_station_evidence,
        named_station_mapping_path=args.named_station_mapping,
        release_multiplier=float(policy.get("release_multiplier", 1.0)),
        solved_network_save_to=args.network,
        build_report_save_to=args.report,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
