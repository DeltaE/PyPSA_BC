"""Build and export an unsolved PyPSA-BC network for validation."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--year", type=int, default=2021)
    parser.add_argument("--capacity-choice", default="investment")
    parser.add_argument("--reservoir-representation", default="A_full_cascade")
    parser.add_argument("--water-policy", default="evidence_based")
    parser.add_argument("--storage-policy", default="zero_pondage")
    parser.add_argument("--external-boundary-policy", default="observed_interchange")
    parser.add_argument(
        "--network-flow-formulation",
        choices=("transport", "dc_load_flow"),
        default="transport",
    )
    parser.add_argument("--release-multiplier", type=float, default=1.0)
    parser.add_argument(
        "--generation-calibration-policy",
        choices=("none", "statcan_monthly_nonhydro"),
        default="none",
    )
    parser.add_argument(
        "--ror-dispatch-policy",
        choices=("curtailable", "fixed_source_profile"),
        default="curtailable",
    )
    parser.add_argument(
        "--hydro-dispatch-cost-policy",
        choices=("source_generic_vom", "uniform_reservoir_vom"),
        default="source_generic_vom",
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
        "--named-station-capacity-policy",
        choices=("none", "bc_hydro_fiscal_2021"),
        default="none",
    )
    parser.add_argument(
        "--named-station-energy-policy",
        choices=("none", "bc_hydro_fiscal_2021_screening_band"),
        default="none",
    )
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
    parser.add_argument("--include-vre-investments", action="store_true")
    args = parser.parse_args()

    from workflow.scripts import build_model

    build_model.main(
        year=args.year,
        capacity_choice=args.capacity_choice,
        include_vre_investments=args.include_vre_investments,
        reservoir_representation=args.reservoir_representation,
        water_policy=args.water_policy,
        storage_policy=args.storage_policy,
        external_boundary_policy=args.external_boundary_policy,
        network_flow_formulation=args.network_flow_formulation,
        release_multiplier=args.release_multiplier,
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
        solved_network_save_to=args.network,
        build_report_save_to=args.report,
        solve_network=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
