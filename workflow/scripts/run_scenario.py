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
        release_multiplier=float(policy.get("release_multiplier", 1.0)),
        solved_network_save_to=args.network,
        build_report_save_to=args.report,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
