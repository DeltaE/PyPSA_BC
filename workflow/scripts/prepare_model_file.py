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
    parser.add_argument("--release-multiplier", type=float, default=1.0)
    parser.add_argument("--include-vre-investments", action="store_true")
    args = parser.parse_args()

    from workflow.scripts import build_model

    build_model.main(
        year=args.year,
        capacity_choice=args.capacity_choice,
        include_vre_investments=args.include_vre_investments,
        reservoir_representation=args.reservoir_representation,
        water_policy=args.water_policy,
        release_multiplier=args.release_multiplier,
        solved_network_save_to=args.network,
        build_report_save_to=args.report,
        solve_network=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
