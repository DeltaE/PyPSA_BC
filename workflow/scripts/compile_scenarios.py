"""Compile reservoir-representation x water-policy cases from workflow YAML."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml

from pypsa_bc.studies.cascade.scenarios import expand_scenarios


def compile_scenarios(config_path: Path, output_directory: Path) -> list[dict]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    cases = expand_scenarios(config)
    output_directory.mkdir(parents=True, exist_ok=True)

    for case in cases:
        case_path = output_directory / f"{case['case_id']}.yaml"
        case_path.write_text(
            yaml.safe_dump(case, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )

    manifest = output_directory / "scenario_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "case_id",
                "reservoir_representation",
                "water_use_policy",
                "storage_policy",
                "runnable",
                "representation_status",
                "policy_constraint_mode",
                "release_multiplier",
                "storage_constraint_mode",
            ],
        )
        writer.writeheader()
        for case in cases:
            writer.writerow(
                {
                    "case_id": case["case_id"],
                    "reservoir_representation": case["reservoir_representation"],
                    "water_use_policy": case["water_use_policy"],
                    "storage_policy": case["storage_policy"],
                    "runnable": case["runnable"],
                    "representation_status": case["representation"].get("status", ""),
                    "policy_constraint_mode": case["policy"].get("constraint_mode", ""),
                    "release_multiplier": case["policy"].get("release_multiplier", 1.0),
                    "storage_constraint_mode": case["storage"].get(
                        "constraint_mode", ""
                    ),
                }
            )
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/workflow.yaml"))
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("results/workflow/scenarios/cases"),
    )
    args = parser.parse_args()
    cases = compile_scenarios(args.config, args.output_directory)
    runnable = sum(bool(case["runnable"]) for case in cases)
    print(f"Compiled {len(cases)} cases ({runnable} currently runnable)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
