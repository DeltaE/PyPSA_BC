"""Perform a read-only preflight and write a machine-readable workflow report."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path

import yaml

from pypsa_bc.studies.cascade.scenarios import expand_scenarios


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/workflow.yaml"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cases = expand_scenarios(config)
    required = [Path(path) for path in config["preflight"]["required_paths"]]
    missing = [path.as_posix() for path in required if not path.exists()]
    report = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "status": "PASS" if not missing else "BLOCKED",
        "python": sys.version,
        "platform": platform.platform(),
        "workflow_config": args.config.as_posix(),
        "scenario_count": len(cases),
        "runnable_scenario_count": sum(bool(case["runnable"]) for case in cases),
        "missing_required_paths": missing,
        "solve_enabled": bool(config.get("execution", {}).get("allow_solve", False)),
        "note": (
            "Preflight is read-only. A PASS does not imply that scientific gates pass "
            "or that production solves are authorized."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Workflow preflight {report['status']}: {args.output}")
    return 0 if not missing else 2


if __name__ == "__main__":
    raise SystemExit(main())
