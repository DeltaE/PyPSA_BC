#!/usr/bin/env python3

import argparse
from pathlib import Path

from pypsa_bc.vis.network_visual_report import build_visual_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a static + interactive visual report for a PyPSA network .nc file."
    )
    parser.add_argument(
        "--network",
        default="results/PyPSA_BC_network_s_2021.nc",
        help="Path to solved PyPSA network NetCDF file.",
    )
    parser.add_argument(
        "--vis-dir",
        default="vis",
        help="Output directory for final HTML report and vis_resources.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_visual_report(network_path=Path(args.network), vis_dir=Path(args.vis_dir))

    print(f"Report written to: {result.report_path}")
    print(f"Plot resources written to: {result.resources_dir}")
    print(f"Artifacts generated: {len(result.artifacts)}")


if __name__ == "__main__":
    main()
