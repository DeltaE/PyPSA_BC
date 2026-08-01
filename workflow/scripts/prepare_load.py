"""Prepare regional hourly load files as the load branch of workflow stage 3."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from workflow.scripts import disaggregate_load

from pypsa_bc.attributes_parser import AttributesParser


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2021)
    args = parser.parse_args()
    cfg = AttributesParser().data_cfg
    source = Path(f"{cfg['data']['load']['bch']}{args.year}.xls")
    if not source.exists():
        raise FileNotFoundError(f"BC Hydro load file not found: {source}")
    raw = pd.read_excel(source)
    hourly = disaggregate_load.fix_hourly_load(raw, args.year)
    disaggregate_load.main(float(hourly["LOAD"].sum()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
