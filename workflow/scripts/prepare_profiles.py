"""
prepare_profiles — stage 3 of the PyPSA-BC workflow (time series).

Builds the hourly profiles that feed the network:
    wind ts            capacity-factor series for existing wind (atlite + GWA)
    solar ts           capacity-factor series for existing solar (atlite)
    ror power          run-of-river power-availability series (atlite + basins)
    reservoir inflows  reservoir inflow series (atlite + basins + WUP calibration)

All four need the ERA5 cutout and the outputs of prepare_assets. Each is also
runnable on its own (create_ext_wind_ts.main(), … / prepare_ror.main()).
"""

import argparse
from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.logger import Pipeline


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main(
    run_wind_ts: bool = True,
    run_solar_ts: bool = True,
    run_ror_power: bool = True,
    run_reservoir_inflows: bool = True,
    force_update: bool = False,
):
    from workflow.scripts import (
        create_ext_wind_ts, create_ext_solar_ts,
        create_ror_ps, create_reservoir_inflows,
    )

    dcfg = AttributesParser().data_cfg
    out = dcfg["output"]

    stage_cfg = {
        "wind ts": {
            "enabled": run_wind_ts,
            "output": Path(out["create_ext_wind_ts"]["fname"]),
            "runner": create_ext_wind_ts.main,
            "artifact": "bc_ext_wind_ts.csv",
            "needs": {
                "cutout": Path(dcfg["data"]["cutout"]),
                "wind assets": Path(out["create_ext_wind_assets"]["fname"]),
            },
        },
        "solar ts": {
            "enabled": run_solar_ts,
            "output": Path(out["create_ext_solar_ts"]["fname"]),
            "runner": create_ext_solar_ts.main,
            "artifact": "bc_ext_solar_ts.csv",
            "needs": {
                "cutout": Path(dcfg["data"]["cutout"]),
                "solar assets": Path(out["create_ext_solar_assets"]["fname"]),
            },
        },
        "ror power": {
            "enabled": run_ror_power,
            "output": Path(out["ror_ps"]["fname"]),
            "runner": create_ror_ps.main,
            "artifact": "bc_ext_ror_ts.csv",
            "needs": {
                "cutout": Path(dcfg["data"]["cutout"]),
                "hydro generation": Path(out["create_hydro_assets"]["hydro_generation"]),
                "na basins": Path(dcfg["basin_files"]["na_file"]),
                "arctic basins": Path(dcfg["basin_files"]["arctic_file"]),
            },
        },
        "reservoir inflows": {
            "enabled": run_reservoir_inflows,
            "output": Path(out["reservoir_inflows"]["fname"]),
            "runner": create_reservoir_inflows.main,
            "artifact": "bc_ext_reservoir_inflows.csv",
            "needs": {
                "cutout": Path(dcfg["data"]["cutout"]),
                "hydro generation": Path(out["create_hydro_assets"]["hydro_generation"]),
                "hydro reservoirs": Path(out["create_hydro_assets"]["hydro_reservoir"]),
                "na basins": Path(dcfg["basin_files"]["na_file"]),
                "arctic basins": Path(dcfg["basin_files"]["arctic_file"]),
                "inflow tables": Path(dcfg["inventory"]["inflow_tables"]),
            },
        },
    }

    # Preflight only for stages that will actually execute.
    missing = {}
    for cfg in stage_cfg.values():
        if not cfg["enabled"]:
            continue
        if cfg["output"].exists() and not force_update:
            continue
        for name, path in cfg["needs"].items():
            if not path.exists():
                missing[name] = str(path)

    if missing:
        raise FileNotFoundError(
            "prepare_profiles: missing required inputs for selected stages. "
            f"Missing: {missing}"
        )

    stages = ["wind ts", "solar ts", "ror power", "reservoir inflows"]
    with Pipeline("Profiles", stages) as pipe:
        with pipe.stage("wind ts"):
            cfg = stage_cfg["wind ts"]
            f = str(cfg["output"])
            if not cfg["enabled"]:
                pipe.deliver(cfg["artifact"], f, "skipped (disabled)")
            elif cfg["output"].exists() and not force_update:
                pipe.deliver(cfg["artifact"], f, f"{_rows(f)} (exists locally)")
            else:
                cfg["runner"]()
                pipe.deliver(cfg["artifact"], f, _rows(f))

        with pipe.stage("solar ts"):
            cfg = stage_cfg["solar ts"]
            f = str(cfg["output"])
            if not cfg["enabled"]:
                pipe.deliver(cfg["artifact"], f, "skipped (disabled)")
            elif cfg["output"].exists() and not force_update:
                pipe.deliver(cfg["artifact"], f, f"{_rows(f)} (exists locally)")
            else:
                cfg["runner"]()
                pipe.deliver(cfg["artifact"], f, _rows(f))

        with pipe.stage("ror power"):
            cfg = stage_cfg["ror power"]
            f = str(cfg["output"])
            if not cfg["enabled"]:
                pipe.deliver(cfg["artifact"], f, "skipped (disabled)")
            elif cfg["output"].exists() and not force_update:
                pipe.deliver(cfg["artifact"], f, f"{_rows(f)} (exists locally)")
            else:
                cfg["runner"]()
                pipe.deliver(cfg["artifact"], f, _rows(f))

        with pipe.stage("reservoir inflows"):
            cfg = stage_cfg["reservoir inflows"]
            f = str(cfg["output"])
            if not cfg["enabled"]:
                pipe.deliver(cfg["artifact"], f, "skipped (disabled)")
            elif cfg["output"].exists() and not force_update:
                pipe.deliver(cfg["artifact"], f, f"{_rows(f)} (exists locally)")
            else:
                cfg["runner"]()
                pipe.deliver(cfg["artifact"], f, _rows(f))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force-update",
        action="store_true",
        help="Rebuild selected profiles even when their output files already exist.",
    )
    args = parser.parse_args()
    main(force_update=args.force_update)
