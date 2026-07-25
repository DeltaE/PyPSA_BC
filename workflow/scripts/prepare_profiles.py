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

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.logger import Pipeline


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main():
    from workflow.scripts import (
        create_ext_wind_ts, create_ext_solar_ts,
        create_ror_ps, create_reservoir_inflows,
    )

    dcfg = AttributesParser().data_cfg
    out = dcfg["output"]

    cutout = Path(dcfg["data"]["cutout"])
    if not cutout.exists():
        raise FileNotFoundError(
            f"prepare_profiles: cutout not found ({cutout}). Build the ERA5 cutout "
            f"first (create_cutout), and run prepare_assets so the asset CSVs exist."
        )

    with Pipeline("Profiles", ["wind ts", "solar ts", "ror power", "reservoir inflows"]) as pipe:
        with pipe.stage("wind ts"):
            create_ext_wind_ts.main()
            f = out["create_ext_wind_ts"]["fname"]
            pipe.deliver("bc_ext_wind_ts.csv", f, _rows(f))

        with pipe.stage("solar ts"):
            create_ext_solar_ts.main()
            f = out["create_ext_solar_ts"]["fname"]
            pipe.deliver("bc_ext_solar_ts.csv", f, _rows(f))

        with pipe.stage("ror power"):
            create_ror_ps.main()
            f = out["ror_ps"]["fname"]
            pipe.deliver("bc_ext_ror_ts.csv", f, _rows(f))

        with pipe.stage("reservoir inflows"):
            create_reservoir_inflows.main()
            f = out["reservoir_inflows"]["fname"]
            pipe.deliver("bc_ext_reservoir_inflows.csv", f, _rows(f))


if __name__ == "__main__":
    main()
