"""
prepare_ror — run-of-river power-availability series.

Computes the RoR power series from the hydro generation assets (hydro_type
'ror'/'ror-water'), the HydroBASINS layers, and the ERA5 cutout. Depends on
prepare_hydro (hydro_generation.csv) and an existing atlite cutout.

NOTE: this produces a time series, so it really belongs to the `prepare_profiles`
stage; kept runnable on its own here for the asset build-out.
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
    from workflow.scripts import create_ror_ps  # heavy deps (atlite/geopandas)

    dcfg = AttributesParser().data_cfg
    ror_out = dcfg["output"]["ror_ps"]["fname"]

    need = {
        "hydro_generation.csv": Path(dcfg["output"]["create_hydro_assets"]["hydro_generation"]),
        "cutout": Path(dcfg["data"]["cutout"]),
        "na_basins": Path(dcfg["basin_files"]["na_file"]),
        "arctic_basins": Path(dcfg["basin_files"]["arctic_file"]),
    }
    missing = {k: str(p) for k, p in need.items() if not p.exists()}
    if missing:
        raise FileNotFoundError(
            f"prepare_ror: missing inputs {missing}. Run prepare_hydro first and "
            f"ensure the atlite cutout + HydroBASINS layers are downloaded."
        )

    with Pipeline("RoR power series", ["ror power series"]) as pipe:
        with pipe.stage("ror power series"):
            create_ror_ps.main()
            pipe.deliver("bc_ext_ror_ts.csv", ror_out, _rows(ror_out))


if __name__ == "__main__":
    main()
