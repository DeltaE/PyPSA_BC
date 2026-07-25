"""
prepare_wind — wind sub-pipeline: fetch CWTD, build the turbine dict, build assets.

One staged run:
    1. fetch CWTD      download canada_turbines.xlsx (skip if present)
    2. turbine dict    resolve models via OEDB (atlite); fall back to the inventory
    3. wind assets     create_ext_wind_assets

Equivalent to running, in order:
    fetch_inputs.main(only=["can_turbines"])
    build_turbine_dict.main()
    create_ext_wind_assets.main()
but under a single progress bar + delivered-files table.
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.data.downloader import download_file
from pypsa_bc.reporting.logger import Pipeline
from workflow.scripts import build_turbine_dict, create_ext_wind_assets


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main(force_download: bool = False):
    aparser = AttributesParser()
    dcfg = aparser.data_cfg
    remote = dcfg.get("remote", {})
    turbine_dict_path = dcfg["data"]["wind"]["turbine_dict"]
    assets_path = dcfg["output"]["create_ext_wind_assets"]["fname"]

    with Pipeline("Wind assets", ["fetch CWTD", "CODERS tables", "turbine dict", "wind assets"]) as pipe:
        with pipe.stage("fetch CWTD"):
            spec = remote["can_turbines"]
            p = download_file(spec["url"], spec["dest"], force=force_download)
            pipe.deliver(Path(p).name, p, f"{Path(p).stat().st_size/1e6:.1f} MB")

        with pipe.stage("CODERS tables"):
            coders = get_coders()   # shared client; caches raw CSVs to disk
            g = coders.load_table("generators", as_gdf=False)
            gg = coders.load_table("generation_generic", as_gdf=False)
            pipe.deliver("generators.csv", coders.table_path("generators"), len(g))
            pipe.deliver("generation_generic.csv", coders.table_path("generation_generic"), len(gg))

        with pipe.stage("turbine dict"):
            mapping, fell_back, unmatched = build_turbine_dict.build(turbine_dict_path)
            note = f"{len(mapping)} models"
            if fell_back:
                note += f", {len(fell_back)} inventory-fallback"
            pipe.deliver("turbine_dict.json", turbine_dict_path, note)

        with pipe.stage("wind assets"):
            create_ext_wind_assets.main()
            pipe.deliver("wind_assets.csv", assets_path, _rows(assets_path))


if __name__ == "__main__":
    main()
