"""
prepare_solar — solar sub-pipeline: ensure CODERS tables, build solar assets.

    CODERS tables   generators + generation_generic (cached via shared client)
    solar assets    create_ext_solar_assets
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.logger import Pipeline
from workflow.scripts import create_ext_solar_assets


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main():
    out = AttributesParser().data_cfg["output"]
    with Pipeline("Solar assets", ["CODERS tables", "solar assets"]) as pipe:
        with pipe.stage("CODERS tables"):
            coders = get_coders()
            g = coders.load_table("generators", as_gdf=False)
            gg = coders.load_table("generation_generic", as_gdf=False)
            pipe.deliver("generators.csv", coders.table_path("generators"), len(g))
            pipe.deliver("generation_generic.csv", coders.table_path("generation_generic"), len(gg))

        with pipe.stage("solar assets"):
            create_ext_solar_assets.main()
            f = out["create_ext_solar_assets"]["fname"]
            pipe.deliver("solar_assets.csv", f, _rows(f))


if __name__ == "__main__":
    main()
