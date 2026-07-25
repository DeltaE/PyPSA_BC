"""
prepare_tpp — thermal power plant sub-pipeline.

    CODERS tables   generators + generation_generic (cached via shared client)
    tpp assets      create_ext_tpp_assets  (maps generators to base-network buses)

Depends on prepare_base_network having written buses.csv.
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.logger import Pipeline
from workflow.scripts import create_ext_tpp_assets


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main():
    dcfg = AttributesParser().data_cfg
    out = dcfg["output"]

    buses_csv = Path(dcfg["output"]["base_network"]) / "buses.csv"
    if not buses_csv.exists():
        raise FileNotFoundError(
            f"{buses_csv} not found — run prepare_base_network first "
            f"(TPP maps generators onto the network buses)."
        )

    with Pipeline("TPP assets", ["CODERS tables", "tpp assets"]) as pipe:
        with pipe.stage("CODERS tables"):
            coders = get_coders()
            g = coders.load_table("generators", as_gdf=False)
            gg = coders.load_table("generation_generic", as_gdf=False)
            pipe.deliver("generators.csv", coders.table_path("generators"), len(g))
            pipe.deliver("generation_generic.csv", coders.table_path("generation_generic"), len(gg))

        with pipe.stage("tpp assets"):
            create_ext_tpp_assets.main()
            f = out["create_ext_tpp_assets"]["fname"]
            pipe.deliver("tpp_assets.csv", f, _rows(f))


if __name__ == "__main__":
    main()
