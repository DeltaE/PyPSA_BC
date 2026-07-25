"""
prepare_hydro — hydro sub-pipeline (existing generation + reservoirs).

    CODERS tables   generators, generation_generic, existing_hydro (via shared client)
    hydro cascade   reconstruct hydro_cascade.csv from existing_hydro + the static
                    inventory (data/static: hydro_topology.csv, hydro_reservoirs.csv)
                    using to_legacy_cascade
    hydro assets    create_hydro_assets  (CODERS + WUP files + reconstructed cascade)

Static inventory (data/static/) is seeded once by
`pypsa_bc.data.bootstrap_hydro_static` from the last-known legacy cascade file.
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data import to_legacy_cascade
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.logger import Pipeline

STATIC_DIR = "data/static"


def _rows(path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def main():
    from workflow.scripts import create_hydro_assets  # imported here (heavy deps)

    dcfg = AttributesParser().data_cfg
    out = dcfg["output"]["create_hydro_assets"]
    hydro_cascade_out = out["hydro_cascade"]          # produced cascade (output, not raw input)

    # Pre-flight: the reconstruction needs the static inventory + the WUP files.
    need = {
        "hydro_topology.csv": Path(STATIC_DIR) / "hydro_topology.csv",
        "hydro_reservoirs.csv": Path(STATIC_DIR) / "hydro_reservoirs.csv",
        "gen_wup": Path(dcfg["inventory"]["gen_wup"]),
        "res_wup": Path(dcfg["inventory"]["res_wup"]),
        "inflow_tables": Path(dcfg["inventory"]["inflow_tables"]),
    }
    missing = {k: str(p) for k, p in need.items() if not p.exists()}
    if missing:
        raise FileNotFoundError(
            f"prepare_hydro: missing inputs {missing}. Seed data/static/ via "
            f"`python -m pypsa_bc.data.bootstrap_hydro_static --legacy <hydro_cascade.csv>` "
            f"and place the WUP files under data/custom/."
        )

    with Pipeline("Hydro assets", ["CODERS tables", "hydro cascade", "hydro assets"]) as pipe:
        with pipe.stage("CODERS tables"):
            coders = get_coders()
            g = coders.load_table("generators", as_gdf=False)
            gg = coders.load_table("generation_generic", as_gdf=False)
            eh = coders.existing_hydro                      # writes existing_hydro.csv
            pipe.deliver("generators.csv", coders.table_path("generators"), len(g))
            pipe.deliver("generation_generic.csv", coders.table_path("generation_generic"), len(gg))
            pipe.deliver("existing_hydro.csv", coders.table_path("existing_hydro"), len(eh))

        with pipe.stage("hydro cascade"):
            existing_hydro_csv = coders.table_path("existing_hydro")
            to_legacy_cascade.to_legacy_cascade(
                str(existing_hydro_csv), STATIC_DIR, hydro_cascade_out)
            pipe.deliver("hydro_cascade.csv", hydro_cascade_out, _rows(hydro_cascade_out))

        with pipe.stage("hydro assets"):
            create_hydro_assets.main()
            pipe.deliver("hydro_generation.csv", out["hydro_generation"], _rows(out["hydro_generation"]))
            pipe.deliver("hydro_reservoirs.csv", out["hydro_reservoir"], _rows(out["hydro_reservoir"]))


if __name__ == "__main__":
    main()
