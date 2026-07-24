"""
prepare_base_network — stage 1 of the PyPSA-BC workflow.

Thin orchestrator: builds the base electrical network (lines -> buses ->
transformers) via the pypsa_bc.network package and writes the CSVs to the
output.base_network directory. Console shows a staged checklist + progress bar
and a delivered-files table; full detail is written to logs/pypsa_bc.log.
"""

from pathlib import Path

import pandas as pd

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.network import buses, lines, transformers
from pypsa_bc.reporting.logger import Pipeline


def main() -> dict[str, pd.DataFrame]:
    aparser = AttributesParser()
    out = Path(aparser.data_cfg["output"]["base_network"])
    out.mkdir(parents=True, exist_ok=True)

    r: dict[str, pd.DataFrame] = {}
    with Pipeline("Base network", ["lines", "buses", "transformers"]) as pipe:
        with pipe.stage("lines"):
            r["lines"], r["line_types"] = lines.get_prepared_lines()
            pipe.deliver("lines.csv", out / "lines.csv", len(r["lines"]))
            pipe.deliver("line_types.csv", out / "line_types.csv", len(r["line_types"]))

        with pipe.stage("buses"):
            r["buses"] = buses.get_prepared_buses(r["lines"])
            pipe.deliver("buses.csv", out / "buses.csv", len(r["buses"]))

        with pipe.stage("transformers"):
            r["transformers"], r["transformer_types"] = \
                transformers.get_prepared_transformers(r["buses"])
            pipe.deliver("transformers.csv", out / "transformers.csv", len(r["transformers"]))
            pipe.deliver("transformer_types.csv", out / "transformer_types.csv",
                         len(r["transformer_types"]))

    # return r


if __name__ == "__main__":
    main()
