"""Build PyPSA buses from the prepared CODERS transmission network."""

from pathlib import Path

import pandas as pd

from pypsa_bc import utils
from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.network import corrections, lines
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger

log = get_logger("buses")

aparser = AttributesParser()
coders_data = get_coders()
REQUIRED_COLUMNS = coders_data.required.get("substations", [])
META = {
    "module": Path(__file__).stem,
    "scenario": aparser.get_scenario,
    "path": aparser.assumption_report_save_to,
}
INTERTIE_SUFFIXES = {"IPT", "INT"}


def check_missing_buses(
    prepared_substations: pd.DataFrame,
    prepared_lines: pd.DataFrame,
) -> None:
    """Warn about CODERS substations that are absent from all line endpoints."""
    sub_codes = prepared_substations["node_code"].str.split("_").str[1].str.strip().unique()
    start_codes = prepared_lines["starting_node_code"].str.split("_").str[1].str.strip()
    end_codes = prepared_lines["ending_node_code"].str.split("_").str[1].str.strip()
    line_codes = set(pd.concat([start_codes, end_codes]).dropna())

    for code in sub_codes:
        if code not in line_codes:
            log.warning("substation code %s is missing from the lines dataset", code)


def create_bus_df(
    prepared_lines: pd.DataFrame,
    prepared_substations: pd.DataFrame,
    prepared_generators: pd.DataFrame,
    node_corrections: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Create unique internal buses referenced by the prepared transmission lines.

    Boundary intertie endpoints are deliberately not invented here. Unresolved
    internal endpoints are also skipped, so validation reports the affected line
    references instead of receiving an invalid nameless bus.
    """
    data: dict[str, dict[str, object]] = {}
    skipped_interties: set[str] = set()
    unresolved_nodes: set[str] = set()
    applied_corrections: set[str] = set()
    correction_register = (
        node_corrections if node_corrections is not None else pd.DataFrame()
    )

    for _, line in prepared_lines.iterrows():
        for node_code in (line["starting_node_code"], line["ending_node_code"]):
            if not isinstance(node_code, str):
                unresolved_nodes.add(str(node_code))
                continue

            if node_code in correction_register.index:
                correction = correction_register.loc[node_code]
                bus_name = f"{line['v_nom']}_{'_'.join(node_code.split('_')[1:])}"
                data.setdefault(
                    bus_name,
                    {
                        "x": correction["longitude"],
                        "y": correction["latitude"],
                        "type": line["type"],
                        "v_nom": line["v_nom"],
                    },
                )
                applied_corrections.add(str(correction["correction_id"]))
                continue

            if node_code.split("_")[-1] in INTERTIE_SUFFIXES:
                skipped_interties.add(node_code)
                continue

            bus_name, bus_x, bus_y = get_bus_name_x_y(
                line,
                node_code,
                prepared_substations,
                prepared_generators,
            )
            if bus_name is None:
                unresolved_nodes.add(node_code)
                continue

            data.setdefault(
                bus_name,
                {"x": bus_x, "y": bus_y, "type": line["type"], "v_nom": line["v_nom"]},
            )

    if skipped_interties:
        log_assumption(
            parameter="intertie boundary bus handling",
            value=f"{len(skipped_interties)} external endpoint(s) excluded",
            unit="endpoints",
            rationale=(
                "Rule 1 currently builds the BC-internal network only; explicit boundary "
                "buses and exchange treatment require a separate modelling decision"
            ),
            source="CODERS transmission line endpoints ending IPT/INT",
            **META,
        )
        log.warning(
            "excluded %d external intertie endpoint(s); validation will retain their line references",
            len(skipped_interties),
        )

    if applied_corrections:
        log_assumption(
            parameter="registered representative bus coordinates",
            value=", ".join(sorted(applied_corrections)),
            unit="correction IDs",
            rationale=(
                "Special non-automated treatments from the Rule 1 correction register; "
                "map-derived points are representative rather than surveyed locations"
            ),
            source="data/validation/base_network/node_corrections.csv",
            **META,
        )

    if unresolved_nodes:
        log.warning(
            "could not resolve %d internal line endpoint(s): %s",
            len(unresolved_nodes),
            ", ".join(sorted(unresolved_nodes)),
        )

    columns = ["name", "x", "y", "type", "v_nom"]
    return (
        pd.DataFrame.from_dict(data, orient="index")
        .rename_axis("name")
        .reset_index()
        .reindex(columns=columns)
    )


def get_bus_name_x_y(
    line: pd.Series,
    node_code: str,
    substations: pd.DataFrame,
    generators: pd.DataFrame,
) -> tuple[str | None, float | None, float | None]:
    """Resolve one line endpoint to a PyPSA bus name and coordinates."""
    bus_name = f"{line['v_nom']}_{'_'.join(node_code.split('_')[1:])}"

    exact_substations = substations[substations["node_code"] == node_code]
    if not exact_substations.empty:
        match = exact_substations.iloc[0]
        return bus_name, match["longitude"], match["latitude"]

    exact_generators = generators[generators["connecting_node_code"] == node_code]
    if not exact_generators.empty:
        match = exact_generators.iloc[0]
        return bus_name, match["longitude"], match["latitude"]

    node_token = node_code.split("_")[1]
    sub_tokens = substations["node_code"].str.split("_").str[1]
    partial_substations = substations[sub_tokens == node_token]
    if not partial_substations.empty:
        match = partial_substations.iloc[0]
        return bus_name, match["longitude"], match["latitude"]

    generator_tokens = generators["connecting_node_code"].str.split("_").str[1]
    partial_generators = generators[generator_tokens == node_token]
    if not partial_generators.empty:
        match = partial_generators.iloc[0]
        return bus_name, match["longitude"], match["latitude"]

    log.warning("no coordinate source found for bus endpoint %s", node_code)
    return None, None, None


def sanitize(substations: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize the configured CODERS substation columns."""
    assert isinstance(substations, pd.DataFrame), "Input must be a pandas DataFrame"
    assert "node_code" in substations.columns, "DataFrame must contain 'node_code' column"

    for column in REQUIRED_COLUMNS:
        assert column in substations.columns, f"DataFrame must contain '{column}' column"
        if pd.api.types.is_string_dtype(substations[column]):
            substations[column] = substations[column].str.replace(r"\s+", "", regex=True)

    return substations


def get_prepared_buses(prepared_lines: pd.DataFrame | None = None) -> pd.DataFrame:
    """Prepare and save the Rule 1 PyPSA bus table."""
    utils.print_update(level=1, message="Preparing base network for PyPSA-BC")

    substations = sanitize(coders_data.substations.copy())
    prepared_lines_df = (
        prepared_lines.copy() if prepared_lines is not None else lines.get_prepared_lines()[0]
    )
    generators = coders_data.generators.copy()
    correction_path = aparser.data_cfg["inventory"]["base_network_node_corrections"]
    node_corrections = corrections.load_node_corrections(correction_path)

    buses = create_bus_df(
        prepared_lines=prepared_lines_df,
        prepared_substations=substations,
        prepared_generators=generators,
        node_corrections=node_corrections,
    )
    check_missing_buses(prepared_substations=substations, prepared_lines=prepared_lines_df)

    output_path = Path(aparser.data_cfg["output"]["base_network"]) / "buses.csv"
    buses.to_csv(output_path, index=False)
    log.debug("buses.csv saved to %s", output_path)
    return buses
