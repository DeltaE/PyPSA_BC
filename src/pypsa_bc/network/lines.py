"""Prepare CODERS transmission lines for the PyPSA base network."""

from pathlib import Path

import pandas as pd

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.network import corrections
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger

log = get_logger("lines")

aparser = AttributesParser()
BASE_NETWORK_CFG = aparser.base_network_cfg
coders_data = get_coders()
REQUIRED_COLUMNS = coders_data.required.get("transmission_lines", [])
META = {
    "module": Path(__file__).stem,
    "scenario": aparser.get_scenario,
    "path": aparser.assumption_report_save_to,
}


def sanitize(lines: pd.DataFrame) -> pd.DataFrame:
    """Normalize ratings and endpoint codes in a transmission-line table."""
    assert isinstance(lines, pd.DataFrame), "Input must be a pandas DataFrame"
    assert "summer_rating_mw" in lines.columns, (
        "DataFrame must contain 'summer_rating_mw' column"
    )

    lines.loc[:, "summer_rating_mw"] = (
        pd.to_numeric(lines["summer_rating_mw"], errors="coerce").fillna(0.0).astype(float)
    )
    for column in ("starting_node_code", "ending_node_code"):
        assert column in lines.columns, f"DataFrame must contain '{column}' column"
        lines[column] = lines[column].str.replace(r"\s+", "", regex=True)
    return lines


def add_line_op_params(lines: pd.DataFrame) -> pd.DataFrame:
    """Add nominal apparent-power capacity to each line."""
    lines["s_nom"] = lines.apply(add_line_s_nom, axis=1)
    log.debug("s_nom calculated for all lines; assumption report updated")
    return lines


def get_line_table() -> pd.DataFrame:
    """Load the configured conductor inventory exported as line types."""
    inventory_path = aparser.data_cfg["inventory"]["line_table"]
    return pd.read_excel(inventory_path)


def add_pypsa_attributes(lines: pd.DataFrame) -> pd.DataFrame:
    """Add PyPSA identifiers and structural attributes to CODERS lines."""
    names: list[str] = []
    line_types: list[str] = []
    buses_0: list[str] = []
    buses_1: list[str] = []
    lengths: list[float] = []
    nominal_voltages: list[float] = []

    for _, line in lines.iterrows():
        voltage = int(line["voltage_kV"])
        voltage_type = f"{voltage}kV"
        names.append(line["starting_node_code"][3:] + line["ending_node_code"][2:])
        line_types.append(voltage_type)
        buses_0.append(f"{voltage}_{'_'.join(line['starting_node_code'].split('_')[1:])}")
        buses_1.append(f"{voltage}_{'_'.join(line['ending_node_code'].split('_')[1:])}")
        lengths.append(line["line_segment_length_km"])
        nominal_voltages.append(line["voltage_kV"])

    lines["name"] = names
    lines["type"] = line_types
    lines["bus0"] = buses_0
    lines["bus1"] = buses_1
    lines["length"] = lengths
    lines["v_nom"] = nominal_voltages
    return lines


def add_line_s_nom(line: pd.Series) -> float:
    """Calculate line capacity from MVA or the MW rating and configured power factor."""
    candidate_columns = ["summer_rating_in_mva", "summer_rating_mw"]
    assert any(column in line.index for column in candidate_columns), (
        "CODERS | table | transmission_lines | Expected at least one rating column "
        f"from {candidate_columns}, but found none in columns: {line.index.tolist()}"
    )

    if "summer_rating_in_mva" in line.index and pd.notna(line["summer_rating_in_mva"]):
        s_nom = float(line["summer_rating_in_mva"])
    else:
        power_factor = BASE_NETWORK_CFG.get("pf", 0.9)
        s_nom = round(float(line["summer_rating_mw"]) / power_factor, 4)

    log_assumption(
        parameter="line s_nom",
        value="summer_rating_mw (ttc_summer)",
        unit="MW",
        rationale=(
            "CODERS Total Transfer Capability (TTC) used as thermal proxy, converted "
            "to MVA using the default power factor "
            f"({BASE_NETWORK_CFG.get('pf', 0.9)}) from base_network.yaml"
        ),
        **META,
    )
    return s_nom


def get_prepared_lines(
    input_lines: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare and save the Rule 1 line and conductor tables."""
    assert input_lines is None or isinstance(input_lines, pd.DataFrame), (
        "Input must be a pandas DataFrame or None (default to load from CODERS)"
    )

    source_lines = coders_data.transmission_lines if input_lines is None else input_lines
    prepared = source_lines.copy()
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in prepared.columns]
    assert not missing_columns, (
        "CODERS | table | transmission_lines | Required columns missing: "
        f"{missing_columns}"
    )

    correction_path = aparser.data_cfg["inventory"]["base_network_line_corrections"]
    correction_register = corrections.load_line_corrections(correction_path)
    prepared = corrections.apply_line_corrections(prepared, correction_register)
    prepared = add_pypsa_attributes(add_line_op_params(sanitize(prepared)))
    line_types = get_line_table()

    output_directory = Path(aparser.data_cfg["output"]["base_network"])
    lines_path = output_directory / "lines.csv"
    line_types_path = output_directory / "line_types.csv"
    prepared.to_csv(lines_path, index=False)
    line_types.to_csv(line_types_path, index=False)
    log.debug("lines.csv saved to %s", lines_path)
    log.debug("line_types.csv saved to %s", line_types_path)

    return prepared, line_types
