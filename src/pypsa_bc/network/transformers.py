"""Create transformers between voltage levels at shared physical nodes."""

from itertools import pairwise
from pathlib import Path

import pandas as pd

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger

log = get_logger("transformers")

aparser = AttributesParser()
BASE_NETWORK_CFG = aparser.base_network_cfg
META = {
    "module": Path(__file__).stem,
    "scenario": aparser.get_scenario,
    "path": aparser.assumption_report_save_to,
}


def create_transformer_df(prepared_buses: pd.DataFrame | None = None) -> pd.DataFrame:
    """Connect adjacent unique voltage levels found at each physical bus."""
    assert prepared_buses is not None, "prepared_buses cannot be None"

    voltages_by_bus: dict[str, set[float]] = {}
    skipped_indices: list[object] = []
    for index, row in prepared_buses.iterrows():
        name = row["name"]
        if not isinstance(name, str):
            skipped_indices.append(index)
            continue
        physical_bus = "_".join(name.split("_")[1:])
        voltages_by_bus.setdefault(physical_bus, set()).add(row["v_nom"])

    if skipped_indices:
        log.warning(
            "skipped %d bus row(s) with non-string names: %s",
            len(skipped_indices),
            skipped_indices,
        )
        log_assumption(
            parameter="transformer build - skipped nameless buses",
            value=f"{len(skipped_indices)} bus row(s) skipped",
            unit="rows",
            rationale="A bus without a valid name cannot form a physical transformer key",
            source="create_transformer_df",
            **META,
        )

    transformers: list[list[object]] = []
    for physical_bus, voltage_set in voltages_by_bus.items():
        voltages = sorted(voltage_set)
        for low_voltage, high_voltage in pairwise(voltages):
            transformer_type = f"{high_voltage}/{low_voltage}"
            transformers.append(
                [
                    f"{physical_bus}_{high_voltage}_{low_voltage}",
                    f"{high_voltage}_{physical_bus}",
                    f"{low_voltage}_{physical_bus}",
                    transformer_type,
                ]
            )

    return pd.DataFrame(transformers, columns=["name", "bus0", "bus1", "type"])


def create_transformer_types_df(transformers: pd.DataFrame) -> pd.DataFrame:
    """Create one standardized parameter row for each unique transformer type."""
    attributes = BASE_NETWORK_CFG.get("transformer", {})
    rows: list[list[object]] = []

    for transformer_type in transformers["type"].drop_duplicates():
        high_voltage, low_voltage = (int(float(value)) for value in transformer_type.split("/"))
        rows.append(
            [
                transformer_type,
                BASE_NETWORK_CFG.get("f_nom", 60),
                attributes.get("s_nom", 2000),
                high_voltage,
                low_voltage,
                attributes.get("vsc", 10),
                attributes.get("vscr", 0.3),
                attributes.get("pfe", 30),
                attributes.get("i0", 0.04),
                attributes.get("phase_shift", 150),
                attributes.get("tap_side", 0),
                attributes.get("tap_neutral", 0),
                attributes.get("tap_min", -9),
                attributes.get("tap_max", 9),
                attributes.get("tap_step", 1.5),
            ]
        )

    columns = [
        "name",
        "f_nom",
        "s_nom",
        "v_nom_0",
        "v_nom_1",
        "vsc",
        "vscr",
        "pfe",
        "i0",
        "phase_shift",
        "tap_side",
        "tap_neutral",
        "tap_min",
        "tap_max",
        "tap_step",
    ]
    result = pd.DataFrame(rows, columns=columns)

    if not result.empty:
        log_assumption(
            parameter="transformer impedance (type parameters)",
            value="vsc=10%, vscr=0.3%, i0=0.04%, pfe=30 kW",
            unit="mixed",
            rationale="Standardized default nameplate values; refine with unit-specific data",
            source="base network preparation",
            **META,
        )
    return result


def get_prepared_transformers(
    prepared_buses: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare and save the Rule 1 transformer and transformer-type tables."""
    assert prepared_buses is not None, "prepared_buses cannot be None"

    transformers = create_transformer_df(prepared_buses=prepared_buses)
    transformer_types = create_transformer_types_df(transformers)

    output_directory = Path(aparser.data_cfg["output"]["base_network"])
    transformers_path = output_directory / "transformers.csv"
    types_path = output_directory / "transformer_types.csv"
    transformers.to_csv(transformers_path, index=False)
    transformer_types.to_csv(types_path, index=False)
    log.debug("transformers.csv saved to %s", transformers_path)
    log.debug("transformer_types.csv saved to %s", types_path)

    return transformers, transformer_types
