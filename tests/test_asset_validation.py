"""Focused unit tests for the Rule 2 asset-validation contract."""

from copy import deepcopy
from pathlib import Path

import pandas as pd
import yaml
from workflow.scripts.validate_assets import validate

CONFIG = yaml.safe_load(
    Path("config/assets_validation.yaml").read_text(encoding="utf-8")
)


def _valid_tables() -> dict[str, pd.DataFrame]:
    common = {"connecting_node_code": "BC_AAA_GSS", "latitude": 49.0, "longitude": -123.0}
    return {
        "wind": pd.DataFrame([{**common, "asset_id": "W1", "generator_capacity": 10.0,
                               "gen_type": "wind_ons", "config_oedb": "E-82/3000"}]),
        "solar": pd.DataFrame([{**common, "asset_id": "S1", "generator_capacity": 5.0,
                                "gen_type": "solar_PV"}]),
        "tpp": pd.DataFrame([{**common, "generation_unit_code": "T1",
                              "generator_capacity": 20.0, "gen_type": "NG_CC"}]),
        "hydro_generation": pd.DataFrame([{
            **common, "asset_id": "H1", "component_id": "HG1", "capacity": 30.0,
            "gen_type": "hydro_daily", "upper_reservoir_id": "R1",
            "lower_reservoir_id": "BC_RIVER_SINK", "spill_lower_reservoir_id": "DEFAULT",
        }]),
        "hydro_reservoirs": pd.DataFrame([{
            "asset_id": "R1", "latitude": 49.1, "longitude": -123.1,
            "min_storage": 0.0, "max_storage": 100.0,
        }]),
        "hydro_cascade": pd.DataFrame([{
            **common, "project_name": "Hydro 1", "gen_node_code": "HG1",
            "gen_type": "Hydro_daily", "install_capacity_in_mw": 30.0,
        }]),
    }


def test_valid_asset_contract_passes() -> None:
    summary, issues = validate(
        _valid_tables(),
        pd.DataFrame([{"name": "230_AAA_GSS"}]),
        {"E82-3.0": "E-82/3000"},
        CONFIG,
    )

    assert summary["status"] == "PASS"
    assert not summary["blocking_checks"]
    assert all(frame.empty for frame in issues.values())


def test_duplicate_identifier_and_unmapped_node_block_rule_2() -> None:
    tables = _valid_tables()
    tables["wind"] = pd.concat([tables["wind"], tables["wind"]], ignore_index=True)
    tables["tpp"].loc[0, "connecting_node_code"] = "BC_UNKNOWN_GSS"

    summary, issues = validate(
        tables,
        pd.DataFrame([{"name": "230_AAA_GSS"}]),
        {"E82-3.0": "E-82/3000"},
        CONFIG,
    )

    assert summary["status"] == "FAIL"
    assert "asset_identifier_integrity" in summary["blocking_checks"]
    assert "base_network_mapping" in summary["blocking_checks"]
    assert len(issues["duplicate_identifiers.csv"]) == 1
    assert len(issues["unmapped_network_nodes.csv"]) == 1


def test_missing_reservoir_storage_is_a_reviewable_warning() -> None:
    tables = deepcopy(_valid_tables())
    tables["hydro_reservoirs"].loc[0, "max_storage"] = None

    summary, issues = validate(
        tables,
        pd.DataFrame([{"name": "230_AAA_GSS"}]),
        {"E82-3.0": "E-82/3000"},
        CONFIG,
    )

    assert summary["status"] == "PASS"
    assert summary["warning_checks"] == ["capacity_evidence_warnings"]
    assert issues["invalid_capacities.csv"].loc[0, "severity"] == "WARNING"
