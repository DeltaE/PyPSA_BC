"""Focused unit tests for Rule 1 base-network preparation components."""

from unittest.mock import patch

import pandas as pd

from pypsa_bc.network import buses, corrections, lines, transformers


def test_bus_builder_skips_external_and_unresolved_endpoints_without_null_bus() -> None:
    prepared_lines = pd.DataFrame(
        [
            {
                "starting_node_code": "BC_AAA_GSS",
                "ending_node_code": "BC_ABBC01_IPT",
                "type": "230kV",
                "v_nom": 230,
            },
            {
                "starting_node_code": "BC_AAA_GSS",
                "ending_node_code": "BC_GST_JCT",
                "type": "230kV",
                "v_nom": 230,
            },
        ]
    )
    substations = pd.DataFrame(
        [{"node_code": "BC_AAA_GSS", "longitude": -123.0, "latitude": 49.0}]
    )
    generators = pd.DataFrame(
        columns=["connecting_node_code", "longitude", "latitude"]
    )

    with patch.object(buses, "log_assumption"):
        result = buses.create_bus_df(prepared_lines, substations, generators)

    assert result["name"].tolist() == ["230_AAA_GSS"]
    assert result["name"].notna().all()


def test_registered_external_endpoint_becomes_representative_bus() -> None:
    prepared_lines = pd.DataFrame(
        [{
            "starting_node_code": "BC_AAA_GSS",
            "ending_node_code": "BC_ABBC01_IPT",
            "type": "138kV",
            "v_nom": 138,
        }]
    )
    substations = pd.DataFrame(
        [{"node_code": "BC_AAA_GSS", "longitude": -115.0, "latitude": 50.0}]
    )
    generators = pd.DataFrame(
        columns=["connecting_node_code", "longitude", "latitude"]
    )
    register = pd.DataFrame(
        [{
            "correction_id": "NODE-TEST",
            "node_code": "BC_ABBC01_IPT",
            "longitude": -114.6,
            "latitude": 49.7,
        }]
    ).set_index("node_code", drop=False)

    with patch.object(buses, "log_assumption"):
        result = buses.create_bus_df(
            prepared_lines,
            substations,
            generators,
            node_corrections=register,
        )

    boundary = result.set_index("name").loc["138_ABBC01_IPT"]
    assert boundary["x"] == -114.6
    assert boundary["y"] == 49.7


def test_line_attributes_use_the_same_bus_naming_contract() -> None:
    source = pd.DataFrame(
        [
            {
                "starting_node_code": "BC_AAA_GSS",
                "ending_node_code": "BC_BBB_DSS",
                "voltage_kV": 230,
                "line_segment_length_km": 12.5,
            }
        ]
    )

    result = lines.add_pypsa_attributes(source.copy())

    assert result.loc[0, "bus0"] == "230_AAA_GSS"
    assert result.loc[0, "bus1"] == "230_BBB_DSS"
    assert result.loc[0, "type"] == "230kV"


def test_line_capacity_prefers_mva_and_falls_back_to_power_factor() -> None:
    with patch.object(lines, "log_assumption"):
        direct = lines.add_line_s_nom(
            pd.Series({"summer_rating_in_mva": 125.0, "summer_rating_mw": 90.0})
        )
        fallback = lines.add_line_s_nom(
            pd.Series({"summer_rating_in_mva": None, "summer_rating_mw": 90.0})
        )

    assert direct == 125.0
    assert fallback == round(90.0 / lines.BASE_NETWORK_CFG.get("pf", 0.9), 4)


def test_transformer_builder_deduplicates_voltage_levels() -> None:
    prepared_buses = pd.DataFrame(
        [
            {"name": "69_AAA_GSS", "v_nom": 69},
            {"name": "69_AAA_GSS", "v_nom": 69},
            {"name": "138_AAA_GSS", "v_nom": 138},
            {"name": "230_AAA_GSS", "v_nom": 230},
        ]
    )

    result = transformers.create_transformer_df(prepared_buses)

    assert len(result) == 2
    assert result["type"].tolist() == ["138/69", "230/138"]
    assert result["name"].is_unique


def test_registered_line_exclusion_and_rating_patch_are_applied() -> None:
    source = pd.DataFrame(
        [
            {"transmission_line_id": 1, "summer_rating_mw": 0.0},
            {"transmission_line_id": 2, "summer_rating_mw": 10.0},
        ]
    )
    register = pd.DataFrame(
        [
            {
                "correction_id": "LINE-A",
                "transmission_line_id": "1",
                "action": "replace_value",
                "field": "summer_rating_mw",
                "original_value": 0.0,
                "replacement_value": 65.0,
            },
            {
                "correction_id": "LINE-B",
                "transmission_line_id": "2",
                "action": "exclude_line",
                "field": "",
                "original_value": None,
                "replacement_value": None,
            },
        ]
    )

    result = corrections.apply_line_corrections(source, register)

    assert result["transmission_line_id"].tolist() == [1]
    assert result.loc[0, "summer_rating_mw"] == 65.0
