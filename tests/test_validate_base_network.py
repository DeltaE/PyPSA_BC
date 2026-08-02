from pathlib import Path

import pandas as pd
import yaml
from workflow.scripts.validate_base_network import validate


def _write_fixture(folder: Path, *, invalid: bool = False) -> None:
    folder.mkdir(parents=True)
    pd.DataFrame(
        [
            {"name": "230_A", "x": -123.0, "y": 49.0, "type": "230kV", "v_nom": 230},
            {"name": "230_B", "x": -122.0, "y": 50.0, "type": "230kV", "v_nom": 230},
        ]
    ).to_csv(folder / "buses.csv", index=False)
    pd.DataFrame(
        [{
            "name": "A_B", "transmission_line_id": 1, "bus0": "230_A",
            "bus1": "MISSING" if invalid else "230_B", "type": "230kV",
            "v_nom": 230, "length": 10.0, "s_nom": 100.0,
        }]
    ).to_csv(folder / "lines.csv", index=False)
    pd.DataFrame([{"code": "230kV"}]).to_csv(folder / "line_types.csv", index=False)
    pd.DataFrame(columns=["name", "bus0", "bus1", "type"]).to_csv(
        folder / "transformers.csv", index=False
    )
    pd.DataFrame(columns=["name", "s_nom", "v_nom_0", "v_nom_1"]).to_csv(
        folder / "transformer_types.csv", index=False
    )


def _config() -> dict:
    path = Path("config/base_network_validation.yaml")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_valid_connected_fixture_passes(tmp_path: Path) -> None:
    network_dir = tmp_path / "network"
    _write_fixture(network_dir)

    summary, issues = validate(network_dir, _config())

    assert summary["status"] == "PASS"
    assert summary["blocking_checks"] == []
    assert issues["missing_line_endpoints.csv"].empty


def test_missing_endpoint_fails_and_is_preserved(tmp_path: Path) -> None:
    network_dir = tmp_path / "network"
    _write_fixture(network_dir, invalid=True)

    summary, issues = validate(network_dir, _config())

    assert summary["status"] == "FAIL"
    assert "line_endpoint_integrity" in summary["blocking_checks"]
    assert len(issues["missing_line_endpoints.csv"]) == 1
    assert bool(issues["missing_line_endpoints.csv"].iloc[0]["missing_bus1"])


def test_exact_approved_island_membership_passes(tmp_path: Path) -> None:
    network_dir = tmp_path / "network"
    _write_fixture(network_dir)
    buses = pd.read_csv(network_dir / "buses.csv")
    buses = pd.concat(
        [
            buses,
            pd.DataFrame([
                {"name": "230_C", "x": -121.0, "y": 51.0, "type": "230kV", "v_nom": 230},
                {"name": "230_D", "x": -120.0, "y": 52.0, "type": "230kV", "v_nom": 230},
                {"name": "230_E", "x": -122.5, "y": 50.5, "type": "230kV", "v_nom": 230},
            ]),
        ],
        ignore_index=True,
    )
    buses.to_csv(network_dir / "buses.csv", index=False)
    lines = pd.read_csv(network_dir / "lines.csv")
    lines = pd.concat(
        [
            lines,
            pd.DataFrame([
                {
                    "name": "C_D", "transmission_line_id": 2, "bus0": "230_C",
                    "bus1": "230_D", "type": "230kV", "v_nom": 230,
                    "length": 5.0, "s_nom": 100.0,
                },
                {
                    "name": "A_E", "transmission_line_id": 3, "bus0": "230_A",
                    "bus1": "230_E", "type": "230kV", "v_nom": 230,
                    "length": 5.0, "s_nom": 100.0,
                },
            ]),
        ],
        ignore_index=True,
    )
    lines.to_csv(network_dir / "lines.csv", index=False)
    policy = {
        "main_component": {"minimum_bus_fraction": 0.5},
        "expected_islands": {"test_island": {"buses": ["230_C", "230_D"]}},
    }

    summary, issues = validate(network_dir, _config(), policy)

    assert summary["status"] == "PASS"
    approved = issues["disconnected_components.csv"].query(
        "classification == 'approved_island'"
    )
    assert approved["policy_id"].tolist() == ["test_island"]


def test_changed_approved_island_membership_fails(tmp_path: Path) -> None:
    network_dir = tmp_path / "network"
    _write_fixture(network_dir)
    policy = {
        "main_component": {"minimum_bus_fraction": 0.5},
        "expected_islands": {"test_island": {"buses": ["230_C", "230_D"]}},
    }

    summary, issues = validate(network_dir, _config(), policy)

    assert summary["status"] == "FAIL"
    assert "connected_component_policy" in summary["blocking_checks"]
    missing = issues["disconnected_components.csv"].query(
        "classification == 'missing_expected_island'"
    )
    assert missing["policy_id"].tolist() == ["test_island"]
