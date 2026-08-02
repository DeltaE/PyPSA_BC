"""Validate prepared PyPSA-BC asset tables and retain QA evidence."""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

ASSET_PATHS = {
    "wind": Path("data/processed_data/wind/existing/bc_ext_wind_assets.csv"),
    "solar": Path("data/processed_data/solar/existing/bc_ext_solar_assets.csv"),
    "tpp": Path("data/processed_data/tpp/existing/bc_ext_tpp_assets.csv"),
    "hydro_generation": Path("data/processed_data/hydro/existing/hydro_generation.csv"),
    "hydro_reservoirs": Path("data/processed_data/hydro/existing/hydro_reservoirs.csv"),
    "hydro_cascade": Path("data/processed_data/hydro/existing/hydro_cascade.csv"),
}
BUS_PATH = Path("data/processed_data/network/buses.csv")
TURBINE_DICT_PATH = Path("data/downloaded_data/wind/turbine_dict.json")

TABLE_SCHEMAS = {
    "wind": ["asset_id", "connecting_node_code", "generator_capacity", "gen_type",
             "latitude", "longitude", "config_oedb"],
    "solar": ["asset_id", "connecting_node_code", "generator_capacity", "gen_type",
              "latitude", "longitude"],
    "tpp": ["generation_unit_code", "connecting_node_code", "generator_capacity",
            "gen_type", "latitude", "longitude"],
    "hydro_generation": ["asset_id", "component_id", "connecting_node_code", "capacity",
                         "gen_type", "latitude", "longitude", "upper_reservoir_id",
                         "lower_reservoir_id", "spill_lower_reservoir_id"],
    "hydro_reservoirs": ["asset_id", "latitude", "longitude", "min_storage", "max_storage"],
    "hydro_cascade": ["project_name", "gen_node_code", "connecting_node_code", "gen_type",
                      "install_capacity_in_mw", "latitude", "longitude"],
}
IDENTIFIER_FIELDS = {
    "wind": "asset_id", "solar": "asset_id", "tpp": "generation_unit_code",
    "hydro_generation": "asset_id", "hydro_reservoirs": "asset_id",
    "hydro_cascade": "project_name",
}
CAPACITY_FIELDS = {
    "wind": "generator_capacity", "solar": "generator_capacity",
    "tpp": "generator_capacity", "hydro_generation": "capacity",
    "hydro_reservoirs": "max_storage", "hydro_cascade": "install_capacity_in_mw",
}
ISSUE_COLUMNS = {
    "schema_issues.csv": ["table", "field", "severity", "issue"],
    "duplicate_identifiers.csv": ["table", "field", "value", "affected_rows", "severity"],
    "invalid_coordinates.csv": ["table", "row", "identifier", "latitude", "longitude",
                                "severity", "issue"],
    "invalid_capacities.csv": ["table", "row", "identifier", "field", "value",
                               "severity", "issue"],
    "invalid_technology_labels.csv": ["table", "row", "identifier", "value", "severity"],
    "unmapped_network_nodes.csv": ["table", "row", "identifier", "connecting_node_code",
                                   "severity", "issue"],
    "wind_turbine_coverage.csv": ["row", "asset_id", "config_oedb", "severity", "issue"],
    "hydro_reference_issues.csv": ["row", "asset_id", "field", "value", "severity", "issue"],
}


def _frame(rows: list[dict[str, Any]], filename: str) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[filename])


def _status(severity: str, count: int) -> str:
    if count == 0:
        return "PASS"
    return "FAIL" if severity == "ERROR" else "WARN" if severity == "WARNING" else "INFO"


def validate(
    tables: dict[str, pd.DataFrame],
    buses: pd.DataFrame,
    turbine_dict: dict[str, str],
    config: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, pd.DataFrame]]:
    """Validate in-memory tables; kept separate from I/O for focused tests."""
    checks: list[dict[str, Any]] = []
    rows = {name: [] for name in ISSUE_COLUMNS}

    def add_check(name: str, severity: str, count: int, denominator: int, message: str) -> None:
        checks.append({
            "check": name, "severity": severity, "status": _status(severity, count),
            "count": int(count), "denominator": int(denominator), "message": message,
        })

    usable: dict[str, bool] = {}
    for table_name, required in TABLE_SCHEMAS.items():
        table = tables.get(table_name, pd.DataFrame())
        missing = sorted(set(required) - set(table.columns))
        usable[table_name] = not missing
        rows["schema_issues.csv"].extend(
            {"table": table_name, "field": field, "severity": "ERROR",
             "issue": "missing_required_column"}
            for field in missing
        )
        add_check(f"{table_name}_required_columns", "ERROR", len(missing), len(required),
                  "Every required asset field must be present.")

    duplicate_count = 0
    for table_name, field in IDENTIFIER_FIELDS.items():
        if not usable[table_name]:
            continue
        values = tables[table_name][field].astype("string").str.strip()
        invalid = values.isna() | values.eq("")
        repeated = values[~invalid].value_counts()
        for value, count in repeated[repeated > 1].items():
            rows["duplicate_identifiers.csv"].append({
                "table": table_name, "field": field, "value": value,
                "affected_rows": int(count), "severity": "ERROR",
            })
            duplicate_count += int(count)
        for row in tables[table_name].index[invalid]:
            rows["schema_issues.csv"].append({
                "table": table_name, "field": f"{field}[{row}]", "severity": "ERROR",
                "issue": "missing_identifier",
            })
            duplicate_count += 1
    add_check("asset_identifier_integrity", "ERROR", duplicate_count,
              sum(len(table) for table in tables.values()),
              "Identifiers must be populated and unique within each table.")

    geo = config["geography"]
    for table_name, table in tables.items():
        if not usable.get(table_name):
            continue
        identifier_field = IDENTIFIER_FIELDS[table_name]
        lat = pd.to_numeric(table["latitude"], errors="coerce")
        lon = pd.to_numeric(table["longitude"], errors="coerce")
        invalid = (~np.isfinite(lat) | ~np.isfinite(lon)
                   | ~lat.between(geo["latitude_min"], geo["latitude_max"])
                   | ~lon.between(geo["longitude_min"], geo["longitude_max"]))
        for row in table.index[invalid]:
            rows["invalid_coordinates.csv"].append({
                "table": table_name, "row": int(row),
                "identifier": table.at[row, identifier_field],
                "latitude": table.at[row, "latitude"], "longitude": table.at[row, "longitude"],
                "severity": "ERROR", "issue": "missing_non_numeric_or_outside_bc_envelope",
            })
    add_check("asset_coordinates", "ERROR", len(rows["invalid_coordinates.csv"]),
              sum(len(t) for t in tables.values()),
              "Coordinates must be numeric and inside the broad modelling envelope.")

    for table_name, field in CAPACITY_FIELDS.items():
        if not usable[table_name]:
            continue
        table = tables[table_name]
        values = pd.to_numeric(table[field], errors="coerce")
        for row in table.index[~np.isfinite(values) | (values <= 0)]:
            missing_reservoir = table_name == "hydro_reservoirs" and not np.isfinite(values.loc[row])
            severity = (config["capacity"]["missing_reservoir_max_storage_severity"]
                        if missing_reservoir else config["capacity"]["nonpositive_severity"])
            rows["invalid_capacities.csv"].append({
                "table": table_name, "row": int(row),
                "identifier": table.at[row, IDENTIFIER_FIELDS[table_name]], "field": field,
                "value": table.at[row, field], "severity": severity.upper(),
                "issue": "missing" if not np.isfinite(values.loc[row]) else "nonpositive",
            })
    capacity_errors = sum(r["severity"] == "ERROR" for r in rows["invalid_capacities.csv"])
    capacity_warnings = sum(r["severity"] == "WARNING" for r in rows["invalid_capacities.csv"])
    add_check("positive_asset_capacities", "ERROR", capacity_errors,
              sum(len(tables[n]) for n in CAPACITY_FIELDS),
              "Generation capacity and known reservoir storage must be positive.")
    add_check("capacity_evidence_warnings", "WARNING", capacity_warnings,
              len(tables["hydro_reservoirs"]),
              "Missing reservoir storage remains visible pending an evidence policy.")

    for table_name, allowed in config["technologies"].items():
        if not usable.get(table_name):
            continue
        table = tables[table_name]
        invalid = ~table["gen_type"].astype(str).isin(set(allowed))
        for row in table.index[invalid]:
            rows["invalid_technology_labels.csv"].append({
                "table": table_name, "row": int(row),
                "identifier": table.at[row, IDENTIFIER_FIELDS[table_name]],
                "value": table.at[row, "gen_type"], "severity": "ERROR",
            })
    add_check("technology_vocabulary", "ERROR", len(rows["invalid_technology_labels.csv"]),
              sum(len(tables[n]) for n in config["technologies"]),
              "Technology labels must use the declared Rule 2 vocabulary.")

    bus_suffixes: set[str] = set()
    if "name" in buses:
        bus_suffixes = {
            name.split("_", 1)[1] for name in buses["name"].dropna().astype(str)
            if "_" in name
        }
    mapping_severity = config["network_mapping"]["severity"].upper()
    for table_name in ("wind", "solar", "tpp", "hydro_generation", "hydro_cascade"):
        if not usable[table_name]:
            continue
        table = tables[table_name]
        for row, code in table["connecting_node_code"].items():
            suffix = str(code).removeprefix("BC_")
            if not suffix or suffix not in bus_suffixes:
                rows["unmapped_network_nodes.csv"].append({
                    "table": table_name, "row": int(row),
                    "identifier": table.at[row, IDENTIFIER_FIELDS[table_name]],
                    "connecting_node_code": code, "severity": mapping_severity,
                    "issue": "no_voltage_specific_bus_for_source_node",
                })
    add_check("base_network_mapping", mapping_severity,
              len(rows["unmapped_network_nodes.csv"]),
              sum(len(tables[n]) for n in ("wind", "solar", "tpp", "hydro_generation", "hydro_cascade")),
              "Every generating asset must resolve to at least one Rule 1 bus.")

    known_turbines = {str(value) for value in turbine_dict.values()}
    if usable["wind"]:
        for row, value in tables["wind"]["config_oedb"].items():
            base_value = str(value).split("*", 1)[0]
            if base_value not in known_turbines:
                rows["wind_turbine_coverage.csv"].append({
                    "row": int(row), "asset_id": tables["wind"].at[row, "asset_id"],
                    "config_oedb": value, "severity": "ERROR",
                    "issue": "configuration_not_resolved_by_turbine_dictionary",
                })
    add_check("wind_turbine_model_coverage", "ERROR",
              len(rows["wind_turbine_coverage.csv"]), len(tables["wind"]),
              "Every wind configuration must resolve through the turbine dictionary.")

    reservoirs = (set(tables["hydro_reservoirs"]["asset_id"].dropna().astype(str))
                  if usable["hydro_reservoirs"] else set())
    hydro_cfg = config["hydro"]
    if usable["hydro_generation"]:
        hydro = tables["hydro_generation"]
        for field in ("upper_reservoir_id", "lower_reservoir_id", "spill_lower_reservoir_id"):
            for row, value in hydro[field].items():
                value = "" if pd.isna(value) else str(value).strip()
                if (not value or value == hydro_cfg["default_reference"]
                        or value.endswith(hydro_cfg["external_sink_suffix"])):
                    continue
                if value not in reservoirs:
                    rows["hydro_reference_issues.csv"].append({
                        "row": int(row), "asset_id": hydro.at[row, "asset_id"], "field": field,
                        "value": value, "severity": "ERROR", "issue": "unknown_reservoir_reference",
                    })
    add_check("hydro_reservoir_references", "ERROR",
              len(rows["hydro_reference_issues.csv"]), len(tables["hydro_generation"]),
              "Internal reservoir references must resolve; named river sinks are external boundaries.")

    issues = {name: _frame(content, name) for name, content in rows.items()}
    status = "FAIL" if any(c["severity"] == "ERROR" and c["count"] for c in checks) else "PASS"
    summary = {
        "validator": "validate_assets", "version": config.get("version", 1),
        "generated_at_utc": datetime.now(UTC).isoformat(), "status": status,
        "confidence": "Needs revision" if status == "FAIL" else "Ready to share",
        "grain": {
            "wind_solar_tpp": "one prepared generating unit per row",
            "hydro_generation": "one prepared hydro generating asset per row",
            "hydro_reservoirs": "one prepared reservoir per row",
            "hydro_cascade": "one source hydro project-unit record per row",
        },
        "counts": {name: len(table) for name, table in tables.items()},
        "checks": checks,
        "blocking_checks": [c["check"] for c in checks if c["severity"] == "ERROR" and c["count"]],
        "warning_checks": [c["check"] for c in checks if c["severity"] == "WARNING" and c["count"]],
        "artifacts": list(ISSUE_COLUMNS),
        "limitations": [
            "The coordinate check uses a broad bounding box, not an administrative polygon.",
            "A bus mapping proves identifier compatibility, not electrical feasibility.",
            "Capacity checks establish structural plausibility, not authoritative nameplate values.",
            "External river sinks are accepted boundaries and are not reservoir assets.",
        ],
    }
    return summary, issues


def _write_report(summary: dict[str, Any], output_dir: Path) -> None:
    checks = "\n".join(
        f"| {c['check']} | {c['severity']} | {c['status']} | {c['count']} | "
        f"{c['denominator']} | {c['message']} |" for c in summary["checks"]
    )
    report = f"""# Asset Validation Report

**Overall assessment:** {summary['confidence']}  
**Gate status:** {summary['status']}  
**Generated (UTC):** {summary['generated_at_utc']}

## Prepared records

{chr(10).join(f'- {name}: {count}' for name, count in summary['counts'].items())}

## Checks performed

| Check | Severity | Status | Findings | Denominator | Interpretation |
|---|---|---:|---:|---:|---|
{checks}

## Blocking checks

{chr(10).join(f'- `{name}`' for name in summary['blocking_checks']) or '- None'}

## Warning checks

{chr(10).join(f'- `{name}`' for name in summary['warning_checks']) or '- None'}

## Evidence artifacts

{chr(10).join(f'- `{name}`' for name in summary['artifacts'])}

## Required caveats

{chr(10).join(f'- {item}' for item in summary['limitations'])}
"""
    (output_dir / "validation_report.md").write_text(report, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    tables = {name: pd.read_csv(path) for name, path in ASSET_PATHS.items()}
    buses = pd.read_csv(BUS_PATH)
    turbine_dict = json.loads(TURBINE_DICT_PATH.read_text(encoding="utf-8"))
    if not isinstance(turbine_dict, dict):
        raise TypeError("turbine_dict.json must contain a JSON object")
    summary, issues = validate(tables, buses, turbine_dict, config)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, frame in issues.items():
        frame.to_csv(args.output_dir / filename, index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8"
    )
    _write_report(summary, args.output_dir)
    print(f"Asset validation: {summary['status']}")
    print(f"Report: {args.output_dir / 'validation_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
