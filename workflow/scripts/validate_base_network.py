"""Validate prepared PyPSA-BC base-network tables and retain QA evidence."""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
import yaml

TABLE_SCHEMAS = {
    "buses": ["name", "x", "y", "type", "v_nom"],
    "lines": ["name", "bus0", "bus1", "type", "v_nom", "length", "s_nom"],
    "line_types": ["code"],
    "transformers": ["name", "bus0", "bus1", "type"],
    "transformer_types": ["name", "s_nom", "v_nom_0", "v_nom_1"],
}

ISSUE_COLUMNS = {
    "invalid_buses.csv": ["row", "name", "x", "y", "v_nom", "severity", "issues"],
    "missing_line_endpoints.csv": [
        "row", "name", "transmission_line_id", "bus0", "bus1",
        "missing_bus0", "missing_bus1", "severity",
    ],
    "self_loops.csv": ["row", "name", "transmission_line_id", "bus0", "bus1", "severity"],
    "disconnected_components.csv": [
        "component_id", "bus_count", "is_largest", "severity", "bus_names",
    ],
    "parameter_outliers.csv": [
        "component", "row", "identifier", "field", "value", "severity", "issue",
    ],
    "duplicate_identifiers.csv": [
        "component", "field", "value", "affected_rows", "severity", "issue",
    ],
    "parallel_lines.csv": [
        "bus_a", "bus_b", "v_nom", "circuit_count", "line_names",
        "transmission_line_ids",
    ],
}


def _severity_status(severity: str, count: int) -> str:
    if count == 0:
        return "PASS"
    if severity == "ERROR":
        return "FAIL"
    return "WARN" if severity == "WARNING" else "INFO"


def _read_tables(network_dir: Path) -> tuple[dict[str, pd.DataFrame], list[dict[str, Any]]]:
    tables: dict[str, pd.DataFrame] = {}
    checks: list[dict[str, Any]] = []
    for name, required_columns in TABLE_SCHEMAS.items():
        path = network_dir / f"{name}.csv"
        if not path.exists():
            tables[name] = pd.DataFrame(columns=required_columns)
            checks.append({
                "check": f"{name}_file_exists",
                "severity": "ERROR",
                "status": "FAIL",
                "count": 1,
                "denominator": 1,
                "message": f"Required file is missing: {path}",
            })
            continue
        table = pd.read_csv(path)
        tables[name] = table
        missing = sorted(set(required_columns) - set(table.columns))
        checks.append({
            "check": f"{name}_required_columns",
            "severity": "ERROR",
            "status": "FAIL" if missing else "PASS",
            "count": len(missing),
            "denominator": len(required_columns),
            "message": "Missing columns: " + ", ".join(missing) if missing else "Required columns present.",
        })
    return tables, checks


def _issue_frame(rows: list[dict[str, Any]], filename: str) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=ISSUE_COLUMNS[filename])


def validate(network_dir: Path, config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, pd.DataFrame]]:
    tables, checks = _read_tables(network_dir)
    buses = tables["buses"]
    lines = tables["lines"]
    line_types = tables["line_types"]
    transformers = tables["transformers"]
    transformer_types = tables["transformer_types"]

    issues = {name: _issue_frame([], name) for name in ISSUE_COLUMNS}

    def add_check(name: str, severity: str, count: int, denominator: int, message: str) -> None:
        checks.append({
            "check": name,
            "severity": severity,
            "status": _severity_status(severity, count),
            "count": int(count),
            "denominator": int(denominator),
            "rate": round(float(count / denominator), 8) if denominator else None,
            "message": message,
        })

    required_available = all(
        set(required).issubset(tables[name].columns)
        for name, required in TABLE_SCHEMAS.items()
    )
    if required_available:
        geo = config["geography"]
        invalid_bus_rows: list[dict[str, Any]] = []
        numeric_x = pd.to_numeric(buses["x"], errors="coerce")
        numeric_y = pd.to_numeric(buses["y"], errors="coerce")
        numeric_voltage = pd.to_numeric(buses["v_nom"], errors="coerce")
        names = buses["name"].astype("string")
        duplicate_names = names.notna() & names.duplicated(keep=False)
        for row in buses.index:
            reasons: list[str] = []
            if pd.isna(names.loc[row]) or not str(names.loc[row]).strip():
                reasons.append("missing_bus_name")
            if bool(duplicate_names.loc[row]):
                reasons.append("duplicate_bus_name")
            if not np.isfinite(numeric_x.loc[row]) or not np.isfinite(numeric_y.loc[row]):
                reasons.append("missing_or_non_numeric_coordinate")
            elif not (
                geo["longitude_min"] <= numeric_x.loc[row] <= geo["longitude_max"]
                and geo["latitude_min"] <= numeric_y.loc[row] <= geo["latitude_max"]
            ):
                reasons.append("outside_broad_bc_bounds")
            if not np.isfinite(numeric_voltage.loc[row]) or numeric_voltage.loc[row] <= 0:
                reasons.append("missing_or_nonpositive_voltage")
            if reasons:
                invalid_bus_rows.append({
                    "row": int(row), "name": buses.at[row, "name"],
                    "x": buses.at[row, "x"], "y": buses.at[row, "y"],
                    "v_nom": buses.at[row, "v_nom"], "severity": "ERROR",
                    "issues": ";".join(reasons),
                })
        issues["invalid_buses.csv"] = _issue_frame(invalid_bus_rows, "invalid_buses.csv")
        add_check("valid_buses", "ERROR", len(invalid_bus_rows), len(buses),
                  "Bus names, coordinates, and nominal voltages must be valid.")

        valid_bus_names = set(names.dropna().astype(str))
        missing_endpoint_rows: list[dict[str, Any]] = []
        self_loop_rows: list[dict[str, Any]] = []
        for row in lines.index:
            bus0, bus1 = str(lines.at[row, "bus0"]), str(lines.at[row, "bus1"])
            missing0, missing1 = bus0 not in valid_bus_names, bus1 not in valid_bus_names
            common = {
                "row": int(row), "name": lines.at[row, "name"],
                "transmission_line_id": lines.at[row, "transmission_line_id"]
                if "transmission_line_id" in lines else None,
                "bus0": bus0, "bus1": bus1,
            }
            if missing0 or missing1:
                missing_endpoint_rows.append({
                    **common, "missing_bus0": missing0, "missing_bus1": missing1,
                    "severity": "ERROR",
                })
            if bus0 == bus1:
                self_loop_rows.append({**common, "severity": "ERROR"})
        issues["missing_line_endpoints.csv"] = _issue_frame(
            missing_endpoint_rows, "missing_line_endpoints.csv"
        )
        issues["self_loops.csv"] = _issue_frame(self_loop_rows, "self_loops.csv")
        add_check("line_endpoint_integrity", "ERROR", len(missing_endpoint_rows), len(lines),
                  "Every line endpoint must reference a prepared bus.")
        add_check("line_self_loops", "ERROR", len(self_loop_rows), len(lines),
                  "Transmission lines may not connect a bus to itself.")

        duplicate_rows: list[dict[str, Any]] = []
        duplicate_specs = [("buses", buses, "name", "ERROR", "Bus names are primary keys.")]
        if "transmission_line_id" in lines:
            severity = config["identifiers"]["duplicate_transmission_line_id_severity"].upper()
            duplicate_specs.append((
                "lines", lines, "transmission_line_id", severity,
                "Source transmission-line IDs must be unique at prepared-line grain.",
            ))
        severity = config["identifiers"]["duplicate_line_name_severity"].upper()
        duplicate_specs.append((
            "lines", lines, "name", severity,
            "Repeated display names may be parallel circuits and require traceable source IDs.",
        ))
        for component, table, field, severity, message in duplicate_specs:
            counts = table[field].dropna().astype(str).value_counts()
            repeated = counts[counts > 1]
            for value, affected in repeated.items():
                duplicate_rows.append({
                    "component": component, "field": field, "value": value,
                    "affected_rows": int(affected), "severity": severity, "issue": message,
                })
            add_check(
                f"{component}_{field}_uniqueness", severity,
                int(repeated.sum()), len(table), message,
            )
        issues["duplicate_identifiers.csv"] = _issue_frame(
            duplicate_rows, "duplicate_identifiers.csv"
        )

        parallel_rows: list[dict[str, Any]] = []
        parallel = lines.copy()
        parallel["bus_a"] = parallel[["bus0", "bus1"]].astype(str).min(axis=1)
        parallel["bus_b"] = parallel[["bus0", "bus1"]].astype(str).max(axis=1)
        for (bus_a, bus_b, voltage), group in parallel.groupby(
            ["bus_a", "bus_b", "v_nom"], dropna=False
        ):
            if len(group) > 1:
                parallel_rows.append({
                    "bus_a": bus_a, "bus_b": bus_b, "v_nom": voltage,
                    "circuit_count": len(group),
                    "line_names": ";".join(group["name"].astype(str)),
                    "transmission_line_ids": ";".join(
                        group["transmission_line_id"].astype(str)
                    ) if "transmission_line_id" in group else "",
                })
        issues["parallel_lines.csv"] = _issue_frame(parallel_rows, "parallel_lines.csv")
        add_check("parallel_line_groups", "INFO", len(parallel_rows), len(lines),
                  "Parallel circuits are retained and reported; they are not automatically errors.")

        outlier_rows: list[dict[str, Any]] = []

        def record_numeric(component: str, table: pd.DataFrame, identifier_field: str,
                           field: str, error_nonpositive: bool = True,
                           warn_below: float | None = None,
                           warn_above: float | None = None) -> None:
            values = pd.to_numeric(table[field], errors="coerce")
            for row, value in values.items():
                severity, reason = None, None
                if not np.isfinite(value):
                    severity, reason = "ERROR", "missing_or_non_numeric"
                elif error_nonpositive and value <= 0:
                    severity, reason = "ERROR", "nonpositive"
                elif warn_below is not None and value < warn_below:
                    severity, reason = "WARNING", f"below_warning_threshold_{warn_below}"
                elif warn_above is not None and value > warn_above:
                    severity, reason = "WARNING", f"above_warning_threshold_{warn_above}"
                if severity:
                    outlier_rows.append({
                        "component": component, "row": int(row),
                        "identifier": table.at[row, identifier_field], "field": field,
                        "value": table.at[row, field], "severity": severity, "issue": reason,
                    })

        length_cfg = config["parameters"]["line_length_km"]
        capacity_cfg = config["parameters"]["line_capacity_mva"]
        record_numeric("lines", lines, "name", "length", True,
                       length_cfg["warn_below"], length_cfg["warn_above"])
        record_numeric("lines", lines, "name", "s_nom", True, None,
                       capacity_cfg["warn_above"])
        record_numeric("buses", buses, "name", "v_nom")
        if not transformer_types.empty:
            record_numeric("transformer_types", transformer_types, "name", "s_nom")

        # lines.type is a voltage label while line_types.csv preserves the
        # conductor inventory. The builder later supplies explicit r/x/b/g and
        # clears lines.type, so these are not a foreign-key relationship.
        expected_voltage_label = (
            pd.to_numeric(lines["v_nom"], errors="coerce")
            .map(lambda value: f"{int(value)}kV" if np.isfinite(value) else "")
        )
        for row in lines.index[lines["type"].astype(str) != expected_voltage_label]:
            outlier_rows.append({
                "component": "lines", "row": int(row), "identifier": lines.at[row, "name"],
                "field": "type", "value": lines.at[row, "type"], "severity": "ERROR",
                "issue": "line_voltage_label_mismatch",
            })
        duplicate_conductor_codes = line_types["code"].dropna().astype(str).duplicated(keep=False)
        for row in line_types.index[duplicate_conductor_codes]:
            outlier_rows.append({
                "component": "line_types", "row": int(row),
                "identifier": line_types.at[row, "code"], "field": "code",
                "value": line_types.at[row, "code"], "severity": "ERROR",
                "issue": "duplicate_conductor_code",
            })
        known_transformer_types = set(transformer_types["name"].dropna().astype(str))
        for row in transformers.index[
            ~transformers["type"].astype(str).isin(known_transformer_types)
        ]:
            outlier_rows.append({
                "component": "transformers", "row": int(row),
                "identifier": transformers.at[row, "name"], "field": "type",
                "value": transformers.at[row, "type"], "severity": "ERROR",
                "issue": "unknown_transformer_type",
            })
        for row in transformers.index:
            for field in ("bus0", "bus1"):
                if str(transformers.at[row, field]) not in valid_bus_names:
                    outlier_rows.append({
                        "component": "transformers", "row": int(row),
                        "identifier": transformers.at[row, "name"], "field": field,
                        "value": transformers.at[row, field], "severity": "ERROR",
                        "issue": "missing_bus_reference",
                    })
        issues["parameter_outliers.csv"] = _issue_frame(
            outlier_rows, "parameter_outliers.csv"
        )
        error_outliers = sum(row["severity"] == "ERROR" for row in outlier_rows)
        warning_outliers = sum(row["severity"] == "WARNING" for row in outlier_rows)
        add_check("electrical_parameter_validity", "ERROR", error_outliers,
                  len(lines) + len(buses) + len(transformer_types),
                  "Required electrical values and type references must be valid.")
        add_check("electrical_parameter_warnings", "WARNING", warning_outliers,
                  len(lines) + len(buses) + len(transformer_types),
                  "Extreme but possible values require review.")

        graph = nx.Graph()
        graph.add_nodes_from(valid_bus_names)
        for table in (lines, transformers):
            for bus0, bus1 in table[["bus0", "bus1"]].dropna().astype(str).itertuples(
                index=False, name=None
            ):
                if bus0 in valid_bus_names and bus1 in valid_bus_names:
                    graph.add_edge(bus0, bus1)
        components = sorted(nx.connected_components(graph), key=len, reverse=True)
        component_rows: list[dict[str, Any]] = []
        for component_id, component in enumerate(components, start=1):
            is_largest = component_id == 1
            component_rows.append({
                "component_id": component_id, "bus_count": len(component),
                "is_largest": is_largest,
                "severity": "INFO" if is_largest else "ERROR",
                "bus_names": ";".join(sorted(component)),
            })
        issues["disconnected_components.csv"] = _issue_frame(
            component_rows, "disconnected_components.csv"
        )
        excess_components = max(0, len(components) - config["connectivity"]["max_components"])
        isolated_count = sum(len(component) == 1 for component in components)
        add_check("connected_components", "ERROR", excess_components, len(components),
                  "The prepared electrical graph must satisfy the configured component limit.")
        isolated_severity = "INFO" if config["connectivity"]["allow_isolated_buses"] else "ERROR"
        add_check("isolated_buses", isolated_severity, isolated_count, len(buses),
                  "Isolated buses cannot participate in network dispatch.")

    status = "FAIL" if any(
        check["severity"] == "ERROR" and check["count"] > 0 for check in checks
    ) else "PASS"
    source_files = [network_dir / f"{name}.csv" for name in TABLE_SCHEMAS]
    existing_mtimes = [path.stat().st_mtime for path in source_files if path.exists()]
    data_as_of = (
        datetime.fromtimestamp(max(existing_mtimes), UTC).isoformat()
        if existing_mtimes else None
    )
    summary = {
        "validator": "validate_base_network",
        "version": config.get("version", 1),
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "data_as_of_utc": data_as_of,
        "network_directory": str(network_dir),
        "status": status,
        "confidence": "Needs revision" if status == "FAIL" else "Ready to share",
        "grain": {
            "buses": "one voltage-specific electrical bus per prepared identifier",
            "lines": "one prepared transmission circuit per row",
            "transformers": "one voltage-level connection per row",
        },
        "counts": {name: len(table) for name, table in tables.items()},
        "checks": checks,
        "blocking_checks": [
            check["check"] for check in checks
            if check["severity"] == "ERROR" and check["count"] > 0
        ],
        "warning_checks": [
            check["check"] for check in checks
            if check["severity"] == "WARNING" and check["count"] > 0
        ],
        "artifacts": list(ISSUE_COLUMNS),
        "limitations": [
            "The coordinate test uses a broad bounding box, not an administrative polygon.",
            "Parallel circuits are reported but retained when source identifiers remain traceable.",
            "This gate validates prepared topology and parameters; it does not validate power-flow feasibility.",
        ],
    }
    return summary, issues


def _write_report(summary: dict[str, Any], output_dir: Path) -> None:
    rows = []
    for check in summary["checks"]:
        rows.append(
            f"| {check['check']} | {check['severity']} | {check['status']} | "
            f"{check['count']} | {check.get('denominator', '')} | {check['message']} |"
        )
    report = f"""# Base Network Validation Report

**Overall assessment:** {summary['confidence']}  
**Gate status:** {summary['status']}  
**Generated (UTC):** {summary['generated_at_utc']}  
**Source data as of (UTC):** {summary['data_as_of_utc']}

## Dataset and grain

- Buses: {summary['counts'].get('buses', 0)}
- Lines: {summary['counts'].get('lines', 0)}
- Line types: {summary['counts'].get('line_types', 0)}
- Transformers: {summary['counts'].get('transformers', 0)}
- Transformer types: {summary['counts'].get('transformer_types', 0)}

## Checks performed

| Check | Severity | Status | Findings | Denominator | Interpretation |
|---|---|---:|---:|---:|---|
{chr(10).join(rows)}

## Blocking checks

{chr(10).join(f'- `{name}`' for name in summary['blocking_checks']) or '- None'}

## Warning checks

{chr(10).join(f'- `{name}`' for name in summary['warning_checks']) or '- None'}

## Evidence artifacts

{chr(10).join(f'- `{name}`' for name in summary['artifacts'])}

## Required caveats

{chr(10).join(f'- {text}' for text in summary['limitations'])}
"""
    (output_dir / "validation_report.md").write_text(report, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    summary, issues = validate(args.network_dir, config)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, frame in issues.items():
        frame.to_csv(args.output_dir / filename, index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8"
    )
    _write_report(summary, args.output_dir)
    print(f"Base network validation: {summary['status']}")
    print(f"Report: {args.output_dir / 'validation_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
