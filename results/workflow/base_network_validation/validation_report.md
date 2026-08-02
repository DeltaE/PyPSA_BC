# Base Network Validation Report

**Overall assessment:** Ready to share  
**Gate status:** PASS  
**Generated (UTC):** 2026-08-02T03:45:40.234747+00:00  
**Source data as of (UTC):** 2026-08-02T03:45:26.090579+00:00

## Dataset and grain

- Buses: 1018
- Lines: 1237
- Line types: 34
- Transformers: 110
- Transformer types: 18

## Checks performed

| Check | Severity | Status | Findings | Denominator | Interpretation |
|---|---|---:|---:|---:|---|
| buses_required_columns | ERROR | PASS | 0 | 5 | Required columns present. |
| lines_required_columns | ERROR | PASS | 0 | 7 | Required columns present. |
| line_types_required_columns | ERROR | PASS | 0 | 1 | Required columns present. |
| transformers_required_columns | ERROR | PASS | 0 | 4 | Required columns present. |
| transformer_types_required_columns | ERROR | PASS | 0 | 4 | Required columns present. |
| valid_buses | ERROR | PASS | 0 | 1018 | Bus names, coordinates, and nominal voltages must be valid. |
| line_endpoint_integrity | ERROR | PASS | 0 | 1237 | Every line endpoint must reference a prepared bus. |
| line_self_loops | ERROR | PASS | 0 | 1237 | Transmission lines may not connect a bus to itself. |
| buses_name_uniqueness | ERROR | PASS | 0 | 1018 | Bus names are primary keys. |
| lines_transmission_line_id_uniqueness | ERROR | PASS | 0 | 1237 | Source transmission-line IDs must be unique at prepared-line grain. |
| lines_name_uniqueness | WARNING | WARN | 348 | 1237 | Repeated display names may be parallel circuits and require traceable source IDs. |
| parallel_line_groups | INFO | INFO | 165 | 1237 | Parallel circuits are retained and reported; they are not automatically errors. |
| electrical_parameter_validity | ERROR | PASS | 0 | 2273 | Required electrical values and type references must be valid. |
| electrical_parameter_warnings | WARNING | WARN | 56 | 2273 | Extreme but possible values require review. |
| connected_component_policy | ERROR | PASS | 0 | 3 | The main component must meet its minimum share and every non-main component must exactly match an approved island membership. |
| isolated_buses | ERROR | PASS | 0 | 1018 | Isolated buses cannot participate in network dispatch. |

## Blocking checks

- None

## Warning checks

- `lines_name_uniqueness`
- `electrical_parameter_warnings`

## Evidence artifacts

- `invalid_buses.csv`
- `missing_line_endpoints.csv`
- `self_loops.csv`
- `disconnected_components.csv`
- `parameter_outliers.csv`
- `duplicate_identifiers.csv`
- `parallel_lines.csv`

## Required caveats

- The coordinate test uses a broad bounding box, not an administrative polygon.
- Registered map-derived coordinates are representative schematic points, not surveyed assets.
- Parallel circuits are reported but retained when source identifiers remain traceable.
- This gate validates prepared topology and parameters; it does not validate power-flow feasibility.
