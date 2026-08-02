# Asset Validation Report

**Overall assessment:** Needs revision  
**Gate status:** FAIL  
**Generated (UTC):** 2026-08-02T18:56:37.151567+00:00

## Prepared records

- wind: 10
- solar: 3
- tpp: 38
- hydro_generation: 118
- hydro_reservoirs: 21
- hydro_cascade: 76

## Checks performed

| Check | Severity | Status | Findings | Denominator | Interpretation |
|---|---|---:|---:|---:|---|
| wind_required_columns | ERROR | PASS | 0 | 7 | Every required asset field must be present. |
| solar_required_columns | ERROR | PASS | 0 | 6 | Every required asset field must be present. |
| tpp_required_columns | ERROR | PASS | 0 | 6 | Every required asset field must be present. |
| hydro_generation_required_columns | ERROR | PASS | 0 | 10 | Every required asset field must be present. |
| hydro_reservoirs_required_columns | ERROR | PASS | 0 | 5 | Every required asset field must be present. |
| hydro_cascade_required_columns | ERROR | PASS | 0 | 7 | Every required asset field must be present. |
| asset_identifier_integrity | ERROR | PASS | 0 | 266 | Identifiers must be populated and unique within each table. |
| asset_coordinates | ERROR | PASS | 0 | 266 | Coordinates must be numeric and inside the broad modelling envelope. |
| positive_asset_capacities | ERROR | PASS | 0 | 266 | Generation capacity and known reservoir storage must be positive. |
| capacity_evidence_warnings | WARNING | WARN | 2 | 21 | Missing reservoir storage remains visible pending an evidence policy. |
| technology_vocabulary | ERROR | PASS | 0 | 245 | Technology labels must use the declared Rule 2 vocabulary. |
| base_network_mapping | ERROR | FAIL | 1 | 245 | Every generating asset must resolve to at least one Rule 1 bus. |
| wind_turbine_model_coverage | ERROR | PASS | 0 | 10 | Every wind configuration must resolve through the turbine dictionary. |
| hydro_reservoir_references | ERROR | PASS | 0 | 118 | Internal reservoir references must resolve; named river sinks are external boundaries. |

## Blocking checks

- `base_network_mapping`

## Warning checks

- `capacity_evidence_warnings`

## Evidence artifacts

- `schema_issues.csv`
- `duplicate_identifiers.csv`
- `invalid_coordinates.csv`
- `invalid_capacities.csv`
- `invalid_technology_labels.csv`
- `unmapped_network_nodes.csv`
- `wind_turbine_coverage.csv`
- `hydro_reference_issues.csv`

## Required caveats

- The coordinate check uses a broad bounding box, not an administrative polygon.
- A bus mapping proves identifier compatibility, not electrical feasibility.
- Capacity checks establish structural plausibility, not authoritative nameplate values.
- External river sinks are accepted boundaries and are not reservoir assets.
