# Assumption Report

Running log of modelling assumptions for PyPSA-BC. One row per
assumption, appended automatically when it is plugged into the model.
Do not edit rows by hand — append via `pypsa_bc.assumptions.log_assumption`.

| Date | Scenario | Module | Parameter | Value | Unit | Rationale | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-07-23 | Base_CNZ_2035 | base_network | line s_nom source | ttc_summer | MW | CODERS TTC used as thermal proxy (pf=1) | transmission_lines.csv 2025OCT28 |
| 2026-07-23 | Base_CNZ_2035 | base_network | reservoir e_initial | 0 | MWh | Cold-start |  |
| 2026-07-23 | Base | base_network | line s_nom | summer_rating_mw (ttc_summer) | MW | CODERS Total Transfer Capability (TTC) used as thermal proxy, converted to MVA using default power factor (0.9) from config.yaml |  |
| 2026-07-23 | Test | prepare_lines | line s_nom | summer_rating_mw (ttc_summer) | MW | CODERS Total Transfer Capability (TTC) used as thermal proxy, converted to MVA using default power factor (0.9) from config.yaml |  |
| 2026-07-23 | Test | lines | line s_nom | summer_rating_mw (ttc_summer) | MW | CODERS Total Transfer Capability (TTC) used as thermal proxy, converted to MVA using default power factor (0.9) from config.yaml |  |
| 2026-07-23 | Test | buses | AB-BC intertie handling | terminated on AB boundary bus; no import path |  | BC-internal scope; interties out of scope for this run | lines end node BC_ABBC03_IPT |
| 2026-07-23 | Test | transformers | transformer build — skipped nameless buses | 1 bus row(s) without a valid name skipped | - | Bus rows with NaN name (e.g. intertie boundary nodes) cannot form a transformer key | create_transformer_df |
| 2026-07-23 | Test | transformers | transformer impedance (type params) | vsc=10%, vscr=0.3%, i0=0.04%, pfe=30kW | mixed | Default nameplate values — refine per unit | base network prep |
| 2026-08-01 | Test | buses | registered representative bus coordinates | NODE-001, NODE-002, NODE-003, NODE-004, NODE-005, NODE-006, NODE-007, NODE-008, NODE-009 | correction IDs | Special non-automated treatments from the Rule 1 correction register; map-derived points are representative rather than surveyed locations | data/validation/base_network/node_corrections.csv |
| 2026-08-01 | Test | transformers | transformer impedance (type parameters) | vsc=10%, vscr=0.3%, i0=0.04%, pfe=30 kW | mixed | Standardized default nameplate values; refine with unit-specific data | base network preparation |
