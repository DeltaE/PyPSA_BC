# Assumption Report

Running log of modelling assumptions for PyPSA-BC. One row per
assumption, appended automatically when it is plugged into the model.
Do not edit rows by hand — append via `pypsa_bc.assumptions.log_assumption`.

| Date | Scenario | Module | Parameter | Value | Unit | Rationale | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-07-23 | Base_CNZ_2035 | base_network | line s_nom source | ttc_summer | MW | CODERS TTC used as thermal proxy (pf=1) | transmission_lines.csv 2025OCT28 |
