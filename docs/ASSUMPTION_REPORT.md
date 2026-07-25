# Assumption Report

Running log of modelling assumptions for PyPSA-BC. One row per
assumption, appended automatically when it is plugged into the model.
Do not edit rows by hand — append via `pypsa_bc.assumptions.log_assumption`.

| Date | Scenario | Module | Parameter | Value | Unit | Rationale | Source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-07-23 | Base_CNZ_2035 | base_network | line s_nom source | ttc_summer | MW | CODERS Total Transfer Capability used as per-line thermal proxy; power factor assumed 1 (MVA≈MW) | transmission_lines.csv (CODERS 2025OCT28) |
| 2026-07-23 | Base_CNZ_2035 | base_network | zero-reactance taps floored | x = 1e-4 | ohm | 10 sub-13 m taps have x=0; floored to avoid singular susceptance | dataQA 2026-07-23 |
| 2026-07-23 | Base_CNZ_2035 | hydro | reservoir e_initial | 0 | MWh | Cold-start; WSC Jan-1-2021 storage estimates not yet integrated | NetworkReport |
| 2026-07-23 | - | hydro | hydro head (storage->energy) | dam height | m | Dam height proxy where rated head unavailable (GRanD adapter) | GRanD/GDW |
