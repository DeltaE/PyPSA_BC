# PyPSA-BC model change and evidence report — 2026-08-01

## Outcome

The calendar-2021 full-horizon reference now solves optimally after three
evidence-driven corrections: model-year asset filtering, fixed source profiles
for run-of-river generation, and monthly Statistics Canada calibration targets
for biomass and non-renewable combustible generation. Seven named hydro-station
capacities are aligned to BC Hydro's fiscal-2021 table.

The corrected model passes all eight computational-integrity checks and the
declared calibration-fit thresholds. This does **not** constitute independent
operational validation: the load chronology and two non-hydro categories are
calibration inputs, while the official station-energy comparison covers a
different fiscal period and remains a screening test.

The hydro-dispatch investigation now identifies the mechanism behind the Seven
Mile and Bridge River screen failures. Their effective CAD 6.88/MWh turbine cost
exceeds the maximum solved local marginal price of approximately CAD 1.968/MWh
in every hour, so the source-cost formulation prefers free spill paths. A solved
±25% fiscal-energy-band sensitivity restores material turbine routing at both
stations, preserves 62.240 TWh of total hydro generation, and raises the model
objective by 4.60%. This is a cross-period calibration sensitivity, not
independent validation.

TTC is excluded from the core reference, calibration, figures, and headline
claims. Observed hourly interchange remains the external boundary. The
`hourly_ttc` implementation is retained only as an optional secondary
sensitivity.

## Implemented model changes

| Change | Implementation | Scientific reason |
|---|---|---|
| 2021 asset-vintage filter | Solar, wind, and thermal records are filtered using commissioning and closure years during model assembly. | Prevents future and retired assets from contaminating the 2021 reference. |
| Run-of-river source policy | For 98 RoR generators, `p_min_pu` equals the source-backed `p_max_pu` profile under `fixed_source_profile`. | Prevents cost-only curtailment from reallocating RoR energy arbitrarily among hydro assets. |
| Monthly non-hydro calibration | Twenty-four monthly equality constraints reproduce Statistics Canada biomass and non-renewable combustible totals. | Represents observed must-run/contract/cogeneration effects missing from a cost-only dispatch, while retaining hourly and plant-level optimization. |
| Named-station capacities | G.M. Shrum, Revelstoke, Mica, Kootenay Canal, Peace Canyon, Seven Mile, and Bridge River are aligned to official fiscal-2021 capacity evidence. | Removes known capacity inconsistencies before comparing station energy. |
| Evidence-role metadata | Network metadata and audit tables distinguish calibration inputs, independent screens, and model diagnostics. | Prevents calibrated agreement from being reported as independent validation. |
| TTC isolation | Core preflight and observed-interchange execution no longer depend on TTC data. | Keeps the study focused on cascade representation and avoids an unvalidated boundary assumption. |
| Explicit hydro dispatch-cost policy | `source_generic_vom` and `uniform_reservoir_vom` are named, audited policies rather than hidden cost edits. | Makes turbine-versus-spill economic ordering reproducible and testable. |
| Named-station annual-energy sensitivity | Six reservoir-link stations can receive explicit annual energy bands; Kootenay Canal remains outside because it is a fixed-profile generator. | Tests whether the hydraulic topology can route water through screened stations without presenting constrained output as validation. |
| Optimal-only solve metadata | Solve summaries now retain UTC timestamps, elapsed time, objective, solver method, and failure details when the optimizer returns. | Strengthens the execution audit trail for future scenario runs. |

The prepared network contains 80 buses, 156 generators, 151 links, 33 stores,
25 loads, and 8,760 hourly snapshots. The build writes separate audit tables for
asset vintage, generation calibration, RoR dispatch, named-station capacity,
external boundary, storage, releases, and network formulation.

## Corrected annual reference

Source: `results/workflow/annual/A_lower_zero_observed_transport_full_horizon.nc`.

| Integrity metric | Result |
|---|---:|
| Solver | HiGHS, optimal |
| Demand | 64.433 TWh |
| Peak demand | 11.679 GW |
| Imports / exports | 0.718 / 6.869 TWh |
| Backstop generation | 0 GWh |
| Maximum hourly electrical-balance residual | 7.28 × 10^-12 MW |
| External-interchange error | 0 MW |
| Release violations | 0 hours |
| Reservoir-bound violations | 0 |
| Cyclic-storage closure residual | 0 |
| Model objective | 419.955 million model monetary units |

The objective is a model quantity, not a validated estimate of observed system
cost. Raw zero-cost transport flows are non-unique; congestion reporting uses
the implemented lexicographic minimum-transfer reconstruction.

## Calendar-2021 calibration fit

| Category | Statistics Canada | PyPSA-BC | Relative difference | Evidence role |
|---|---:|---:|---:|---|
| Total generation | 71.985 TWh | 70.585 TWh | -1.95% | Calibration-dataset residual alignment |
| Hydro | 63.802 TWh | 62.240 TWh | -2.45% | Calibration-dataset residual alignment |
| Wind | 1.999 TWh | 2.162 TWh | +8.19% | Calibration-dataset residual alignment |
| Biomass | 4.066 TWh | 4.066 TWh | 0.00% | Enforced monthly calibration target |
| Non-renewable combustible | 2.115 TWh | 2.115 TWh | 0.00% | Enforced monthly calibration target |
| Solar | 0.0043 TWh | 0.0018 TWh | -58.9% | Diagnostic; very small and boundary-sensitive |

Load energy and peak are each 0.95% below the BC Hydro series and the hourly
correlation is 1.0 because that series supplies the modeled chronology. These
are input-fidelity checks only.

## Named-station screen

All seven official capacity values match after calibration. Calendar-2021 model
energy is within the ±25% screening band of fiscal-2021 BC Hydro energy for G.M.
Shrum, Revelstoke, Mica, and Peace Canyon. Kootenay Canal is 34.8% lower. Seven
Mile and Bridge River are approximately 100% lower and require investigation.
The water-routing diagnostic shows that only 0.34% of Seven Mile's routed water
and none of Bridge River's routed water use the turbine path; the remainder uses
the modeled spill/release route. The economic-ordering audit explains this
result: the two stations' CAD 6.88/MWh turbine cost exceeds the solved local
marginal price in every hour. Kootenay Canal is a fixed-profile generator, so its
gap instead points to the source availability profile or period boundary.

This is not a like-for-like annual validation because the modeled and observed
periods differ. The new sensitivity shows that the existing topology can route
water through both turbine paths, so a missing hydraulic edge is not the primary
cause. It does not establish which dispatch policy is empirically correct. Do not
adopt fiscal-year station-energy bands as the primary calendar-year reference
without period-matched evidence.

## Hydro dispatch-policy sensitivity

The source-cost policy remains the unconstrained reference. A uniform CAD
1.97/MWh reservoir-hydro cost policy was prepared, but controlled simplex and
interior-point attempts produced no exported optimal annual solution. This is a
computational performance result and does not prove infeasibility.

The accepted diagnostic retained source-generic costs and imposed ±25% annual
energy bands around six BC Hydro fiscal comparators. It solved optimally and
passed all eight integrity checks. Seven Mile increased from 11.6 to 2,279.2 GWh
and Bridge River from 0 to 1,664.3 GWh. G.M. Shrum decreased from 14,686.3 to
12,608.6 GWh and Mica from 8,356.7 to 6,501.8 GWh. Total hydro generation and the
provincial generation mix were unchanged. Kootenay Canal remained at 1,713.2 GWh
and still failed the cross-period screen.

The result establishes dispatch responsiveness under an imposed calibration
band. It cannot validate station energy against the same evidence and does not
clear the independent-validation gate.

## Reproducibility map

| Purpose | Path |
|---|---|
| Generation calibration engine | `src/pypsa_bc/studies/generation_calibration.py` |
| Named-station alignment and comparison | `src/pypsa_bc/studies/station_validation.py` |
| Hydro dispatch-cost policy | `src/pypsa_bc/studies/hydro_dispatch_policy.py` |
| Model assembly integration | `workflow/scripts/build_model.py` |
| Solve callback integration | `src/pypsa_bc/studies/cascade/constraints.py` |
| Reproducible BC Hydro annual-report fetch | `workflow/scripts/fetch_bc_hydro_annual_validation.py` |
| Reproducible station-table extraction | `workflow/scripts/extract_bc_hydro_supply_validation.py` |
| Provincial calibration audit | `workflow/scripts/studies/chapter2/audit_baseline_calibration.py` |
| Named-station audit | `workflow/scripts/studies/chapter2/audit_named_station_validation.py` |
| Workflow configuration | `config/workflow.yaml` |
| Prepared-network audits | `results/workflow/base_model/model_preparation_*.csv` |
| Integrity audit and figures | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/` |
| Calibration tables and figures | `results/workflow/calibration/annual_reference_2021/` |
| Station screen and figure | `results/workflow/calibration/named_stations_2021/` |
| Economic-ordering audit and figure | `results/workflow/sensitivities/hydro_dispatch_economics/` |
| Station-band solved network and audits | `results/workflow/sensitivities/station_energy_band_25pct/` |
| Hydro sensitivity comparison | `results/workflow/sensitivities/hydro_dispatch_comparison/` |
| Publication Methods, Results, and decision schematic | `studies/chapter2/hydro_dispatch_sensitivity_methods_results_2026-08-01.md` |
| Combined stop/go decision | `results/workflow/readiness/publication_readiness_report.md` |

## Publication decision

The reference is computationally ready for controlled diagnostics, but it is not
ready for a claim that PyPSA-BC reproduces realized station operation. The Seven
Mile and Bridge River mechanism is diagnosed and bounded by a solved sensitivity;
the remaining gate is empirical selection or justification of the hydro dispatch
policy, resolution of Kootenay Canal's source profile, and held-out operational
validation. A matched A/B/C comparison and the full scenario matrix remain on
hold until those requirements are met.
