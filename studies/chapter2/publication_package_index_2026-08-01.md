# PyPSA-BC Chapter 2 publication package index

**Package date:** 2026-08-01  
**Reference year:** calendar 2021  
**Overall decision:** **HOLD for publication claims and A/B/C production**  
**Reason:** computational integrity and calibration fit pass, but held-out
operational validation is not established.

## What is complete

| Layer | Evidence | Decision |
|---|---|---|
| Input and evidence audit | Gate 1 manifests, hashes, units, topology, source roles, and assumption register | PASS with documented evidence limitations |
| Hydraulic physics | Seven solved analytic cases plus B/C conservation, release-envelope, and storage-envelope audits | PASS |
| Annual reference | Optimal 8,760-hour station-resolved network | PASS all eight integrity checks |
| Provincial calibration | Load fidelity and Statistics Canada generation comparisons | PASS calibration fit; not independent validation |
| Named-station capacity | Seven BC Hydro capacities aligned and audited | PASS as calibrated capacity layer |
| Named-station energy | Fiscal-2021 energy versus calendar-2021 model | REVIEW for Kootenay Canal, Seven Mile, and Bridge River in the unconstrained reference |
| Hydro economic mechanism | Effective turbine VOM versus solved hourly local marginal price | Complete; source-cost ordering explains Seven Mile and Bridge River spill preference |
| Hydro station-band sensitivity | Optimal full-year solution with six ±25% station bands | PASS computationally; calibration sensitivity only |
| Publication figures | Static PNG/SVG economic and station-sensitivity figures with source tables and captions | Complete and visually inspected |
| Manuscript Methods/Results | Claim-safe hydro sensitivity text and schematic | Complete for diagnostic Results or supplement |
| TTC | Optional code path only | Deferred; excluded from core results and gates |
| Regression verification | 98 tests plus 6 subtests | PASS on 2026-08-01 |

## Authoritative model and result artifacts

| Artifact | Path | Use |
|---|---|---|
| Prepared reference network | `results/workflow/base_model/network_prepared.nc` | Audited input model |
| Solved reference network | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon.nc` | Primary computational reference |
| Reference solve summary | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_summary.json` | Solver contract |
| Reference integrity summary | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/diagnostic_summary.json` | Eight-check result gate |
| Provincial calibration summary | `results/workflow/calibration/annual_reference_2021/calibration_summary.json` | Calibration and evidence roles |
| Named-station reference summary | `results/workflow/calibration/named_stations_2021/named_station_summary.json` | Capacity and energy screen |
| Publication-readiness summary | `results/workflow/readiness/publication_readiness_summary.json` | Combined stop/go decision |
| Solved station-band network | `results/workflow/sensitivities/station_energy_band_25pct/network_solved.nc` | Cross-period calibration sensitivity |
| Station-band integrity summary | `results/workflow/sensitivities/station_energy_band_25pct/diagnostics/diagnostic_summary.json` | Sensitivity integrity gate |
| Hydro sensitivity comparison | `results/workflow/sensitivities/hydro_dispatch_comparison/hydro_dispatch_sensitivity_summary.json` | Comparative quantitative result |
| Economic-ordering summary | `results/workflow/sensitivities/hydro_dispatch_economics/hydro_dispatch_economics_summary.json` | Mechanism diagnostic |

## Reproducible code map

| Function | Source |
|---|---|
| Cascade construction and physics | `src/pypsa_bc/studies/cascade/` |
| Hydro dispatch-cost policies | `src/pypsa_bc/studies/hydro_dispatch_policy.py` |
| Station capacity and energy policies | `src/pypsa_bc/studies/station_validation.py` |
| Generation calibration | `src/pypsa_bc/studies/generation_calibration.py` |
| External boundary | `src/pypsa_bc/studies/external_boundary.py` |
| Model preparation | `workflow/scripts/prepare_model_file.py` |
| Optimal-only annual solver | `workflow/scripts/solve_model.py` |
| Reference integrity audit | `workflow/scripts/studies/chapter2/audit_full_horizon_results.py` |
| Calibration audit | `workflow/scripts/studies/chapter2/audit_baseline_calibration.py` |
| Named-station audit | `workflow/scripts/studies/chapter2/audit_named_station_validation.py` |
| Economic-ordering audit | `workflow/scripts/studies/chapter2/audit_hydro_dispatch_economics.py` |
| Hydro sensitivity comparison | `workflow/scripts/studies/chapter2/compare_hydro_dispatch_sensitivities.py` |
| Combined readiness gate | `workflow/scripts/studies/chapter2/audit_publication_readiness.py` |
| Workflow orchestration | `workflow/Snakefile` and `workflow/rules/common.smk` |

## Schematics and manuscript documentation

| Document | Contents |
|---|---|
| `hydraulic_representation_methods.md` | A/B/C definitions, invariants, and aggregation schematic |
| `hydro_dispatch_sensitivity_methods_results_2026-08-01.md` | Hydro policy Methods, Results, decision schematic, captions, and limitations |
| `model_change_report_2026-08-01.md` | Implemented corrections, audited results, and publication decision |
| `baseline_calibration_decision_2026-08-01.md` | Calibration thresholds, evidence roles, and claim boundary |
| `production_pilot_status_2026-08-01.md` | Solver, transport-flow, annual-result, and sensitivity status |
| `publication_execution_tracker.md` | Strategy, gates, remaining work, and timeline |
| `publication_figure_chart_map_2026-08-01.md` | Figure contracts and visual QA |

## Key numerical results

| Result | Reference | Station-band sensitivity |
|---|---:|---:|
| Objective (million model monetary units) | 419.955 | 439.254 |
| Objective change | — | +4.60% |
| Total hydro generation (TWh) | 62.240 | 62.240 |
| Seven Mile generation (GWh) | 11.6 | 2,279.2 |
| Seven Mile turbine capture | 0.34% | 67.4% |
| Bridge River generation (GWh) | 0 | 1,664.3 |
| Bridge River turbine capture | 0% | 64.1% |
| Emergency backstop (GWh) | 0 | 0 |
| Annual integrity gate | PASS | PASS |

The station-band values cannot be described as independently validated because
the sensitivity uses the fiscal station comparator as a constraint.

## Required work before the publication gate can pass

1. Obtain period-matched, held-out station generation, reservoir, or release
   observations and use them to evaluate hydro dispatch policy.
2. Resolve the Kootenay Canal fixed-profile discrepancy.
3. Close or formally bound remaining facility-level release and operating-rule
   evidence gaps.
4. Re-run the combined readiness gate with independent validation evidence.
5. Only after a `GO` decision, execute one matched A/B/C comparison and require
   the same integrity, calibration, and held-out validation checks.
6. Expand to the scenario matrix only if the matched comparison is interpretable.

## Verification record

| Verification | Result |
|---|---|
| Full Python regression suite | PASS — 98 tests and 6 subtests |
| Ruff checks on changed modelling, audit, solve, and test modules | PASS |
| Workflow preflight | PASS — 30 scenarios, 24 runnable, no missing required paths, solving disabled |
| Snakemake parser compilation | PASS — `workflow/Snakefile` compiles with 39 rules; `workflow/rules/common.smk` compiles |
| Generated JSON parsing | PASS for attempt register, economic summary, and sensitivity summary |
| Static hydro figures | PASS after original-resolution visual inspection |
| Combined publication-readiness gate | HOLD — named-station energy screen and independent validation remain blocking |
| Publication-package completeness audit | PASS — 29/29 required artifacts found; figure format and resolution checks pass |

The installed Snakemake CLI timed out after announcing the workflow-specific
profile when asked to list rules. The independent Snakemake parser compilation
passes, and the underlying workflow Python targets pass their tests. Resolve the
CLI/profile startup issue before relying on unattended orchestration on this
Windows host.

Machine-readable and reader-facing audit outputs are stored at
`results/workflow/readiness/publication_package_audit.json` and
`results/workflow/readiness/publication_package_audit.md`. The reproducible audit
is implemented in
`workflow/scripts/studies/chapter2/audit_publication_package.py`.

## Claim that is currently supportable

PyPSA-BC provides a reproducible, internally balanced calendar-2021 reference
and shows that generic hydro cost ordering can materially alter station-level
water routing without changing provincial hydro energy. The current evidence
does not establish that modeled station dispatch reproduces realized BC Hydro
operation, that modeled corridor stress represents observed congestion, or that
the objective represents observed system cost.
