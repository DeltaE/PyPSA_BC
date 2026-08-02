# Chapter 2 publication workspace

This directory contains the reproducible evidence for:
**How Hydropower Cascade Aggregation Affects Operational Feasibility Estimates
in a Hydro-Dominant Power System**.

## Current decision

- **Gate 1: PASS with two claim-boundary warnings** - no input-gate errors
  remain. Seven release rules and Waneta pondage are covered by predeclared
  scenario envelopes rather than observed/legal values.
- **Gate 2: PASS** - all seven required solved physics cases pass with PyPSA 1.2.4
  and HiGHS.
- **Gate 2B: PASS** - all 12 basin-treatment conservation audits pass for the
  implemented B two-store and C single-bucket transformations.
- **Gate 2C: PASS** - all 42 release-bound checks pass across seven stations,
  two bounds, and three hydraulic representations.
- **Gate 2D: PASS** - all six Waneta pondage checks pass across two bounds and
  three hydraulic representations.
- Travel-time coverage is 20/20 as a declared structural scenario, not as
  station-calibrated evidence.
- Minimum-release execution coverage is 20/20: 13 stations have observed or
  not-applicable decisions and seven have scenario-class bounds. The latter are
  not estimates of legal or observed operations.
- The calendar-2021 reference is optimal and passes all eight annual integrity
  checks. Calibration fit passes, but independent operational validation is not
  established.
- The Seven Mile and Bridge River spill mechanism is diagnosed and a solved
  station-energy-band sensitivity passes computational checks. Because that
  sensitivity uses cross-period fiscal evidence, the publication gate remains
  **HOLD**.

## Start here

1. `publication_package_index_2026-08-01.md` - authoritative code, result,
   figure, schematic, documentation, and readiness map.
2. `publication_execution_tracker.md` - current scientific decision, strategy,
   dated 22-week execution map, and stop/go rules.
3. `../chapter2_publication_strategy_and_timeline_report.md` - full publication
   strategy and modelling plan.
4. `gate1/gate1_report.md` - latest input/evidence stop-go decision.
5. `gate2/gate2_report.md` - solved physical-verification results.
6. `hydraulic_representation_methods.md` - publishable A/B/C definitions,
   schematic, invariants, and claim boundary.
7. `gate2/representation_verification.md` - actual-input B/C conservation audit.
8. `release_uncertainty_methods.md` and
   `gate2/release_uncertainty_verification.md` - bound definitions, claim rules,
   and actual-input verification.
9. `storage_uncertainty_methods.md` and
   `gate2/storage_uncertainty_verification.md` - Waneta pondage bounds and gate.
10. `inputs/cascade_evidence_register.csv` - station-level parameter evidence.
11. `../../data/validation/chapter2/source_manifest.csv` - downloaded source
   hashes.
12. `../../data/validation/chapter2/constraint_source_hits.md` - page-level
   candidate passages for manual evidence approval.
13. `bc_hydro_data_integration_report_2026-07-31.md` - source-by-source decisions,
   implemented routes and schedules, and the remaining evidence gaps.
14. `../../data/validation/chapter2/storage_representation_decisions.csv` -
    facility-level storage decisions and unresolved owner-evidence requirements.
15. `hydro_dispatch_sensitivity_methods_results_2026-08-01.md` - publishable
    Methods and Results text, decision schematic, captions, and claim boundary.
16. `publication_figure_chart_map_2026-08-01.md` - figure contracts, evidence
    roles, output paths, and visual QA status.
17. `model_change_report_2026-08-01.md` - complete implemented-change and result
    record for the current reference and hydro sensitivity.

## Reproduce the evidence layer

Run from the repository root in the locked model environment:

```powershell
$env:PYTHONPATH='src;.'
python workflow/scripts/studies/chapter2/fetch_validation_sources.py
python workflow/scripts/studies/chapter2/report_validation_sources.py
python workflow/scripts/studies/chapter2/build_cascade_evidence_register.py
python workflow/scripts/studies/chapter2/resolve_cascade_evidence.py
python workflow/scripts/studies/chapter2/build_release_schedules.py
python workflow/scripts/studies/chapter2/apply_cascade_evidence.py
python workflow/scripts/studies/chapter2/audit_inputs.py --protocol config/studies/chapter2_cascade.yaml
python workflow/scripts/studies/chapter2/verify_physics.py
python workflow/scripts/studies/chapter2/audit_hydraulic_representations.py
python -m workflow.scripts.studies.chapter2.audit_release_uncertainty
python -m workflow.scripts.studies.chapter2.audit_storage_uncertainty
python -m unittest discover -s tests/studies/chapter2 -v
```

The focused regression suite contains 45 tests. The workflow compiles 30
A/B/C × water-policy × storage-policy cases; 24 are runnable. The six
evidence-only cases remain blocked because seven exact facility rules are
unresolved. Production solving remains disabled pending review of both
uncertainty protocols.

## Scheduled route-specific releases

`inputs/release_schedules.csv` is the production input for seasonal or otherwise
time-varying minimum releases that apply to a particular hydraulic route. It
uses one row per route and network snapshot:

| Column | Meaning |
|---|---|
| `link` | Exact PyPSA Link name of the constrained river, bypass, outlet, or spill route |
| `snapshot` | Network timestamp; every snapshot must occur exactly once for each link |
| `minimum_release_m3_per_hour` | Finite, non-negative minimum route flow |
| `schedule_id` | Stable identifier joining the profile to an approved operating rule |
| `source` | Evidence-registry source identifier |
| `evidence_class` | `observed`, `derived`, `calibrated`, `imputed`, or `scenario` |

The loader rejects missing hours, extra hours, duplicate link/timestamp rows,
unknown links, missing source metadata, unresolved evidence classes, and
negative values. It does not infer a route from a station name. The production
table contains complete 2021 schedules for Seton River, Alouette River, Elk
Falls Canyon, and Cranberry Creek. Rules that depend on reservoir state, tides,
operator-selected event timing, or prior dispatch remain explicitly deferred
in `release_rule_definitions.yaml`.
