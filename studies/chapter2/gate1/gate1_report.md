# Chapter 2 Gate 1 — Input, Unit, and Topology Integrity

**Gate status:** PASS
**Protocol:** 0.1.0 (2026-07-29)
**Generated:** 2026-08-01T23:15:14.830155+00:00

## Decision

Gate 1 passes. The frozen inputs can proceed to physical verification.

## Evidence snapshot

- Inputs declared: 23
- Cascade groups represented: 6
- Cascade stations represented: 20
- Errors: 0
- Warnings: 2
- Unit checks passing: 29/29
- Builder checks passing: 14/14

## Publication-critical parameter coverage

| Parameter | Available | Required station records |
|---|---:|---:|
| max_discharge_m3_per_hour | 20 | 20 |
| minimum_release_m3_per_hour | 20 | 20 |
| spill_capacity_m3_per_hour | 19 | 20 |
| travel_time_h | 20 | 20 |
| usable_storage_m3 | 20 | 20 |
| water_to_power_efficiency_mwh_per_m3 | 20 | 20 |

### Minimum-release coverage by evidence class

- not_applicable: 6 station(s)
- observed: 7 station(s)
- scenario: 7 station(s)

Scenario-class records are structural lower/upper sensitivities; they are not observed, calibrated, or legal release rules.

## Open issues

| Severity | Check | Subject | Finding | Required action |
|---|---|---|---|---|
| warning | storage_uncertainty_claim_boundary | 1 cascade station(s) | Controllable storage is bracketed from 0 to 24 hours of mean inflow, not measured volume. | Report core conclusions only when robust across both storage bounds. |
| warning | release_uncertainty_claim_boundary | 7 cascade stations | Release rules are covered by structural sensitivity bounds, not observed/legal rules. | Report core conclusions only when their direction and material interpretation are robust across both bounds. |

## Reproduction

Run from the repository root in the locked model environment:

```powershell
$env:PYTHONPATH='src'
python workflow/scripts/studies/chapter2/audit_inputs.py --protocol config/studies/chapter2_cascade.yaml
```

The CSV and JSON files beside this report are the machine-readable evidence.
