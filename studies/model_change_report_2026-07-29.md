# PyPSA-BC model change report

**Reporting date:** 2026-07-31 (updated from the 2026-07-29 baseline)  
**Scope:** publication-readiness work for the Chapter 2 hydropower-cascade study  
**Current decision:** Gate 1 PASS (no errors, two warnings); Gate 2 PASS (7/7);
Gate 2B PASS (12/12 basin-treatment checks); Gate 2C PASS (42/42
station-treatment checks); Gate 2D PASS (6/6 representation-case checks)

## Executive summary

The model now represents the Bridge River system with separate turbine and dam-
release routes, restores La Joie as an upstream generating reservoir, treats
Duncan as a documented non-generating control reservoir, and maintains water in
m3/h consistently through both reservoir and water-tracking run-of-river
components. The inflow build was corrected so monthly calibration occurs after
UTC-to-local horizon alignment; every modelled source-backed series now
reconciles exactly to its monthly target.

These changes materially improve hydraulic fidelity and auditability. They do
not make the study publication-ready by themselves. Minimum-release execution
coverage is 20/20: 13 observed/not-applicable records and seven predeclared
scenario bounds. The latter remain evidence gaps and cannot be presented as
observed or legal rules. One additional warning concerns Waneta's omitted
controllable storage. Walden North is now resolved from the accepted Bridge
River WUP as a no-storage project. Walter Hardman's documented
459,350 m³ live headpond storage is now represented explicitly.

## Changes to physical topology

### Bridge River and Seton routes

The previous topology could not distinguish water passing through a generating
station from water released at a dam into a different river reach. A
source-labelled route override table now supplies separate turbine and spill
destinations.

| Station | Turbine destination | Dam/spill destination |
|---|---|---|
| La Joie | Carpenter Lake / `BC_BR1_RES` | Carpenter Lake / `BC_BR1_RES` |
| Bridge River 1+2 aggregate | Seton Lake / `BC_SON_RES` | Lower Bridge River terminal sink |
| Walden | Seton Lake / `BC_SON_RES` | Cayoosh Creek terminal sink |
| Seton | Fraser River terminal sink | Seton River terminal sink |
| Alouette | Stave Lake / `BC_SFL_RES` | Alouette River terminal sink |
| John Hart | Lower Campbell River terminal sink | Elk Falls Canyon terminal sink |
| Walter Hardman | Arrow Lakes / `BC_ALH_RES` | Cranberry Creek terminal sink |
| Whatshan | Arrow Lakes / `BC_ALH_RES` | Whatshan River terminal sink |

Bridge River 1 and 2 remain electrically aggregated, but their maximum turbine
flows now add (65 + 95 = 160 m3/s) while the shared Terzaghi dam-release
capacity is represented once, not double-counted.

### La Joie restoration

The station/reservoir alias mismatch (`BC_LAJ_DFS` versus `BC_LAJ_GSS`) was
resolved in the reproducible asset builder. La Joie is now:

- included as a reservoir generator;
- supplied with a 12-month source-backed inflow profile;
- connected hydraulically to Carpenter Lake;
- serialized into `hydro_reservoirs.pickle`; and
- retained through the full hydro formatting step without a skipped-asset
  record.

### Terminal routing

Explicit non-cyclic terminal buses and Stores now receive water leaving the
modelled cascade. Separate terminal components are created when turbine and
spill routes leave by different river paths. This prevents artificial disposal
constraints and prevents end-of-horizon water from reappearing through a cyclic
terminal condition.

### Duncan role

Duncan is registered as a `non_generating_control_reservoir`, supported by the
official BC Hydro facility description that no power generation occurs at
Duncan Dam. Its source inflow table is retained for provenance but is not
treated as a missing generating-station input by Gate 1.

## Changes to hydraulic units and inflows

### Local-time calibration

ERA5 samples begin at 07:00 UTC, corresponding to midnight Pacific Standard
Time for the model horizon. Previously, monthly calibration was performed
before replacing the UTC index with the local model index. Boundary hours could
therefore move into an adjacent month after calibration.

The pipeline now:

1. verifies that cutout and model horizons contain the same number of hours;
2. aligns samples to the local model horizon;
3. calibrates each local calendar month to its source mean; and
4. rejects a length mismatch instead of silently relabelling it.

All modelled monthly source profiles now reconcile to numerical precision.

### Water-tracking run-of-river stations

The prior `ror-water` builder placed an energy-availability series on a bus
whose carrier was Water, then converted it back to water at the turbine Link.
For Walden, the model now uses its observed monthly-calibrated inflow directly
in m3/h. The turbine Link:

- consumes m3/h;
- produces electricity using MW/(m3/h);
- returns the same water quantity downstream; and
- uses a spill Link with unit water efficiency.

For water-tracking stations without observed flow tables, the fallback converts
the existing energy profile to m3/h before component construction. The internal
water-bus contract is therefore consistent in both paths.

## Operating constraints

The model now supports:

- constant combined turbine-plus-spill minimum release;
- evidence-labelled hourly schedules on an explicitly named hydraulic route;
- strict rejection of duplicate, missing, extra, negative, or unknown schedule
  rows; and
- preservation of schedule metadata through NetCDF export/import.

The reproducible schedule builder now generates 35,040 rows for four routes and
8,760 snapshots: Seton River (5 m3/s), Alouette River (1.52-3 m3/s), Elk Falls
Canyon (4-7 m3/s), and Cranberry Creek (0.1 m3/s). Optional pulses, combined-
route targets, and state/intertemporal rules remain explicitly deferred rather
than assigned arbitrary dates or constants.

Site C now carries the federal 390 m3/s minimum. Whatshan, Stave Falls, and
Seven Mile have evidence-backed zero values only for the minimum-release
parameter; their other reservoir or coordination rules remain separate.

## Matched hydraulic aggregation treatments

The production builder now implements the study's two counterfactual hydraulic
representations. `B_two_store` pools non-terminal levels into an upper bucket
and terminal levels into a lower bucket. `C_single_bucket` pools all modeled
storage and inflow in each basin into one bucket. Both transformations begin
from Configuration A and preserve total controllable storage, every hourly
inflow series, station turbine water capacities, water-to-power efficiencies,
electrical output buses, and explicit external environmental-release routes.

An actual-input verifier applies B and C to all six modeled cascade groups. All
12 basin-treatment rows pass storage, inflow, turbine-capacity, electrical-bus,
external-route, obsolete-link, placeholder-destination, and water-self-loop
checks. The workflow now compiles 15 A/B/C × water-policy cases: 12 sensitivity
cases are runnable and three evidence-only cases remain blocked. Solve
authorization remains disabled pending scientific review of the uncertainty
release and storage uncertainty protocols.

The same audit exposed and fixed a stale Waneta route: `DEFAULT Water Bus` had
been treated as a real reservoir because the token appeared in the global
reservoir-ID list. The formatter now removes blank and placeholder reservoir
identifiers before resolving downstream routes. Waneta turbine and spill water
now enter `Seven Mile Water Bus`, backed by the non-cyclic Seven Mile terminal
Store.

## Evidence and audit changes

- Added source-labelled hydraulic route overrides.
- Added a reservoir-role decision table.
- Added the official BC Hydro Columbia facility page to the validation source
  registry and SHA-256 manifest.
- Added route and role evidence to the cascade evidence register.
- Extended Gate 1 to distinguish evidence-backed non-generating reservoirs
  from unresolved orphan reservoirs.
- Retained a strict error when a publication-critical minimum-release rule is
  unresolved.

## Verification status

| Check | Result |
|---|---|
| Focused Chapter 2 regression suite | 45 passing tests |
| Gate 2 analytic physics | PASS, 7/7 required cases |
| Gate 2B matched hydraulic representations | PASS, 12/12 basin-treatment checks |
| Gate 2C release uncertainty | PASS, 42/42 station-treatment checks |
| Gate 2D storage uncertainty | PASS, 6/6 representation-case checks |
| Gate 1 unit checks | 29/29 pass |
| Gate 1 builder checks | 14/14 pass |
| Gate 1 overall | PASS: 0 errors, 2 warnings |

The regression suite covers route splitting, shared dam-release capacity,
non-generating reservoir roles, local-time inflow alignment, water-unit
conversion, analytic delay and mass balance, route schedules, terminal
boundaries, NetCDF schedule persistence, matched B/C storage aggregation,
placeholder-route rejection, and preservation of electrical and external water
routes.

## Remaining publication blockers

1. Replace the seven release-rule sensitivities with exact rules if primary
   evidence becomes available; until then, retain the lower/upper claim boundary.
2. Resolve Waneta active storage or predeclare and test a defensible storage
   envelope. Walden North is resolved as having no storage.
3. Review the 24 runnable cases, then authorize the full scenario matrix and
   independent water-balance/result validation. The evidence-only cases remain blocked.

## Reproduction order

From the repository root in the project environment:

```powershell
$env:PYTHONPATH='src;.'
python -m workflow.scripts.create_hydro_assets
python -m workflow.scripts.create_reservoir_inflows
python -m workflow.scripts.enrich_format_hydro
python -m workflow.scripts.studies.chapter2.fetch_validation_sources
python -m workflow.scripts.studies.chapter2.build_cascade_evidence_register
python -m workflow.scripts.studies.chapter2.resolve_cascade_evidence
python -m workflow.scripts.studies.chapter2.build_release_schedules
python -m workflow.scripts.studies.chapter2.audit_inputs
python -m workflow.scripts.studies.chapter2.verify_physics
python -m workflow.scripts.studies.chapter2.audit_hydraulic_representations
python -m pytest tests/studies/chapter2 -q
```

The Gate 1 command is expected to exit non-zero while the stated
publication-critical error remains.
