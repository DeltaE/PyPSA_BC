# Chapter 2 model-production status — 2026-08-01

## Decision

The 2021 `A_full_cascade × release_uncertainty_lower_bound × zero_pondage`
base case has a complete, optimal, full-horizon solution and passes the annual
computational-integrity gate. It is the current reference result for continued
calibration and scenario development. Publication readiness remains on hold until
external calibration, representation comparisons, and uncertainty results are
complete.

The corrected provincial-generation calibration passes all declared fit
thresholds. Hydro is 2.45% below the Statistics Canada total; biomass and
non-renewable combustible output reproduce 4.07 TWh and 2.11 TWh because their
monthly totals are enforced calibration inputs. This is not independent
validation. The official named-station capacity check passes, while the
cross-period energy screen flags Kootenay Canal, Seven Mile, and Bridge River.
The complete decision and correction sequence is recorded in
`studies/chapter2/baseline_calibration_decision_2026-08-01.md`.

The Seven Mile and Bridge River mechanism is now diagnosed. Their effective
CAD 6.88/MWh turbine cost exceeds the maximum solved local marginal price in
every hour, producing a spill-path preference in the source-cost reference. A
solved ±25% station-energy-band sensitivity restores material turbine routing,
preserves total hydro generation, passes all eight integrity checks, and raises
the model objective by 4.60%. This result bounds dispatch-policy uncertainty but
does not independently validate station operation.

The primary external-boundary policy is observed BC Hydro interchange. Hourly
Total Transfer Capability (TTC) is deprioritized and excluded from the core solve,
result gate, figures, and headline claims. TTC may later be tested as a secondary
boundary sensitivity; it is not an internal corridor rating.

## Current reference case

| Contract | Setting |
|---|---|
| Weather, load, and boundary year | 2021 |
| Spatial resolution | Regional district |
| Hydro representation | A — station-resolved full cascade |
| Water-policy case | Release-uncertainty lower bound |
| Short-duration storage case | Zero pondage |
| External boundary | Observed BC Hydro interchange |
| Internal electrical formulation | Capacity-constrained transport |
| KVL | Not enforced; aggregated equivalent admittances are unvalidated |
| Transport-flow reporting | Lexicographic minimum-transfer postprocessing |
| Optimization horizon | 8,760 hours, single full-year horizon |
| Solver | HiGHS dual simplex |

## Annual result and integrity gate

The solver returned `status=ok` and `termination=optimal`; only then was the
network exported. Independent diagnostics give:

| Check | Result |
|---|---|
| Complete chronology | PASS — 8,760 consecutive hourly snapshots |
| Annual period | 2021-01-01 00:00 to 2021-12-31 23:00 |
| Provincial demand | 64.43 TWh |
| Peak provincial load | 11.68 GW |
| Observed imports imposed | 0.72 TWh |
| Observed exports imposed | 6.87 TWh |
| External-boundary fidelity | PASS — 0 MW maximum error |
| Electrical balance | PASS — 7.28e-12 MW maximum residual |
| Emergency backstop | PASS — 0 MWh |
| Release constraints | PASS — 0 violation hours |
| Reservoir bounds | PASS — 0 violations |
| Cyclic annual closure | PASS — 0 capacity-normalized residual at exported precision |
| Model objective | 419,954,818 model monetary units |

The objective is reported only as a model diagnostic. It has not yet been
validated as an observed BC system-cost estimate and must not be described as
such in a manuscript.

## Transport-flow identifiability and congestion reporting

The first annual flow plot exposed a zero-cost circulation defect in reporting.
Lossless bidirectional transport links can satisfy the same nodal injections with
arbitrary cycle flows, so raw solver flows are not unique. In the raw solution,
many corridors appeared saturated for most of the year; those values were not a
valid basis for congestion claims.

The reporting workflow now performs a second, hourly linear program after the
cost-optimal solve. It holds each node's solved transport divergence fixed and
minimizes the sum of absolute corridor transfers. This is a lexicographic
minimum-transfer reconstruction: generation, load, trade, water routing,
reservoir states, capacity limits, and the primary objective are unchanged.

The reconstruction removed 405.71 TWh of absolute raw-solver transfer attributable
to non-unique circulation. Maximum reconstruction error was 5.01e-12 MW. Eight
modeled corridors reached at least 99% of their capacity in one or more hours after
reconstruction. These are diagnostic stress signals, not yet publishable congestion
findings, because internal corridor capacities remain proxy data requiring
source-level validation.

An attempted alternative with two directional links and a small in-objective
flow penalty passed a 336-hour pilot but did not achieve a feasible annual optimum
within the 1,200-second bound. No network was exported. It is retained as failed
method evidence, and the lexicographic postprocessor is the primary reporting
method because it preserves the proven annual optimum exactly.

## Defects resolved before the annual result

1. Three wind assets had 8,760 missing availability values. CRS-aware raster
   sampling now produces finite profiles and reproduces the CODERS annual-energy
   targets for all ten wind assets.
2. Observed interchange contains imports and exports, while the original external
   generator bounds prohibited negative dispatch. External-market generators are
   now bidirectional and fixed setpoints are independently checked against bounds.
3. The solve adapter previously allowed warning-status exports. It now exports only
   `status=ok` plus `termination=optimal`; failures produce no solved network.
4. Terminal water sinks previously used `1e15 m3` placeholders. Each is now bounded
   by total modeled horizon inflow plus a 1% numerical margin.
5. Water is internally scaled to million-m3 units, with inflow, storage, turbine
   efficiency, costs, and release constraints transformed together.
6. Biomass and natural-gas inventory stores now have finite non-binding bounds and
   use a 1,000-MWh-fuel internal unit, removing excessive numerical ranges without
   changing electrical output or total cost.
7. Enforcing DC load flow on unvalidated regional equivalent admittances created
   artificial loop-flow constraints and 17.85% backstop use in the 336-hour pilot.
   The primary model now uses a capacity-constrained transport formulation; DC load
   flow is retained only as a documented sensitivity.
8. Raw transport-flow circulation is now removed through reproducible
   lexicographic minimum-transfer postprocessing before utilization is reported.
9. Model-year filtering now excludes generation components outside their source
   inventory start/closure years. For the 2021 build this removes two 2025 solar
   projects and two thermal facilities recorded as closed in 2010. The corrected
   annual reference is optimal and passes all eight integrity checks.
10. Observed-interchange construction and verification no longer load TTC data.
    TTC workbooks are dependencies only when the optional `hourly_ttc` policy is
    explicitly selected.

## Rolling-horizon evidence

A two-window, 336-hour transport pilot solved optimally, transferred reservoir
states exactly, reproduced observed interchange, used no backstop generation, and
had no release or storage-bound violations. Longer rolling sequences became
infeasible because finite look-ahead and terminal storage policy distorted water
value near window boundaries. Rolling horizon is therefore not the primary annual
method at present. The successful full-horizon solution supersedes it for the base
case; rolling horizon remains a method-development branch requiring a stronger
terminal-value or target-trajectory formulation.

The original DC-load-flow failure evidence is preserved in
`results/workflow/pilots/rolling_A_lower_zero_observed_dc_sensitivity/`. The
successful transport pilot is in
`results/workflow/pilots/rolling_A_lower_zero_observed/`.

## Hydro dispatch-policy sensitivity

Two explicit hydro-cost policies are implemented: `source_generic_vom` and
`uniform_reservoir_vom`. The first remains the unconstrained reference. The
uniform CAD 1.97/MWh policy produced no exported optimal annual solution in the
controlled simplex or interior-point attempts; this is not an infeasibility
finding.

The accepted station-band diagnostic retained source-generic costs and imposed
±25% fiscal-energy bands at six reservoir links. It solved optimally. Seven Mile
generation increased from 11.6 to 2,279.2 GWh and Bridge River from 0 to 1,664.3
GWh, while G.M. Shrum and Mica decreased. Total hydro remained 62.240 TWh.
Kootenay Canal was not constrained and remained 34.8% below its fiscal
comparator. The full execution record, Methods and Results text, decision
schematic, tables, and figures are in
`studies/chapter2/hydro_dispatch_sensitivity_methods_results_2026-08-01.md`.

## Evidence locations

- Solved annual network:
  `results/workflow/annual/A_lower_zero_observed_transport_full_horizon.nc`
- Solver contract:
  `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_summary.json`
- Annual audit, tables, SVG figures, and PNG previews:
  `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/`
- Prepared input network and build audits:
  `results/workflow/base_model/`
- Failed directional-regularization attempt:
  `results/workflow/annual/A_lower_zero_observed_transport_regularized_full_horizon_summary.json`
- Baseline calibration audit, tables, and figures:
  `results/workflow/calibration/annual_reference_2021/`
- Pinned Statistics Canada generation evidence and manifest:
  `data/validation/statcan/`
- Hydro economic-ordering audit and figure:
  `results/workflow/sensitivities/hydro_dispatch_economics/`
- Station-band solved network and complete audits:
  `results/workflow/sensitivities/station_energy_band_25pct/`
- Hydro dispatch comparison table, summary, and figure:
  `results/workflow/sensitivities/hydro_dispatch_comparison/`

## Remaining publication gates

1. Obtain period-matched station generation, release, or reservoir evidence and
   use it to select or reject the source-cost and station-band dispatch policies.
2. Resolve the Kootenay Canal fixed-profile discrepancy and repeat the station
   screen on a compatible period.
3. Resolve or bound the remaining facility-level release-rule evidence gaps and
   document each scenario-labelled assumption.
4. Execute the implemented reservoir-representation, water-policy, and storage-policy
   matrix only after one comparison case reproduces this base-case integrity gate.
5. Quantify representation error, operational differences, robustness, and
   uncertainty rather than presenting a single deterministic dispatch.
6. Validate the highest-stress internal corridor capacity proxies before making
   congestion claims; otherwise retain them as model diagnostics and limitations.
7. Produce manuscript figures from audited tables only, with code, units, data
   provenance, and claim boundaries attached.

## Claim boundary

The annual reference case is computationally valid and internally consistent. It
is not yet a calibrated representation of realized BC Hydro operations, a causal
scenario result, or a publishable estimate of system cost or congestion. TTC is
outside the current primary claim path and remains deferred.
