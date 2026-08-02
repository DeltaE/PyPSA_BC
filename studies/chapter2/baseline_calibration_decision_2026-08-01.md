# Baseline calibration decision — 2026-08-01

## Decision

The corrected optimal 2021 annual reference passes all computational-integrity
and declared calibration-fit checks. This is **calibration fidelity, not
independent validation**: monthly biomass and combustible totals are imposed as
calibration constraints and the load chronology is an input. Full scenario
production remains on hold pending period-matched hydro-operation evidence,
resolution of the Kootenay Canal source-profile discrepancy, and held-out
operational validation. The Seven Mile and Bridge River mechanism has been
diagnosed and bounded by a solved sensitivity, but not independently validated.

Observed BC Hydro interchange remains the primary external-boundary condition.
TTC is deferred to a separately labelled secondary sensitivity and must not be
used in the core solve, calibration gate, principal figures, or headline claims.

## Reproducible evidence

The independent generation benchmark is Statistics Canada Table
25-10-0015-01, calendar 2021, British Columbia, all producer classes. The fetch
and filter procedure is implemented in
`workflow/scripts/fetch_statcan_generation_validation.py`. It retains the
official bulk archive, table metadata, a filtered monthly table, retrieval
metadata, source URLs, and SHA-256 hashes under `data/validation/statcan/`.

The comparison is reproduced by
`workflow/scripts/studies/chapter2/audit_baseline_calibration.py`. Audited
tables, machine-readable decisions, SVG figures, and PNG previews are written to
`results/workflow/calibration/annual_reference_2021/`.

## Observed–model comparison

| Quantity | Observed | Modeled | Relative error | Gate |
|---|---:|---:|---:|---|
| BC Hydro gross telemetered load | 65.05 TWh | 64.43 TWh | -0.95% | PASS, input-fidelity only |
| Peak BC Hydro load | 11.79 GW | 11.68 GW | -0.95% | PASS, input-fidelity only |
| Total provincial generation | 71.99 TWh | 70.58 TWh | -1.95% | PASS, calibration alignment |
| Hydro generation | 63.80 TWh | 62.24 TWh | -2.45% | PASS, calibration alignment |
| Wind generation | 2.00 TWh | 2.16 TWh | +8.19% | PASS |
| Biomass generation | 4.07 TWh | 4.07 TWh | 0.00% | PASS, enforced input |
| Non-renewable combustible generation | 2.11 TWh | 2.11 TWh | 0.00% | PASS, enforced input |
| Solar generation | 0.004 TWh | 0.002 TWh | -58.9% | Diagnostic |

The load chronology correlation is 1.000 because that observed series supplies
the modeled chronology. It is therefore an input-alignment check, not independent
validation. Hydro monthly shape and annual magnitude now pass. Generation-mix
total-variation distance is 0.0046. Because Statistics Canada is now explicitly
used for calibration, these values cannot be relabelled as held-out validation.

## Scientific interpretation

The apparently close total-generation result conceals a structural mix error.
In a cost-only dispatch with ample modeled water, hydro displaces biomass and
gas. Real production can reflect cogeneration heat demand, energy contracts,
minimum stable operation, outages, fuel and environmental constraints, and
owner decisions that are not represented in the current objective and
constraints. Reducing hydro inflow solely to match the observed total would
hide rather than resolve this error.

The Statistics Canada benchmark also has a broader accounting boundary: it
covers all producer classes in the province, while PyPSA-BC is organized around
a balancing-authority-oriented asset and load representation. Category-level
constraints must not be imposed until the asset/accounting boundary is
reconciled.

The inventory audit found a concrete vintage defect: two 15 MW solar projects
with 2025 start years were present in the 2021 input dictionaries. Two thermal
facilities with recorded 2010 closure years were also present. Model assembly now
filters solar, wind, and thermal component dictionaries using source inventory
start and closure years and writes an asset-vintage audit beside the build report.
For 2021 this retains one of three solar components, all ten wind components, and
35 of 37 thermal components. The corrected annual reference is optimal and again
passes all eight computational-integrity checks. It removes the artificial solar
surplus but does not resolve the material biomass/combustible dispatch gap.

## Hydro dispatch-policy finding

The reference assigns Seven Mile and Bridge River an effective turbine cost of
CAD 6.88/MWh. Their maximum solved local marginal price is approximately CAD
1.968/MWh, so price never exceeds turbine cost. This economic ordering explains
why the optimization uses spill routes at those facilities; it does not prove
that the generic cost represents actual station operations.

A ±25% fiscal-energy-band sensitivity solved optimally, passed all eight annual
integrity checks, increased Seven Mile generation from 11.6 to 2,279.2 GWh, and
increased Bridge River generation from 0 to 1,664.3 GWh. Total hydro remained
62.240 TWh, while the model objective increased by 4.60%. Because the sensitivity
uses the fiscal station values as constraints, those values cannot independently
validate its output. Kootenay Canal remained 34.8% below its fiscal comparator.

## Frozen calibration gates

The following thresholds were set before representation comparisons:

| Check | Threshold |
|---|---:|
| Annual load magnitude | absolute relative error ≤ 2% |
| Peak load magnitude | absolute relative error ≤ 2% |
| Hourly load chronology | Pearson `r ≥ 0.99`, peak position ≤ 1 h |
| Total generation magnitude | absolute relative error ≤ 5% |
| Hydro generation magnitude | absolute relative error ≤ 5% |
| Hydro monthly shape | Pearson `r ≥ 0.90` |
| Wind generation magnitude | absolute relative error ≤ 10% |
| Material non-hydro operation | modeled dispatch required when observed share > 1% |
| Generation mix | total-variation distance ≤ 0.05 |

Passing these checks does not prove station-level validity; it only authorizes
the matched A/B/C representation comparison to proceed to the next validation
layer.

## Model-development sequence

1. Produce a facility crosswalk from every modeled generator to the Statistics
   Canada producer class and the 2021 operating boundary. Record commissioning,
   retirement, ownership, fuel, cogeneration role, and inclusion decision.
2. Reconcile the provincial load, internal generation, imports, exports, losses,
   and any behind-the-meter or non-BA production in one annual accounting table.
3. Rebuild the prepared and solved reference with the implemented 2021 asset-
   vintage filter, then verify that the solar discrepancy and every integrity
   check respond as expected.
4. For biomass and combustible plants, encode source-labelled operating drivers.
   Prefer facility or facility-group monthly energy profiles where evidence
   exists. Otherwise use transparent monthly calibration bands as a baseline
   scenario, with an unconstrained economic-dispatch sensitivity.
5. Keep hydro inflows source-calibrated. Let non-hydro correction alter hydro
   dispatch endogenously, then re-evaluate hydro magnitude, monthly shape,
   storage, releases, spill, and cyclic closure.
6. Rerun the annual solve and require both the eight computational-integrity
   checks and the frozen calibration gate to pass.
7. Obtain a period-matched station or water-operation dataset, or predeclare an
   independently justified hydro-dispatch policy, before selecting the primary
   reference from the source-cost and station-band formulations.
8. Resolve the Kootenay Canal fixed-profile discrepancy and repeat the held-out
   station screen.
9. Only then run the first matched A/B/C representation comparison. Do not run
   the full scenario matrix until that paired comparison reproduces all gates.

## Claim boundary

At present, PyPSA-BC supports claims about a reproducible, internally balanced
2021 reference formulation and exposes a diagnosed baseline calibration error.
It does not yet support claims that modeled resource dispatch reproduces realized
BC operation, that internal corridors are empirically congested, or that the
model objective estimates observed system cost. TTC remains outside this claim
path until explicitly reactivated as a secondary sensitivity.

The hydro policy sensitivity supports one narrower claim: the existing cascade
topology can route material water through Seven Mile and Bridge River, and their
reference spill preference arises from the modeled economic ordering. It does
not identify the empirically correct station dispatch.
