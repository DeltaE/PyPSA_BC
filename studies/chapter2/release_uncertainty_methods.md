# Release-rule uncertainty methods

## Purpose and claim boundary

Seven cascade stations lack an exact, implementation-ready release rule for the
2021 study horizon. The model therefore uses two deterministic sensitivity
policies. These policies bracket structural uncertainty; they do not estimate
legal compliance, observed dispatch, or the probability distribution of flows.
The evidence-only policy remains blocked.

The lower and upper policies preserve every verified static minimum and all four
source-labelled hourly route schedules. They change only the seven unresolved
stations. We will report a core aggregation result only when its direction and
material interpretation remain stable across both bounds.

## Predeclared bounds

| Station | Scope | Lower bound | Upper bound | Interpretation |
|---|---|---:|---:|---|
| La Joie | Combined turbine plus spill | 5.7 m³/s | 18.4 m³/s | Constants bracket the documented elevation-dependent rule |
| Ruskin | Combined turbine plus spill | 0 m³/s | 100 m³/s | Upper case applies the WUP peaking floor continuously and is deliberately conservative |
| Arrow Lakes | Combined turbine plus spill | 0 | 10% of mean natural inflow | Structural stress range; owner/legal rule unresolved |
| Bridge 1/Terzaghi | Lower Bridge River spill route only | 0 | 10% of mean natural inflow | Never applied to the turbine diversion into Seton Lake |
| Mica | Combined turbine plus spill | 0 | 10% of mean natural inflow | Structural stress range; study-year Treaty schedule unavailable |
| Waneta | Combined turbine plus spill | 0 | 10% of mean natural inflow | Structural stress range; owner/regulatory rule unresolved |
| Walden North | Cayoosh Creek spill route only | 0 | 10% of mean natural inflow | Never applied to the Seton Lake diversion |

The 10% values are not calibrated estimates. They provide a consistent stress
test where no defensible facility-specific upper rule is available.

## Mathematical implementation

For an absolute bound \(q_s\) in m³/s, the hourly constraint uses

\[
q_h = 3600 q_s.
\]

For an inflow-fraction bound \(f\), the hourly constraint uses

\[
q_h = f\,\overline{I},
\]

where \(\overline{I}\) is the station reservoir's mean 2021 natural inflow in
m³/h. Combined-station rules constrain the sum of turbine and spill discharge.
Named-route rules constrain only the registered spill or environmental-release
link. This distinction prevents a river-flow requirement from being applied to
a power diversion.

## Verification

Gate 2C reconstructs the prepared hydro network for `A_full_cascade`,
`B_two_store`, and `C_single_bucket`, then applies both bounds. All 42
station-level treatments pass. The audit verifies the numerical minimum, route
specificity, scenario evidence label, and preservation of the four
evidence-backed hourly route schedules. Unit and regression tests also verify
unit conversion, mean-inflow scaling, lower/upper ordering, route isolation,
and schedule merging.

Reproduce the audit from the repository root:

```powershell
$env:PYTHONPATH='src;.'
python -m workflow.scripts.studies.chapter2.audit_release_uncertainty
```

Machine-readable outputs are stored in
`studies/chapter2/gate2/release_uncertainty_verification.csv` and `.json`.
