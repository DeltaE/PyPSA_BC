# Waneta controllable-storage uncertainty methods

## Evidence and decision

BC Hydro's Waneta transaction application reports a normal maximum forebay
elevation of 462.6 m, a minimum operating elevation of 457.8 m, and hydraulic
balancing with Seven Mile. It does not report active storage volume. Fisheries
and Oceans Canada classifies Waneta as run-of-river, describes small daily
volume and level fluctuations, and states that the reservoir's basic dimensions
were not recorded. These sources support a low-storage interpretation but do not
support a numeric active-volume estimate.

The study therefore brackets Waneta controllable pondage with two deterministic
scenarios:

- `zero_pondage`: no controllable Store, matching the prior run-of-river model;
- `upper_pondage`: 24 hours of 2021 mean natural inflow.

The upper case is a conservative structural stress bound. It is not a measured,
calibrated, or probabilistic storage estimate.

## Implementation

For storage case \(c\), controllable volume is

\[
S_c = h_c\,\overline{I}_{WAN},
\]

where \(h_c\) is 0 or 24 h and \(\overline{I}_{WAN}\) is the mean hourly Waneta
natural inflow in m³/h. With the prepared 2021 input,
\(\overline{I}_{WAN}=999{,}042.36\) m³/h and the upper bound is
23,977,016.71 m³. The upper case adds a cyclic Store at the Waneta water bus
before applying the A/B/C hydraulic transformation. Aggregation therefore
preserves the added volume while changing only hydraulic resolution.

## Verification and claim rule

Gate 2D checks both storage cases across `A_full_cascade`, `B_two_store`, and
`C_single_bucket`. All six representation-case checks pass. The verifier
confirms the storage increment within a 1×10⁻⁶ relative tolerance, unchanged
hourly Waneta inflow, scenario evidence labels, and protocol metadata.

Core aggregation conclusions must retain their direction and material
interpretation across both storage cases. If they do not, Waneta active-volume
evidence becomes a prerequisite for publication rather than a limitation.

Reproduce the gate from the repository root:

```powershell
$env:PYTHONPATH='src;.'
python -m workflow.scripts.studies.chapter2.audit_storage_uncertainty
```

Machine-readable outputs are stored in
`studies/chapter2/gate2/storage_uncertainty_verification.csv` and `.json`.
