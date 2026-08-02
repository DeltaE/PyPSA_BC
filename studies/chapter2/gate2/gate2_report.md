# Chapter 2 Gate 2 — Physical verification

**Gate status:** PASS
**PyPSA:** 1.2.4
**Solver:** HiGHS
**Tolerance:** 1e-07

## Decision

All required analytic cases pass. The hydraulic mechanics are ready for integration testing with evidence-resolved BC cascade data.

## Results

| Test | Status | Observed | Expected | Unit | Detail |
|---|---|---:|---:|---|---|
| two_reservoir_analytic_test | pass | 5.0 | 5.0 | MWh per one-hour snapshot | Ten water units at efficiency 0.5 yield five electrical units. |
| two_reservoir_analytic_test | pass | 10.0 | 10.0 | m3/h | All turbine discharge reaches the lower reservoir bus. |
| three_reservoir_delay_test | pass | 1.0 | 1.0 | snapshot hour | The second link withdraws water one hour after the first release. |
| three_reservoir_delay_test | pass | 3.0 | 3.0 | snapshot hour | Cumulative one-hour plus two-hour delay delivers water at hour three. |
| minimum_release_test | pass | 5.0 | 5.0 | m3/h | Combined turbine and spill release meets the minimum in every snapshot. |
| scheduled_route_release_test | pass | 0.0 | 0.0 | m3/h | The named river route follows its time-varying minimum profile. |
| scheduled_route_release_test | pass | 0.0 | 0.0 | m3/h | The route-specific schedule does not impose release on the turbine path. |
| turbine_spill_mass_balance_test | pass | 8.0 | 8.0 | m3/h | Four electrical units at 0.5 conversion require eight water units. |
| turbine_spill_mass_balance_test | pass | 2.0 | 2.0 | m3/h | The remaining two water units spill. |
| turbine_spill_mass_balance_test | pass | 10.0 | 10.0 | m3/h | Turbine discharge plus spill equals upstream inflow. |
| delay_boundary_test | pass | 0.0 | 0.0 | m3 | A non-cyclic tail release has no in-horizon downstream delivery. |
| delay_boundary_test | pass | 10.0 | 10.0 | m3/h | A cyclic two-hour delay wraps the final-hour release to hour one. |
| terminal_storage_test | pass | 10.0 | 10.0 | m3 | A non-cyclic terminal Store accumulates all terminal release. |

Machine-readable evidence: `analytic_verification.csv` and `summary.json`.
