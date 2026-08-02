# Chapter 2 publication execution tracker

**Status date:** 2026-08-01 (optimal annual reference and independent audit)  
**Target submission window:** 2027-01-25 to 2027-02-12  
**Primary path:** rigorous specialist energy-systems journal  
**Stretch path:** broader, higher-selectivity journal if validation and decision
consequences are strong  
**Planning assumption:** one lead researcher at 0.7-0.8 full-time equivalent,
with short domain-expert reviews at scientific gates

## 1. Executive decision

The study has moved from planning to controlled implementation. The research
question, estimands, A/B/C configurations, 81-case core scenario matrix, unit
contract, evidence classes, and stop/go rules are frozen in
`config/studies/chapter2_cascade.yaml`.

Gate 1 **passes with no errors and two warnings**. Travel-time coverage
is 20/20, but every delay is explicitly classified as a scenario because the
prior hourly BC study found no station-specific travel-time observations.
Minimum-release execution coverage is 20/20: 13 observed/not-applicable records
and seven scenario-class lower/upper bounds. Seton River, Alouette River, Elk
Falls Canyon, and Cranberry Creek retain source-labelled hourly schedules. The
seven bounds bracket unresolved state-dependent, treaty, or owner-specific
rules; they do not resolve the underlying evidence gaps.
Gate 2 passes all seven required analytic physics cases. Gate 2B passes all
12 actual-input basin-treatment conservation audits. Gate 2C passes all 42
station-treatment checks across both release bounds and all representations.
Gate 2D passes all six Waneta pondage checks. The workflow compiles 30 cases;
24 are runnable and six evidence-only cases remain blocked. The full-horizon
2021 A/lower-bound/zero-pondage reference case is optimal and passes eight
independent annual integrity checks. The full production matrix remains opt-in
until calibration and protocol review are complete. Observed BC Hydro interchange
is the primary external boundary; TTC is deferred to a secondary sensitivity.

The corrected calendar-2021 reference passes the declared calibration-fit
thresholds: total generation is 1.95% low, hydro 2.45% low, and wind 8.19% high.
Monthly biomass and combustible totals match by construction under explicit,
source-labelled calibration constraints. This agreement is not independent
validation. The official named-station capacity gate passes, but the cross-period
energy screen flags Kootenay Canal, Seven Mile, and Bridge River for review;
the economic-ordering audit now explains the Seven Mile and Bridge River spill
preference. A solved ±25% station-energy-band sensitivity restores material
turbine routing at both stations, leaves total hydro unchanged, raises the model
objective by 4.60%, and passes all eight integrity checks. It is a cross-period
calibration sensitivity rather than independent validation. Kootenay Canal and
the absence of held-out operational evidence remain blocking issues.

The paper remains publishable if it makes a narrow claim: hydraulic aggregation
can change operational-feasibility estimates under otherwise matched model and
scenario assumptions. It must describe PyPSA-BC as **delay-capable**, not
delay-calibrated, and it must distinguish observed constraints from scenarios,
approximations, and exclusions.

## 2. Evidence-backed status

| Item | Evidence | Status |
|---|---|---|
| Research question and paired estimands | frozen protocol | Complete |
| A/B/C representation contracts and production builders | protocol plus `src/pypsa_bc/studies/cascade/representation.py` | Complete |
| Full-horizon annual reference | optimal solved network plus independent audit | PASS — 8,760 h, zero backstop, exact observed interchange, release/bounds/closure/balance pass |
| Primary transport-flow reporting | hourly lexicographic minimum-transfer reconstruction | Complete; raw zero-cost circulating flows excluded from congestion diagnostics |
| Baseline generation calibration | official Statistics Canada 2021 monthly generation by type | PASS for calibration fidelity; independent validation not established |
| Named-station validation layer | official BC Hydro fiscal-2021 capacity and energy | Capacity PASS; cross-period energy REVIEW for Kootenay Canal, Seven Mile, and Bridge River |
| Hydro economic-ordering audit | effective turbine VOM versus 8,760 solved local marginal prices | Complete; CAD 6.88/MWh Seven Mile and Bridge River cost exceeds maximum local price in all hours |
| Named-station energy-band sensitivity | optimal full-year network plus integrity, calibration, and station audits | PASS computationally; Seven Mile 2,279.2 GWh, Bridge River 1,664.3 GWh, total hydro unchanged, objective +4.60%; calibration sensitivity only |
| Uniform reservoir-VOM sensitivity | controlled HiGHS simplex and IPM attempts | No exported optimal solution; computational outcome only, not infeasibility evidence |
| Hydro publication methods and figures | scripted economic audit, comparison, chart map, and claim-safe text | Complete for diagnostic/supplement use |
| Calibration evidence pipeline | pinned source archive, hashes, filtered table, scripted audit and figures | Complete and reproducible |
| 2021 asset-vintage filtering | solar, wind, and thermal start/closure years checked during build | Implemented; corrected annual rerun is optimal and passes all eight integrity checks |
| TTC execution dependency | observed-interchange path reads actual tie-line data only | Removed from core; retained only for optional `hourly_ttc` sensitivity |
| Scenario matrix | 30 compiled cases; 24 runnable | Frozen, not yet run beyond the annual reference |
| Unit contract | m3/h water, m3 storage, MW electricity | Complete |
| Input hashes and topology audit | `gate1/` artifacts | Complete |
| Cascade evidence register | `inputs/cascade_evidence_register.csv` | 20/20 delay scenarios; exact minimum-release evidence remains incomplete |
| Hydro inflow unit correction | `workflow/scripts/enrich_format_hydro.py` | Complete |
| Local-time inflow calibration | `src/pypsa_bc/hydro.py` | Complete; all 13 modelled source-backed series reconcile exactly by month |
| La Joie inflow and route repair | processed inflows and hydro dictionaries | Complete |
| Walden water-unit correction | `workflow/scripts/enrich_format_hydro.py` | Complete; observed m3/h inflow now feeds the water bus |
| Duncan reservoir role | official facility evidence and role decision table | Complete; non-generating control reservoir |
| Walter Hardman storage | WUP plan pages 4-5 and storage-decision table | Complete; 459,350 m³ live headpond Store modeled |
| Walden North storage | accepted Bridge River WUP, section 2.2.4.1 | Resolved; explicit no-storage project |
| Waneta storage | BC Hydro/DFO sources plus storage protocol | 0-24 h mean-inflow scenario envelope; not a measured volume |
| Delay assignment | turbine and spill water paths | Complete |
| Constant combined turbine-plus-spill release | Linopy callback | Complete |
| Scheduled route-specific release engine | Linopy callback plus strict CSV loader | Complete; four BC route profiles populated for all 8,760 hours |
| Terminal water and delay boundaries | solved analytic cases | Complete |
| Gate 2 verification | 7 required cases, PyPSA 1.2.4, HiGHS | PASS |
| Gate 2B representation verification | 12 basin-treatment checks across B/C | PASS |
| Gate 2C release-envelope verification | 42 station-treatment checks across A/B/C | PASS |
| Gate 2D storage-envelope verification | 6 representation-case checks across A/B/C | PASS |
| Scenario catalogue | 3 representations × 5 water policies × 2 storage policies | 30 compiled; 24 runnable; evidence-only cases blocked |
| Project regression suite | 98 tests plus 6 subtests | Passing on 2026-08-01 |
| Clean editable installation | Python 3.12 workspace environment | Passing after root README repair |
| Primary-source package | 17/17 retained or downloaded files/pages; 17/17 SHA-256 hashes recorded | Complete |
| Page-level source search | `constraint_source_hits.csv` and `.md` | Complete; manual interpretation required |

## 3. Gate 1 decision

### Remaining scientific limitation

Seven stations lack an exact observed/legal minimum-release representation.
They now have a predeclared lower/upper structural sensitivity, so this gap no
longer produces a Gate 1 error. It still limits inference and blocks the
evidence-only case.
This count combines three distinct scientific tasks:

1. Resolve genuinely unknown rules from primary sources or facility owners.
2. Instantiate verified seasonal schedules on approved facility routes.
3. Represent state- or route-dependent rules on the correct hydraulic path.

A zero in the legacy CODERS/WUP feature table does not prove that no constraint
exists. A station can have:

- no minimum-release condition;
- a constant combined release;
- a seasonal schedule;
- a reservoir-level- or tide-dependent rule;
- a constraint on a bypass/spill route rather than the turbine route; or
- an owner-specific rule outside the BC Hydro WUP corpus.

### Source findings that change the model design

| Facility/system | Source-backed finding | Required modelling decision |
|---|---|---|
| G.M. Shrum | Peace WUP specifies no additional water-management constraint | Mark no constant minimum; retain other reservoir constraints separately |
| Peace Canyon | Peace WUP specifies 283.2 m3/s minimum release | Use constant combined release |
| Strathcona and Ladore | Campbell WUP specifies no minimum discharge requirement | Mark not applicable |
| John Hart | Seasonal and route-specific lower-river/spill rules | Add scheduled, route-specific link constraints |
| Revelstoke | Year-round 5 kcfs (141.58 m3/s; commonly reported as 142 m3/s), with outage qualification | Use a constant rule and document the exception |
| La Joie | Minimum discharge varies with Downton reservoir elevation | Restore the missing generating edge and implement a state-dependent rule or a declared approximation |
| Terzaghi/Bridge | Lower Bridge release uses a separate route and annual/seasonal flow logic | Do not apply the rule to Bridge turbine diversion into Seton Lake |
| Seton | Target release follows a schedule to Seton River | Use a separate downstream-river route |
| Alouette | Low-level outlet base flow plus a seasonal surface release | Split Alouette River release from diversion to Stave |
| Ruskin | Tailwater, tide, and block-loading conditions affect flow | Implement state dependence or a predeclared conservative approximation |
| Seven Mile/Waneta | Coordinated, state-dependent discharge distribution; upstream minimum block can be zero | Do not invent a constant station minimum |

### Remaining warning

| Warning | Required closure |
|---|---|
| Three water-passing stations lack a Store | Confirm zero controllable storage or add storage evidence |

## 4. Scientific strategy

### Primary estimand

For each matched scenario, estimate:

- Configuration B minus Configuration A; and
- Configuration C minus Configuration A.

Report effects for unserved energy, loss-of-load hours, peak shortfall, VRE
curtailment, corridor congestion, operating cost, hydro generation, spill,
reservoir range, and net imports. Do not report p-values without a probability
model. Do not call deterministic backstop dispatch expected unserved energy.

### Claim hierarchy

1. **Core claim:** aggregation changes operational-feasibility estimates in this
   BC case under matched assumptions.
2. **Mechanism claim:** changes arise through altered water timing, routing,
   storage coordination, and binding release constraints.
3. **Robustness claim:** the sign and magnitude of the effect across demand,
   VRE, hydrology, delay, terminal-storage, and release-rule sensitivities.
4. **Generalization boundary:** results do not establish a universal bias for all
   hydro systems or probabilistic reliability.

### Required validation

Use independent observations not used to calibrate inflows:

- provincial/system load;
- station or facility-group generation;
- reservoir levels or storage;
- river releases/flows where public;
- annual/seasonal energy and water-balance checks.

Freeze thresholds before viewing comparative A/B/C results. If independent
validation is weak, narrow the claim and target a specialist journal.

### Minimum publication figure set

1. BC network and cascade-system map with explicit hydraulic routes.
2. A/B/C schematic showing exactly what aggregation removes.
3. Validation small multiples for load, generation, storage/level, and flow.
4. Paired effect plot across the 27 matched scenario families.
5. Stress-event chronology linking reservoir state, releases, congestion, and
   shortfall.
6. Sensitivity plot for 0, 1, 3, 6, and 12 hour routing.
7. Sensitivity plot for release-rule approximation and terminal storage.
8. Water and electricity balance diagnostics in the supplement.

## 5. Dated execution map

| Week | Dates | Primary work | Exit criterion |
|---:|---|---|---|
| 1 | Aug 3-9 | Approve station-by-station rule crosswalk; resolve owner/facility names | Every station has a rule type, source, locator, and decision |
| 2 | Aug 10-16 | Restore La Joie; classify Duncan; split diversion, river, turbine, and spill routes | No unexplained missing edge or misdirected water |
| 3 | Aug 17-23 | Implement time-varying and route-specific schedules | John Hart, Alouette, Bridge/Seton, and similar rules have explicit paths |
| 4 | Aug 24-30 | Implement state-dependent rules or predeclare justified approximations | Gate 1 has no publication-critical unknown |
| 5 | Aug 31-Sep 6 | Add analytic tests for scheduled, route-specific, and state-dependent releases | Gate 2 remains PASS |
| 6 | Sep 7-13 | Implement Configuration B | Hydraulic and electrical equivalence report passes |
| 7 | Sep 14-20 | Implement Configuration C and one-bus sensitivity | Gate 3 passes |
| 8 | Sep 21-27 | Prepare independent validation observations | Calibration and validation roles are separated |
| 9 | Sep 28-Oct 4 | Validate Configuration A and diagnose mismatches | Gate 4 passes or claim is narrowed |
| 10 | Oct 5-11 | Freeze solver, assumptions, manifests; run A/B/C pilots | Pilot runs reproduce and balances close |
| 11 | Oct 12-18 | Run first half of the 81 core solves | Automated run QA passes |
| 12 | Oct 19-25 | Complete the 81 core solves | Matched result cube complete |
| 13 | Oct 26-Nov 1 | Run delay, release, terminal-storage, and inflow sensitivities | Critical assumption bounds quantified |
| 14 | Nov 2-8 | Run bus aggregation and selected ablations | Hydraulic and electrical effects separated |
| 15 | Nov 9-15 | Compute paired effects and event attribution | Frozen analysis tables produced by scripts |
| 16 | Nov 16-22 | Produce publication figures and supplement | Figure QA and source crosswalk pass |
| 17 | Nov 23-29 | Replace manuscript placeholders; write Methods and Results | Complete 7,000-word internal manuscript |
| 18 | Nov 30-Dec 6 | Internal scientific review | No unsupported claim in the crosswalk |
| 19 | Dec 7-20 | Independent domain review and response | Hydraulic interpretation approved or qualified |
| 20 | Dec 21-Jan 10 | Format journal package, supplement, code/data statement, and cover letter | Gate 6 checklist passes |
| 21 | Jan 11-24 | Clean-environment reproducibility rerun and final proof | Tables and figures reproduce from archive |
| 22 | Jan 25-Feb 12 | Submit | Immutable submission archive created |

## 6. Time and effort estimate

| Workstream | Lead effort | Reviewer effort | Main risk |
|---|---:|---:|---|
| Gate 1 evidence and topology closure | 10-15 days | 6-10 hours | Owner-specific rules and route interpretation |
| Scheduled/state-dependent constraints and verification | 10-15 days | 4-8 hours | Model redesign for route/state dependence |
| B/C aggregation and equivalence | 8-12 days | 3-5 hours | Confounded electrical/hydraulic aggregation |
| Independent validation and Gate 4 | 10-15 days | 8-12 hours | Data access and calibration leakage |
| Core runs and sensitivities | 10-15 days plus compute | 2-4 hours | Failed solves and compute limits |
| Analysis, figures, manuscript, supplement | 18-25 days | 15-25 hours | Review cycles and reproducibility defects |

The January-February 2027 window is feasible if Gate 1 passes by August 30 and
validation observations are ready by September 21. A delay in either milestone
moves submission to March-May 2027.

## 7. Stop/go rules

- No core scenario production before Gate 1 and Gate 2 pass.
- No calibration target may also serve as independent validation.
- No manual manuscript number or figure; every result must resolve to a
  versioned artifact and run manifest.
- No aggregation comparison unless capacity, inflow, storage, hydraulic routes,
  electrical geography, and terminal conditions are matched or declared as an
  ablation.
- No claim of calibrated travel delays.
- No constant minimum-flow parameter for a seasonal, route-specific, or
  state-dependent rule unless the approximation and sensitivity are explicit.
- If Gate 4 is weak, narrow the claim and journal target.

## 8. Immediate next actions

1. Obtain calendar-2021 station generation, reservoir, or release observations
   that were not used to construct the hydro sensitivity.
2. Audit and correct the Kootenay Canal fixed source profile against a
   period-compatible station record.
3. Obtain the operative Terzaghi variance/order and the study-year Mica Treaty
   operating schedule.
4. Obtain owner/regulatory evidence for Arrow Lakes and Waneta.
5. Approve a Downton stage-storage mapping and implement the La Joie piecewise
   release rule.
6. Choose and predeclare an exact mixed-integer or conservative sensitivity
   formulation for Ruskin tailwater and block loading.
7. Test both predeclared Waneta storage bounds in every primary comparison.
   Walden North is evidence-backed as a no-storage project. Walter Hardman's documented
   459,350 m³ live headpond storage is now modeled explicitly.
8. Expand Gate 2 with solved state-dependent release cases.
9. Run the first matched A/B/C comparison only after the held-out station layer
   and Kootenay Canal review pass.

## 9. Evidence locations

- Protocol: `config/studies/chapter2_cascade.yaml`
- Gate 1 decision: `studies/chapter2/gate1/gate1_report.md`
- Gate 2 decision: `studies/chapter2/gate2/gate2_report.md`
- Gate 2B decision: `studies/chapter2/gate2/representation_verification.md`
- A/B/C methods and schematic: `studies/chapter2/hydraulic_representation_methods.md`
- Hydro dispatch Methods, Results, and decision schematic:
  `studies/chapter2/hydro_dispatch_sensitivity_methods_results_2026-08-01.md`
- Publication figure chart map:
  `studies/chapter2/publication_figure_chart_map_2026-08-01.md`
- Hydro economic-ordering outputs:
  `results/workflow/sensitivities/hydro_dispatch_economics/`
- Hydro station-band sensitivity outputs:
  `results/workflow/sensitivities/station_energy_band_25pct/`
- Evidence register: `studies/chapter2/inputs/cascade_evidence_register.csv`
- Source registry: `data/validation/chapter2/source_registry.yaml`
- Download/hash manifest: `data/validation/chapter2/source_manifest.csv`
- Page-level search report: `data/validation/chapter2/constraint_source_hits.md`
- Evidence decisions: `data/validation/chapter2/cascade_parameter_decisions.yaml`
- Full publication plan: `studies/chapter2_publication_strategy_and_timeline_report.md`
