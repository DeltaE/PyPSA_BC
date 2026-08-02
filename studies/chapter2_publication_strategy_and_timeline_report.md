# Publication Strategy, Scientific Quality Plan, and Delivery Timeline

## Hydropower Cascade Representation and Operational Feasibility in BC-PyPSA

**Prepared:** 28 July 2026  
**Project directory:** `E:\CoWork\PROJECTS\PyPSA`  
**Primary manuscript:** `E:\CoWork\PROJECTS\QE\CH_02\2026 07 01 - Chapter2_Cascade_Draft_EL_TN_EL.docx`  
**Planning horizon:** 20 working weeks for a strong submission, subject to validation-data access  
**Status:** implementation and evidence closure; analytical comparison results remain to be produced

> **Evidence update - 29 July 2026:** Gate 2 now passes all six required
> solved physics cases. Gate 1 has improved to one error and four warnings.
> Travel time is covered for 19/19 stations as an explicit structural scenario,
> not as calibrated station evidence. Constant or explicitly absent minimum
> releases are resolved for 5/19 stations; four additional facilities have
> verified seasonal or state-dependent rules awaiting route-aware
> implementation. The defensible submission window is now 25 January to
> 12 February 2027, with March-May 2027 as contingency. The live scientific
> decision and dated 22-week schedule are maintained in
> `studies/chapter2/publication_execution_tracker.md`.

---

## Executive summary

The proposed research is scientifically viable and potentially publishable because it addresses a precise and unresolved question: how much does aggregating a hydropower cascade change estimates of hourly operational feasibility when demand, variable renewable energy, transmission, hydrology, and all other system assumptions remain fixed? British Columbia provides a demanding case because large regulated cascades supply most provincial electricity while electrification and renewable integration increase the value of hourly flexibility.

The existing research and modelling plans cover the major elements required for a rigorous paper:

- a station-resolved reference cascade;
- controlled two-store and single-bucket simplifications;
- explicit water, storage, turbine, spill, and electrical accounting;
- operational travel-time and minimum-flow constraints;
- matched demand, renewable, and hydrological scenarios;
- independent structural, numerical, operational, and system validation;
- mechanism ablations and sensitivity tests;
- reproducible run manifests and common metric definitions;
- publication-quality maps, event plots, effect-size figures, and supplementary evidence.

The strategy is sufficiently mapped, and the first implementation gates now
produce machine-readable evidence. The model assigns delay to downstream water
paths and enforces a constant combined turbine-plus-spill minimum release. Six
solved analytic physics cases pass. The remaining implementation risk is
facility-specific: several accepted Water Use Plans define seasonal,
route-specific, reservoir-level-dependent, tide-dependent, or block-loading
rules that cannot be reduced honestly to one constant station value. The second
critical risk is access to credible independent reservoir, flow, and generation
records. Without independent operational validation, the study can demonstrate
structural sensitivity but should not claim measured real-world bias.

For one researcher working mainly on this project, the current evidence supports
a **22-week path to submission**, targeting 25 January to 12 February 2027.
March-May 2027 is the contingency window if state-dependent rule
implementation or validation-data access slips. A broader, higher-selectivity
submission would still require multiple weather years, stronger uncertainty
treatment, and evidence that cascade representation changes a planning
decision.

The recommended baseline is a 22-week programme with six formal scientific gates:

1. input and unit integrity;
2. verified cascade physics;
3. configuration equivalence;
4. independent validation;
5. decision-relevant and robust comparative results;
6. clean-environment reproducibility and submission readiness.

The research should not advance automatically from one gate to the next. Failed gates require correction or a narrower claim.

---

## 1. Research objective and publication claim

### 1.1 Primary research question

How much does hydropower cascade aggregation alter estimates of operational feasibility in a hydro-dominant electricity system under electrification, variable renewable energy integration, and hydrological stress?

### 1.2 Primary scientific estimand

For each matched scenario \(s\), estimate the difference between the simplified representation and the station-resolved reference:

\[
\Delta Y_{C,s}=Y_{\text{single bucket},s}-Y_{\text{full cascade},s}
\]

and:

\[
\Delta Y_{B,s}=Y_{\text{two store},s}-Y_{\text{full cascade},s}.
\]

The experiment must hold demand, renewable portfolios, hydrology, conventional generation, transmission ratings, trade assumptions, costs, time horizon, and solver settings constant within each matched scenario.

### 1.3 Primary outcomes

The paper should predeclare four primary operational outcomes:

1. unserved energy in MWh per year;
2. loss-of-load hours and maximum hourly shortfall in MW;
3. wind and solar curtailment in GWh and percentage points;
4. transmission congestion by corridor, measured as corridor-hours and unique congested hours.

Secondary outcomes should include:

- total operating cost;
- generation by cascade and station;
- turbine utilization;
- spill;
- reservoir storage range;
- imports and exports;
- binding minimum-flow and storage constraints;
- solver runtime and peak memory.

Predeclaring primary outcomes reduces selective reporting. The final manuscript should report all planned primary outcomes, including null or contradictory findings.

### 1.4 Claim hierarchy

The study must distinguish three levels of evidence:

**Structural effect:** Configurations A, B, and C produce different results under matched inputs.

**Mechanism:** exported constraint and event data show why the configurations diverge, for example delayed downstream arrival or binding minimum flow.

**Accuracy or bias:** independent operational observations support Configuration A more strongly than the simplified alternatives.

Only the third level justifies describing the difference as “aggregation bias.” If operational observations remain incomplete, the paper should use “structural sensitivity to hydropower representation.”

---

## 2. Publication-readiness assessment

### 2.1 Strengths

The research has six clear strengths.

First, it asks a falsifiable question. The experiment can show a material difference, a negligible difference, or a difference limited to specific stress conditions. Each outcome contributes useful evidence.

Second, it uses a controlled comparison. Many regional model papers compare scenarios that change several assumptions simultaneously. The proposed A/B/C design can isolate hydropower representation if the aggregation-equivalence checks succeed.

Third, the formulation connects hydrology and power-system operation. Water availability, storage, delayed transfer, turbine conversion, spill, transmission, and demand interact in a single chronological optimization.

Fourth, the case is policy relevant. BC depends heavily on regulated hydropower and faces increasing electricity demand, renewable integration, and transmission investment pressure.

Fifth, the contribution can travel beyond BC. The reusable outputs are the topology schema, constraint layer, validation ladder, controlled comparison, quality checks, and decision rules.

Sixth, the project can release open code and reproducible derived inputs. Reproducibility would distinguish the paper from studies that depend entirely on proprietary production-cost models.

### 2.2 Current limitations

The present research package has four material limitations.

1. Operational travel delay and minimum-flow constraints are planned but not yet implemented in the active network workflow.
2. Independent hourly reservoir, discharge, spill, and facility-generation observations may be incomplete.
3. The current manuscript contains analytical placeholders rather than completed results.
4. A single weather or inflow year would not establish robustness to hydrological conditions.

These limitations are manageable, but they define the project’s critical path.

### 2.3 Journal positioning

The journal choice should follow the final contribution rather than precede it.

| Final contribution | Best-fit journal class | Positioning |
|---|---|---|
| Strong BC case study with transparent model comparison | Energy 360; Energy Reports | Broad energy-system and early-career contribution |
| New cascade formulation with detailed power-system verification | Electric Power Systems Research; IJEPES | Engineering, optimization, and reliability |
| Representation changes strategic planning conclusions | Energy Strategy Reviews; Smart Energy | Energy-system modelling and planning |
| Multi-year validation plus generalizable and decision-changing results | Applied Energy | High-impact applied energy-system contribution |
| Renewable integration and hydropower flexibility dominate the paper | Renewable Energy Focus; Renewable Energy | Renewable-system integration |

Energy 360 is a legitimate and realistic initial target. It has a broad energy scope and explicitly supports early-career researchers. The study could aim higher if validation covers both major cascades and the results demonstrate planning consequences rather than small dispatch differences.

---

## 3. Scientific quality framework

### 3.1 Verification versus validation

The study must treat verification and validation as different activities.

**Verification** asks whether the code solves the equations intended by the researcher. Synthetic cascade tests, water-balance residuals, constraint inspection, and A/B/C equivalence checks provide verification.

**Validation** asks whether modelled quantities resemble independent observations. Facility generation, reservoir levels, provincial load, trade, and other system records provide validation.

A model can be verified but poorly validated. Conversely, a model can match one aggregate statistic while violating water conservation. Publication requires both.

### 3.2 Minimum verification evidence

The following evidence is mandatory:

- all cascade topology edges resolve to valid upstream and downstream components;
- no unexplained cycles, self-loops, orphan Stores, or duplicated physical reservoirs;
- flow, storage, and electrical units pass explicit unit tests;
- two- and three-reservoir synthetic systems reproduce known delayed-flow solutions;
- minimum flow binds correctly in a controlled low-demand test;
- turbine and spill pathways conserve water;
- full-network hourly and annual water-balance residuals remain within a predeclared tolerance;
- electrical balance residuals remain within solver tolerance;
- terminal storage behaves according to the selected cyclic, observed, or water-value condition;
- aggregated configurations preserve intended totals.

Numerical tolerance should be scaled to the problem. A normalized residual close to the solver feasibility tolerance is preferable to a fixed absolute threshold that becomes meaningless for large reservoirs. The report should publish both absolute and normalized residuals.

### 3.3 Configuration-equivalence evidence

Before solving, the workflow should compare A, B, and C for:

- total natural inflow by hour and year;
- total usable storage;
- installed turbine capacity;
- maximum discharge;
- water-to-energy conversion assumptions;
- initial and terminal storage treatment;
- electrical capacity by bus;
- renewable, load, trade, and transmission inputs.

Intended totals should match to floating-point tolerance. Any planned difference must appear in a configuration table. The single energy bucket should initially retain the original generating-station electrical buses. A separate one-bus sensitivity may test electrical aggregation, but it must not be part of the primary hydraulic comparison.

### 3.4 Minimum validation coverage

The minimum defensible validation package should include:

- provincial load magnitude, peak, timing, monthly energy, and load-duration shape;
- calibrated inflow checks for every reservoir with monthly reference statistics;
- independent generation or storage observations for both Peace and Columbia;
- at least five major facilities or facility groups, where data permit;
- validation assets representing a substantial share of cascaded hydro capacity;
- annual system energy and trade checks;
- clear separation between calibration and independent validation data.

A stronger paper would validate more than half of cascaded generation capacity and include storage trajectories for the major head reservoirs. Coverage should be reported directly rather than hidden inside average error metrics.

### 3.5 Suggested validation metrics

| Quantity | Core metrics | Interpretation |
|---|---|---|
| Provincial load | annual bias, monthly bias, peak error, RMSE, MAE, correlation, duration-curve error | Input and time-alignment quality |
| Natural inflow | monthly bias, seasonal timing, KGE or NSE if independent daily data exist | Hydrological timing and magnitude |
| Facility generation | annual/seasonal energy bias, NRMSE, correlation, ramp distribution, peak coincidence | Operational performance |
| Reservoir storage or level | normalized RMSE, range bias, timing of seasonal maximum/minimum | Storage dynamics |
| Water balance | maximum and p99 residual, annual closure | Physical verification |
| Transmission | corridor loading patterns, known interface comparison, congestion distribution | Network behavior |

No universal error threshold guarantees validity. Tentative targets may include annual facility-generation bias below 5–10% and correct seasonal timing, but the final acceptance criteria must reflect observation quality and must be set before inspecting A/B/C outcomes.

### 3.6 Robustness evidence

The minimum robustness matrix should cross:

- three hydro representations;
- three demand conditions;
- three VRE portfolios;
- three hydrological conditions.

This produces 81 core solves. If computation limits the full factorial design, the final retained matrix and selection logic must be documented before comparative analysis.

Secondary mechanism ablations should include:

- Configuration A without travel delay;
- Configuration A without minimum flow;
- Configuration A with instantaneous coupling;
- Configuration A with relaxed seasonal storage bounds;
- alternative terminal condition.

Sensitivity tests should examine:

- low/core/high travel time;
- missing downstream local inflow;
- water-to-power conversion efficiency;
- initial storage;
- environmental-flow coverage;
- spill limits;
- value of lost load;
- solver and tolerance.

Deterministic scenarios are not independent statistical observations. The paper should emphasize paired effect sizes, ranges, consistency of sign, and stress-condition dependence. It should avoid p-values unless a defensible stochastic sample or probabilistic ensemble is introduced.

### 3.7 Reliability terminology

Backstop dispatch in one deterministic year measures unserved energy. “Expected energy not served” should be reserved for a probability-weighted ensemble covering uncertain outages, hydrology, demand, or forecast conditions. Precise terminology will prevent a predictable reviewer criticism.

---

## 4. Modelling strategy

### 4.1 Configuration A: full cascade

Configuration A should represent:

- one Store per physical reservoir;
- incremental natural inflow assigned without upstream double counting;
- one turbine pathway per station or documented aggregate station;
- electrical output at the station’s original network bus;
- downstream water output from the same turbine discharge;
- parallel spill transfer;
- explicit integer-hour travel time;
- minimum combined turbine-and-spill release where supported;
- physical storage bounds and evidence-supported seasonal rules;
- terminal storage treatment.

Travel delay should be implemented through explicit linear sending and receiving constraints rather than an assumed software attribute until the installed PyPSA version and serialization behavior are confirmed.

The first delayed hours require a boundary convention. The preferred method is a warm-up period longer than the maximum travel delay, followed by extraction of the analysis year. Cyclic wrapping should be a sensitivity because it moves end-of-year releases into the first hours.

### 4.2 Configuration B: two-store approximation

Configuration B should aggregate each major basin to upper and lower Stores while preserving:

- total storage;
- total natural inflow;
- installed turbine capacity;
- original electrical output buses where feasible;
- one upper-to-lower water transfer.

It should remove station-level delays and reach-specific operating rules unless a valid aggregate form exists. Its purpose is practical: determine whether limited disaggregation recovers most of Configuration A at lower data and computational cost.

### 4.3 Configuration C: single energy bucket

Configuration C should contain one basin-level Store but preserve station electrical delivery locations in the primary comparison. Stations may draw from the common Store through individual output Links. This removes hydraulic ordering without moving generation around the electrical network.

A separate `C_bucket_one_bus` sensitivity can test the compound effect of hydraulic and electrical aggregation.

### 4.4 Scenario design

The core experiment should use matched combinations:

| Factor | Minimum levels | Required definition |
|---|---|---|
| Demand | baseline, medium, high electrification | annual TWh, peak GW, spatial distribution |
| VRE | low, medium, high | wind/solar MW, available TWh, share of load |
| Hydrology | dry, median, wet | historical years or documented chronological transformations |
| Hydro model | A, B, C | exact component and constraint switches |

The study should use chronological hydrological conditions rather than uniformly multiplying every hour whenever suitable historical years exist. Uniform scalars preserve neither drought timing nor seasonal correlation.

### 4.5 Event-based causal analysis

Annual metrics alone cannot demonstrate mechanism. The workflow should select events using predeclared rules:

1. maximum hourly A–C unserved-energy difference;
2. window with maximum cumulative curtailment difference;
3. corridor-hour with maximum congestion redistribution.

For each event, export:

- demand and residual load;
- renewable availability and dispatch;
- station generation;
- turbine discharge and spill;
- upstream and downstream storage;
- delayed receiving flow;
- minimum-flow bounds;
- imports and exports;
- line loading;
- relevant dual values when reliable.

The Results section should trace at least one event from upstream storage through river transfer to electricity-system outcome.

---

## 5. Work packages and deliverables

### WP0 — Study governance and freeze

**Objective:** establish a controlled and auditable research baseline.

**Tasks:**

- freeze the manuscript hypothesis and primary outcomes;
- create the study configuration;
- record package and solver versions;
- hash model inputs;
- define run naming and immutable outputs;
- create an issue and decision log.

**Deliverables:**

- study protocol;
- scenario configuration;
- environment record;
- input manifest;
- initial risk register.

**Effort:** approximately one week.

### WP1 — Input, unit, and topology audit

**Objective:** establish trustworthy cascade inputs.

**Tasks:**

- reconcile generation, reservoir, cascade, and inflow IDs;
- verify units and conversions;
- detect missing, duplicated, and excluded assets;
- build the canonical directed cascade graph;
- classify every parameter as observed, derived, calibrated, imputed, scenario, or unknown.

**Deliverables:**

- input-audit table;
- topology edge table;
- parameter-coverage matrix;
- exclusions table;
- unit tests.

**Effort:** one to two weeks.

### WP2 — Physics implementation and synthetic verification

**Objective:** implement operational travel delay, minimum flow, and boundary treatment.

**Tasks:**

- construct two- and three-reservoir analytic networks;
- implement delayed sending/receiving constraints;
- implement minimum combined release;
- verify spill and turbine water accounting;
- test warm-up and terminal storage;
- export independent balance diagnostics.

**Deliverables:**

- constraint module;
- synthetic test suite;
- expected-solution fixtures;
- verification report.

**Effort:** three to four weeks.

### WP3 — Alternative hydro configurations

**Objective:** build scientifically controlled B and C configurations.

**Tasks:**

- define reservoir-to-aggregate mappings;
- build upper/lower Stores;
- build the common bucket;
- preserve electrical geography;
- produce aggregation-equivalence reports;
- test pilot networks.

**Deliverables:**

- A/B/C builders;
- mapping tables;
- equivalence report;
- pilot matched triplet.

**Effort:** two to three weeks.

### WP4 — Validation-data preparation

**Objective:** assemble independent evidence.

**Tasks:**

- catalogue public and restricted sources;
- fetch public records reproducibly;
- request restricted records early;
- harmonize identifiers, time zones, and units;
- document every correction;
- separate calibration and validation data.

**Deliverables:**

- validation-data manifest;
- fetch scripts;
- patch log;
- station coverage matrix;
- clean validation tables.

**Effort:** three to six weeks, partly parallel with WP2 and WP3.

### WP5 — Baseline and operational validation

**Objective:** determine whether Configuration A provides a defensible reference.

**Tasks:**

- validate load;
- validate inflow;
- compare facility generation;
- compare storage/elevation;
- check annual energy and trade;
- explain mismatches and revise only through documented changes.

**Deliverables:**

- validation metric tables;
- facility plots;
- storage plots;
- acceptance-gate report.

**Effort:** two to four weeks.

### WP6 — Scenario execution

**Objective:** solve the locked matched experiment.

**Tasks:**

- benchmark runtime and memory;
- run baseline A/B/C triplets;
- run the core scenario matrix;
- preserve failed runs and solver logs;
- verify every solved network automatically.

**Deliverables:**

- immutable solved networks;
- run manifests;
- solve-quality table;
- failed-run log.

**Effort:** two to three weeks of wall-clock time after workflow stabilization.

### WP7 — Ablations and sensitivity

**Objective:** identify mechanisms and test robustness.

**Tasks:**

- run one-at-a-time ablations;
- run hydrological, efficiency, storage, delay, and policy sensitivities;
- test at least one alternate solver or tolerance configuration;
- measure computational cost.

**Deliverables:**

- ablation metrics;
- sensitivity metrics;
- robustness summary;
- accuracy–complexity table.

**Effort:** two to three weeks.

### WP8 — Analysis and visualization

**Objective:** translate solved networks into reproducible evidence.

**Tasks:**

- calculate primary and secondary metrics from explicit component registries;
- calculate paired effects;
- select event windows deterministically;
- produce static vector and interactive figures;
- generate complete manuscript and supplementary tables.

**Deliverables:**

- tidy metric tables;
- event exports;
- manuscript figures;
- supplementary figures;
- figure-source data.

**Effort:** two to three weeks.

### WP9 — Manuscript completion and submission

**Objective:** convert the analytical package into a journal submission.

**Tasks:**

- replace analytical placeholders;
- rewrite Abstract last with quantitative results;
- separate Results from interpretation;
- revise title according to evidence;
- complete data/code availability and contribution statements;
- conduct internal scientific and prose review;
- prepare journal-specific files and cover letter.

**Deliverables:**

- clean manuscript;
- supplementary information;
- data and code archive;
- cover letter;
- submission checklist.

**Effort:** three to four weeks, overlapping the final analysis.

---

## 6. Baseline 20-week schedule

### 6.1 Week-by-week plan

| Week | Primary work | Parallel work | Exit criterion |
|---:|---|---|---|
| 1 | Freeze research protocol, outcomes, scenario schema, environment | Start validation-data catalogue and requests | Approved protocol and manifest format |
| 2 | Complete input, unit, and topology audit | Clean available validation records | No unresolved core topology or unit defects |
| 3 | Build two-reservoir analytic test | Continue validation acquisition | Exact no-delay water-balance solution |
| 4 | Build delayed three-reservoir test | Map facility and reservoir IDs | Delayed arrival reproduced exactly |
| 5 | Implement operational delay in PyPSA/Linopy | Prepare initial/storage observations | Full-network delay constraints compile and solve |
| 6 | Implement minimum flow, warm-up, and terminal treatment | Prepare validation metric scripts | Physics verification gate passes |
| 7 | Build Configuration B | Load/inflow validation | B totals reconcile with A |
| 8 | Build Configuration C and electrical-location control | Generation/storage validation | A/B/C equivalence and pilot triplet pass |
| 9 | Diagnose pilot event and finalize metric registry | Complete validation coverage matrix | Mechanism observable or hypothesis revised |
| 10 | Validate Configuration A across available facilities | Lock acceptance criteria | Reference validation report complete |
| 11 | Resolve documented validation issues | Freeze demand/VRE/hydrology definitions | Validation gate passes or claims narrowed |
| 12 | Benchmark runtime and run baseline triplets | Prepare automated quality checks | Stable run orchestration |
| 13 | Run first half of core matrix | Inspect solve diagnostics | Complete matched triplets only |
| 14 | Run second half of core matrix | Begin metric aggregation | Core scenario matrix complete |
| 15 | Run mechanism ablations | Generate event candidates | Main mechanisms attributed |
| 16 | Run hydrology and parameter sensitivities | Calculate runtime/memory metrics | Robustness package complete |
| 17 | Finalize paired metrics and causal event traces | Draft Results structure | Reproducible primary result table |
| 18 | Produce maps, plots, tables, and supplementary evidence | Internal figure review | Publication-quality figure set |
| 19 | Complete Results, Discussion, Abstract, and Conclusion | Reproducibility audit | Complete internal manuscript |
| 20 | Scientific review, journal formatting, cover letter, archive | Final checklist | Submission-ready package |

### 6.2 Timeline map

| Work package | 1–2 | 3–4 | 5–6 | 7–8 | 9–10 | 11–12 | 13–14 | 15–16 | 17–18 | 19–20 |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| WP0 Governance | ● |  |  |  |  |  |  |  |  | ● |
| WP1 Audit | ● |  |  |  |  |  |  |  |  |  |
| WP2 Physics |  | ● | ● |  |  |  |  |  |  |  |
| WP3 A/B/C builders |  |  |  | ● |  |  |  |  |  |  |
| WP4 Validation data | ● | ● | ● | ● | ● |  |  |  |  |  |
| WP5 Validation |  |  |  |  | ● | ● |  |  |  |  |
| WP6 Core runs |  |  |  |  |  | ● | ● |  |  |  |
| WP7 Robustness |  |  |  |  |  |  |  | ● |  |  |
| WP8 Analysis/figures |  |  |  |  |  |  | ● | ● | ● |  |
| WP9 Manuscript |  |  |  |  |  |  |  |  | ● | ● |

### 6.3 Alternative schedules

**Accelerated 12–16 week submission:** possible if validation data are already accessible, synthetic implementation succeeds quickly, and the study targets Energy 360 with three hydrological conditions rather than an extensive multi-year ensemble.

**Strong 20–26 week submission:** recommended for IJEPES, Electric Power Systems Research, Energy Strategy Reviews, or a stronger Energy 360 paper. Includes broader operational validation and full ablations.

**Six-to-nine-month stretch programme:** suitable for Applied Energy. Adds multiple independent weather years, stronger uncertainty, broader facility validation, and a capacity/planning consequence.

---

## 7. Critical path and dependencies

The critical path is:

```text
unit/topology audit
    → verified delayed cascade
    → equivalent A/B/C configurations
    → validated Configuration A
    → locked scenario matrix
    → core solves
    → ablations and sensitivities
    → paired effects and event evidence
    → completed manuscript
```

Validation-data acquisition should begin in Week 1 because external access can delay the entire programme. Public-data preparation can proceed in parallel with model development. Restricted-data delay should not stop synthetic verification or configuration construction.

The full scenario schedule cannot be estimated confidently until one 8,760-hour matched triplet is benchmarked. If one solve requires one to three CPU-hours, 100–120 core and sensitivity solves require roughly 100–360 CPU-hours before failed-run repetition. Parallel execution may reduce wall-clock time, but workflow debugging and storage management remain necessary.

---

## 8. Scientific gates and decision rules

### Gate 1 — Input integrity

**Pass when:**

- unit tests pass;
- cascade topology resolves;
- exclusions are declared;
- parameter evidence classes are complete;
- no unexplained duplicate reservoir or inflow exists.

**Failure response:** correct inputs before modifying optimization behavior.

### Gate 2 — Physics verification

**Pass when:**

- analytic delayed-flow tests pass;
- minimum-flow behavior is correct;
- turbine and spill conserve water;
- warm-up and terminal conventions are verified;
- full-network residuals remain within tolerance.

**Failure response:** do not build the comparative scenario matrix.

### Gate 3 — Configuration equivalence

**Pass when:**

- A/B/C preserve all intended aggregate totals;
- electrical geography remains controlled;
- one baseline matched triplet solves;
- quality checks pass for every configuration.

**Failure response:** correct mappings or redefine the estimand.

### Gate 4 — Reference validation

**Pass when:**

- load and inflow transformations are validated;
- operational evidence covers both major cascades where possible;
- mismatches are understood and documented;
- Configuration A is defensible as a reference.

**Failure response:** narrow the claim to structural sensitivity or obtain additional evidence.

### Gate 5 — Scientific contribution

**Pass for an operational-effect paper when:**

- A–C differences are measurable;
- event traces demonstrate mechanism;
- results persist under relevant sensitivities;
- at least one planning or operational interpretation changes.

**Alternative outcome:** if differences remain negligible, publish a domain-of-validity result showing when aggregation is adequate. A well-verified null finding can still be valuable.

### Gate 6 — Submission readiness

**Pass when:**

- every planned primary outcome is reported;
- manuscript values match versioned tables;
- figures have archived source data;
- code and public inputs are reproducible;
- limitations match the evidence;
- no analytical placeholder remains.

---

## 9. Risk register

| Risk | Probability | Impact | Mitigation | Possible schedule effect |
|---|---|---|---|---|
| Facility generation/storage data unavailable | High | High | Request early; use public WUP records; publish coverage matrix; narrow claims | +4–12 weeks |
| Water-unit inconsistency | Medium | High | Establish explicit unit contract and analytic tests before refactoring | +1–3 weeks |
| Delay constraints difficult to integrate | Medium | High | Use small Linopy examples; isolate variables; verify serialization | +2–4 weeks |
| A/B/C totals do not reconcile | Medium | High | Automated equivalence report; preserve station electrical buses | +1–2 weeks |
| Full-year solves too slow | Medium | Medium | Benchmark pilot; parallel runs; limit sensitivities by design | +1–4 weeks |
| Some scenarios fail numerically | Medium | Medium | preserve logs; use staged stress levels; solver/tolerance checks | +1–3 weeks |
| Aggregation effect is negligible | Medium | Medium | reposition as adequacy/domain-of-validity result | no necessary delay |
| Aggregation effect depends on one extreme event | Medium | Medium | event analysis and multi-year hydrology | +2–6 weeks |
| Validation data used for calibration | Medium | High | preassign calibration and validation roles; held-out periods | +1–2 weeks |
| Electrical aggregation confounds hydraulic effect | Medium | High | retain original electrical buses in Configuration C | +1 week |
| Manuscript overclaims “bias” | Medium | High | use evidence hierarchy and neutral title until validation | no schedule effect |
| Input or code changes during study | Medium | High | immutable manifests, hashes, and run directories | +1 week setup |

---

## 10. Staffing and effort

### 10.1 Solo-researcher model

A focused 20-week programme represents approximately 700–900 research hours, including coding, data preparation, analysis, figures, and manuscript revision. The upper range applies when validation records require manual extraction or identifiers require extensive reconciliation.

### 10.2 Recommended support

Even if one researcher performs most work, three external reviews would improve quality:

- a hydrology or reservoir-operations reviewer for flows, storage, and operating-rule interpretation;
- a power-system optimization reviewer for PyPSA/Linopy constraints and reliability terminology;
- a domain or policy reviewer for BC-specific interpretation and planning implications.

A co-author or research assistant responsible for validation-data preparation could shorten the calendar schedule by three to five weeks.

### 10.3 Meeting rhythm

Recommended governance:

- weekly 30-minute progress and blocker review;
- gate review at Weeks 2, 6, 8, 11, 16, and 20;
- decision log updated whenever a parameter, constraint, or scenario definition changes;
- manuscript result table generated only from versioned metrics.

---

## 11. Reproducibility and data governance

Every run should record:

- run ID;
- source revision or checksum;
- environment and package versions;
- solver and options;
- input checksums;
- hydro configuration;
- demand, VRE, and hydrology levels;
- active constraint switches;
- solve status;
- objective and runtime;
- network path;
- quality-check results.

Public validation data should be fetched through scripts and stored under `data/validation/<provider>`. Each provider folder should contain a README with URL, access date, license, coverage, checksum, transformation, and limitations.

Corrections to observations must be reproducible. A patch log should identify the row, field, original value, replacement, transformation rule, and reason.

The publication archive should contain:

- configuration files;
- input manifest;
- environment lock file;
- scenario manifest;
- tidy metric tables;
- figure-source data;
- synthetic tests;
- public derived inputs;
- instructions for reproducing at least one complete matched scenario.

Restricted operational data should not be redistributed. The archive should state how qualified researchers may request access and should provide public or synthetic substitutes for testing the workflow.

---

## 12. Publication figure and table plan

### 12.1 Main figures

1. BC network and hydropower cascade map.
2. PyPSA cascade component and constraint schematic.
3. Validation dashboard for load, inflow, generation, and storage.
4. Paired A/B/C effect plot for primary outcomes.
5. Demand–VRE feasibility heat map by hydrology and configuration.
6. Congestion maps for A, C, and C–A.
7. Event timeline tracing hydraulic constraint to electrical outcome.
8. Cascade storage and spill heat map.
9. Accuracy-versus-complexity comparison.
10. Sensitivity interval or tornado plot.

### 12.2 Main tables

1. Model and data sources with evidence classes.
2. Operating-rule and validation coverage.
3. Exact A/B/C structural differences.
4. Scenario definitions.
5. Validation performance.
6. Primary paired outcomes.
7. Mechanism ablations.
8. Sensitivity and computational results.

### 12.3 Figure standards

- use BC Albers EPSG:3005 for analytical maps;
- retain full geographic extent and remove axes where appropriate;
- preserve common scales across compared panels;
- include units and sample coverage;
- use colorblind-safe palettes;
- export SVG or PDF plus at least 300-dpi PNG;
- archive figure-source tables;
- retain interactive figures for exploration but use static archival figures in the paper.

---

## 13. Manuscript completion strategy

Write the final manuscript in this order:

1. finalize Methods from the implemented workflow;
2. write Results from versioned tables;
3. write Discussion from demonstrated mechanisms and limitations;
4. shorten and refocus the Introduction;
5. write Conclusions;
6. write the Abstract last.

The Results should lead with findings rather than figure descriptions. It should state absolute Configuration A values before reporting A–C differences. The Discussion should open with the answer and then explain mechanism, prior-work context, planning implications, transferability, and limitations.

Use a neutral title until validation justifies stronger language:

> How Hydropower Cascade Aggregation Affects Operational Feasibility Estimates: A PyPSA Study of British Columbia

If independent validation confirms systematic optimism in the simplified model, the title may use “operational bias.”

The final manuscript should contain approximately 7,000 words excluding references and supplementary material. Tables and figures should carry detail that would otherwise burden the main text.

---

## 14. Recommended immediate actions

The first ten actions are:

1. approve the primary research question and outcomes;
2. create the study YAML and run-manifest schema;
3. freeze the current processed inputs;
4. complete the water-unit audit;
5. build the canonical cascade topology table;
6. request operational validation data;
7. construct the two-reservoir analytic test;
8. implement delayed water transfer;
9. implement minimum combined release;
10. benchmark one full-year baseline network.

The project should reassess journal ambition after Gate 4. If validation and effects are strong, target IJEPES, Energy Strategy Reviews, or Applied Energy. If the study remains a transparent but data-limited structural case study, Energy 360 remains a credible target.

---

## 15. Final recommendation

Proceed with the study. The scientific question, controlled comparison, modelling architecture, and publication strategy are sufficiently developed. Adopt the 20-week schedule as the working baseline and begin validation-data acquisition immediately.

The project’s success should not be measured by whether Configuration C performs worse than Configuration A. Success means producing a verified and reproducible estimate of the effect of aggregation, explaining the mechanism, defining where the simplification is adequate, and limiting claims to the strength of the validation evidence.

Under favourable data access, the project can produce a strong submission within four to five months. Restricted-data delays or an Applied Energy-level validation programme can extend the work to six to nine months.

---

## Appendix A — Deliverable checklist

- [ ] Approved study protocol
- [ ] Input and parameter evidence manifest
- [ ] Canonical cascade topology
- [ ] Unit-contract tests
- [ ] Synthetic delayed-cascade tests
- [ ] Minimum-flow tests
- [ ] A/B/C configuration builders
- [ ] Aggregation-equivalence report
- [ ] Validation-data manifest and patch log
- [ ] Configuration A validation report
- [ ] Locked scenario matrix
- [ ] Core solved networks and manifests
- [ ] Ablation and sensitivity runs
- [ ] Tidy metric tables
- [ ] Event datasets
- [ ] Publication figures and source data
- [ ] Supplementary information
- [ ] Completed manuscript without placeholders
- [ ] Reproducibility archive
- [ ] Journal-specific cover letter and checklist

## Appendix B — Publication decision matrix

| Evidence at completion | Recommended claim | Recommended target |
|---|---|---|
| Verified model, limited operational validation | Structural sensitivity to hydro representation | Energy 360; Energy Reports |
| Good facility validation and clear operational effects | Aggregation changes operational feasibility estimates | EPSR; IJEPES |
| Strong multi-year validation and planning implications | Model representation changes system-planning conclusions | Energy Strategy Reviews; Smart Energy |
| Broad validation, uncertainty, generalization, and decision consequence | Generalizable methodological contribution | Applied Energy |

## Appendix C — Reference journal pages

- Energy 360: https://www.sciencedirect.com/journal/energy-360
- Applied Energy: https://www.sciencedirect.com/journal/applied-energy
- Energy Strategy Reviews: https://www.sciencedirect.com/journal/energy-strategy-reviews
- International Journal of Electrical Power & Energy Systems: https://www.sciencedirect.com/journal/international-journal-of-electrical-power-and-energy-systems
- Electric Power Systems Research: https://www.sciencedirect.com/journal/electric-power-systems-research
- Smart Energy: https://www.sciencedirect.com/journal/smart-energy
