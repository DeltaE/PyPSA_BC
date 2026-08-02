# Chapter 2 — PyPSA Cascade Analysis and Modelling Plan

**Study:** Operational bias from aggregated hydropower reservoirs in British Columbia  
**Companion manuscript:** `E:\CoWork\PROJECTS\QE\CH_02\2026 07 01 - Chapter2_Cascade_Draft_EL_TN_EL.docx`  
**Repository:** `E:\CoWork\PROJECTS\PyPSA`  
**Status:** implementation plan based on a source and workflow audit; no manuscript result placeholder should be replaced until the corresponding acceptance checks pass.

## 1. Decision question and estimand

The study asks whether simplifying BC's linked reservoirs changes an operational conclusion, not whether a more detailed model merely produces different numbers.

The primary estimand for outcome \(Y\) in matched scenario \(s\) is:

\[
\Delta Y_{C,s}=Y_{\text{single bucket},s}-Y_{\text{full cascade},s}
\]

with an equivalent expression for the two-store configuration. Primary outcomes are:

1. expected/unserved energy (MWh/year), loss-of-load hours, and peak shortfall (MW);
2. wind and solar curtailment (GWh and percentage points of available VRE);
3. transmission congestion (corridor-hours, unique congested hours, and location).

Secondary outcomes are operating cost, imports/exports, hydro generation, spill, reservoir range, binding operating constraints, runtime, and memory.

Three claims must remain separate:

- **Structural effect:** A and C differ under matched inputs.
- **Mechanism:** a traced constraint, water transfer, or electrical bottleneck explains the difference.
- **Accuracy/bias:** observed operation supports A more strongly than C.

Without independent operational observations, report “structural sensitivity to cascade representation,” not measured real-world error.

## 2. Current implementation audit

### 2.1 Existing assets and workflows

The current repository already provides much of the required base:

| Concern | Current implementation | Main files |
|---|---|---|
| Hydro inventory | Generates asset, reservoir, and cascade tables; classifies `ror`, `ror-water`, `reservoir`, and `reservoir-impute` | `workflow/scripts/create_hydro_assets.py` |
| Basin inflow | Intersects HydroBASINS units, constructs incremental cascade catchments, routes runoff, and calibrates monthly inflow | `src/pypsa_bc/hydro.py`; `workflow/scripts/create_reservoir_inflows.py` |
| Run-of-river | Converts runoff shape to availability and calibrates annual energy | `src/pypsa_bc/hydro.py`; `workflow/scripts/create_ror_ps.py` |
| Reservoir components | Creates water buses, Stores, inflow Generators, store/release Links, multi-output discharge Links, spill Links, and terminal water Stores | `workflow/scripts/enrich_format_hydro.py` |
| Network assembly | Imports the electrical network and adds hydro dictionaries before clustering and solving | `workflow/scripts/build_model.py` |
| Sequential workflow | Runs preparation, validation-load ingestion, load disaggregation, and model build | `workflow/run_pypsa_bc.py` |
| Input/report visuals | Maps, basin plots, cascade schematic, Sankey, load validation, and network report | `workflow/scripts/build_input_visuals.py`; `src/pypsa_bc/vis/network_visual_report.py` |
| Validation fetching | Reproducible download location for public validation data | `workflow/scripts/fetch_validation_data.py`; `data/validation` |

Audited processed inputs currently contain:

- 118 hydro generation rows and approximately 17.73 GW;
- 21 reservoir rows;
- 16 `reservoir` or `reservoir-impute` generating assets;
- 102 `ror` or `ror-water` assets;
- 8,760 hourly reservoir-inflow rows for 16 reservoir IDs;
- cascade labels for Peace, Mica/Columbia, Bridge, Stave, Campbell, Seven Mile, and other/default assets.

These are input-table counts, not guaranteed solved-network counts. The analysis must report both and reconcile them.

### 2.2 Gaps that must be resolved before experiments

1. **No operational station-to-station delay is currently enforced.** The routing delay in `calculate_plant_inflows()` applies to catchment runoff before optimization. It is not a constraint on controlled turbine release or spill between reservoir Stores.
2. **`min_water_discharge` is carried in the hydro table but is not applied in `get_reservoir_dict()`.** A dedicated minimum-flow constraint layer is needed.
3. **Storage and flow units need a formal contract.** Current comments mix m³/s, m³/h, MW, and MWh-equivalent terms. The inflow calibration formula appears to return an hourly volume even though one docstring calls it cms. Unit tests must establish the actual convention before refactoring.
4. **The manuscript proposal overstates implementation status.** Delay and minimum-flow features are planned additions, not current baseline features.
5. **Configuration B and C builders do not exist.**
6. **There is no immutable study run manifest or matched-scenario result table.**
7. **Water-balance validation is not independent of PyPSA's own equations.**
8. **Evidence coverage is uneven.** Some downstream reservoirs receive zero imputed local inflow; facility flow, travel-time, rule-curve, and operational data may be unavailable.
9. **A two-mechanism confound is possible in Configuration C.** If the common bucket is connected to one representative electrical bus, the experiment changes hydraulic and network geography together.
10. **Current skipping behavior can hide missing hydro assets.** `logs/hydro_skipped_assets.csv` is useful, but a study run should fail when a required cascade asset is skipped unless the exclusion is declared.

## 3. Reproducible study architecture

Do not place study logic directly into general plotting scripts. Keep the reusable physical formulation in `src`, thin entry points in `workflow/scripts`, configuration in YAML, and all derived outputs in run-specific folders.

Proposed structure:

```text
config/
  studies/
    chapter2_cascade.yaml
src/pypsa_bc/studies/cascade/
  __init__.py
  schema.py
  units.py
  topology.py
  configurations.py
  constraints.py
  scenarios.py
  validation.py
  metrics.py
  events.py
  plotting.py
workflow/scripts/studies/chapter2/
  audit_inputs.py
  build_configurations.py
  run_scenarios.py
  validate_runs.py
  compute_metrics.py
  make_figures.py
  build_supplement.py
tests/studies/chapter2/
  test_units.py
  test_topology.py
  test_two_reservoir.py
  test_three_reservoir_delay.py
  test_aggregation_equivalence.py
  test_metrics.py
studies/chapter2/
  manifests/
  networks/
  metrics/
  validation/
  figures/
  tables/
  logs/
```

The current project uses a Python workflow runner rather than Snakemake. Add a study orchestrator such as:

```powershell
python -m workflow.scripts.studies.chapter2.run_scenarios `
  --config config/studies/chapter2_cascade.yaml `
  --stage all
```

The `--stage` options should include `audit`, `build`, `solve`, `validate`, `metrics`, `figures`, and `all`. Every stage must be restartable and must read upstream artifacts rather than recomputing them implicitly.

## 4. Configuration contract

`config/studies/chapter2_cascade.yaml` should specify:

```yaml
study_id: chapter2_cascade
reference_year: 2021
timezone: America/Vancouver
snapshot:
  start: 2021-01-01 00:00
  end: 2021-12-31 23:00
hydro:
  unit_system: water_volume_hourly
  configurations: [A_full, B_two_store, C_bucket]
  delays:
    source: data/inventory/custom/cascade_travel_time.csv
    boundary: warmup
  minimum_flow:
    source: data/inventory/custom/hydro_minimum_flow.csv
    missing: sensitivity_only
  local_inflow_imputation:
    core: zero
    sensitivities: [zero, upstream_ratio, catchment_area]
  terminal_condition: cyclic
scenarios:
  demand: [baseline, electrification_mid, electrification_high]
  vre: [low, medium, high]
  hydrology: [dry, median, wet]
solver:
  name: highs
  options: {}
metrics:
  congestion_threshold_pu: 1.0
  lost_load_carrier: unserved_energy
```

Every input field needs a unit, provenance source, missing-data rule, and evidence class:

- `observed`: direct public or licensed record;
- `derived`: deterministic transformation of observed data;
- `calibrated`: adjusted to a known statistic;
- `imputed`: filled from another asset or generic assumption;
- `scenario`: deliberately varied;
- `unknown`: prohibited in core runs.

## 5. Implementation phases

### Phase 0 — Freeze and inventory the baseline

Script: `audit_inputs.py`

Tasks:

1. Hash or record the modification time and size of every input.
2. Validate required columns in hydro generation, reservoirs, cascade topology, inflows, load, VRE, and network files.
3. Report duplicates, missing IDs, nonnumeric capacities, negative values, out-of-BC coordinates, orphan reservoirs, self-loops, disconnected water components, and invalid downstream references.
4. Reconcile generation rows with reservoir IDs and inflow columns.
5. Export `input_audit.csv`, `topology_edges.csv`, `parameter_coverage.csv`, and `excluded_assets.csv`.
6. Fail if a core cascade asset is skipped without a declared exclusion.

BC coordinate validation should use the BC polygon in EPSG:3005 rather than a longitude/latitude bounding box. Coordinates outside BC can remain only when explicitly representing a cross-border terminal; they must not be plotted as BC assets.

Acceptance:

- all core topology edges resolve;
- no duplicate physical Store is instantiated under different names;
- totals by cascade reconcile with the source tables;
- all exclusions have reasons and manuscript implications.

### Phase 1 — Establish and test the unit system

Module: `units.py`; tests: `test_units.py`

Recommended internal contract for one-hour snapshots:

- natural inflow, turbine discharge, and spill: m³/h;
- Store energy: m³;
- water-link nominal power: m³/h;
- water-to-power efficiency: MWh/m³;
- electrical output: MW over a one-hour snapshot;
- travel time: integer hours.

Refactor only after tests document current behavior. Add explicit conversion helpers:

```python
cms_to_m3_per_hour(q_cms)
m3_per_hour_to_cms(q)
rated_water_to_power_efficiency(capacity_mw, rated_discharge_cms)
```

Check `mean_inflow_calibration()` empirically: each modeled month's mean should equal the target monthly flow after converting to the declared unit. Remove misleading docstrings and avoid compensating scale factors whose purpose is not explicit.

Acceptance:

- a 1 m³/s constant series becomes exactly 3,600 m³/h;
- a one-hour 3,600 m³ inflow increases Store state by 3,600 m³ absent outflow;
- rated discharge times conversion efficiency equals rated MW;
- unit checks pass for every cascade asset.

### Phase 2 — Build a canonical directed cascade graph

Module: `topology.py`

Create one canonical table:

| field | meaning |
|---|---|
| `cascade_id` | stable cascade name |
| `station_id` | physical generating station |
| `upstream_reservoir_id` | source Store |
| `downstream_reservoir_id` | destination Store or terminal |
| `electric_bus` | BC electrical output location |
| `capacity_mw` | installed power |
| `q_max_m3h` | maximum turbine flow |
| `spill_max_m3h` | maximum spill |
| `q_min_m3h` | minimum combined release, nullable |
| `delay_h` | station-to-station lag |
| `efficiency_mwh_per_m3` | fixed conversion |
| `evidence_*` | evidence class by parameter |

Use a directed-acyclic-graph check within each cascade. Parallel plants that share reservoirs should remain parallel edges unless the aggregation is documented. Terminal nodes must be explicit.

Acceptance:

- topological ordering is deterministic;
- no cycle or self-loop exists;
- every path ends at a terminal;
- upstream and downstream maps reproduce the cascade Sankey;
- graph totals match the processed inventory.

### Phase 3 — Implement Configuration A

Modules: `configurations.py`, `constraints.py`

Retain the existing multi-output Link concept but construct the configuration through one tested builder. The builder must return both the PyPSA components and a mapping table from physical IDs to component names.

Delay implementation:

1. Create separate sending and receiving water-flow variables if necessary.
2. Add a Linopy constraint `receive[t] == send[t-delay]`.
3. Define boundary behavior. Preferred: run a warm-up horizon long enough to cover maximum delay, then retain the analysis year. Alternative: fix pre-horizon arrivals from observed/base dispatch. Cyclic wrap is allowed only as a documented sensitivity because it moves end-of-year water into the first hours.
4. Apply the same delay to turbine discharge and spill unless evidence supports different routing.

Minimum flow:

```text
turbine_discharge[g,t] + spill[g,t] >= minimum_release[g,t]
```

When obligations apply at a reach or gauge, build the appropriate sum rather than copying the requirement to every station.

Storage/rule curves:

- physical storage limits for all reservoirs with defensible volume;
- monthly operating bands only for covered reservoirs;
- cyclic terminal state in the core;
- observed initial state or a bounded initialization in validation runs.

Acceptance:

- two- and three-reservoir analytic tests pass;
- receiving flow shifts by exactly the declared delay;
- total water residual stays below the predeclared tolerance;
- minimum-flow constraints bind correctly in a low-demand synthetic case;
- disabling the constraint reproduces the unconstrained reference.

### Phase 4 — Implement Configurations B and C

Module: `configurations.py`

Configuration B:

- map physical reservoirs to `upper` and `lower` basin stores;
- sum storage volume without double counting;
- aggregate incremental natural inflow to the correct store;
- retain station electrical output locations where possible;
- use one upper-to-lower water transfer;
- remove station-level delays and reach-level minimum flows unless a valid aggregate constraint exists.

Configuration C:

- create one basin Store;
- sum storage and natural inflow;
- preserve electrical geography by allowing station output Links at their original buses to draw from the shared Store;
- sum or retain station capacities according to this mapping;
- remove internal hydraulic ordering and delays.

Do not connect all basin generation to one bus in the primary experiment. That would combine reservoir aggregation with transmission aggregation. It may be a separate sensitivity labeled `C_bucket_one_bus`.

Before solving, write `aggregation_equivalence.csv` comparing A, B, and C for:

- total storage;
- annual and hourly aggregate natural inflow;
- installed MW;
- maximum turbine flow;
- conversion-weighted annual energy potential;
- initial/terminal water convention;
- electrical capacity by bus.

Acceptance: exact equality where intended; explained tolerance where conversion requires weighting.

### Phase 5 — Scenario generation and solves

Modules: `scenarios.py`; script: `run_scenarios.py`

Core factors:

| Factor | Suggested levels | Required reporting |
|---|---|---|
| Hydro representation | A, B, C | exact components and switches |
| Demand | baseline, medium, high electrification | annual TWh, peak GW, peak hour, spatial changes |
| VRE | low, medium, high | wind/solar MW, available TWh, share of load |
| Hydrology | dry, median, wet | source years or scaling method, annual basin inflow |

If all levels are used, the core is 81 solves. Run A/B/C for a matched scenario consecutively or under one manifest group so partial triples are obvious.

Secondary ablations:

- A without operational delay;
- A without minimum flow;
- A with instantaneous hydraulic coupling;
- A without seasonal bounds;
- A with a noncyclic or terminal-value condition.

Sensitivities:

- travel-time low/core/high;
- local downstream inflow zero/area-based/ratio-based;
- efficiency or head low/core/high;
- observed versus optimized initial storage;
- partial versus full environmental-flow coverage;
- value of lost load;
- spill limit assumptions;
- solver and tolerance.

Each run writes:

```text
manifest.json
network.nc
solver.log
quality_checks.json
component_mapping.parquet
```

Do not overwrite successful networks. Use a deterministic `run_id` built from factor levels plus a short input/config hash.

### Phase 6 — Validation

Modules: `validation.py`; script: `validate_runs.py`

Validation ladder:

1. **Schema/topology:** IDs, graph, coordinates, components, capacities.
2. **Synthetic physics:** exact water accounting, delay, spill, minimum flow, terminal behavior.
3. **Input time series:** time zone, gaps, leap day, load totals, inflow calibration, VRE bounds.
4. **Hydro operation:** station generation, storage/elevation, spill/release where available.
5. **System operation:** annual energy balance, imports/exports, technology totals, load, known interfaces.

Suggested metrics:

- load: annual bias, monthly bias, peak error, RMSE, MAE, correlation, load-duration error;
- inflow: monthly bias, seasonal timing, KGE/NSE if using independent daily data;
- hydro generation: energy bias, NRMSE, correlation, ramp distribution, peak coincidence;
- storage: normalized level RMSE, seasonal maxima/minima timing, range bias;
- water balance: maximum and p99 absolute residual plus annual closure;
- system: energy by carrier, trade bias, and qualitative/quantitative corridor evidence.

Validation data acquisition:

- add every public URL to `config/data.yaml` under a study-specific validation section;
- fetch through `workflow/scripts/fetch_validation_data.py` or a new idempotent extension;
- store raw files under `data/validation/<provider>/`;
- write `data/validation/README.md` with URL, access date, license, checksum, transformation, and known limitations;
- never silently patch observations. A patch script must write a change log containing row, field, old value, new value, rule, and reason.

Predeclare acceptance thresholds after inspecting measurement resolution but before viewing A–C comparative outcomes.

### Phase 7 — Metrics and event attribution

Module: `metrics.py`; script: `compute_metrics.py`

Write one tidy `scenario_metrics.parquet` with one row per run and fields for all outcomes and quality flags. Write:

- `corridor_metrics.parquet`;
- `hydro_station_metrics.parquet`;
- `reservoir_metrics.parquet`;
- `hourly_events.parquet`;
- `paired_effects.parquet`.

Metrics must be computed from explicit carrier/component mappings, not name substring matching. Maintain a registry for:

- unserved-energy generators;
- wind and solar generators;
- hydro turbine and spill Links;
- reservoir Stores;
- transmission corridors.

Event selection must be deterministic:

1. hour of maximum A–C unserved-energy difference;
2. multi-hour window with maximum cumulative curtailment difference;
3. corridor-hour with maximum congestion redistribution.

`events.py` should export all causal series for a window: demand, net load, generation, imports, storage, turbine flow, spill, delayed arrival, minimum-flow bound, line loading, and shadow prices where reliable.

### Phase 8 — Figures and tables

Module: `plotting.py`; script: `make_figures.py`

Every figure reads only the tidy tables or declared network files and writes both a vector/static publication file and the underlying data.

Core manuscript visuals:

1. **BC system and cascade map:** BC Albers EPSG:3005, full province extent, axes removed, hydro markers sized by MW, reservoirs distinct from generators, and Peace/Columbia insets.
2. **Cascade component schematic:** one reach showing Store, inflow, turbine, electricity, delayed water, spill, and minimum flow.
3. **Validation dashboard:** observed versus modeled load, inflow, generation, and storage with units and sample coverage.
4. **Paired effect plot:** A/B/C values connected within each scenario for EENS, curtailment, and congestion.
5. **Feasibility heat map:** demand × VRE, faceted by hydrology and configuration.
6. **Congestion maps:** A, C, and C–A difference with identical extent and scale.
7. **Event timeline:** common x-axis for net load, hydro output, flow, storage, spill, and line loading.
8. **Cascade storage heat map:** reservoirs ordered upstream to downstream.
9. **Accuracy–complexity plot:** operational error versus runtime/component count.
10. **Sensitivity interval/tornado plot:** changes in the A–C effect.

Supplementary visuals:

- complete cascade Sankey plus one Sankey per cascade;
- travel-time and evidence-quality maps;
- load-duration and residual-load curves;
- inflow hydrographs and monthly calibration;
- facility generation duration curves;
- reservoir rule-curve plots;
- line loading duration curves;
- shadow-price distributions;
- solve runtime/memory diagnostics;
- missingness and validation-coverage matrix.

Publication rules:

- static maps at least 300 dpi and also SVG/PDF;
- no cropped province or legends;
- common scales for compared panels;
- colorblind-safe palette;
- units in every axis/legend;
- no dual axis unless unavoidable and explicitly explained;
- observed data styled consistently across figures;
- interactive Plotly versions may accompany the report but do not replace archival static figures.

Core tables:

1. data sources and evidence classes;
2. model/framework literature comparison;
3. operating-rule coverage;
4. A/B/C structural equivalence;
5. scenario definitions;
6. validation performance;
7. primary paired outcomes;
8. sensitivity and ablation summary;
9. exclusions and unresolved limitations.

## 6. Manuscript placeholder crosswalk

| Manuscript placeholder | Required artifact |
|---|---|
| R0, R5 | `scenario_metrics.parquet`, run manifest summary |
| R1 | input audit and component mapping |
| R2 | `quality_checks.json`, synthetic test report |
| R3 | load/inflow/RoR validation tables |
| R4 | station/reservoir validation coverage and metrics |
| R6 | paired EENS table and scarcity event export |
| R7 | wind/solar curtailment effects |
| R8 | corridor metrics and difference maps |
| R9 | deterministic event trace |
| R10 | B recovery and runtime/component metrics |
| R11 | ablation paired effects |
| R12 | sensitivity paired effects |
| C1–C2 | final effect table plus validation qualification |

Add a script that checks the manuscript source for `[RESULT PLACEHOLDER`, `[MODEL PLACEHOLDER`, and related tags. The publication build should fail if unresolved analytical placeholders remain in a submission-mode document.

## 7. Quality gates

### Gate A — Input readiness

- unit contract passes;
- topology graph passes;
- parameter coverage published;
- all exclusions declared;
- A/B/C totals reconcile.

### Gate B — Physics readiness

- analytic tests pass;
- max water residual is within tolerance;
- delay and boundary behavior verified;
- minimum-flow and terminal constraints verified.

### Gate C — Solve readiness

- matched scenario triplets solve;
- objective and balance diagnostics pass;
- run manifests are complete;
- no output overwrite.

### Gate D — Validation readiness

- calibration and validation roles separated;
- load and hydro performance reported;
- missing observed coverage visible;
- claims calibrated to evidence.

### Gate E — Publication readiness

- all primary metrics reproducible from tidy tables;
- every comparison uses matched inputs;
- mechanism claims trace to event data or ablations;
- figures have source data and consistent scales;
- manuscript placeholders resolved or explicitly retained as planned-analysis callouts;
- code release and data availability statement point to an immutable version.

## 8. Recommended implementation order

1. Freeze inputs and produce the topology/parameter audit.
2. Resolve the water unit contract and add unit tests.
3. Build two- and three-reservoir analytic networks.
4. Implement Configuration A delay and minimum flow.
5. Independently verify water balance from solved exports.
6. Build and equivalence-test B and C.
7. Run one baseline matched triplet and inspect event mechanics.
8. Obtain/fetch validation data and complete the validation matrix.
9. Lock the scenario matrix and acceptance thresholds.
10. Run core scenarios, then ablations and sensitivities.
11. Compute tidy metrics and generate publication figures.
12. Replace manuscript placeholders only from versioned artifacts.

## 9. Minimum publishable analysis

If operational data remain limited, the minimum defensible paper is:

- fully tested A/B/C physical formulations;
- transparent parameter-evidence and missing-data tables;
- matched scenario effects;
- independent load validation and water-conservation validation;
- at least one observed hydro-generation or storage comparison;
- event-based mechanism tracing;
- strong qualification that the reported difference is a model-structural sensitivity.

The stronger paper adds multi-facility hourly validation, observed reservoir trajectories, defensible travel times and minimum-flow coverage, one-at-a-time mechanism ablations, and a held-out hydrological year.

## 10. Immediate risks to address

1. Do not state that PyPSA currently applies a native delay attribute until the exact installed version and serialized component schema are verified.
2. Do not describe calibrated monthly inflow as independent validation.
3. Do not use a single-bus bucket in the primary comparison.
4. Do not infer environmental flows directly from generic `min_water_discharge` values without source and temporal semantics.
5. Do not interpret perfect-foresight dispatch as operator behavior.
6. Do not compare unmatched renewable capacity, load, transmission, or hydrology across A/B/C.
7. Do not report congestion from axes or maps alone; calculate it from rated flow with a declared threshold and corridor mapping.
8. Do not allow failed or skipped hydro assets to disappear from the study without an exclusion table.
9. Do not publish only annual totals; include scarcity-event traces that demonstrate mechanism.
10. Do not replace analytical placeholders manually from ad hoc notebook output; use the versioned metric tables.
