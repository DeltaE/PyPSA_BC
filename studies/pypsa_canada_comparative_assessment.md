# PyPSA-Canada and PyPSA-BC: Comparative Assessment and Strategic Development Plan

**Assessment date:** 2026-07-29  
**Purpose:** Identify developments from PyPSA-Canada that can strengthen PyPSA-BC without weakening the scientific design of the Chapter 2 cascade-aggregation study.

## Executive verdict

PyPSA-Canada and PyPSA-BC are complementary rather than competing models.

- **PyPSA-Canada is currently stronger as a general modelling framework.** Its main advantages are the Snakemake workflow, command-line interface, scenario configuration, separation of planning and dispatch, representative-day options, rolling-horizon dispatch, policy constraints, run logging, continuous integration, and national/interprovincial scope.
- **PyPSA-BC is currently stronger as a scientific model of British Columbia hydro operation.** It has more detailed BC electrical geography, station-resolved generation, explicit water stores and water-transfer paths, cascade connectivity, turbine and spill releases, travel-delay mechanics, terminal-water handling, evidence classes, source hashes, and formal input/physics gates.
- **The best strategy is not to replace PyPSA-BC with the national model.** PyPSA-BC should adopt the software-engineering and scenario-management strengths of the PyPSA-Canada framework, while retaining its own hydro representation. The national model should be used as an external planning-scale benchmark and as the deliberately aggregated comparison configuration in the paper.

The strongest publication contribution is therefore:

> National-scale models necessarily simplify reservoir hydro. PyPSA-BC can quantify when those simplifications materially change estimates of adequacy, curtailment, congestion, imports, spill, and operational feasibility in a hydro-dominant system.

### Clarification: staged execution versus workflow orchestration

PyPSA-BC's staged S1-S5 design is not intrinsically weaker than Snakemake. The stages provide a sensible scientific decomposition from source-data preparation through asset construction, profile preparation, network assembly, solution, and visualization. The current limitation lies in **machine-enforced orchestration**, not in the staged concept or the BC-specific model formulation.

At present, stage order and artifact freshness depend more heavily on the operator. The workflow does not yet expose every dependency through a central directed acyclic graph (DAG), automatically rebuild all affected downstream products after an upstream change, or consistently bind each result to an immutable configuration and input manifest. PyPSA-Canada's advantage is therefore execution engineering: dependency resolution, incremental rebuilding, standardized scenario invocation, logging, and failure recovery. Snakemake does not establish scientific validity; incorrect assumptions or incomplete output declarations remain incorrect inside a Snakemake workflow.

The appropriate development path is to retain the PyPSA-BC stages and wrap them with explicit workflow contracts. Each production stage should declare its inputs, outputs, configuration, validation conditions, log, and provenance record. Notebooks should remain available for exploration and review, while production scripts create the authoritative artifacts consumed by later stages.

## Material reviewed

This assessment distinguishes two NRCan repositories:

1. The [PyPSA-Canada framework](https://github.com/NRCan/pypsa-canada), a reusable scenario, planning, dispatch, and workflow package.
2. The [PyPSA-Canada National model](https://github.com/NRCan/pypsa-canada-national), an aggregated model of the bulk electricity systems of Canada's ten provinces.

The local repository snapshots inspected were:

| Repository | Inspected commit | Commit date | Observable maturity |
|---|---:|---:|---|
| `NRCan/pypsa-canada` | `c8de2986897f08fe43b4a190a08a6a0ce372c361` | 2026-07-29 | 148 commits shown by GitHub; Python package, documentation, tests, CI, and workflow |
| `NRCan/pypsa-canada-national` | `df7156ad5a600ebdccf5840688be6332b2ca2484` | 2026-07-23 | 29 commits shown by GitHub; seven sequential model-building notebooks plus visualization scripts |

The PyPSA-Canada framework describes scenario-based modelling, planning and dispatch, representative days, configurable spatial scales, cost optimization, custom constraints, and a Snakemake workflow. The national repository describes an API-assisted, notebook-based build for a spatially aggregated capacity-expansion and planning model. These stated purposes are important: they are not claiming station-level reservoir operations.

For PyPSA-BC, this review inspected the package and workflow code, the Chapter 2 protocol, its study tests, and the current Gate 1 and Gate 2 outputs. The current evidence, updated 2026-08-01, is:

- 45 focused Chapter 2 tests pass;
- Gate 2 physical verification passes seven of seven analytic cases;
- Gate 1 input/topology/evidence audit passes with two explicit claim-boundary warnings;
- the matched representation audit passes 12/12 rows, release uncertainty passes 42/42 checks, and storage uncertainty passes 6/6 checks;
- 20 cascade stations in six cascade groups are covered by the audited representation;
- 30 cases are declared across representation, release-policy, and Waneta-storage axes, of which 24 are runnable and six evidence-only reference cases remain deliberately blocked;
- travel-time values remain scenarios rather than measured station-specific observations.

BC Hydro is also a major public operational-data provider for this comparison.
Its website publishes historical hourly balancing-authority load, directional
hourly transfer capability with Alberta and the United States, net actual flow,
ACE reports, a hydrometeorologic station network, reservoir-level and selected
discharge pages, IPP supply lists, transmission diagrams, and engineering
studies. These sources strengthen PyPSA-BC's provincial calibration and
external-boundary validation. They do not supply a complete station-level
dispatch record or a fully machine-readable current network. The resource-level
claim boundaries and retrieval rules are catalogued in
`data/validation/bc_hydro/public_data_catalog.csv`.

## Comparative picture

Scores are relative to the current intended use of each model, not universal quality grades. A low national-model score for cascade physics is expected for a planning model; a low PyPSA-BC score for national policy analysis reflects its narrower scope.

| Scientific or technical dimension | PyPSA-Canada framework / national model | PyPSA-BC | Comparative finding |
|---|---:|---:|---|
| Reproducible workflow orchestration | 5/5 | 3/5 | PyPSA-Canada has a modular Snakemake DAG, CLI, logs, copied configuration, timestamps, benchmarks, and crash artifacts. PyPSA-BC has a sound staged design and many reproducible scripts, but stage dependencies, artifact freshness, and selective rebuilding are not yet enforced by a central workflow engine. |
| Scenario and configuration management | 5/5 | 4/5 | Canada cleanly separates defaults, scenarios, planning, dispatch, solvers, policies, and post-processing. PyPSA-BC now compiles 30 immutable case definitions across three scientific axes; 24 are runnable and six evidence-only cases are blocked by design. |
| Long-term capacity expansion | 5/5 | 2/5 | Canada supports multiple investment periods, lifetimes, discounting, and extendable assets. PyPSA-BC is currently oriented toward operational comparison. |
| Operational dispatch workflow | 4/5 | 3/5 | Canada has planning-to-dispatch conversion and rolling-horizon controls. PyPSA-BC has detailed operational components but not yet an equally standardized dispatch runner. |
| National and interprovincial context | 5/5 | 2/5 | Canada represents all ten provinces and trade interfaces. PyPSA-BC resolves BC internally but needs carefully defined external boundary conditions. |
| BC spatial and asset detail | 2/5 | 5/5 | The national model clusters facilities and corridors. PyPSA-BC retains substantially more detailed substations, lines, generators, and BC geography. |
| Reservoir storage realism | 1/5 | 4/5 | The national model aggregates storage by cluster and assumes duration. PyPSA-BC uses explicit reservoir stores and volumes, although some parameters still lack sufficient evidence. |
| Cascade connectivity and water conservation | 1/5 | 5/5 mechanics; 3/5 evidence | The national model does not route water station-to-station. PyPSA-BC explicitly represents turbine, spill, downstream routing, and terminal water. |
| Travel-time and routing mechanics | 1/5 | 4/5 mechanics; 2/5 calibration | PyPSA-BC has verified integer-delay mechanics, but station-specific delays remain scenario assumptions. |
| Environmental/operational release rules | 1/5 | 3/5 | PyPSA-BC supports combined turbine-plus-spill minimum release, but schedule-, route-, and state-dependent rules remain incomplete. |
| Hydrological data breadth | 4/5 | 4/5 | Canada uses MERIT/GRFR discharge and national sources. PyPSA-BC uses HydroBASINS, BC-specific inventories, inflows, WUP evidence, and facility data. |
| Independent hydro validation | 2/5 | 2/5 currently | Canada's visible validation is primarily aggregate energy matching. PyPSA-BC has stronger validation architecture, but independent operational validation is not yet complete. |
| Policy constraints | 5/5 | 2/5 | Canada includes reserve margin, emissions and policy hooks. These are useful additions only where they support the BC research question. |
| Software tests and integration controls | 5/5 | 3/5 | Canada has package tests, CI coverage, and pre-commit checks. PyPSA-BC has valuable physics/audit tests but no equivalent visible repository-level CI workflow. |
| Scientific evidence traceability | 3/5 | 5/5 | PyPSA-BC's evidence classes, decisions, source registry, manifests, hashes, and publication gates are unusually strong and should be preserved. |
| Publication-specific falsifiability | 3/5 | 5/5 design; production pending | PyPSA-BC predeclares configurations, matched comparisons, outcomes, and failure gates. The 24 runnable production solves and independent validation are still outstanding. |

## What PyPSA-Canada does especially well

### 1. Workflow architecture

The framework divides the model into explicit rules for loading the base network, adding components and loads, selecting snapshots, adding representative days, modifying components, solving planning, creating a dispatch network, solving dispatch, post-processing, mapping corridors, summarizing results, and finalizing the run.

This is the largest immediately transferable improvement. It would strengthen the execution of the existing PyPSA-BC stages by making the workflow:

- restartable from the last valid step;
- easier to audit and reproduce;
- safer for the declared multi-axis experiment;
- explicit about the dependency between data, model, solve, validation, figures, and tables;
- able to preserve the exact configuration and logs of every publication result.

### 2. Planning and dispatch are separate products

The framework treats capacity planning and operational dispatch as different workflows. It can build a dispatch network from planning results and offers dispatch horizons and overlaps. This separation is scientifically useful because investment decisions and operational feasibility should not be conflated.

For the Chapter 2 paper, the primary configurations should use **fixed, matched electrical capacity** and full hourly operational chronology. Capacity expansion can be a separate extension or sensitivity, not the central comparison.

### 3. Mature scenario controls

The default configuration contains:

- multiple investment periods and a discount rate;
- alternative solver profiles and random seeds;
- configurable load shedding and unit commitment;
- planning reserve-margin constraints;
- emissions and policy constraints;
- component-modification and custom-constraint hooks;
- planning and dispatch post-processing.

PyPSA-BC's Chapter 2 protocol already has a scientifically stronger matched scenario matrix, but it needs an execution layer that turns every declared case into an immutable run specification.

### 4. Chronology reduction methods

The framework implements multiple representative-day approaches, including methods that consider variable renewables and hydro, k-medoids, CARPE-DIEM, peak/average selection, and all days.

This is valuable for future planning studies, but it must be used cautiously in a cascade paper. Representative days can destroy inter-day storage continuity, seasonal water budgets, drought persistence, and travel-delay chronology. The publication reference should remain 8,760 sequential hourly snapshots. Representative days can be tested only as a **temporal-aggregation sensitivity**, with its induced error reported.

### 5. National context and policy analysis

The national model provides a route to compare BC against a consistent Canada-wide representation:

- interprovincial corridors and U.S. trade;
- clustered regional load and generation;
- national capacity-expansion assumptions;
- federal emissions and clean-electricity constraints;
- projected end-use demand.

This is particularly useful for boundary-condition sensitivity and for explaining why a BC-specific operational model is needed.

### 6. Engineering quality controls

The framework includes tests, a coverage pipeline, pre-commit checks, documentation, a changelog, contribution guidance, and run finalization. Successful runs retain configuration and logs; failed runs collect crash artifacts. These are strong practices for a transparent scientific codebase.

## Where PyPSA-BC is scientifically stronger

### 1. It represents water, not only hydro energy

The national model represents aggregated hydro as `StorageUnit` objects. Its hydro notebook states that:

- reservoir hydro is modelled on a monthly timescale;
- provincial inflows are scaled to historical aggregate generation;
- reservoirs are aggregated by cluster;
- daily storage is assumed as 24 times rated output;
- monthly storage is assumed as 31 days times rated output;
- per-unit storage and inflows are not represented;
- Québec reservoir duration is multiplied by nine to compensate for under-sizing.

Those assumptions can be reasonable screening approximations at national planning scale, but they cannot enforce physical station-to-station water balance.

PyPSA-BC instead has explicit water buses/stores and distinct turbine, spill, routing, run-of-river water, and terminal paths. Its analytic tests verify immediate electricity production, delayed downstream water, minimum combined release, turbine-plus-spill mass balance, delay boundaries, and terminal storage behavior.

### 2. It can isolate cascade-aggregation bias

The Chapter 2 design has:

- Configuration A: station-resolved full cascade;
- Configuration B: upper/lower stores;
- Configuration C: one basin bucket;
- a one-bus sensitivity that also changes electrical geography;
- matched demand, VRE, and hydrology cases;
- predeclared adequacy, curtailment, congestion, storage, spill, import, cost, and runtime outcomes.

That controlled hierarchy is precisely what is needed to determine whether an aggregated national-style hydro representation gives an optimistic or pessimistic estimate of operational feasibility.

### 3. Evidence is treated as a model variable

PyPSA-BC distinguishes observed, derived, calibrated, imputed, scenario, not-applicable, and unknown parameters. It also records source registries, manifests, hashes, and parameter decisions. This is stronger scientific practice than silently filling missing reservoir or release data.

### 4. BC operational geography is preserved

PyPSA-BC's detailed lines, substations, generators, hydro stations, and maps can expose congestion and locational adequacy effects hidden by provincial or coarse regional aggregation. This is central to the paper: aggregation can modify both the hydraulic and electrical feasible regions.

## Weaknesses and cautions in PyPSA-Canada

These limitations should not be read as defects outside the model's intended planning scope.

1. **Hydro is an energy-budget approximation, not a cascade model.** It has no visible explicit downstream water routing, station travel time, route capacity, environmental release schedule, or reservoir elevation/head relationship.
2. **Storage duration includes broad assumptions.** Fixed 24-hour and 744-hour durations and the Québec factor-of-nine adjustment require sensitivity analysis if used for hydro conclusions.
3. **Calibration and validation are not independent at plant level.** Run-of-river conversion is calibrated to annual station production, while reservoir inflows are scaled to provincial historical energy. Matching the data used for scaling is not independent validation.
4. **The national build is notebook-driven.** Notebooks are transparent for exploration, but sequential manual execution is less reproducible than a tested data-build DAG.
5. **Some required data are access-controlled.** CODERS requires an account and API key, which complicates full third-party reproduction unless raw-data redistribution and cached derived products are addressed.
6. **Chronology reduction is risky for long-duration hydro.** A representative-day planning solution should not be treated as evidence of hourly cascade feasibility.
7. **Loose/fast solver profiles need care.** Any load shedding, relaxed feasibility tolerances, or approximate solver termination must be reported and subjected to a strict rerun before publication.
8. **No visible station-level operational validation suite exists in the national repository.** Aggregate generation agreement alone cannot validate cascade constraints.

## Weaknesses and cautions in PyPSA-BC

1. **Passing gates do not complete validation.** Gate 1 now passes with two claim-boundary warnings and the Gate 2 verification suites pass, but independent observed-operation validation and production results remain outstanding.
2. **Key operating evidence remains bounded rather than fully observed.** Release and Waneta-storage uncertainty are now explicit scenario axes; travel times remain structural scenarios, not calibrated observations.
3. **Important hydro rules are not yet encoded.** Route-specific, calendar/schedule-dependent, and state-dependent release rules are open work items.
4. **The production experiment is orchestrated but not yet executed.** The workflow compiles 30 cases, blocks six evidence-only reference cases, and exposes 24 runnable configurations. The remaining publication task is to solve, validate, summarize, and archive those cases with complete run manifests.
5. **Planning and dispatch are not cleanly separated.** The current model is strongest for operational analysis, but future expansion studies require a dedicated planning formulation.
6. **External boundary conditions need a benchmark.** BC imports, exports, and neighbouring-system flexibility could dominate adequacy results if treated as unconstrained backstops.
7. **Independent validation is still outstanding.** The validation architecture is strong, but observed station/reservoir/flow series must be separated into calibration and hold-out data.
8. **Repository-level engineering is behind the NRCan framework.** CI, coverage reporting, pre-commit enforcement, packaged documentation, and run/crash artifact handling should be added.

## Strategic developments to adopt

### Priority 1 — adopt before the production experiment

| Development | How to implement in PyPSA-BC | Scientific benefit | Indicative effort |
|---|---|---|---:|
| Snakemake experiment DAG | Add rules for fetch, prepare, audit, build A/B/C, solve, validate, summarize, plot, and manuscript tables | Exact dependency graph, restartability, and auditable paper build | 4–7 days |
| Immutable run manifest | Save merged YAML, Git/source state where available, input hashes, environment package versions, solver/options, timestamps, status, and runtime with every case | Makes every figure and table traceable | 2–3 days |
| Strict solver profiles | Define `screening`, `production`, and `verification` profiles; production must record feasibility/optimality gaps and rerun flagged cases | Prevents false feasibility or cost differences from solver settings | 1–2 days |
| Standard scenario expansion | Retain the compiled 3 representations × 5 release policies × 2 Waneta-storage policies; solve the 24 runnable cases and preserve the six blocked evidence cases | Eliminates manual case drift while keeping evidence gaps visible | Implemented; execution pending |
| Planning/dispatch separation | Keep Chapter 2 capacities fixed; place any capacity-expansion solve upstream and feed its assets to the operational model | Preserves causal interpretation of aggregation effects | 2–4 days |
| CI and pre-commit | Run the 45 focused tests, Gate 1, all Gate 2 verification suites, formatting, and static checks on each change | Protects cascade physics during development | 1–2 days |
| Failure artifact collection | Preserve failed configuration, logs, solver file/status, and input hashes | Makes infeasibility diagnosable and reportable | 1–2 days |

### Priority 2 — adopt after the core hydro evidence gaps

| Development | Recommended adaptation | Scientific use |
|---|---|---|
| Rolling-horizon dispatch | Add state handoff for water stores and delay queues; test several horizon/overlap lengths against the full-year perfect-foresight reference | Quantify perfect-foresight bias and enable operational studies |
| National boundary benchmark | Extract comparable BC load, generation, storage, and transfer aggregates from PyPSA-Canada National | Bound import/export and neighbouring-system assumptions |
| Spatial clustering sensitivity | Use the national model's Leiden/census-region concept to construct a deliberately coarser electrical network | Separate hydraulic aggregation error from electrical aggregation error |
| Standard policy hooks | Add emissions cap, reserve margin, and expansion constraints as modular optional constraints | Supports future decarbonization studies without contaminating the core cascade experiment |
| Structured summaries | Emit case-level Parquet/CSV tables for energy balance, adequacy, hydro balance, storage, congestion, cost, and runtime | Enables automatic manuscript tables and uncertainty plots |

### Priority 3 — useful future extensions

- multi-investment-period capacity expansion;
- federal policy scenarios;
- Canada-wide linked studies;
- API adapters with cache manifests and redistribution-safe derived data;
- IDEA or other model-exchange export;
- representative-day planning studies after hydro chronology error has been quantified.

## What should not be copied

1. Do not replace explicit cascade water stores and routes with a `StorageUnit` in Configuration A.
2. Do not set reservoir energy capacity from a uniform `max_hours` assumption without evidence and sensitivity bounds.
3. Do not scale storage by an unexplained correction factor in the publication reference case.
4. Do not use representative days as the primary chronology for cascade-feasibility claims.
5. Do not validate a model only against the same aggregate energy series used for calibration.
6. Do not allow unreported load shedding to turn physical infeasibility into apparently feasible results.
7. Do not treat a notebook run as a complete reproducibility package unless its inputs, environment, cell order, outputs, and hashes are controlled.
8. Do not import every national policy feature before the cascade paper; added model scope can obscure the central causal comparison.

## Recommended cross-model benchmark

The models should not be expected to produce identical networks. Instead, construct an aligned BC benchmark for a common historical year using quantities that both models can represent:

| Benchmark family | Common comparisons |
|---|---|
| Demand | annual TWh, hourly peak GW, seasonal profile, regional/provincial allocation |
| Installed capacity | MW by hydro, wind, solar, biomass, gas, and other thermal |
| Generation | annual and monthly TWh by carrier |
| Hydro | run-of-river vs reservoir capacity and generation; aggregate storage-energy assumption; spill where available |
| Trade | annual imports/exports, peak interface flow, interface capacity |
| Transmission | total corridor transfer capacity; major BC internal bottlenecks where mappings are possible |
| Adequacy | unserved energy, loss-of-load hours, peak shortfall, and assumed value of lost load/load-shedding cost |
| Emissions and costs | consistent units, currency year, fuel costs, carbon price, and scope |

Every difference must be classified as:

- data-source difference;
- spatial aggregation difference;
- temporal aggregation difference;
- component-formulation difference;
- boundary-condition difference;
- solver/configuration difference.

This classification prevents an apparent “model disagreement” from being incorrectly attributed to cascade aggregation.

## Publication strategy

### Use PyPSA-Canada as the external planning benchmark

The manuscript can state that an official Canadian planning framework provides national context but uses an intentionally aggregated hydro representation. PyPSA-BC then tests the operational consequences of progressively aggregating a physically explicit BC cascade.

The comparison should be methodological, not adversarial. The paper should not claim that PyPSA-Canada is “wrong”; it should demonstrate the conditions under which a planning approximation is or is not adequate for particular operational metrics.

### Add a national-style configuration as a transparent sensitivity

Configuration C should be parameterized to resemble the national model as closely as scientifically defensible:

- one aggregate storage bucket per basin or BC region;
- a documented energy-capacity duration;
- monthly inflow energy;
- no station travel time;
- aggregated turbine capacity;
- the same electrical buses and boundary conditions as Configuration A wherever possible.

This creates a direct bridge from the paper's experimental design to Canadian planning practice. A second sensitivity can apply coarse electrical clustering. Comparing these separately avoids confounding hydraulic and electrical aggregation.

### Frame the novelty precisely

The likely novelty is not “another PyPSA model of Canada.” It is:

1. a station- and route-aware operational representation of major BC hydro cascades;
2. a controlled hierarchy of aggregation configurations;
3. matched stress scenarios;
4. explicit input and physics gates;
5. quantified error in adequacy, curtailment, congestion, imports, spill, and storage behavior caused by aggregation;
6. practical guidance on the resolution required for different planning questions.

### Claims that are supportable only after completion

Do not claim that cascade aggregation materially biases a metric until:

- Gate 1 passes;
- scheduled/route/state-dependent rules are implemented or explicitly bounded;
- configurations A, B, and C conserve matched water and energy budgets;
- production solver tolerances are demonstrated not to drive the differences;
- independent validation is completed;
- all 24 runnable matched cases are solved and the six evidence-only exclusions remain transparently documented;
- sensitivity to perfect foresight and external trade assumptions is reported.

## Recommended decision

Proceed with PyPSA-BC as the publication model. Do not migrate the study into PyPSA-Canada National, and do not discard the existing S1-S5 staged architecture.

Adopt the PyPSA-Canada **framework practices** in this order:

1. run manifests and strict solver profiles;
2. Snakemake experiment orchestration;
3. CI, coverage, and crash artifacts;
4. planning/dispatch separation;
5. structured result summaries;
6. national boundary and spatial-aggregation benchmarks;
7. rolling-horizon and representative-period sensitivities.

Retain PyPSA-BC's **hydro mechanics and evidence system** as the scientific core. The immediate modelling priority remains resolving Gate 1, encoding the outstanding operating rules, and obtaining independent validation—not adding broad national planning scope.

## Source links

- [NRCan PyPSA-Canada framework](https://github.com/NRCan/pypsa-canada)
- [NRCan PyPSA-Canada National model](https://github.com/NRCan/pypsa-canada-national)
- [PyPSA-Canada default scenario configuration](https://github.com/NRCan/pypsa-canada/blob/main/example/config/default_config.yaml)
- [PyPSA-Canada workflow](https://github.com/NRCan/pypsa-canada/tree/main/pypsa_canada/workflow)
- [PyPSA-Canada tests](https://github.com/NRCan/pypsa-canada/tree/main/tests)
- [PyPSA-Canada National hydro notebook](https://github.com/NRCan/pypsa-canada-national/blob/main/3a%20-%20Hydro%20Modeling.ipynb)
- [PyPSA-Canada National network-clustering notebook](https://github.com/NRCan/pypsa-canada-national/blob/main/2%20-%20Network%20Clustering.ipynb)
- [PyPSA project](https://github.com/PyPSA/PyPSA)
