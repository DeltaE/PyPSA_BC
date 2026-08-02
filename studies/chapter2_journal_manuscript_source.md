# Quantifying Operational Bias from Aggregated Hydropower Reservoirs: A PyPSA Framework for British Columbia

**Running title:** Hydropower cascade aggregation and operational feasibility

**Authors:** [AUTHOR PLACEHOLDER: insert author names, affiliations, ORCID identifiers, and corresponding-author email]

## Highlights

- A controlled three-configuration experiment isolates the effect of hydropower spatial aggregation.
- A station-resolved PyPSA formulation conserves water through turbine, spill, and delayed downstream flows.
- Operational bias is evaluated using unserved energy, renewable curtailment, and transmission congestion.
- British Columbia provides a hydro-dominant test system with large Peace and Columbia River cascades.
- The workflow separates verified model outputs from assumptions necessitated by unavailable operating data.

## Abstract

Hydropower reservoirs are a principal source of flexibility in low-carbon electricity systems, but large-scale planning models commonly represent a river basin as one equivalent energy store. This abstraction removes the ordering of generating stations, the travel time of water, and operating restrictions that determine whether stored water can be converted into electricity at a particular location and hour. We develop a controlled experiment to quantify the operational bias introduced by reservoir aggregation in British Columbia, Canada, a hydro-dominant system facing simultaneous electrification and variable renewable energy growth. The experiment is implemented in BC-PyPSA, an hourly nodal dispatch model, and compares three hydropower representations while holding demand, renewable portfolios, transmission topology, costs, and hydrological input constant: a station-resolved cascade with delayed water transfer and minimum-flow constraints; a two-store basin representation; and a single energy-bucket representation. Each representation is evaluated across matched demand, renewable, and hydrological stress cases. Primary outcomes are annual unserved energy, wind and solar curtailment, and corridor congestion hours; secondary outcomes include peak shortfall, spill, reservoir drawdown, dispatch cost, and computational burden. Validation follows a hierarchy from component topology and hourly water balance to facility generation, seasonal storage, provincial energy, and transmission-flow patterns.

**[RESULT PLACEHOLDER R0—ABSTRACT: Report the number of modeled reservoirs, generating stations, cascade groups, scenarios, and successful runs. Insert the principal effect size for Configuration C relative to Configuration A for unserved energy, curtailment, and congestion, each with units and uncertainty or scenario range. State whether Configuration B recovered most of the full-cascade result. Do not use “significant” unless a defined statistical test is reported.]**

The study provides a reproducible test of when an equivalent-reservoir abstraction is adequate and when it creates planning-relevant optimism about hydro flexibility. The framework is intended for transfer to other cascade-dominated systems where physical topology is known but operating records are incomplete.

**Keywords:** cascaded hydropower; reservoir aggregation; power-system dispatch; PyPSA; operational feasibility; variable renewable energy; British Columbia; model validation

## 1. Introduction

Electricity systems with large regulated reservoirs are frequently described as flexible. Reservoirs can shift generation from low-value to high-value hours, provide rapid balancing, and reduce the curtailment associated with wind and solar generation. That flexibility is real, but it is not unlimited and it is not located in a single abstract store. Water is held in specific reservoirs, released through specific turbines, and passed through a sequence of downstream stations. Its availability is conditioned by river topology, storage limits, travel time, turbine capacity, environmental releases, contractual obligations, and the electricity network that connects generation to demand. A model can therefore conserve annual hydro energy while overstating the amount of electricity that the system can deliver at a critical hour.

This distinction becomes consequential as electrification changes both the magnitude and timing of electricity demand. Space heating, industrial conversion, transport electrification, and new large loads can raise winter peaks more quickly than annual energy requirements. At the same time, wind and solar generation introduce periods of surplus and scarcity that increase the value of temporal shifting [1–6]. Hydro-dominant systems such as British Columbia, Québec, and Norway are often expected to absorb these changes with existing reservoirs. Planning conclusions about transmission reinforcement, renewable siting, firm capacity, trade, and resource adequacy then depend on how accurately the dispatch model represents the physical conversion of stored water into electricity.

Most large-scale open electricity models simplify hydropower because detailed operational data are difficult to obtain. A common approach aggregates all reservoirs and generating stations in a basin into one equivalent energy store. Natural inflow is converted to energy and added to this store; generation removes energy subject to an aggregate power limit. This representation is computationally efficient and may be adequate for estimating annual energy. It does not, however, preserve the hydraulic sequence of a cascade. Water released at an upstream station may take several hours to arrive at the next reservoir. The same cubic metre may produce electricity at multiple stations, but only in the correct order and within each station's discharge limit. A required environmental release may force water through a station when electricity has little value, while a reservoir rule curve may prevent the model from using water that appears available in an unconstrained aggregate store.

Four simplifications are especially relevant to operational feasibility. First, aggregation removes inter-reservoir coupling, allowing the optimizer to use basin energy without demonstrating that the required water reaches the relevant plant. Second, instantaneous transfer eliminates travel time and can make downstream generation available too early during scarcity events. Third, fixed or omitted operating rules remove lower bounds on release and seasonal bounds on storage. Fourth, an unconstrained terminal condition allows a finite-horizon optimization to consume water that should remain valuable after the modeled year. Head-dependent efficiency creates an additional nonlinear relationship between storage level, discharge, and power, but it can be held at a common linearized value across experimental configurations so that the effect of topology is identifiable.

The direction of the expected bias is intuitive but its magnitude is not. Removing constraints enlarges the mathematical feasible region, so a simplified model cannot have a higher minimum operating cost than an otherwise identical detailed model. Yet the resulting differences in unserved energy, curtailment, and congestion depend on the coincidence of demand, inflow, renewable availability, storage state, and transmission bottlenecks. Aggregation may have little effect in an energy-abundant year and a large effect during a short winter scarcity event. It may shift rather than simply increase congestion by changing where and when hydro is dispatched. A priori reasoning therefore cannot establish whether the abstraction changes a planning conclusion.

British Columbia (BC) is an informative test case. Hydropower supplies approximately 90% of provincial electricity in a typical year [8], and much of the controllable output is concentrated in large regulated cascades. On the Peace River, water stored in Williston Reservoir is released through G.M. Shrum and then passes Peace Canyon and Site C. On the Columbia River system, large upstream storage and generating facilities interact with downstream plants and cross-border obligations. Smaller cascades, including Bridge, Stave, Campbell, and Seven Mile, provide additional examples of linked storage and generation. The province is also experiencing rapid demand growth and major transmission investment pressure [9–11]. These characteristics make BC both a policy-relevant case and a stringent test of hydro representation.

Existing work establishes the system value of hydropower flexibility [2,4–6] and identifies missing hydropower detail in planning and reliability models [7]. Regional PyPSA applications generally use aggregated hydrological inputs rather than station-resolved cascades [12–16,20]. Other open frameworks can represent linked reservoirs, and previous BC studies have developed production-cost and capacity-expansion models with varying levels of hydro detail [17–19,21–23]. Mauel [19] provides the closest BC precedent by representing multiple cascade units and comparing facility output with observed records. However, the literature does not provide a controlled, within-model estimate of the bias caused by replacing the same cascade with an equivalent reservoir while all other inputs remain fixed.

This study asks: **How much does hydropower aggregation alter estimates of operational feasibility in a hydro-dominant electricity system under electrification and variable renewable energy integration?** We answer this question through three matched configurations. Configuration A resolves reservoirs and generating stations, conserves water through turbine and spill pathways, represents station-to-station travel time, applies minimum-flow constraints where defensible data exist, and enforces a terminal storage condition. Configuration B represents each major basin using an upper and lower store with limited coupling. Configuration C collapses each basin into one energy bucket without station-to-station coupling. We compare the configurations across a common scenario matrix using unserved energy, curtailment, congestion, cost, spill, and storage indicators.

The contribution is methodological rather than a claim that maximum physical detail is always preferable. First, we define an auditable PyPSA cascade formulation in which water and electricity balances can be tested separately. Second, we use a controlled experimental design that attributes differences to hydropower representation instead of conflating them with alternative demand or investment assumptions. Third, we establish a validation hierarchy suited to data-limited hydro systems and distinguish parameters supported by public records from assumptions. Fourth, we report both operational outcomes and computational cost, allowing users to judge whether a two-store approximation captures the decision-relevant behavior of a full cascade.

The manuscript is organized as follows. Section 2 reviews hydropower representation and defines the unresolved research gap. Section 3 describes the BC-PyPSA case, cascade equations, experimental configurations, scenarios, validation, and outcome metrics. Section 4 provides a results structure with explicit analytical placeholders. Section 5 sets out the interpretation, transferability, and limitations that should be assessed after the experiment is complete. Section 6 concludes.

## 2. Background and research contribution

### 2.1 Hydropower flexibility and the importance of chronology

The value of reservoir hydropower arises from its ability to move generation through time. Aubin et al. [2] describe how Québec's hydro fleet can shift from energy supply toward balancing as variable renewable production grows. Arvesen et al. [4] show that flexible Norwegian hydropower and transmission can jointly manage renewable variability, while HydroWIRES studies emphasize the contribution of hydropower to firm capacity and integration of variable resources [5]. System-level optimization studies likewise find that reservoir flexibility can reduce cost and facilitate wind deployment [6]. These benefits depend on chronological operation: the state of storage after one hour constrains every subsequent hour.

A cascade adds spatial chronology to this temporal dependence. Let an upstream release at hour \(t\) arrive downstream at hour \(t+\tau\). The downstream station cannot use that water before the lag has elapsed, and an upstream decision changes the feasible downstream schedule. Multiple tributaries and stations create a directed network of such relationships. An equivalent-reservoir formulation preserves a basin-wide energy budget but normally removes the directed water network. It can therefore allocate generation across stations and hours without reproducing the physical sequence that creates that energy.

The distinction between energy adequacy and deliverability is central. A system may possess enough annual water to serve annual demand but still experience hourly shortfalls because discharge, turbine, travel-time, or transmission constraints bind simultaneously. Conversely, an environmental release may increase generation or spill during low-demand hours and reduce storage available for later peaks. Chronological dispatch at hourly resolution is consequently required to measure the operational manifestations of aggregation, even when long-term capacity-expansion decisions are outside scope.

### 2.2 Representation in open power-system models

PyPSA provides linear optimization components for generators, stores, storage units, links, and networks [24]. Regional implementations such as PyPSA-Eur and PyPSA-Earth have enabled transparent, reproducible analysis at continental scale [12–14]. Their hydropower inputs are necessarily simplified to match broad spatial coverage and public data. Related applications for Brazil, Vietnam, and Canada demonstrate the versatility of the framework but do not establish a validated station-resolved, delay-aware cascade for a hydro-dominant jurisdiction [15,16,20].

Some non-PyPSA frameworks provide richer hydro structures. openTEPES represents multi-reservoir networks through directed relationships [17]. OSeMOSYS applications have represented station sequences for long-term planning, although representative time slices cannot reproduce all chronological water transfers [21]. SILVER and COPPER support detailed Canadian production-cost and capacity-planning studies [18,22,23]. Mauel's extension of SILVER represents 25 units in 12 BC cascade groups and validates selected stations against available generation records [19]. That work shows both the feasibility and the data difficulty of detailed BC hydro modeling.

The relevant gap is not simply that one software package lacks a feature. Three evidence gaps remain across the reviewed literature. No published experiment compares station-resolved and aggregated reservoirs in the same hourly nodal model with identical non-hydro inputs. No published PyPSA application validates a delay-aware cascade at jurisdiction scale in a near-fully hydro system. Finally, no study reports the error in operational-feasibility indicators caused by aggregation as a function of renewable, demand, and hydrological stress.

**[TABLE 1 ABOUT HERE: Cascade-hydropower capabilities across the reviewed open-source frameworks. Verify each entry against the cited software version and publication before submission. Use “documented,” “partial,” “not documented,” and “not applicable” rather than checkmarks when evidence is ambiguous.]**

### 2.3 Experimental logic

Model-detail studies can become circular if the most detailed model is assumed to be true. We instead treat Configuration A as the reference representation and validate as many of its intermediate quantities as the data allow. The experiment then measures the difference created by controlled simplification. A difference from Configuration A is an aggregation effect, not automatically an error relative to the real system. Only comparison with observed operation supports an accuracy claim.

Three hypotheses guide the analysis. **H1:** relative to the full cascade, the single energy bucket understates unserved energy during stressed hours because it makes basin-wide stored energy instantaneously deliverable. **H2:** the energy bucket understates renewable curtailment when physically constrained hydro cannot reduce or shift output as freely as the aggregate representation assumes. **H3:** aggregation changes the location and duration of transmission congestion by redistributing hydro dispatch among electrical buses. A fourth, practical hypothesis is that the two-store representation captures a large share of the full-cascade effect at lower data and computational cost.

These hypotheses are directional only where the feasible-set relationship is clear. System outcomes can remain unchanged if the omitted hydro constraints never bind. Congestion on an individual corridor can increase or decrease because dispatch is spatially redistributed. Results are therefore reported as effect sizes and scenario distributions, not solely as hypothesis confirmations.

## 3. Methods

### 3.1 BC-PyPSA system boundary

BC-PyPSA is an hourly nodal dispatch model of the British Columbia electricity system built using PyPSA [24]. The present study uses a one-year chronological horizon and preserves the model's existing electrical network, load disaggregation, generator fleet, trade interfaces, and renewable profiles. The current processed hydro inventory contains 118 generating assets totaling approximately 17.73 GW, including 16 assets classified as reservoir or reservoir-imputed and 102 assets classified as run-of-river or water-linked run-of-river. The processed reservoir table contains 21 reservoirs across the Peace, Mica/Columbia, Bridge, Stave, Campbell, Seven Mile, and Kootenay groups. These inventory counts describe the audited input files and may differ from the final solved network after asset aggregation, filtering, or skipped-asset checks.

The base year is **[MODEL PLACEHOLDER M1: state the final weather, load, and operating year]**. The temporal resolution is one hour, with **[M2: number of snapshots after leap-day and time-zone handling]** snapshots. All time series are converted to one documented time zone before alignment. Snapshot weights equal one hour in the core experiment. The electrical network contains **[M3: final bus, line, transformer, and controllable-link counts]**, and load is disaggregated to **[M4: spatial unit and number of load buses]** using **[M5: data source and weighting method]**. Trade is represented as **[M6: fixed exchange, price-responsive link, or bounded generator/load]**. Backstop generators with a high marginal cost preserve mathematical feasibility and measure unserved energy.

The study is an operational comparison, not a capacity-expansion optimization. Renewable portfolios, conventional capacity, and transmission ratings are fixed within each matched scenario. This choice prevents the optimizer from compensating for a simplified hydro representation by building a different system. A separate planning extension may use the validated representation to estimate investment consequences, but that extension is outside the core causal comparison.

### 3.2 Input data and provenance

The hydro topology is assembled from the processed BC-PyPSA hydro generation, reservoir, and cascade tables. Each generating asset has an electrical connection, installed capacity, hydraulic upstream and downstream identifiers where known, discharge limits, spill limits, and a cascade group. Water Use Plan records and manually curated reservoir files provide storage bounds and monthly inflow statistics for a subset of reservoirs [25]. HydroBASINS catchments and weather-derived runoff are used to construct the temporal shape of natural inflow. For reservoirs with public monthly inflow statistics, the modeled hourly runoff series is scaled within each month to match the reported mean. For identified downstream reservoirs without a separate inflow series, local inflow is currently imputed as zero; this assumption is exposed as a validation flag and sensitivity case rather than silently treated as observed.

Run-of-river availability is derived from upstream catchment runoff and calibrated to facility annual average energy subject to installed capacity. The existing preprocessing estimates basin-to-plant timing from HydroBASINS main-channel distance and an assumed flow speed. This meteorological routing shapes natural inflow before optimization. It is distinct from the operational station-to-station delay introduced below: the first routes runoff from catchments to an asset, while the second delays controlled turbine discharge and spill between modeled reservoirs.

Demand uses the BC-PyPSA provincial load profile and spatial disaggregation workflow. Provincial totals are validated against BC Hydro balancing-authority load where available. Wind and solar profiles are derived from the same weather year and mapped to proposed or existing resource locations. Generator technical and cost parameters are taken from the project's processed asset tables and generic technology file. Every scenario records checksums or modification times for these inputs, the active configuration, solver, package versions, and code revision.

**[DATA PLACEHOLDER D1: insert a complete data-source table with source organization, variable, native resolution, modeled transformation, license/access condition, coverage, and validation role.]**

**[DATA PLACEHOLDER D2: report missingness and imputation counts for storage, minimum discharge, maximum discharge, spill capacity, travel time, annual energy, and inflow statistics. Identify all assets excluded from the solved model and give reasons.]**

### 3.3 Station-resolved cascade formulation

The cascade is represented as a directed graph \(G=(R,E)\), where \(R\) is the set of reservoirs and terminal outlets and \(E\) is the set of turbine and spill pathways. Each physical reservoir is a PyPSA Store connected to a water bus. Natural inflow is a fixed time-series injection. A generating station is a multi-output Link whose input is water from the upstream station bus, whose first output is electricity at the associated electrical bus, and whose second output is the same discharged water at the downstream water bus. A parallel spill Link transfers water without electricity generation. A terminal water Store with sufficiently large capacity receives water leaving the modeled cascade and prevents an artificial disposal constraint.

For reservoir \(r\) and hour \(t\), the water balance is

\[
S_{r,t}=S_{r,t-1}+I_{r,t}
+\sum_{g\in U_r}Q_{g,t-\tau_g}
+\sum_{k\in U_r}W_{k,t-\tau_k}
-\sum_{g\in D_r}Q_{g,t}
-\sum_{k\in D_r}W_{k,t},
\]

where \(S_{r,t}\) is stored water, \(I_{r,t}\) is local natural inflow, \(Q_{g,t}\) is turbine discharge, \(W_{k,t}\) is spill, \(U_r\) denotes upstream pathways terminating at \(r\), \(D_r\) denotes pathways leaving \(r\), and \(\tau\) is an integer travel delay in hours. All flow and storage quantities use one internally consistent water unit. If the implementation uses cubic metres per hour as link power and cubic metres as Store energy, one hourly snapshot makes the dimensional correspondence explicit. For non-hourly weighted snapshots, the balance must include snapshot duration.

Electric output from station \(g\) is

\[
P_{g,t}=\eta_g Q_{g,t},
\]

where \(\eta_g\) is a fixed water-to-electricity conversion in MWh per cubic metre for an hourly time step. It is calculated from rated electrical capacity and rated discharge when both are available. This approach linearizes the effects of head and turbine efficiency. The same \(\eta_g\) is used in all three configurations to isolate the effect of cascade structure. A sensitivity case perturbs \(\eta_g\) within a documented range or substitutes head-specific values where reservoir level and head data are defensible.

Discharge and storage satisfy

\[
0 \le Q_{g,t}\le \overline{Q}_g,\qquad
0 \le W_{k,t}\le \overline{W}_k,\qquad
\underline{S}_{r,t}\le S_{r,t}\le \overline{S}_{r,t}.
\]

The lower and upper storage bounds may vary by month when defensible rule-curve data exist. Where only physical extremes are known, the model uses those extremes and labels the result as a physical-bound rather than operating-policy representation.

The current BC-PyPSA workflow already creates water buses, reservoir Stores, inflow injections, multi-output discharge Links, spill Links, and cyclic storage. The controlled study adds station-to-station delays and operating-rule constraints in a dedicated scenario layer. Delay is implemented through explicit linear constraints that equate the receiving flow at time \(t\) with the sending flow at \(t-\tau\), with boundary treatment documented for the first \(\tau\) hours. This avoids relying on a software-version-specific component attribute and allows water-balance tests to inspect sending and receiving series directly.

**[IMPLEMENTATION PLACEHOLDER I1: insert the final PyPSA and Linopy versions, the exact variable names used in the delayed-flow constraints, and the boundary convention. State whether pre-horizon flows are taken from a warm-up run, wrapped cyclically, or fixed from observations.]**

### 3.4 Operating constraints

Minimum environmental flow is imposed at each station or controlled river reach for which a defensible requirement is available:

\[
Q_{g,t}+W_{g,t}\ge \underline{F}_{g,t}.
\]

The right-hand side may be constant, monthly, or time-varying. Because a legal release may be defined at a downstream gauge rather than a single facility, each constraint is mapped to the combination of model pathways that controls the relevant reach. Units, averaging periods, exceptions, and seasonal windows are retained from the source rather than converted to an hourly minimum without justification.

Maximum change in turbine discharge is represented as

\[
-\Delta_g^- \le Q_{g,t}-Q_{g,t-1}\le \Delta_g^+.
\]

Ramping constraints are included in the core model only where source data are available and clearly refer to a compatible temporal measure. Generic technology ramp rates are not treated as observed facility water-operation limits. Where facility data are unavailable, ramping is a sensitivity test and not part of the reference claim.

Monthly reservoir constraints take the form \(\underline{S}_{r,m(t)} \le S_{r,t}\le \overline{S}_{r,m(t)}\). Treaty or license obligations are represented only when they can be translated into a documented minimum flow, target release, or storage band over a specified interval. The 1964 Columbia River Treaty remains distinct from the 2024 agreement-in-principle; the final model must not represent a proposed term as a binding current rule [26,27].

At the annual boundary, each reservoir uses a cyclic terminal condition,

\[
S_{r,T}=S_{r,0}.
\]

This removes the incentive to empty reservoirs at the end of the modeled year. It does not identify the correct initial storage level because a cyclic optimization may choose that level endogenously. The preferred validation run fixes or bounds initial storage using observed data and includes a warm-up period; the controlled comparison uses an identical terminal convention across configurations. A terminal water-value sensitivity is included because very large reservoirs can carry value across years [28].

**[TABLE 2 ABOUT HERE: Operating constraints, legal or technical source, temporal form, assets covered, transformation, uncertainty class, and inclusion status. Separate “core,” “sensitivity,” and “not implemented.”]**

### 3.5 Hydropower configurations

Configuration A is the station-resolved reference. It contains one Store for each modeled reservoir, one turbine pathway for each modeled generating station, explicit downstream water transfer, travel delay, station discharge and spill limits, cyclic terminal storage, and minimum-flow or reservoir constraints supported by the evidence hierarchy. Natural inflow is assigned to the incremental catchment between adjacent reservoirs so that upstream runoff is not counted again downstream.

Configuration B is a two-store approximation. Each major basin has an upper and lower Store. Reservoir energy or water volume, inflow, and turbine capacity are aggregated to these two locations using a published mapping. The stores are coupled by one transfer pathway. Station-level travel times and reach-level minimum flows are removed unless a basin-wide rule can be represented without creating a different total obligation. Configuration B tests whether limited spatial disaggregation captures decision-relevant behavior with lower data requirements.

Configuration C is one equivalent energy bucket per basin. Aggregate inflow is converted to energy using a capacity-weighted or energy-weighted conversion factor; aggregate storage and turbine capacity equal the mapped basin totals after duplicate checks. The bucket has no internal station sequence or travel delay. Its electrical output is connected using **[MODEL PLACEHOLDER M7: specify whether capacity and energy are distributed among existing electrical buses or assigned to one representative bus]**. This choice is critical: connecting all basin output to a single bus changes both hydraulic and electrical aggregation. The preferred design preserves station electrical locations using dispatch Links drawing from the common bucket, so that Configuration C removes hydraulic structure without erasing the network geography.

All configurations conserve the same aggregate natural inflow, initial/terminal water condition, storage volume or equivalent energy, installed generating capacity, and water-to-energy assumptions to the maximum extent possible. An automated equivalence report compares these totals before solving. If aggregation changes any total, the scenario fails validation rather than proceeding silently.

**[TABLE 3 ABOUT HERE: Exact features of Configurations A–C, including stores, turbine pathways, spill pathways, delays, environmental flows, rule curves, terminal condition, electrical interconnection treatment, and aggregation formulas.]**

The A–B–C comparison changes several mechanisms at once. To attribute the Configuration A–C difference, we also run one-at-a-time ablations: A without travel delay; A without minimum-flow constraints; A with all station stores but instantaneous coupling; and A with relaxed seasonal bounds. These runs are secondary because their number can be large, but they are necessary if the manuscript claims which mechanism drives the aggregate effect.

### 3.6 Scenario design

The core experiment crosses hydropower representation with demand, renewable, and hydrological conditions. Demand includes **[SCENARIO PLACEHOLDER S1: define baseline, medium-electrification, and high-electrification cases, including annual energy and peak GW]**. Renewable portfolios include **[S2: define low, medium, and high wind/solar cases using installed MW, annual available TWh, and share of load]**. Hydrology includes **[S3: define median, dry, and wet years or inflow scalars, preserving chronology]**. A full three-by-three-by-three design produces 27 non-hydro stress combinations and 81 core solves across A, B, and C. If data or computation limits the design, the final manuscript reports the retained combinations and the reason for omission.

The comparison uses matched pairs: scenario \(s\) has identical non-hydro inputs in A, B, and C. This pairing supports direct effect sizes \(\Delta Y_{B,s}=Y_{B,s}-Y_{A,s}\) and \(\Delta Y_{C,s}=Y_{C,s}-Y_{A,s}\). The distribution across stress cases shows whether aggregation bias is concentrated in particular conditions. No probabilistic interpretation is attached to the scenarios unless explicit probabilities are assigned.

Sensitivity tests address assumptions not fully supported by public data: travel-time values, flow speed used in natural-inflow routing, missing local inflow at downstream reservoirs, conversion efficiency, initial storage, terminal condition, environmental-flow coverage, spill limits, and the value of lost load. Solver tolerances and a second solver are checked for a representative subset. Each sensitivity changes one parameter family while retaining matched A–C pairs.

### 3.7 Optimization and reproducibility

For each scenario, the model minimizes total variable operating cost, start and shutdown cost where active, trade cost, and the penalty assigned to unserved energy, subject to electricity balance, transmission, generator, water-balance, storage, and operating-rule constraints. The objective does not assign an arbitrary reward to water transfer. Spill receives either zero cost or a small documented tie-breaking penalty that is too small to alter economic dispatch relative to genuine costs.

The workflow uses configuration files rather than hard-coded scenario names. A scenario manifest stores the run identifier; Git commit or source checksum; package and solver versions; input hashes; hydropower configuration; demand, renewable, and hydrology case; constraint switches; solve status; objective; runtime; peak memory if measured; and output path. Outputs are written to immutable run-specific directories. A post-processing script reads the manifest and solved NetCDF files, applies common metric definitions, and writes tidy Parquet or CSV tables used by every figure and manuscript result.

Numerical quality checks include solver termination status, primal and dual feasibility where available, the maximum hourly electricity-balance residual, maximum reservoir water-balance residual, negative storage or flow beyond tolerance, capacity violations, and terminal-condition residual. A run is excluded from comparative figures only according to predeclared failure rules; failed runs remain visible in the manifest.

### 3.8 Validation framework

Validation is hierarchical because no single observed dataset covers the full water-power network. Level 1 tests structure. Every nonterminal turbine and spill pathway must lead to a valid downstream water bus; every reservoir must be reachable from a natural inflow or an upstream reservoir; cascade graphs must be acyclic except for explicit annual storage boundary treatment; and every hydro electrical output must connect to a valid BC bus. Duplicate IDs, self-loops, orphan stores, out-of-province coordinates, and nonpositive capacities are reported.

Level 2 tests conservation and dimensions. Synthetic two- and three-reservoir networks have analytically known solutions and verify delayed arrival, mass conservation, spill, minimum flow, terminal storage, and first-hour boundary behavior. In full runs, the hourly water-balance residual is calculated independently from exported dispatch. Annual water entering each cascade must equal terminal outflow plus storage change within tolerance. Aggregated configurations must conserve the corresponding total energy or water.

Level 3 validates inputs. Modeled provincial load is compared with BC Hydro balancing-authority observations using annual energy, peak demand, peak timestamp, monthly energy, load-duration curves, root-mean-square error, mean absolute error, bias, and correlation. Runoff-derived inflow is compared with public monthly inflow statistics used for calibration and, where independent data exist, with held-out records. Annual run-of-river generation potential is checked against facility average annual energy without reusing the same quantity as an independent validation claim.

Level 4 validates hydro operation. Facility or plant-group generation is compared with observed hourly, daily, or monthly generation for the stations and years that can be obtained. Metrics include energy bias, normalized root-mean-square error, correlation, Kling–Gupta efficiency or Nash–Sutcliffe efficiency where hydrological interpretation is appropriate, ramp distributions, seasonal timing, and peak coincidence. Reservoir level or storage trajectories are compared with observed records after converting both series to a common normalized range if volume-elevation curves are unavailable. Environmental-flow constraints are checked at the temporal resolution of the source obligation.

Level 5 validates system outcomes. Modeled provincial energy balance, imports and exports, technology generation, peak net load, and transmission-flow patterns are compared with independent reports. Because public line-flow observations may be unavailable, congestion validation may be qualitative or limited to known interfaces. Such evidence is labeled corroboration rather than calibration.

Data are divided into calibration and validation roles before model tuning. Monthly inflow statistics used to scale a runoff series cannot subsequently be cited as independent evidence that the same monthly mean is correct. If only one year of operational data is available, event-based cross-validation with held-out months is used where feasible; otherwise the limitation is stated directly.

**[VALIDATION PLACEHOLDER V1: provide a station-by-station validation coverage matrix and identify which observed series were unavailable.]**

**[VALIDATION PLACEHOLDER V2: insert acceptance thresholds established before reviewing comparative results. At minimum include water-balance residual, annual load bias, monthly inflow bias, hydro energy bias, storage timing, and solve quality.]**

### 3.9 Outcome metrics and statistical analysis

Unserved energy is the annual electricity supplied by designated backstop generators:

\[
EENS_s=\sum_t\sum_{g\in B}p_{g,t,s}w_t,
\]

where \(B\) is the backstop set and \(w_t\) is snapshot duration. We also report loss-of-load hours, maximum hourly shortfall, and the timestamp and location of the maximum shortfall. Results are insensitive to backstop marginal cost only if that cost remains above all genuine supply options; this is tested.

Available variable renewable energy is \(\sum p^{\mathrm{avail}}_{g,t}\), calculated from capacity, availability, and snapshot duration. Curtailment is available energy minus nonnegative dispatched energy. The annual curtailment ratio divides curtailed energy by available energy, and is reported separately for wind and solar as well as jointly. Storage charging or negative dispatch conventions are handled explicitly so they cannot produce negative curtailment.

Line loading is absolute active power divided by the applicable thermal rating. A corridor is congested when loading is at least **[METRIC PLACEHOLDER K1: define 1.0 p.u. or a justified operational threshold]**. We report congested hours by corridor, energy transmitted while constrained, maximum loading, and shadow price if the dual is reliable. Parallel circuits are either retained individually or aggregated using one consistent corridor mapping.

Additional metrics include total system operating cost; hydro generation by station and cascade; turbine utilization; spill volume; storage range and terminal residual; renewable utilization; net imports; emissions if active; and solver runtime. Computational burden is summarized because a representation that changes no decision-relevant outcome but multiplies solve time may not be justified.

**[TABLE 4 ABOUT HERE: Outcome definitions, units, model signals, aggregation rules, and paired comparison measures. Replace provisional backstop and congestion conventions with the final component registry and threshold.]**

For each scenario, paired effects compare B and C with A. Across scenarios, results are shown as medians, ranges, and distributions rather than as repeated independent observations. A factorial regression or analysis of variance may estimate interactions between representation and stress factors, but scenario cases are deterministic. Confidence intervals are used only for uncertainty ensembles or bootstrap resampling with a defensible sampling unit, not as decoration around deterministic runs.

## 4. Results

### 4.1 Model assembly and validation

**[RESULT PLACEHOLDER R1—ASSEMBLY: Report final counts of electrical buses, lines, loads, generators, hydro reservoirs, turbine pathways, spill pathways, cascade groups, and time steps. Report processed-input counts separately from instantiated-network counts. List exclusions and duplicates resolved.]**

**[RESULT PLACEHOLDER R2—CONSERVATION: Report the maximum and 99th-percentile hourly water-balance residual by configuration, electricity-balance residual, terminal storage residual, and any constraint violations. State tolerance and units. Include the synthetic-network test results.]**

**[RESULT PLACEHOLDER R3—INPUT VALIDATION: Report annual modeled and observed provincial load, peak GW and timestamp, monthly load errors, inflow calibration performance, and run-of-river energy checks. Distinguish calibration matches from independent validation.]**

**[RESULT PLACEHOLDER R4—OPERATIONAL VALIDATION: For each observed facility or reservoir, report temporal coverage, energy bias, RMSE or normalized RMSE, correlation, seasonal timing, and storage-level performance. Explain mismatches before presenting A–C effects.]**

**[FIGURE 1 ABOUT HERE: BC map showing modeled electrical network, load centers, reservoirs, generating stations, and named cascades. Use BC Albers for analysis, remove axes, and size hydro markers by installed capacity.]**

**[FIGURE 2 ABOUT HERE: Auditable component diagram for one cascade reach, showing natural inflow, reservoir Store, turbine discharge, electricity output, delayed downstream water, spill, and minimum-flow constraint.]**

**[FIGURE 3 ABOUT HERE: Validation dashboard with observed versus modeled provincial load, facility generation, inflow, and reservoir storage. Do not combine quantities with incompatible units on one axis.]**

### 4.2 Operational effects of hydropower representation

**[RESULT PLACEHOLDER R5—PRIMARY TABLE: For every core scenario and configuration, report EENS in MWh, loss-of-load hours, curtailment in GWh and percent, total congestion hours, operating cost, hydro generation, spill, and runtime. Include solve status.]**

The primary comparison should begin with absolute outcomes in Configuration A before reporting differences. This prevents a large percentage change from obscuring a small baseline and shows whether the detailed reference itself is operationally acceptable. For scenario \(s\), the aggregation effect is \(Y_{C,s}-Y_{A,s}\); the two-store effect is \(Y_{B,s}-Y_{A,s}\).

**[RESULT PLACEHOLDER R6—UNSERVED ENERGY: State whether Configuration C produced lower, equal, or higher EENS than A. Report the median paired difference, maximum difference, and the demand–VRE–hydrology case in which it occurred. Identify the hours, buses, and cascades associated with the largest hidden shortfall. Report results in MWh and as a fraction of annual load.]**

**[RESULT PLACEHOLDER R7—CURTAILMENT: Report wind, solar, and total renewable curtailment for A–C. Give percentage-point differences, not only relative percentages. Identify whether changes arose from minimum releases, constrained storage, delayed downstream generation, transmission limits, or a combination.]**

**[RESULT PLACEHOLDER R8—CONGESTION: Report the number of corridors whose congestion changed materially, the five largest corridor-level changes, and whether aggregation shifted congestion spatially. Distinguish a change in total congested corridor-hours from a change in the number of unique congested hours.]**

**[FIGURE 4 ABOUT HERE: Paired effect-size plot for EENS, curtailment ratio, and congestion hours across all matched scenarios. Use one panel per metric and connect A, B, and C within scenarios.]**

**[FIGURE 5 ABOUT HERE: Feasibility frontier or heat map across demand and renewable stress, faceted by hydrology and hydropower representation. Color should encode EENS or loss-of-load hours using a common scale.]**

**[FIGURE 6 ABOUT HERE: BC congestion maps for A and C plus a difference map. Use the same projection, bounds, line-width scale, and legend in all panels.]**

### 4.3 Event-based explanation

Annual metrics do not reveal how aggregation changes operation. The final analysis therefore selects two or three events using predeclared criteria: the maximum Configuration C–A shortfall difference, the largest curtailment difference, and the largest congestion redistribution. For each event, a common time window shows demand, renewable availability, system net load, station generation, turbine discharge, spill, storage, and delayed arrival through the relevant cascade.

**[RESULT PLACEHOLDER R9—EVENT MECHANISM: Explain one event by tracing water from upstream storage through each affected station. State which constraint first diverged between A and C, how long the divergence persisted, and how it propagated to generation, imports, curtailment, congestion, or load shedding. Verify the explanation from exported constraint and dispatch series rather than inference from annual totals.]**

**[FIGURE 7 ABOUT HERE: Multi-panel event timeline for the Peace or Columbia cascade. Align all panels on a shared hourly axis and mark travel-time offsets and binding constraints.]**

**[FIGURE 8 ABOUT HERE: Reservoir-by-reservoir storage and spill heat map for the same event, ordered upstream to downstream.]**

### 4.4 Value of the two-store approximation and mechanism ablations

Configuration B is useful only if its loss of detail can be quantified. Define recovery for metric \(Y\) as

\[
Recovery_{B,s}=1-\frac{|Y_{B,s}-Y_{A,s}|}{|Y_{C,s}-Y_{A,s}|},
\]

when the denominator exceeds a predeclared small threshold. Values near one indicate that B recovers the A outcome; values below zero indicate that B is farther from A than C. This statistic is reported alongside absolute differences because it becomes unstable when A and C are nearly identical.

**[RESULT PLACEHOLDER R10—TWO-STORE: Report recovery distributions for each primary metric, scenario classes in which B failed, and the reduction in component count, runtime, and memory relative to A.]**

**[RESULT PLACEHOLDER R11—ABLATIONS: Attribute the A–C difference among delay, minimum flow, store aggregation, rule curves, and terminal treatment using one-at-a-time runs. Report interactions if the sum of individual effects does not reproduce the combined effect.]**

**[FIGURE 9 ABOUT HERE: Accuracy–complexity plot with component count or runtime on the x-axis and normalized operational error on the y-axis for A, B, C, and ablations.]**

### 4.5 Sensitivity and robustness

**[RESULT PLACEHOLDER R12—SENSITIVITY: Report the effect of travel-time uncertainty, downstream local-inflow imputation, efficiency, initial storage, terminal condition, environmental-flow coverage, spill limits, value of lost load, and solver choice. Identify assumptions that can reverse the sign or materially change the magnitude of the primary effects.]**

**[FIGURE 10 ABOUT HERE: Sensitivity tornado or interval plot showing the change in the A–C effect for EENS, curtailment, and congestion. Use absolute units and show the core case as a reference line.]**

## 5. Discussion

### 5.1 Interpreting aggregation bias

The central interpretation must be proportional to the validation evidence. If Configuration C produces less unserved energy than A during matched stress cases, the result shows that hydraulic aggregation creates additional modeled flexibility. If the relevant A behavior is also supported by observed reservoir or station operation, the difference can be described as an optimistic bias. If operational validation is weak, the correct wording is that the result is a structural sensitivity: a realistic constraint set changes the feasibility estimate, but the detailed model has not been proven to reproduce all real operations.

**[DISCUSSION PLACEHOLDER P1: Insert the evidence-backed interpretation of the primary effect. Link each causal claim to an event plot, binding constraint, or ablation result.]**

Aggregation can change outcomes through several mechanisms. A common bucket pools water that is physically separated among reservoirs. Instantaneous transfer lets downstream capacity respond before upstream water arrives. Removing minimum flow can retain water during low-value hours, increasing its apparent availability during later scarcity. Aggregating electrical delivery at one location can bypass a transmission bottleneck. These mechanisms should not be conflated. In particular, an energy bucket connected to one representative bus tests a compound hydraulic-and-electrical simplification, while a common bucket feeding the original station buses more cleanly isolates hydraulic aggregation.

The system may also be insensitive to cascade detail over broad operating ranges. When inflow and storage are abundant, transmission is unconstrained, and turbine limits do not bind, all configurations can produce similar dispatch. This is not evidence that cascade physics are unimportant; it identifies the conditions under which aggregation is adequate for the selected outputs. A useful model-detail recommendation therefore states a domain of validity, such as a range of net-load stress, reservoir state, and renewable penetration, rather than declaring one representation universally sufficient.

### 5.2 Implications for planning

Operational bias matters when it changes a decision. Understated EENS can reduce the apparent need for firm capacity, demand response, imports, or transmission. Understated curtailment can exaggerate the utilization and economic value of remote renewable resources. Redistributed congestion can change the ranking of network reinforcements. The core experiment fixes investment to identify these operational effects; a subsequent capacity-expansion experiment would be required to quantify changes in optimal investment.

**[DISCUSSION PLACEHOLDER P2: Identify any scenario in which A, B, and C imply different planning classifications—for example, acceptable versus unacceptable reliability, a changed priority corridor, or a changed renewable-utilization threshold. State the criterion used for the classification.]**

If Configuration B consistently approaches A while reducing data and computation, it offers a practical intermediate representation. Such a result would be valuable in regions where reservoir locations, capacity, and broad basin topology are public but station-level rules are not. If B fails specifically during short scarcity events, a hybrid approach may retain explicit detail only for cascades whose travel time, storage distribution, or electrical location is operationally consequential.

### 5.3 Transferability

The component mapping and validation logic are transferable beyond BC. Required inputs are a directed reservoir–station graph, storage bounds, turbine and spill capacities, natural inflow, water-to-energy conversion, station electrical buses, and travel-time estimates. Operating rules can be added progressively and tagged by evidence quality. The method is particularly relevant to systems where the same water is reused through multiple large stations and where reservoir geography interacts with a constrained transmission network.

Transfer does not imply that BC parameters can be reused elsewhere. Travel time, operating rules, head variation, market structure, and water obligations are jurisdiction specific. The reusable contribution is the schema, constraint layer, comparison design, tests, and reporting format. A new application should first reproduce synthetic conservation tests and then validate against local observations before interpreting an A–C difference as real-world bias.

### 5.4 Limitations

The first limitation is data access. Public Water Use Plans provide valuable storage and flow information for selected facilities, but they do not constitute a complete hourly operational dataset. Facility generation, reservoir elevation, spill, travel time, and binding operating rules may be unavailable or inconsistent across years. The final paper must publish a coverage table and avoid presenting imputed parameters as observations.

Second, power conversion is linearized at a fixed head and efficiency. Real power output varies with forebay and tailwater elevation, turbine loading, and unit availability. Holding conversion constant across configurations supports identification of topology effects but may misstate absolute generation. A piecewise-linear head model is a useful extension after the cascade and validation pipeline are stable.

Third, deterministic perfect foresight can overstate operational flexibility in every configuration. Operators do not know future inflow, load, wind, and solar without error. Because forecast uncertainty interacts with reservoir positioning, aggregation bias under rolling-horizon or stochastic dispatch may differ from the perfect-foresight result. The present study should therefore be interpreted as a structural comparison under common information, not a complete simulation of operator behavior.

Fourth, an annual cyclic condition prevents terminal depletion but may not represent multi-year water value. Large reservoirs can carry energy across years, and drought sequences may make a one-year return condition restrictive or optimistic depending on the selected year. Multi-year or terminal-water-value sensitivities bound this issue.

Fifth, representation error outside hydropower remains. Public transmission data require network aggregation, load is spatially disaggregated, plant outages may be stylized, and trade may be simplified. These assumptions are held fixed across A–C and therefore do not mechanically create the paired difference, but they can change whether a hydro constraint binds and thus modify effect magnitude.

Finally, Configuration A is a reference, not a digital twin. The manuscript should reserve causal language for mechanisms demonstrated by constraints and event traces, and reserve accuracy claims for quantities compared against independent observations. Transparent uncertainty labels are more scientifically useful than an unsupported appearance of precision.

## 6. Conclusions

This study establishes a reproducible design for measuring how hydropower aggregation changes operational-feasibility estimates in a hydro-dominant system. It compares a station-resolved, delay-aware cascade with two-store and single-bucket alternatives while preserving non-hydro inputs. The formulation keeps water and electricity accounting explicit, incorporates operating constraints according to evidence quality, and validates the model from component topology through observed system behavior.

**[CONCLUSION PLACEHOLDER C1: Insert three quantitative conclusions: the largest or median aggregation effect on EENS, curtailment, and congestion; the stress condition in which it emerged; and the performance of the two-store approximation. Include units.]**

**[CONCLUSION PLACEHOLDER C2: State the practical model-detail recommendation and its domain of validity. If validation remains incomplete, describe the outcome as structural sensitivity rather than measured real-world error.]**

The broader lesson is that annual hydro energy is not synonymous with hourly deliverable capacity. When water must pass through ordered stations, delayed river reaches, operating constraints, and a spatially constrained grid, a basin-wide energy bucket may grant flexibility that the physical system does not possess. Quantifying that difference—and the conditions under which it matters—is necessary before simplified hydro models are used to support high-stakes electrification and renewable-integration decisions.

## Data and code availability

BC-PyPSA source code and public workflow materials are available at https://github.com/DeltaE/PyPSA_BC. **[AVAILABILITY PLACEHOLDER A1: insert the archived release DOI or commit hash used for the paper, exact commands to reproduce the scenario matrix, public input-data manifest, and locations of derived metric tables.]** Some facility operating records may be subject to provider restrictions. **[A2: identify restricted datasets, the process for requesting access, and the public or synthetic substitutes included in the reproducibility package.]**

## Author contributions

**[AUTHOR PLACEHOLDER: use the CRediT taxonomy, including conceptualization, methodology, software, validation, formal analysis, data curation, visualization, writing—original draft, writing—review and editing, supervision, and funding acquisition.]**

## Competing interests

The authors declare **[PLACEHOLDER: no competing interests / specify competing interests]**.

## Acknowledgements

**[ACKNOWLEDGEMENT PLACEHOLDER: funding, data providers, model contributors, and computing resources. Do not imply endorsement by BC Hydro or another organization without permission.]**

## References

[1] Wolinetz M. British Columbia Electrification Impacts Study. Vancouver: Navius Research Inc.; 2020.

[2] Aubin V, Blais M, Anjos MF. How Québec can support the energy transition of northeastern North America. The Electricity Journal 2021;34:106972. https://doi.org/10.1016/j.tej.2021.106972.

[3] Fingrid, Energinet, Statnett, Svenska Kraftnät. Nordic Grid Development Perspective 2023. 2023.

[4] Arvesen A, Hansen OM, Harby A, et al. On the potential role of flexible Norwegian hydropower in managing challenging renewable energy variability events in Europe. IOP Conference Series: Earth and Environmental Science 2025;1442:012003. https://doi.org/10.1088/1755-1315/1442/1/012003.

[5] HydroWIRES, Stark G, Brinkman G. The Role of Hydropower Flexibility in Integrating Renewables in a Low-Carbon Grid. 2023. https://doi.org/10.2172/2000741.

[6] Wang Y, Levin T, Kwon J, et al. The value of hydropower flexibility for electricity system decarbonization. Energy Reports 2025;13:2711–2721. https://doi.org/10.1016/j.egyr.2025.02.019.

[7] Mitra B, Datta S, Kincic S, et al. Gaps in Representations of Hydropower Generation in Steady-State and Dynamic Models. 2023. https://doi.org/10.48550/arXiv.2311.03277.

[8] BC Hydro. BC Hydro quick facts 2024–25. Vancouver, BC: BC Hydro; 2025.

[9] BC Hydro. BC Hydro making it easier and more affordable to connect new homes for people. 2024. https://www.bchydro.com/news/press_centre/news_releases/2024/distribution-extension-policy-updates.html.

[10] BC Hydro. Our plan for B.C.'s future. https://www.bchydro.com/toolbar/about/strategies-plans-regulatory/supply-operations/clean-energy-action-plan.html.

[11] Government of British Columbia. New legislation powers economy with clean energy, North Coast Transmission Line. BC Gov News; 2025. https://news.gov.bc.ca/releases/2025ECS0044-001032.

[12] Hörsch J, Hofmann F, Schlachtberger D, et al. PyPSA-Eur: An open optimisation model of the European transmission system. Energy Strategy Reviews 2018;22:207–215. https://doi.org/10.1016/j.esr.2018.08.012.

[13] Parzen M, Abdel-Khalek H, Fedotova E, et al. PyPSA-Earth: A new global open energy system optimization model demonstrated in Africa. Applied Energy 2023;341:121096. https://doi.org/10.1016/j.apenergy.2023.121096.

[14] Abdel-Khalek H, Schumm L, Jalbout E, et al. PyPSA-Earth sector-coupled: A global open-source multi-energy system model showcased for hydrogen applications in countries of the Global South. Applied Energy 2025;383:125316. https://doi.org/10.1016/j.apenergy.2025.125316.

[15] Deng Y, Cao K-K, Hu W, et al. Harmonized and Open Energy Dataset for Modeling a Highly Renewable Brazilian Power System. Scientific Data 2023;10:103. https://doi.org/10.1038/s41597-023-01992-9.

[16] [REFERENCE PLACEHOLDER: complete and verify the PyPSA-VN conference citation and DOI.]

[17] Ramos A, Alvarez EF, Lumbreras S. OpenTEPES: Open-source Transmission and Generation Expansion Planning. SoftwareX 2022;18:101070. https://doi.org/10.1016/j.softx.2022.101070.

[18] Miri M, Saffari M, Arjmand R, et al. Integrated models in action: Analyzing flexibility in the Canadian power system toward a zero-emission future. Energy 2022;261:125181. https://doi.org/10.1016/j.energy.2022.125181.

[19] Mauel J. Impact of Cascaded Hydro Operational Constraints on Power System Flexibility Requirements for Variable Renewable Energy Integration. 2023. **[REFERENCE PLACEHOLDER: insert thesis institution, degree, repository URL, and access date.]**

[20] Natural Resources Canada. PyPSA-Canada. 2026. **[REFERENCE PLACEHOLDER: insert version, repository DOI or URL, and access date.]**

[21] Fejzić E, Niet T, Wade C, et al. Aligning the Western Balkans power sectors with the European Green Deal. Environmental Research Communications 2024;6:115008. https://doi.org/10.1088/2515-7620/ad8ca4.

[22] McPherson M, Karney B. A scenario-based approach to designing electricity grids with high variable renewable energy penetrations in Ontario, Canada: Development and application of the SILVER model. Energy 2017;138:185–196. https://doi.org/10.1016/j.energy.2017.07.027.

[23] Arjmand R, McPherson M. Canada's electricity system transition under alternative policy scenarios. Energy Policy 2022;163:112844. https://doi.org/10.1016/j.enpol.2022.112844.

[24] Brown T, Hörsch J, Schlachtberger D. PyPSA: Python for Power System Analysis. Journal of Open Research Software 2018;6:4. https://doi.org/10.5334/jors.188.

[25] BC Hydro. Water Use Plans. https://www.bchydro.com/toolbar/about/sustainability/environmental_responsibility/water-use-plans.html.

[26] Government of British Columbia. Columbia River Treaty. 2024. https://engage.gov.bc.ca/columbiarivertreaty/agreement-in-principle/.

[27] U.S. Department of State. Summary of the Agreement in Principle to Modernize the Columbia River Treaty Regime. https://2021-2025.state.gov/summary-of-the-agreement-in-principle-to-modernize-the-columbia-river-treaty-regime/.

[28] Pereira MVF, Pinto LMVG. Multi-stage stochastic optimization applied to energy planning. Mathematical Programming 1991;52:359–375. https://doi.org/10.1007/BF01582895.

[29] IEA PVPS. Recommended Practices for Wind/PV Integration Studies. 2024. **[REFERENCE PLACEHOLDER: verify report title, authors, report number, and publication year against the cited PDF.]**

[30] Bird L, Lew D, Milligan M, et al. Wind and solar energy curtailment: A review of international experience. Renewable and Sustainable Energy Reviews 2016;65:577–586. https://doi.org/10.1016/j.rser.2016.06.082.

## Supplementary material plan

**Supplement S1:** Full input-data dictionary and provenance manifest.

**Supplement S2:** Reservoir–station–downstream topology table with coordinates, capacities, delays, and evidence classes.

**Supplement S3:** Mathematical formulation and PyPSA/Linopy variable mapping.

**Supplement S4:** Synthetic cascade test cases and expected solutions.

**Supplement S5:** Full scenario manifest, solver diagnostics, and failed-run log.

**Supplement S6:** Validation metrics and plots for every facility with observations.

**Supplement S7:** Complete results table for every scenario and configuration.

**Supplement S8:** Sensitivity and one-at-a-time ablation results.
