# Cascade Reservoir Representation and Operational Feasibility in Hydro-dominant System

## Introduction

Electrification is changing how much electricity hydro-dominant power systems must deliver and when they must deliver it. Such systems include British Columbia, Québec, and Norway, where regulated reservoirs supply most of the annual electricity. Across these systems, peak demand is rising, seasonal swings are widening, and load is concentrating in fewer locations [1–3]. Planners assess whether a grid can meet this changing demand by running dispatch models that schedule generator output hour by hour. The accuracy of this assessment depends on how well the model represents the physical structure of the reservoirs that supply most of the power. A dispatch model that aggregates these reservoirs into a single equivalent store dissolves that physical structure. The grid can then appear able to meet demand it cannot deliver. A planner relying on such a model could assess a grid as adequate when, under electrification, it would shed load. This chapter quantifies that gap, and how it grows as variable renewable energy (VRE) takes on more of the system's balancing.

Hydropower reservoirs give hydro-dominant grids much of their flexibility. They hold water back when generation exceeds demand. They release it when demand exceeds generation. In doing so they absorb the variability that wind and solar add to the system. This balancing role is well documented in several studies [2,4–6]. Aubin et al. [2] showed that Québec’s hydro fleet shifted from baseload supply toward system balancing as variable renewables grew. Arvesen et al. [4] reported that hydro flexibility and cross-border transmission together reduce curtailment. The same balancing role is the organizing premise of the U.S. HydroWIRES program [5], where dispatchable hydro flexibility is estimated to provide roughly 24 Gigawatts (GW) of firm capacity and to support about 140 GW of variable renewables. Wang et al. [6] showed that reservoir flexibility lowers system cost while allowing more wind. In a hydro-dominant system, the available flexibility is set by the basin’s hydrology and by its operating rules. Operating rules include the reservoir-level limits, minimum-discharge requirements, and seasonal release schedules that govern when water can be stored and released. In operational dispatch models such as PyPSA (Python for Power System Analysis), the available flexibility is set instead by how the reservoirs are represented. PyPSA is an open-source framework for optimizing generator, storage, and network dispatch.

Most reservoirs in hydro-dominant systems sit in cascades, and the cascade is where a dispatch model's representation tends to depart from the physical system. A cascade is a sequence of reservoirs and generating stations along a river, where water released from one station becomes the inflow to the next after a travel-time delay. The most common simplification collapses a basin's reservoirs and stations into a single equivalent reservoir. Such single equivalent reservoir tracks its contents as stored energy rather than stored water, and schedules releases without regard to the order of stations along the river. This chapter refers to this practice as the single-reservoir aggregation, consistent with what is elsewhere called an equivalent-reservoir or lumped-reservoir representation. Such aggregation drops four features of a real cascade as follows:

Inter-reservoir coupling and the travel-time delay between stations disappear. Water released from one dam becomes inflow to the next after a defined delay, which constrains the basin’s whole generation schedule. Ignoring it allows dispatch solutions that the physical system could not realize.

Head-dependent power generation efficiency is replaced by a fixed conversion factor. Such factor biases estimate of storage value and generation cost.

Minimum environmental flow requirements are omitted. Omitting these requirements removes a legal lower bound on discharge that constrains dispatch regardless of cost.

End-of-horizon water valuation per reservoir is not set. Such setup lets the dispatch model’s optimizer drain reservoirs as the planning horizon ends.

Each of these omissions makes the model treat the cascade as more flexible than it really is. Mitra et al. [7], writing on steady-state and dynamic reliability models, document two of these effects directly.  Neglecting inter-reservoir interdependencies and assuming constant turbine efficiency both distort a model's estimate of available capacity. Omitting minimum environmental flow or end-of-horizon storage valuation produces the same direction of bias by construction — each removes a limit that would otherwise reduce deliverable generation. Their combined effect on system feasibility has not been quantified for any near-100%-hydro jurisdiction. Electrification makes this gap pressing. It is pushing peak demand higher, increasing how fast seasonal ramping must occur, and concentrating new demand in specific regions of the hydro-dominant grids. Under these conditions a single-reservoir aggregation treats more generation as dispatchable than the physical cascade can deliver. This chapter asks how much that overstatement matters. It measures how far estimates of grid performance move when a cascade is replaced by a single-reservoir aggregation.

British Columbia is selected as a test case because it combines near-100% hydro supply with large, regulated cascades under active electrification pressure. Approximately 90% of provincial electricity supply is hydro-based [8], delivered through two regulated cascades large enough to be analytically tractable at province scale, yet complex enough to expose the limits of basin-level aggregation. Figure 1 shows a simplified overview of the cascades. The Peace River system runs the G.M. Shrum, Peace Canyon, and John Horgan (site-C) stations in series below Williston Reservoir. The Columbia system links the Mica, Revelstoke, and Arrow Lakes plants. This hydro base is now under active electrification pressure: BC Hydro's $36 billion ten-year capital plan [9,10] responds to rising peak demand and variable renewable energy (VRE) integration pressure under the CleanBC electrification mandate [11]. Getting the cascade representation right therefore has direct bearing on how that capital is allocated.

Figure : Major hydropower cascades in BC.

The remainder of this chapter is organized as follows. Section 2.2 reviews how existing open-source models represent cascade hydropower and identifies what they leave out. Section 2.3 describes the methods of cascade representation and its operating-rule constraints. Section 2.4 reports the expected results, and Section 2.5 discusses their implications for planning in hydro-dominant systems.

## Literature Review

No published open-source operational modelling framework represents a hydropower cascade in full. Most omit inter-reservoir coupling, travel-time delay, and minimum-flow limits, and instead drive each basin with a single aggregated inflow. This review uses six capabilities to distinguish a full cascade representation from the aggregated one: multi-reservoir hydraulic coupling, travel-time delay, head-dependent generation efficiency, minimum environmental flow constraints, end-of-horizon water valuation, and published validation in a hydro-dominant system. Table 1 applies these six capabilities to the frameworks most comparable to BC-PyPSA.

Table : Cascade hydro capability across eight open-source operational modelling frameworks. 
✓ = implemented; ~ = partially implemented or linearized; ✗ = absent.

The PyPSA regional models share a limitation that matters for this chapter: none represents a hydropower cascade. The five PyPSA-based regional implementations (PyPSA-Eur, PyPSA-Earth, PyPSA-Brazil, PyPSA-VN, and PyPSA-Canada), all suppress inter-reservoir coupling, travel-time delay, and minimum-flow constraints. These models rely on aggregated basin-level inflow. PyPSA-Eur, PyPSA-Earth, and PyPSA-Canada documentation explicitly acknowledges this limitation. This is a limitation shared across every PyPSA regional implementation reviewed, not an artifact of any single model. Mitra et al. [7] reviewed hydropower representations across steady-state and dynamic models in several jurisdictions and find the same suppression of cascade coupling, head-dependent efficiency, and environmental flow constraints. The omission is a systemic feature of large-scale modelling practice, not a PyPSA-specific oversight.

A few open-source frameworks do represent a cascade, but none within the PyPSA ecosystem and none validated as an operational dispatch model for a near-100%-hydro system. openTEPES [17] represents multi-reservoir networks as a directed graph with arc-flow constraints, yet no study has applied it at jurisdiction scale to a near-100%-hydro system, and no bias estimate exists. Fejzić et al. [21] model the Drina River Basin, Western Balkan Region, cascade station by station in OSeMOSYS, but as a long-term capacity-expansion model built on non-sequential representative time steps. That structure cannot carry travel-time delay or chronological storage dynamics, and the system studied is fossil fuel dominated rather than hydro-dominated.

British Columbia has the most developed cascade prior art, yet none of it answers the question this chapter raises. McPherson and Karney (2017) developed SILVER [22] , a standalone production cost and economic dispatch model for provincial-scale hourly electricity analysis, not PyPSA-based. Arjmand and McPherson (2022) developed COPPER [23], a capacity expansion model implemented in Python using the Pyomo library and the CPLEX solver, also not PyPSA-based. Miri et al. [18] coupled these tools into COPPER-SILVER, linking capacity planning with hourly dispatch to quantify how curtailment ratios, hydro storage utilization, and dispatch cost shift as VRE penetration grows across the BC-Alberta corridor. Mauel [19] is the exception and the closest prior work. Mauel extended SILVER with cascaded-hydro representation, represented the BC system as 25 units in 12 cascade groups, and compared its output against BC Hydro generation records for ten facilities. Mauel did not do three things this chapter requires. The model is not built in PyPSA. It couples reservoirs within the same hour, with travel time assumed under one hour for every river because BC delay data were unavailable, so it carries no explicit delay. It never sets the cascade against its aggregated form, so the feasibility bias of aggregation stays unmeasured.

Three gaps remain, and although British Columbia is the case study, each applies to other hydro-dominated systems. The first gap is the absence of a cascade representation in the PyPSA ecosystem. The PyPSA regional models in Table 1 aggregate each basin into a single inflow, and the three open-source frameworks that do represent a cascade, openTEPES, the OSeMOSYS Drina model, and the SILVER extension, are all built outside PyPSA. The second gap is the absence of a validated, delay-aware cascade in a near-100%-hydro jurisdiction. The one validated BC cascade model [19] treats inter-reservoir transport as instantaneous and matches BC Hydro records only coarsely at the facility level. The third gap is the absence of any estimate of how much single-reservoir aggregation overstates operational feasibility. No study has compared a cascade against its aggregated counterpart in the same model to measure that bias. This chapter addresses all three. It integrates a delay-aware cascade into BC-PyPSA, validates it against BC Hydro operations, and measures the feasibility bias using unserved energy, curtailment, and congestion.

## Methods

BC-PyPSA is an hourly nodal dispatch model built on the open-source PyPSA framework [24]. It minimizes total operating cost subject to generation, storage, and transmission constraints over an 8,760-hour annual horizon. The network has 28 nodes aligned with BC Hydro's bulk transmission corridors and RESource siting categories. This resolution matches BC Hydro's transmission planning zones and is partly comparable to BC Hydro’s integrated resource planning. Finer disaggregation is not possible, because sub-station topology at hourly resolution is not public.

Hydropower in BC-PyPSA is represented as coupled reservoir cascades, not as basin-level energy buckets. This study advances that representation and then tests how much the added detail matters. It adds an explicit travel-time delay to the transfer between reservoirs, using the native PyPSA v1.2.0 delay attribute, where the base model currently couples stations within the same hour. It adds the operating rules that govern real reservoir releases, beginning with a minimum environmental flow, as explicit constraints. It then measures the value of this detail by comparing the full cascade against two deliberately simplified versions of the same model. All other model components stay unchanged.

### Cascade Reservoir Formulation

A cascade is built from two kinds of PyPSA component. A Store holds a commodity and tracks its level over time. Here each Store is a reservoir, and its level is the water in storage, expressed in megawatt-hours at the reservoir's mean head. A Link moves a commodity between two buses and can convert it on the way. Here Links carry water between reservoirs and stations and convert water into electricity.

Each reservoir obeys a water-balance constraint. Storage at the end of an hour equals storage at the start, plus inflow, minus the water released through the station and the water spilled. Natural and upstream inflow enters through a per-reservoir Inflow Generator. Every term in this balance is linear, so the cascade keeps the model a linear program.

Figure : PyPSA component mapping blueprint for cascaded units

Each cascade unit couples a reservoir to the station immediately downstream through a delayed water transfer. Figure 2 shows the component schematic of one cascade unit. This structure repeats along each river. Each generating station is modelled as four Links. A Release Link moves water from the reservoir Store to the station bus, and a Store Link returns unused water to the Store. A Discharge Link is the multi-port component. Its first output delivers electricity to the grid bus, with a conversion factor equal to mean turbine efficiency times mean hydraulic head, in megawatt-hours per cubic metre. Its second output passes the discharged water to the downstream reservoir Store. A Spill Link carries water released without generating power to that same downstream Store.

The model adds a travel-time delay that the base model lacks. The travel-time delay attribute is configured using the native delay attribute introduced in PyPSA v1.2.0. This attribute lags the arriving water by a set number of hours. It also adds a minimum environmental flow, a lower bound on the combined turbine and spill discharge at each station, to enforce the legally required release. Neither feature is in the base network, and both are part of the work proposed here

The contribution of this section is not new mathematics. It is the integration of an explicit travel-time delay and enforced operating constraints into a PyPSA cascade at jurisdiction scale, in a reusable form. The results will be validated against BC Hydro operations with available reservoir records.

### Operating-Rule Constraint Layer

Operating rules enter the model as bounds on the Store and Link variables defined above. Four rules apply to any regulated hydro-dominant system. Table 2 lists the four rules, data sources for BC and the status of the data.

Table : Operating Rules

The minimum environmental flow is a lower bound on the combined turbine and spillway flow at each station, in absolute discharge units, which guarantees the legally required release. The maximum ramping rate limits how much a station's flow can change between two consecutive hours. The reservoir-level limits, also called a refill constraint, is a lower bound on reservoir level that varies by month and reserves water for later in the year. Treaty and licence obligations are mandatory minimum releases over defined contractual windows.

The constraint structure defined by these four rules (Table 2) is generic and applies to any regulated cascade. The values that instantiate them for British Columbia are case-study parameters, and their availability is uneven. The constraint structure itself does not change for another jurisdiction. Completing the BC values is part of the calibration phase, not a settled input.

### Terminal Storage Condition

Without a terminal storage condition, a finite-horizon dispatch model has no incentive to retain water past the last hour. A model prevents this by constraining end-of-year storage to equal its starting level. Each reservoir Store uses PyPSA's cyclic state-of-charge option (e_cyclic), which forces the storage level in the final hour to equal the level in the first hour. The starting level itself is set by the optimization rather than fixed to an observed value. This conserves water across the annual cycle. Without such a condition, a finite-horizon linear dispatch empties its reservoirs in the final hours, because stored water has no value once the horizon ends [28]. The cyclic condition removes that artefact. It differs from the seasonal rule curves of Section 2.3.2, which bind in every hour, because it binds only at the final hour.

The terminal condition is adequate only when storage does not carry value across years; it is a standard PyPSA feature, documented here for completeness rather than as a contribution. A single annual return is inadequate when storage carries over across several years, for example during multi-year droughts or in systems with very large reservoirs. In those cases, a piecewise-linear water-value function can price stored water at the horizon end. That extension needs multi-year stochastic inflow data and is outside the scope of this study.

### Outcome Metrics and Operational Feasibility

The comparison quantifies operational feasibility using three metrics recommended by the IEA Wind and IEA-PVPS expert-group report on wind and solar integration studies [29]. These metrics are unserved energy, curtailment, and congestion. Unserved energy is the demand that cannot be met without shedding load, measured in megawatt-hours per year. Curtailment is the share of available wind and solar output that the network cannot absorb in a given hour Congestion is the number of hours per year in which a transmission corridor reaches its rated transfer limit [30].

Table : Outcome metrics: definition, unit, model signal, and comparison measure

All three metrics describe demand and network outcomes that the model produces as its hourly solution. A dispatch model reproduces them only as well as it represents the constraints that produce them. These constraints are generation limits, storage balances, network flows, and reserve requirements. A model solution is mathematically feasible when it satisfies these constraints. The BC-PyPSA model is allowed to shed load at a high penalty cost rather than being forbidden to shed it. A solution can therefore satisfy every constraint and still leave demand unserved. Such a solution is mathematically feasible but not operationally acceptable. A single-reservoir aggregation can reach a low-cost feasible solution that assumes more deliverable generation than the physical cascade could provide, which hides shortfalls the real system would face. Measuring that hidden shortfall is the central aim of this chapter.

### Three-Configuration Comparison Study

The value of the cascade representation is measured by comparing three versions of the model that differ only in how they represent hydropower. Each version simplifies the representation by one further step, so the comparison isolates the effect of that simplification. Table 4 lists what changes across configurations A, B, and C: spatial resolution, inter-reservoir coupling, travel-time delay, and the minimum-flow constraint.

The three configuration set-up isolates the effect of cascade detail on feasibility outcomes.  Configuration A is the full cascade of Sections 2.3.1 - 2.3.3 and serves as the reference. Configuration C is the single-reservoir aggregation, one Store per basin with no inter-station coupling, matching the representation used in prior BC modelling. Configuration B lies between them. It keeps two Stores per basin, an upper and a lower, with limited coupling between them, but it removes the travel-time delay, the minimum-flow constraint, and the storage detail within each river reach.

Table : Structural differences in hydro representation across Configurations A, B, and C.

All three configurations run on identical VRE portfolios, demand profiles, and network topologies, so any difference in outcomes comes from the hydro representation alone. The configurations are compared on three differences. Unserved-energy error is the annual difference in megawatt-hours between a simplified configuration and Configuration A. Curtailment bias is the difference in the annual curtailment ratio, in percentage points. Congestion-hour deviation is the difference in the annual count of hours in which a transmission corridor sits at its rated limit. Configuration B tests whether a coarse two-store split recovers most of the accuracy of the full cascade at lower data cost. The stress tests span at least three VRE penetration levels, which shows whether the bias grows as renewable share rises.

## Expected Results

The three-configuration comparison is designed to show how much of a hydro-dominant grid's apparent flexibility disappears once travel-time delay and minimum-flow limits are enforced. First, a quantified feasibility-bias estimate: the difference in unserved energy, curtailment, and transmission congestion between Configuration A and Configuration C. The estimation is evaluated across VRE penetration levels. Configuration C, a single aggregate store per basin with no inter-station coupling, is the coarsest point in this comparison. This configuration informs how much accuracy hydro flexibility loses when reduced to its simplest possible form; it is a more aggressive simplification than the existing operational BC-PyPSA model. The existing model retains multi-link station-level coupling but does not represent travel-time delay and operating-policy constraints. The existing model has already identified congestion concentrated on the Fraser-Thompson-Okanagan corridor and the broader north-south corridor linking northern generation to southern load centres. Figure 3 and Figure 4 illustrate this congestion pattern in the existing single-representation BC-PyPSA model. Figure 3 maps the transmission network and shows where is the congestion, with line colour indicating the number of hours each corridor operates above 90 percent of its rated capacity. Figure 4 shows loading duration curves for the five most congested corridors, several of which remain above that threshold for nearly the entire 8,760-hour year.

Figure : Transmission line loading across the BC network in the existing single-representation model (preliminary results). Line colour indicates the number of hours each corridor operated at or above 90% of its rated capacity, out of 8,760 hours simulated. Node colour indicates each regional district's dominant generation technology. Hydro-dominant nodes connect to the grid through discharge links rather than the AC lines shown here, so their apparent low connectivity does not indicate isolation.

Figure : Top five congested transmission corridors under the existing operational hydro representation (preliminary results of an aggressive electrification scenario)

The two preliminary findings show that congestion is not dispersed. It concentrates on a small number of corridors loaded near their limit almost continuously. These results come from the model's current hydro representation, prior to the comparison proposed in this chapter. The Configuration A, B, and C runs  will test whether this concentrated congestion pattern holds, weakens, or shifts under coarser representation of the cascade.

Together, these outputs are expected to answer the central research question: how much operational bias the energy-bucket simplification introduces. They also supply the feasibility signals used in the soft-linking architecture described in Chapter 3.

## Discussion

The comparison is expected to show that the energy-bucket abstraction overstates operational flexibility, with the overstatement growing with VRE penetration. The bias is expected to appear strongly in curtailment: basin aggregation allows the model to schedule generation across stations in ways that are physically impossible when upstream storage is constrained, so enforcing inter-reservoir coupling raises curtailment estimates. If Configuration B recovers most of the accuracy gain of Configuration A, it provides a practical middle tier for systems lacking station-level operational data. This result is transferable to any hydro-dominant jurisdiction where basin geometry and mean flow records are publicly available but hourly station-level data are not.

The existing regional PyPSA models that underpin continental energy planning, rely on basin-level hydro aggregation. If the bias quantified here is large enough to affect investment conclusions in cascade-dominated regions, it motivates a systematic reassessment of hydro representation in those models. This chapter does not claim to resolve that question; it establishes the bias magnitude for one well-documented cascade system and provides the reusable implementation needed to extend the test to others.

For BC planning, a material cascade-versus-bucket bias would mean basin-level tools understate the flexibility cost of high-VRE scenarios. That, in turn, would overstate the competitiveness of remote wind and solar against near-load alternatives. Two limits bound these conclusions: a modeling assumption, and a data dependency. First, head-dependent efficiency is linearized at mean reservoir level; the resulting generation error varies across the operating range, is documented in the validation experiment, and is not corrected in the comparison runs, so the experimental contrast remains attributable to coupling structure alone. Second, the validation experiment depends on BC Hydro providing hourly reservoir and generation records. If that access is not confirmed, the comparison remains a structural contrast without an observed-operations benchmark. The broader contribution is integrating cascade physics into the PyPSA ecosystem in a reusable form.

An energy-bucket model does not know that water takes hours to travel between dams — and that missing hour is exactly where deliverable capacity disappears. Closing that gap matters for BC because it separates a capital plan built on the grid's real flexibility from one built on flexibility a simplified model only assumes it has. It matters beyond BC because the same aggregation shortcut sits inside every PyPSA regional model reviewed in this chapter — which is why the cascade built here is designed to travel with it.

## References

[1]	Wolinetz M. British Columbia Electrification Impacts Study. Vancouver: Navius Research Inc.; 2020.

[2]	Aubin V, Blais M, Anjos MF. How Québec can support the energy transition of northeastern North America. Electr J 2021;34:106972. https://doi.org/10.1016/j.tej.2021.106972.

[3]	Nordic Grid Development Perspective 2023. FINGRID, ENERGINET, Statnett, SVENSKA KRAFTNAT; 2023.

[4]	Arvesen A, Hansen OM, Harby A, et al. On the potential role of flexible Norwegian hydropower in managing challenging renewable energy variability events in Europe. IOP Conf Ser Earth Environ Sci 2025;1442:012003. https://doi.org/10.1088/1755-1315/1442/1/012003.

[5]	Hydropower and Water Innovation for a Resilient Electricity System (HydroWIRES), Stark G, Brinkman G. The Role of Hydropower Flexibility in Integrating Renewables in a Low-Carbon Grid. 2023. https://doi.org/10.2172/2000741.

[6]	Wang Y, Levin T, Kwon J, et al. The value of hydropower flexibility for electricity system decarbonization. Energy Rep 2025;13:2711–21. https://doi.org/10.1016/j.egyr.2025.02.019.

[7]	Mitra B, Datta S, Kincic S, et al. Gaps in Representations of Hydropower Generation in Steady-State and Dynamic Models 2023. https://doi.org/10.48550/arXiv.2311.03277.

[8]	BC Hydro quick facts 2024-25. Vancouver, BC: BC Hydro; 2025.

[9]	BC Hydro making it easier and more affordable to connect new homes for people 2024. https://www.bchydro.com/news/press_centre/news_releases/2024/distribution-extension-policy-updates.html (accessed October 28, 2025).

[10]	Our plan for B.C.’s future n.d. https://www.bchydro.com/toolbar/about/strategies-plans-regulatory/supply-operations/clean-energy-action-plan.html (accessed October 28, 2025).

[11]	Solutions E and C. New legislation powers economy with clean energy, North Coast Transmission Line. BC Gov News 2025. https://news.gov.bc.ca/releases/2025ECS0044-001032 (accessed November 30, 2025).

[12]	Hörsch J, Hofmann F, Schlachtberger D, et al. PyPSA-Eur: An open optimisation model of the European transmission system. Energy Strategy Rev 2018;22:207–15. https://doi.org/10.1016/j.esr.2018.08.012.

[13]	Parzen M, Abdel-Khalek H, Fedotova E, et al. PyPSA-Earth. A new global open energy system optimization model demonstrated in Africa. Appl Energy 2023;341:121096. https://doi.org/10.1016/j.apenergy.2023.121096.

[14]	Abdel-Khalek H, Schumm L, Jalbout E, et al. PyPSA-Earth sector-coupled: A global open-source multi-energy system model showcased for hydrogen applications in countries of the Global South. Appl Energy 2025;383:125316. https://doi.org/10.1016/j.apenergy.2025.125316.

[15]	Deng Y, Cao K-K, Hu W, et al. Harmonized and Open Energy Dataset for Modeling a Highly Renewable Brazilian Power System. Sci Data 2023;10:103. https://doi.org/10.1038/s41597-023-01992-9.

[16]	PyPSA-VN: An open model of the Vietnamese electricity system | IEEE Conference Publication | IEEE Xplore n.d. https://ieeexplore.ieee.org/document/9303096 (accessed June 15, 2026).

[17]	Ramos A, Alvarez EF, Lumbreras S. OpenTEPES: Open-source Transmission and Generation Expansion Planning. SoftwareX 2022;18:101070. https://doi.org/10.1016/j.softx.2022.101070.

[18]	Miri M, Saffari M, Arjmand R, et al. Integrated models in action: Analyzing flexibility in the Canadian power system toward a zero-emission future. Energy 2022;261:125181. https://doi.org/10.1016/j.energy.2022.125181.

[19]	Mauel J. Impact of Cascaded Hydro Operational Constraints on Power System Flexibility Requirements for Variable Renewable Energy Integration 2023.

[20]	NRCan/pypsa-canada 2026.

[21]	Fejzić E, Niet T, Wade C, et al. Aligning the Western Balkans power sectors with the European Green Deal. Environ Res Commun 2024;6:115008. https://doi.org/10.1088/2515-7620/ad8ca4.

[22]	McPherson M, Karney B. A scenario based approach to designing electricity grids with high variable renewable energy penetrations in Ontario, Canada: Development and application of the SILVER model. Energy 2017;138:185–96. https://doi.org/10.1016/j.energy.2017.07.027.

[23]	Arjmand R, McPherson M. Canada’s electricity system transition under alternative policy scenarios. Energy Policy 2022;163:112844. https://doi.org/10.1016/j.enpol.2022.112844.

[24]	Brown T, Hörsch J, Schlachtberger D. PyPSA: Python for Power System Analysis. J Open Res Softw 2018;6:4. https://doi.org/10.5334/jors.188.

[25]	Water Use Plan n.d. https://www.bchydro.com/toolbar/about/sustainability/environmental_responsibility/water-use-plans.html (accessed June 15, 2026).

[26]	Columbia River Treaty 2024. https://engage.gov.bc.ca/columbiarivertreaty/agreement-in-principle/ (accessed June 15, 2026).

[27]	Summary of the Agreement in Principle to Modernize the Columbia River Treaty Regime. U S Dep State n.d. https://2021-2025.state.gov/summary-of-the-agreement-in-principle-to-modernize-the-columbia-river-treaty-regime/ (accessed June 15, 2026).

[28]	Pereira MVF, Pinto LMVG. Multi-stage stochastic optimization applied to energy planning. Math Program 1991;52:359–75. https://doi.org/10.1007/BF01582895.

[29]	IEA-PVPS-T14-28-2024-REPORT-Wind-PV-Integration n.d. https://iea-pvps.org/wp-content/uploads/2025/01/IEA-PVPS-T14-28-2024-REPORT-Wind-PV-Integration.pdf (accessed October 28, 2025).

[30]	Bird L, Lew D, Milligan M, et al. Wind and solar energy curtailment: A review of international experience. Renew Sustain Energy Rev 2016;65:577–86. https://doi.org/10.1016/j.rser.2016.06.082.

## Extracted table 1

Models | Multi-res coupling | Travel delay | Head-dep gen | Min-flow | End-of-horizon value | Validated (hydro-dom.) | Notes
PyPSA-Eur [12] | ✗ | ✗ | ✗ | ~ | ✗ | ✗ | Basin inflow aggregated; cascade limitation acknowledged in docs
PyPSA-Earth [13,14] | ✗ | ✗ | ✗ | ~ | ✗ | ✗ | Same core architecture as PyPSA-Eur
PyPSA-Brazil [15] | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | Aggregated inflow; no cascade coupling
PyPSA-VN [16] | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | Aggregated inflow; no cascade coupling
openTEPES [17] | ✓ | ~ | ~ | ✓ | ✓ | ✗ | Native cascade-basin topology; no jurisdiction-scale validation; no bias estimate
SILVER [18,19] | ✓ | ✗ | ~ | ✓ | ~ | ✓ | SILVER extension, not PyPSA; instantaneous inter-reservoir coupling; 25 units / 12 cascade groups; validated against BC Hydro WUP records.
PyPSA-Canada [20] | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | NRCan/CanmetENERGY PyPSA-based framework; applied to Atlantic Loop and Saskatchewan scenarios; aggregated hydro representation inherited from PyPSA core; no cascade topology documented
BC-PyPSA  / (this work) | ✓ | ✓ | ~ | ✓ | ✓ | Planned | Full cascade, head at mean; delay via extra functionality, validation pending BC Hydro data

## Extracted table 2

Operating rule | BC data source | Status
Minimum environmental flow | Water Use Plans and supporting technical reports | Available for a subset of facilities; not yet compiled across the full cascade
Maximum ramping rate | BC Hydro operating specifications | May not be publicly available; a data challenge
Reservoir-level limits | Water Use Plans (e.g., Walter Hardman, Whatshan [25]) | Available for specific reservoirs
Treaty / licence obligations | Columbia River Treaty, 1964 instrument [26,27] | 1964 treaty in force; the 2024 Canada–US agreement is a negotiating framework, not a ratified amendment

## Extracted table 3

Metric | What it measures | Unit | How the model reports it | Comparison measure (vs Configuration A)
Unserved energy | Demand that cannot be met without shedding load | MWh per year | Nonzero dispatch on the $1,000/MWh backstop generators | Unserved-energy error: annual MWh difference
Curtailment | Share of available wind and solar output the network cannot absorb in an hour | % of available VRE | Available VRE output not dispatched | Curtailment bias: difference in the annual curtailment ratio, in percentage points
Congestion | Hours in which a transmission corridor is at its rated transfer limit | hours per year | Corridor flow at its rated limit | Congestion-hour deviation: difference in annual congestion hours

## Extracted table 4

Modelling Component | Configuration A  / (full cascade, reference) | Configuration B* / (two-store) | Configuration C (energy-bucket)
Spatial resolution | One Store per reservoir, one multi-port Link per generating station | Two Stores per basin: upper-basin, lower-basin | One Store per basin
Inter-reservoir coupling | Full, station-by-station cascade | Limited, upper-basin to lower-basin only | None
Station-level travel-time delay | Yes | No | No
Minimum-flow constraint | Yes (combined turbine and spillway flow, Section 2.3.2) | No | No
*Note: Configuration B simplifies travel-time delay, the minimum-flow constraint, and storage granularity together rather than one at a time. Its results therefore show the effect of partial disaggregation, not which single mechanism drives any difference from Configuration A. A one-at-a-time design, with one intermediate configuration per mechanism, would separate them at the cost of more model runs and is outside the scope of this study. | *Note: Configuration B simplifies travel-time delay, the minimum-flow constraint, and storage granularity together rather than one at a time. Its results therefore show the effect of partial disaggregation, not which single mechanism drives any difference from Configuration A. A one-at-a-time design, with one intermediate configuration per mechanism, would separate them at the cost of more model runs and is outside the scope of this study. | *Note: Configuration B simplifies travel-time delay, the minimum-flow constraint, and storage granularity together rather than one at a time. Its results therefore show the effect of partial disaggregation, not which single mechanism drives any difference from Configuration A. A one-at-a-time design, with one intermediate configuration per mechanism, would separate them at the cost of more model runs and is outside the scope of this study. | *Note: Configuration B simplifies travel-time delay, the minimum-flow constraint, and storage granularity together rather than one at a time. Its results therefore show the effect of partial disaggregation, not which single mechanism drives any difference from Configuration A. A one-at-a-time design, with one intermediate configuration per mechanism, would separate them at the cost of more model runs and is outside the scope of this study.
