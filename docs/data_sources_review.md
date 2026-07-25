# Benchmarking Open Canadian Energy-System Datasets Against CODERS

## TL;DR
- For a BC-focused OSeMOSYS/PyPSA/CLEWs workflow, CODERS is not the best available source in any *single* raw-data category — but it is the best pre-integrated, model-ready Canadian package. In 5 of 6 categories a better open primary source exists on at least one criterion (spatial/temporal resolution, recency, or license), so the correct framing is **pipeline-vs-frozen-product**: build a reproducible RESource pipeline on ERA5/RDRS(CaSR) and utility primary data, and use CODERS as a validation baseline and gap-filler.
- The single biggest wins over CODERS are: (1) VRE capacity factors — replace the frozen MERRA-2 0.5°×0.625° grid with ECCC's RDRS/CaSR (~10 km, hourly, OGL-equivalent license) or ERA5 (0.25°), both redistributable; and (2) hourly demand — pull directly from BC Hydro's Balancing Authority load feed, AESO, IESO and EIA-930, which are more recent and continuously updated than CODERS' 2018–2022/23 harmonized tables.
- The categories where CODERS remains hardest to beat are its **node-resolved transmission network** (voltage, reactance, seasonal ratings) and its **hydro renewal/greenfield cost compilations from regulatory filings** — no open dataset matches these, so here CODERS wins and should be retained, with OSM/OpenInfraMap as a geometry cross-check only.

## Key Findings

| Category | Best open alternative | Wins on | Loses on | Verdict for BC workflow |
|---|---|---|---|---|
| VRE capacity factors | RDRS/CaSR (~10 km) + ERA5 (0.25°) via atlite/RESource | Resolution, recency, license (redistributable) | Requires processing; not model-ready | **Better alternative exists** — build pipeline; validate vs CODERS |
| Hydropower | HYDAT streamflow + StatCan + atlite runoff | Spatial/temporal granularity, daily data | No plant-level generation; needs modeling | **Better raw data**; CODERS' flat monthly CF is weak but convenient |
| Hourly demand | BC Hydro BA load, AESO, IESO, EIA-930 | Recency (near-real-time), continuity | Not harmonized; licensing varies | **Better alternative exists** for recency; CODERS good for harmonized history |
| Technology costs | NREL ATB, CER EF2023, Danish Catalogue | Update cadence, documentation | Not Canada-specific (ATB/DEA) | **Mixed** — combine; CODERS unique for hydro renewal costs |
| Transmission topology | OSM/OpenInfraMap (ODbL) | Geometry, openness | No electrical parameters/ratings | **CODERS wins** — retain it |
| Generator inventory | Global Energy Monitor, WRI GPPD, StatCan/CER | GEM recency; StatCan authority | GPPD stale (2022); coverage gaps | **Mixed** — CODERS competitive; supplement with GEM |

## Details

### Category 1 — VRE resource / capacity factor data

**CODERS baseline:** VRE capacity factors are computed with the CoRES tool from ERA5, but delivered on the **legacy MERRA-2 grid (~3,000 cells at 0.5° lat × 0.625° lon)**, hourly 8760-h series. This is a coarse, frozen product — the spatial grid is a downgrade from the native ERA5 resolution the data was derived from.

**The alternatives:**
- **ERA5 (native 0.25°, ~31 km, hourly, 1940–present, ~5-day latency):** Global, continuously updated, distributed via the Copernicus Climate Data Store under a license permitting redistribution with attribution. This is the de-facto standard for energy-system modeling and is the default in atlite/PyPSA-Earth.
- **ECCC RDRS / Canadian Surface Reanalysis (CaSR), ~10 km, hourly:** The highest-resolution reanalysis focused on Canada. v2.1 covers 1980–2018; **final v3.1 has been publicly available since April 2025** at the same ~10 km resolution, and the successor v3.2 extends coverage to 1980–2024. Distributed via CaSPAr (free, requires a Globus account) and ECCC's collaboration server in NetCDF. **License: the ECCC Data Servers End-use Licence (v2.1, Sept 2022), functionally equivalent to the Open Government Licence – Canada**, which grants a "worldwide, royalty-free, perpetual, non-exclusive licence to use the Information, including for commercial purposes" with the right to "copy, modify, publish, translate, adapt, distribute or otherwise use," conditioned only on attribution. This is materially more permissive for redistributing derived model inputs than CODERS' partial EULA.
- **Global Wind Atlas 3 (250 m, DTU/World Bank, CC BY 4.0):** Microscale wind climatology downscaled from ERA5 via WAsP; validated against tall-mast measurements with per-country validation reports. Best for *long-run-average* siting/bias-correction, **not** hourly time series (it is a climatology, not a timeseries).
- **Global Solar Atlas 2.0 (Solargis/World Bank, CC BY 4.0, ~250 m / nominally ~1 km, covers 60°N–55°S):** Validated against 228 ground stations. Note the 60°N latitude cap excludes far-northern Canada; still fine for populated BC.
- **NREL WIND Toolkit / NSRDB:** NSRDB's "Americas" domain extends to **60°N** (bordered by latitudes to 60°N), covering "much of Canada"; the WIND Toolkit / Wind Resource Database now explicitly includes Canada. Latitude limits and US-centric design make these secondary for BC (much of BC sits near/above the useful edge).
- **Renewables.ninja (MERRA-2 0.5°×0.625° or SARAH; CC BY-NC 4.0):** Same coarse grid as CODERS' delivery grid, and its **non-commercial license (CC BY-NC 4.0)** is more restrictive than ERA5/RDRS for redistribution of derived inputs.
- **Pipelines:** atlite and PyPSA-Earth build ERA5 (or SARAH-2) "cutouts" (recommended 30 km × 30 km, hourly) and apply conversion functions to produce technology-specific CF time series with land-use exclusions from Copernicus/protectedplanet — a fully reproducible, scriptable chain. RESource fits this same pattern.

**Validation evidence:** Multiple peer-reviewed studies establish that **ERA5 outperforms MERRA-2** for wind-power simulation. Olauson (2018, *Renewable Energy*, "ERA5: The new champion of wind power modelling?") found that "ERA5 performs better than MERRA-2 in all analysed aspects; correlations are higher, mean absolute and root mean square errors are in average around 20% lower," validated on aggregated output in five countries and 1,051 individual Swedish turbines. Gruber et al. (2022, *Energy*) concluded: "(i) ERA5 outperforms MERRA-2 in terms of the assessed error measures. (ii) Bias-correction with GWA2 does not improve simulation quality substantially, while bias-correction with GWA3 is detrimental." For Canada specifically, ClimateData.ca reports the CaSR developers note an RMSE of ~2 m/s between reanalysis and observed wind, and several studies find ERA5 underestimates extreme wind speeds in Canada. A BC-specific inter-comparison (*Atmosphere-Ocean*, 2024) evaluated RDRS and ERA5-Land against stations in the Skeena and Nechako watersheds.

**Verdict — BETTER ALTERNATIVE EXISTS.** For BC, RDRS/CaSR (~10 km, hourly, OGL-equivalent, updatable) is the superior validated basis, with ERA5 as the global fallback and GWA3/GSA for bias-correction/siting. The right comparison is **RESource-on-ERA5/RDRS pipeline vs CODERS' frozen MERRA-2-grid product**: the pipeline wins on resolution, recency, license, and reproducibility; CODERS wins only on being immediately model-ready.

### Category 2 — Hydropower

**CODERS baseline:** Monthly-resolution hydro CFs held constant per hour, derived from StatCan monthly hydraulic generation ÷ fleet installed capacity, for **2013 and 2021 only**. This is a crude approximation that erases sub-monthly variability and reservoir operations.

**The alternatives:**
- **Water Survey of Canada / HYDAT:** Per Environment and Climate Change Canada, the network comprises "over 2700 active and 5080 discontinued hydrometric monitoring stations across Canada," with real-time data from over 1,900 stations. HYDAT provides **daily** mean flow (some sub-daily/real-time), is updated quarterly, and is accessible via WaterOffice, the MSC GeoMet OGC API, and the `tidyhydat` R package. Governed by the Open Government Licence – Canada. This is a vastly richer temporal and spatial basis than CODERS' flat monthly CF — but it is *streamflow, not generation*, so it requires a runoff-to-generation model.
- **atlite/PyPSA-Earth hydro module:** converts ERA5 runoff into per-plant inflow, then rescales to match annual generation totals (e.g., from IRENA/StatCan) — a reproducible way to turn HYDAT/ERA5 into hourly hydro availability.
- **BC Hydro regulatory filings** (reservoir operations, water licences) provide plant/reservoir-level operational context but are not a clean open dataset.
- StatCan monthly generation tables remain the authoritative generation series (same source CODERS uses) but at monthly resolution.

**Verdict — BETTER RAW DATA EXISTS, but not model-ready.** HYDAT + atlite runoff rescaling offers far better spatial/temporal granularity than CODERS' 2013/2021 flat monthly CF. For a BC hydro-dominated system this matters enormously (BC's system is storage-hydro-driven). Recommendation: build a HYDAT/ERA5-runoff hydro pipeline for BC's major basins; retain CODERS/StatCan annual-monthly totals for calibration. This is the category where CODERS' approach is weakest.

### Category 3 — Hourly electricity demand

**CODERS baseline:** harmonized provincial hourly demand 2018–2022 (2023 for BC). Convenient and multi-province, but ends in 2022/23 and is provincial-only (no nodal/zonal resolution).

**The alternatives:**
- **BC Hydro Balancing Authority load data:** gross telemetered load in **hourly and daily** formats, updated daily (previous-day), with historical BC load and hourly path TTC. This is the authoritative BC series and is continuously current — beating CODERS on recency.
- **AESO** (Alberta): hourly Alberta Internal Load and metered volumes, published with historical archives; some detailed data via data request.
- **IESO** (Ontario): hourly Ontario/market demand back to 1994, plus FSA-level (forward sortation area) consumption — offering *sub-provincial* granularity CODERS lacks.
- **Hydro-Québec, SaskPower, Manitoba Hydro, Maritime utilities:** publish demand with varying granularity/recency.
- **EIA-930 (via PUDL/Catalyst Cooperative):** hourly demand, generation-by-source, and **interchange including international exchanges with Canada** (BC Hydro, AESO, etc. appear as adjacent BAs). Harmonized, machine-readable, continuously updated — the closest thing to an "Open Power System Data for North America." PUDL packages it under open terms.

**Licensing:** Utility data licensing varies (BC Hydro and AESO terms are custom/utility-specific; IESO data is broadly public; EIA-930 is US-Government public domain, and PUDL redistributes under open license).

**Verdict — BETTER ALTERNATIVE EXISTS for recency and granularity.** For BC, pull BC Hydro BA load directly (current, authoritative) and cross-check imports/exports via EIA-930. CODERS remains valuable as a *pre-harmonized, multi-province historical* set for 2018–2022/23, but for anything beyond 2023 or sub-provincial detail the primary utility feeds win.

### Category 4 — Technology cost and performance

**CODERS baseline:** a `generation_generic` cost table plus **unique hydro renewal/greenfield project cost compilations drawn from regulatory filings** — the latter genuinely hard to find elsewhere.

**The alternatives:**
- **NREL Annual Technology Baseline (ATB) 2024:** the gold standard for transparency. The 2024 Electricity ATB is NREL's **10th annual update**, providing CAPEX, O&M, capacity factors and LCOE "both at present and with projections through 2050," including a "Market + Policies case [that] includes some effects of the Inflation Reduction Act" and an "R&D Only case." It is delivered as an Excel workbook plus database-friendly CSVs on the Open Energy Data Initiative data lake, updated annually. Free/public but US-centric (US$ cost basis, US financing), and it carries a disclaimer/technical-limitations agreement rather than a clean CC license.
- **CER Canada's Energy Future 2023 (EF2023):** Canada-specific technology assumptions (Appendix 2), all-scenario data appendices downloadable from Open Government machine-readable files; costs in $2022 CAD. The most authoritative *Canadian* cost/assumption reference.
- **Danish Energy Agency Technology Catalogues:** highly detailed, standardized, regularly updated techno-economic data sheets across nine sectors; free public resource, widely used in open modeling.
- **EIA capital cost studies and Lazard LCOE:** useful cross-checks; Lazard is a report (annual) rather than a structured dataset.

**Verdict — MIXED / COMBINE.** No single open source dominates. For a BC OSeMOSYS/PyPSA model: use **CER EF2023 as the Canadian anchor**, **NREL ATB and DEA catalogues for technology trajectories and parameters not in EF2023**, and **retain CODERS specifically for its hydro renewal/greenfield cost compilations**, which are its genuine unique contribution here. Flag currency/financing-basis conversions explicitly.

### Category 5 — Transmission network topology and electrical parameters

**CODERS baseline:** a node-resolved Canadian network with **voltage, reactance, and summer/winter thermal ratings**, plus interties/interface capacities. This is rare and high-value.

**The alternatives:**
- **OpenStreetMap / OpenInfraMap (ODbL):** rich geometry — voltage tags, substations, HVDC converters — with good coverage of high-voltage lines in developed countries. But OSM has four documented fundamental limitations for power modeling: **(1) no electrical parameters (impedance, ratings), (2) missing parallel circuits, (3) incomplete voltage tags (~15–30% of lines untagged in the US), and (4) no demand or cost data.** Extraction via earth-osm/PyPSA-Earth or Infrageomatics exports.
- **Academic pipelines** (e.g., the OSM-to-OPF pipeline literature) synthesize missing impedances/ratings from voltage and geometry — a method, not a dataset.
- **Provincial utility GIS portals** (BC Hydro, etc.) publish some infrastructure GIS but rarely with electrical parameters openly.
- **COPPER/SILVER (McPherson/SESIT) publications:** these are zonal/interface transmission representations (pipe-and-bubble), consistent with CODERS' intertie data, not a competing nodal-parameter dataset.

**Verdict — CODERS WINS.** Nothing open matches CODERS' node-resolved network *with electrical parameters and seasonal ratings*. OSM/OpenInfraMap can cross-check geometry and topology (ODbL, freely redistributable), but cannot supply reactance or thermal ratings without synthesis. **Retain CODERS** for network parameters; use OSM for geometry validation and gap-filling. This is the strongest argument for keeping CODERS in the stack.

### Category 6 — Generator / asset inventories

**CODERS baseline:** a generators table (excludes <1 MW facilities), regularly updated by the Energy Modelling Hub.

**The alternatives:**
- **WRI Global Power Plant Database (CC BY 4.0):** v1.3.0 contains "nearly 35000 power plants in 167 countries, representing about 72% of the world's capacity" (grid-scale ≥1 MW). However, the WRI GitHub repo states verbatim: "This project is not currently maintained by WRI. There are no planned updates as of this time (early 2022). The last version of this database is version 1.3.0." — so it is stale for a fast-moving fleet.
- **Global Energy Monitor trackers (wind/solar/hydro/gas):** actively maintained, granular, project-level with status and dates; widely used and more current than GPPD, though licensing is custom (research/attribution terms).
- **StatCan and CER:** authoritative Canadian generation/capacity statistics; CER publishes generator-level data. Provincial regulators (e.g., BC) publish facility lists.

**Verdict — MIXED / COMPETITIVE.** CODERS is competitive and Canada-tailored (with sensible <1 MW cutoff). Supplement with **Global Energy Monitor** for recent additions/retirements and **StatCan/CER** for authoritative capacity. GPPD is too stale to rely on. No clear single winner over CODERS, but GEM improves recency.

## Recommendations

**Stage 1 — Adopt CODERS as the integrated baseline and validation reference (now).** It is the only pre-harmonized, model-ready Canadian package and is unmatched for transmission electrical parameters and hydro-cost compilations. Use it to bootstrap BCNexus and to validate any pipeline outputs.

**Stage 2 — Replace the two weakest CODERS categories with reproducible pipelines (immediate priority).**
1. **VRE CFs:** Build RESource on **ERA5 (0.25°) and RDRS/CaSR (~10 km)** via atlite. Bias-correct with GWA3/GSA for siting. This resolves the resolution, recency, license, and reproducibility deficits of CODERS' frozen MERRA-2-grid product simultaneously. *Benchmark that would change this:* if a validated Canadian bias-corrected CF product is published openly, prefer it.
2. **Hydro:** Build a HYDAT (daily streamflow) + ERA5-runoff pipeline rescaled to StatCan annual generation for BC's major basins, replacing CODERS' 2013/2021 flat monthly CF.

**Stage 3 — Swap demand to live primary feeds.** Use **BC Hydro Balancing Authority load** as the BC series and **EIA-930/PUDL** for interties; keep CODERS' harmonized tables only for the 2018–2022/23 historical multi-province window.

**Stage 4 — Layer costs.** CER EF2023 (Canadian anchor, $2022 CAD) + NREL ATB + DEA catalogues; retain CODERS hydro renewal/greenfield costs. Document all currency/financing conversions.

**Stage 5 — Keep CODERS for transmission and inventory; cross-check.** Use OSM/OpenInfraMap (ODbL) to validate network geometry; supplement generator lists with Global Energy Monitor.

**Licensing decision rule:** For any *redistributed derived model input*, prefer sources under CC BY 4.0 (GWA, GSA, WRI GPPD), OGL-Canada/ECCC End-use Licence (RDRS, HYDAT, StatCan/CER), or public-domain (EIA-930). Treat CODERS' EULA-encumbered fields as "compute-with, don't-redistribute" unless SESIT/EMH grants explicit terms, and isolate them so your published pipeline outputs remain redistributable.

## Caveats
- **CODERS specifics** (MERRA-2 grid cell count, hydro CF years, demand vintages, table contents) are taken from the user's October 2025 API tutorial description; I could not independently re-verify every field against the live CODERS API, which sits behind a SESIT/EMH access portal.
- **CaSR v3.1 exact end-year** is not explicitly stated on ECCC's primary overview page (only "available since April 2025," same ~10 km resolution); v3.2 covers 1980–2024. Confirm the v3.1 window on the ECCC dataset-specifics page before committing.
- **RDRS wind heights:** 10 m and 40 m winds and downward shortwave flux are confirmed variables; 80 m hub-height wind availability was not confirmed and should be checked directly — hub-height extrapolation may be required.
- **Reanalysis validation for Canada** draws substantially on European/US validation studies (Olauson 2018; Gruber et al. 2022) plus a smaller set of Canada-specific comparisons; a Canada-wide validation of reanalysis CFs against observed provincial wind/solar generation remains thin, which is itself an argument for local bias-correction.
- Utility and GEM **licensing terms are custom** and can restrict redistribution; verify before publishing derived products.
- Some retrieved figures (e.g., OSM voltage-tag completeness) are US-based estimates; Canadian completeness may differ and should be spot-checked for BC.