# BC Hydro data integration for Chapter 2

**Evidence date:** 2026-07-31  
**Scope:** operating constraints, hydraulic routes, and publication-readiness  
**Historical decision at evidence-integration stage:** useful primary evidence
was integrated; Gate 1 then remained **FAIL**. A later predeclared release-rule
uncertainty envelope moved Gate 1 to **PASS with two warnings**; see
`release_uncertainty_methods.md` and `gate1/gate1_report.md`.
because seven station rules still lack a complete, current representation.

## Answer first

BC Hydro's Water Use Plan (WUP) library closes several important model gaps.
It provides accepted operating rules, route descriptions, reservoir limits, and
monitoring context that generic generator inventories do not contain. The model
now uses these documents to distinguish turbine discharge from ecological dam
release and to generate four source-labelled hourly route schedules.

The new evidence improves minimum-release coverage from 5/20 in the original
audit to 13/20. It does not justify running the publication scenario matrix yet.
La Joie and Ruskin need endogenous state/intertemporal formulations, while five
facilities need a current order, treaty schedule, or non-BC-Hydro owner source.

The BC Hydro website also provides operational datasets that are distinct from
the WUP evidence package. These resources should now be used as the primary
system-level calibration and validation layer. Their availability does not mean
that station-level generation, complete reservoir volumes, or all operative
constraints are public; each dataset is assigned a bounded evidentiary role.

## Public operational-data inventory

The machine-readable catalog is
`data/validation/bc_hydro/public_data_catalog.csv`. It separates observations,
provisional observations, forecasts, and documentary evidence so that a public
web value is not automatically treated as a calibrated model parameter.

| BC Hydro resource | Direct use in PyPSA-BC | Required scientific check |
|---|---|---|
| [Historical transmission data](https://www.bchydro.com/energy-in-bc/operations/transmission/transmission-system/balancing-authority-load-data/historical-transmission-data.html) | Hourly provincial load from 2001 onward; hourly directional TTC for BC-AB and BC-US | Load: annual TWh, peak GW/timestamp, monthly energy, duration curve and residuals. TTC: apply by direction and hour; never substitute it for internal branch ratings. |
| [Historical net actual flow](https://www.bchydro.com/energy-in-bc/operations/transmission/transmission-system/actual-flow-data/historical-data.html) | Hourly interchange validation from 2007 onward | Confirm sign convention, imports/exports, annual energy, directional peaks and coincidence with scarcity hours. |
| [Balancing-authority and ACE data](https://www.bchydro.com/energy-in-bc/operations/transmission/transmission-system/balancing-authority-load-data.html) | Current load checks and balancing-variability context | ACE is not a direct plant-dispatch or reserve-procurement series; use it only for a separately declared balancing sensitivity. |
| [Hydrometeorologic stations](https://www.bchydro.com/energy-in-bc/operations/transmission-reservoir-data/hydrometeorologic-data.html) | Crosswalk more than 200 BC Hydro station IDs to WSC/snow stations; diagnose recent inflow, level and weather conditions | BC Hydro explicitly labels the posted four-day automated feed as not quality-controlled. Historical, quality-assured WSC series should be the validation backbone. |
| [Reservoir levels](https://www.bchydro.com/energy-in-bc/operations/transmission-reservoir-data/previous-reservoir-elevations.html) and [discharges](https://www.bchydro.com/energy-in-bc/operations/transmission-reservoir-data/reservoir-discharges.html) | Validate storage-state timing, operating ranges and selected release events | Convert elevation to volume only with source-backed curves; distinguish actual observations from forecasts and short live windows. |
| [IPP supply lists](https://www.bchydro.com/work-with-us/selling-clean-energy/meeting-energy-needs.html) | Reconcile CODERS assets against the dated in-operation list and resource class | Contracted energy is not hourly output or dependable capacity; record unmatched and ambiguous facilities. |
| [Transmission maps](https://www.bchydro.com/energy-in-bc/operations/transmission/transmission-system/maps.html) and [engineering data](https://www.bchydro.com/energy-in-bc/operations/transmission/transmission-system/engineering-studies-data.html) | Validate bulk corridors, voltage classes, topology, selected ratings and interface assumptions | Public maps are schematics rather than complete GIS; detailed GIS access is directed through ICIS. Preserve uncertainty where ratings or geometry are unavailable. |

The 2021 load workbook, four directional TTC workbooks, and net actual-flow
workbook are pinned under `data/validation/bc_hydro`. Each downloaded file
has 8,760 hourly records after its workbook header is interpreted, uses the
expected binary XLS format, and has a URL, retrieval timestamp, byte count, and
SHA-256 digest in `data/validation/manifest.json`. These files permit an observed
external-boundary benchmark for the same year as the model. Reservoir
validation should next use a station crosswalk to longer historical Water
Survey of Canada records rather than scraping only the short live BC Hydro
windows.

The reproducible descriptive summary in
`data/validation/bc_hydro/summary_2021.json` reports 65.054 TWh gross
balancing-authority load, a peak of 11,792 MW at hour ending 18 on 27 December,
and no missing hourly load records. The publisher's sign convention is now
implemented explicitly: positive actual interchange is BC export and negative
is BC import. The directional TTC series also contain
8,760 valid observations each. Their annual ranges are 0-850 MW (BC to Alberta),
0-1,000 MW (Alberta to BC), 1,100-3,150 MW (BC to U.S.), and 400-2,500 MW (U.S.
to BC). These are observed path limits, not constant interconnector capacities;
hourly TTC is deferred to a secondary boundary sensitivity. The primary
production path uses observed interchange and does not enforce TTC. TTC remains
available for later boundary robustness work and must not be interpreted as an
internal branch rating.

## Sources added to the reproducible package

| Source | Model use | Scientific interpretation |
|---|---|---|
| [Walter Hardman WUP](https://www.bchydro.com/content/dam/hydro/medialib/internet/documents/environment/pdf/wup_walter_hardman_wup.pdf) | Cranberry Creek route and 0.1 m3/s rule | Route-specific; turbine discharge to Arrow Lakes cannot satisfy it |
| [Whatshan WUP](https://www.bchydro.com/content/dam/hydro/medialib/internet/documents/environment/pdf/wup_whatshan_water_use_plan_pdf.pdf) | Whatshan River route; no positive dam minimum | Resolves the minimum-release parameter, but reservoir-elevation limits remain separate |
| [Stave River WUP](https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/corporate/environment-sustainability/water-use-planning/lower-mainland/wup-stave-river-water-use-plan.pdf) | Ruskin tailwater and block-loading rules; Stave Falls classification | Ruskin requires an endogenous intertemporal formulation, not a constant scalar |
| [Bridge River WUP](https://www.bchydro.com/content/dam/hydro/medialib/internet/documents/planning_regulatory/wup/lower_mainland/2011q2/bridge_river_wup_rev.pdf) | Seton River, Lower Bridge River, Cayoosh Creek, and Fraser River route separation | The time-limited Terzaghi trial cannot be presented as a settled current rule |
| [2026 Lower Bridge River monitoring report](https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/corporate/environment-sustainability/water-use-planning/lower-mainland/brgmon-03-yr12-2026-02-26.pdf) | Current adaptive-management context | Confirms that monitoring continues to inform the long-term strategy |
| [Campbell River WUP](https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/corporate/environment-sustainability/water-use-planning/vancouver-island/campbell-river/campbell-river-water-use-plan.pdf) | Elk Falls Canyon route and seasonal minima | Canyon release is distinct from combined Lower Campbell River flow |
| [Alouette WUP](https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/corporate/environment-sustainability/water-use-planning/lower-mainland/wup-alouette-water-use-plan.pdf) | Alouette River outlet route and seasonal floor | Generation is diverted to Stave Lake and cannot satisfy the river rule |
| [Seven Mile WUP](https://www.bchydro.com/content/dam/hydro/medialib/internet/documents/environment/pdf/wup_seven_mile_project_water_use_plan_december_2006.pdf) | Seven Mile minimum-release classification | Zero-discharge blocks are permitted; coordination with Boundary and Waneta remains a separate rule |
| [Columbia River WUP](https://www.bchydro.com/content/dam/hydro/medialib/internet/documents/environment/pdf/wup_columbia_water_use_plan_revised_for_acceptance_by_th.pdf) | Mica–Revelstoke topology and treaty-operation context | Mica is schedule-governed; no unsupported zero was inserted |
| [Site C operating range](https://www.sitecproject.com/downstream-effects) and [federal Decision Statement](https://iaac-aeic.gc.ca/050/documents/p63919/164236E.pdf) | 390 m3/s combined minimum release | The regulatory decision supplies the binding value; BC Hydro corroborates the normal range |

All 15 registered documents/pages are stored under
`data/validation/chapter2/sources`. The generated manifest records the URL,
retrieval time, byte count, and SHA-256 digest for each file.

## Model changes supported by these sources

### Hydraulic routes

| Station | Turbine destination | Non-power release destination |
|---|---|---|
| Seton | Fraser River | Seton River |
| Alouette | Stave Lake | Alouette River |
| John Hart | Lower Campbell River | Elk Falls Canyon |
| Walter Hardman | Arrow Lakes Reservoir | Cranberry Creek |
| Whatshan | Arrow Lakes Reservoir | Whatshan River |

The route table also retains the earlier La Joie, Bridge 1+2, and Walden route
corrections. The builder creates explicit terminal water sinks when a river
leaves the modelled cascade.

### Implemented release rules

The production schedule contains 35,040 rows: four routes times 8,760 hourly
snapshots for 2021.

| Route | Implemented reference rule | Range |
|---|---|---:|
| Seton Dam to Seton River | continuous route minimum | 5.0 m3/s |
| Alouette Dam to Alouette River | 1.52 m3/s floor; 3.0 m3/s from Apr 15 to Jun 14 | 1.52–3.0 m3/s |
| John Hart spillway to Elk Falls Canyon | 7.0 m3/s from Apr 1 to Apr 15; 4.0 m3/s otherwise | 4.0–7.0 m3/s |
| Walter Hardman diversion dam to Cranberry Creek | historical-period floor | 0.1 m3/s |

The Walter Hardman WUP requires all available natural inflow when inflow falls
below 0.1 m3/s. The WUP's cited historical record has winter lows of about
0.5–2.0 m3/s, so the fixed 0.1 m3/s floor is valid only for hydrologies that do
not cross that threshold. Future stress hydrologies must use an inflow-
conditional formulation.

Optional pulse events at Alouette and John Hart are not silently assigned to
arbitrary dates. They remain declared event sensitivities. The combined Lower
Campbell River target also remains separate because turbine and canyon flows
can jointly satisfy it.

### Constant or absent rules

- Site C now has a 390 m3/s combined minimum-release constraint.
- Whatshan now has an evidence-backed zero for the minimum-release parameter.
- Stave Falls now has an evidence-backed zero for a separate station minimum;
  downstream ecological requirements remain assigned to Ruskin.
- Seven Mile has an evidence-backed zero for the minimum-release parameter,
  while its Boundary–Seven Mile–Waneta coordination rule remains open.

## Remaining publication-critical gaps

| Station | Why it remains open | Required closure |
|---|---|---|
| La Joie | Required release varies with Downton Lake elevation | Obtain/approve a stage–storage mapping and implement the piecewise state rule or predeclare a tested approximation |
| Ruskin | Tailwater, tide, prior block level, outage state, and peaking duration affect the rule | Implement a mixed-integer or explicitly conservative block-loading formulation and sensitivity |
| Bridge 1+2 / Terzaghi | The 2011–2015 trial schedule is not a settled current rule; later variances altered high-flow operation | Obtain the operative order/variance or model a clearly labelled trial/variance sensitivity |
| Mica | Treaty and non-treaty schedules govern releases | Obtain the study-year Detailed Operating Plan/request series or define a treaty-schedule scenario |
| Arrow Lakes Generating Station | The station is not a BC Hydro-owned generating asset even though the adjacent dam is a Treaty project | Obtain Columbia Power/operating-party evidence and clarify whether the model conflates the generator and Hugh Keenleyside dam |
| Waneta | BC Hydro documents the forebay operating range; DFO classifies it as run-of-river but reports that basic reservoir dimensions were not recorded | Resolve active volume or approve a predeclared storage sensitivity; release-rule evidence remains separate |
| Walden North | The accepted Bridge WUP explicitly states that the project has no storage | Zero-controllable-storage representation resolved; owner-specific release rule remains in the sensitivity envelope |

These seven cases are not equivalent. La Joie and Ruskin have primary evidence
but need new optimization formulations; the other five need current or owner-
specific evidence before the reference case can be frozen.

## Reproduction

Run from the project root:

```powershell
$env:PYTHONPATH='src;.'
python -m workflow.scripts.studies.chapter2.fetch_validation_sources
python -m workflow.scripts.fetch_validation_data --year 2021 --dataset all
python -m workflow.scripts.summarize_bc_hydro_validation --year 2021
python -m workflow.scripts.studies.chapter2.report_validation_sources
python -m workflow.scripts.studies.chapter2.resolve_cascade_evidence
python -m workflow.scripts.studies.chapter2.build_release_schedules
python -m workflow.scripts.create_hydro_assets
python -m workflow.scripts.enrich_format_hydro
python -m workflow.scripts.studies.chapter2.audit_inputs --protocol config/studies/chapter2_cascade.yaml
python -m unittest discover -s tests/studies/chapter2 -v
```

At this historical evidence-integration checkpoint, Gate 1 failed with one
error and one warning; minimum-release coverage was 13/20. The later
release-uncertainty integration supersedes those status counts: Gate 1 now
passes with two claim-boundary warnings, Gate 2C passes 42/42 release checks,
Gate 2D passes 6/6 storage checks, and the focused suite contains 45 passing
tests. Gate 2 remains 7/7, and the matched B/C
representation audit passes 12/12 basin-treatment rows.

The validation archive now contains 17/17 registered sources. The reproducible
fetch added BC Hydro's Waneta transaction application and the DFO hydro-impact
report; all 17 manifest rows contain a SHA-256 digest.

## Manuscript claim boundary

The manuscript may state that the model uses accepted, source-labelled WUP and
regulatory rules for the implemented facilities. It may not state that all BC
hydropower operations are calibrated or fully reproduced. Travel time remains a
structural sensitivity, five operating rules still need current/owner evidence,
and two verified rules still need endogenous formulations. Gate 1 now passes
with two explicit claim-boundary warnings. Core comparative claims remain
embargoed until the 24 runnable cases are solved, independently validated, and
shown to be robust across the declared release and Waneta-storage uncertainty
bounds.
