# Hydro Cascade Modelling — Region-Agnostic Data Contract

**Author:** Md Eliasinul Islam
**Purpose:** Specify the minimum set of inputs a new hydro-dominant region must supply to reuse the BC-PyPSA cascade preprocessing pipeline without modifying source code.
**Scope:** Defines interfaces only. The processing logic (basin routing, cascade subtraction, flow-time delay, monthly calibration) is region-independent and is not restated here.

---

## 1. Design principle

The pipeline separates *physics* (transferable algorithms) from *jurisdiction* (region-specific data and regulation). A region is fully described by three input tables plus one configuration file and, where applicable, one small flow-policy adapter. If a region produces conforming inputs, the existing pipeline runs unchanged. The test of conformance is operational: a new region should require **no edits to `hydro.py`** — only new data files and one config.

The water-policy step (Water Use Plans in BC) is the single component that is not portable as data alone, because regulatory flow records differ in format across jurisdictions. It is isolated behind an adapter interface (Section 5).

---

## 2. Table A — Generation/turbine sites

One row per generating unit (turbines may later be aggregated to an asset by the pipeline).

| Field | Type | Units | Required | Notes |
|-------|------|-------|----------|-------|
| `component_id` | string | — | yes | Unique unit identifier |
| `asset_id` | string | — | yes | Generating station the unit belongs to; many-to-one |
| `lat` | float | deg | yes | Used to locate the containing basin |
| `lon` | float | deg | yes | — |
| `capacity` | float | MW | yes | Rated unit capacity; summed to asset |
| `annual_avg_energy` | float | GWh/yr | yes | Used to normalise run-of-river power series |
| `hydro_type` | enum | — | yes | `reservoir`, `reservoir-impute`, or `ror` |
| `cascade_order` | int | — | for cascades | Topological rank within a cascade (headwater = lowest) |
| `upper_reservoir_id` | string | — | for cascades | Immediately upstream reservoir `asset_id` |
| `lower_reservoir_id` | string | — | for cascades | Immediately downstream reservoir `asset_id`; empty if terminal |

The cascade topology is defined **entirely** by `upper_reservoir_id` / `lower_reservoir_id`. The pipeline reconstructs the release/spill chain from these two fields; no separate adjacency file is needed.

---

## 3. Table B — Reservoir sites

One row per reservoir (storage) asset.

| Field | Type | Units | Required | Notes |
|-------|------|-------|----------|-------|
| `asset_id` | string | — | yes | Matches `upper/lower_reservoir_id` in Table A |
| `lat`, `lon` | float | deg | yes | Reservoir outlet location |
| `min_storage`, `max_storage` | float | m³ or MWh | optional | Live storage bounds; require head assumption if volumetric |
| `min_level`, `max_level` | float | m | optional | Forebay elevation range |
| `cascade_order` | int | — | for cascades | As in Table A |

If storage bounds are omitted, the network defaults to a single `e_nom` per reservoir (the current BC build behaviour).

---

## 4. Table C — Regulatory flow statistics (per reservoir)

One file per reservoir with a regulatory flow target. Drives the simplified water-policy calibration.

| Field | Type | Units | Required | Notes |
|-------|------|-------|----------|-------|
| `Month` | enum | — | yes | Full month name (`January` … `December`) |
| `Mean Monthly Inflow` | float | m³/s (cms) | yes | Regulatory or historical monthly mean used to calibrate inflow |

Reservoirs of type `reservoir-impute` need not supply this file; the pipeline imputes a zero local-inflow series for downstream reservoirs lacking statistics (documented as a limitation, not a silent default).

---

## 5. Geospatial and reanalysis inputs

| Input | Source | Region action |
|-------|--------|---------------|
| Basin polygons + topology (`HYBAS_ID`, `NEXT_DOWN`, `DIST_MAIN`) | HydroBASINS / HydroSHEDS | Download the continental tile covering the region |
| Runoff + surface height | ERA5 via Atlite cutout | Define the bounding box from site coordinates |

Both are global products; only the spatial extent changes between regions. No schema change is required.

---

## 6. Region configuration (single YAML)

Every BC-specific constant currently in code becomes a declared value. Illustrative keys:

```yaml
region: BC
basins_file:   data/hydrobasins/hybas_na_lev08.shp
cutout:        data/cutouts/bc_era5.nc
flowspeed:     1.0          # m/s, basin-to-plant flow celerity
inflow_tables: data/wup/    # directory of Table C files, one per reservoir
calibration:   mean_inflow_calibrate   # or mean_inflow_constant
env_flow_pu:   0.10         # minimum environmental flow as fraction of p_nom
spill_cost:    1.0          # $/MWh penalty on spill links
reservoir_classes:          # replaces hardcoded RESERVOIR_DISTRICTS
  - Peace River
  - Columbia-Shuswap
```

A new region is then `region_<name>.yaml`. No `.py` file is edited.

---

## 7. Water-policy adapter (the one non-portable step)

Water Use Plans are a BC regulatory artifact; other jurisdictions hold equivalent records in different formats (e.g. Norway's manøvreringsreglement, Quebec licensing, US FERC minimum-flow orders). The *concept* — a regulatory seasonal flow constraint — is universal; the *source format* is not.

The pipeline therefore depends on an abstract interface, not on the WUP file format:

```python
def load_flow_policy(reservoir_id: str) -> pd.DataFrame:
    """Return a frame with columns ['Month', 'Mean Monthly Inflow'] (cms)."""
```

Each region implements one small loader that maps its native records onto Table C. The calibration routine (`mean_inflow_calibration`) consumes the output and is blind to provenance. This confines the irreducibly local logic to roughly thirty lines per region.

---

## 8. Conformance checklist

A region is blueprint-ready when:

1. Tables A–C validate against the schemas above (types, enums, key integrity: every `upper/lower_reservoir_id` resolves to a Table B `asset_id`).
2. A `region_<name>.yaml` supplies all former code constants.
3. A `load_flow_policy` adapter exists if regulatory flow targets are used.
4. The pipeline runs end-to-end producing `inflow_*.csv` **without edits to `hydro.py`**.

Meeting items 1–4 substantiates the narrow, defensible claim: a region-agnostic data contract and pipeline, demonstrated on BC, with an adapter for jurisdiction-specific water policy — rather than an unvalidated general-framework claim.
