# Hydro Data Pipeline — building the cascade data from scratch

**Author:** Md Eliasinul Islam
**Last updated:** 2026-07-23

This is the single, canonical procedure for rebuilding the hydro cascade data
after CODERS withdrew its `hydro_cascade` endpoint. There is **one** active
pipeline; the earlier single-file variant is archived under
`_archive/hydro_pipelineA/`.

---

## The idea in one sentence

The cascade data splits into **static physics** (reservoir storage, discharge,
cascade order — properties of concrete structures, effectively constant) and
**volatile fleet attributes** (capacity, energy, commissioning years — which
change every CODERS release). We pin the static part as version-controlled
inventory and re-join the volatile part from CODERS on each release.

---

## Inputs you need

| Input | Where from | Used in |
|-------|-----------|---------|
| `hydro_cascade.csv` (legacy, withdrawn) | your last saved copy | Stage 1 seed |
| `generators` / `existing_hydro` (live) | CODERS API, `province=BC` | Stage 2 |
| GRanD export + `reservoir→GRAND_ID` crosswalk | GRanD/GDW + hand-made | Stage 1 (optional storage upgrade) |

---

## Data flow

```
STAGE 1  — build static inventory (ONCE, then curate)
────────────────────────────────────────────────────
  legacy hydro_cascade.csv
        │   bootstrap_hydro_static.py         (+ optional GRanD swap)
        ▼
  data/static/hydro_topology.csv      facility → reservoir → cascade_group, cascade_order
  data/static/hydro_reservoirs.csv    reservoir → storage, discharge, sill, spill, provenance


STAGE 2  — rebuild cascade table (EVERY CODERS release)
────────────────────────────────────────────────────
  data/static/hydro_topology.csv    ┐
  data/static/hydro_reservoirs.csv  ├─ build_hydro_cascade.py ──► rebuilt hydro_cascade table
  generators (live CODERS)          ┘


STAGE 3  — consume
────────────────────────────────────────────────────
  rebuilt hydro_cascade table  ──►  create_hydro_assets.main()
```

---

## Step by step

### Stage 1 — seed the static inventory (run once)

```bash
python codebase/pypsa_bc/data/bootstrap_hydro_static.py \
    --legacy path/to/hydro_cascade.csv \
    --out-dir data/static
```

Produces `data/static/hydro_topology.csv` and `data/static/hydro_reservoirs.csv`.
Every reservoir row is stamped `method="legacy_CODERS_hydro_cascade"`,
`verified=False` — a seed, not yet independently checked.

**Optional — upgrade storage from GRanD** (recommended, replaces the frozen
volumes with a primary source). Maintain a two-column crosswalk
`reservoir_grand_crosswalk.csv` (`reservoir`, `grand_id`), then:

```bash
python codebase/pypsa_bc/data/bootstrap_hydro_static.py \
    --legacy path/to/hydro_cascade.csv \
    --out-dir data/static \
    --grand path/to/grand_export.csv reservoir_grand_crosswalk.csv
```

Matched reservoirs flip to `method="GRanD"`, `verified=True`; unmatched keep the
legacy freeze. Repeat as you verify more reservoirs — the inventory improves
incrementally, one reservoir at a time.

### Stage 2 — rebuild the cascade table (every CODERS release)

```bash
python codebase/build_hydro_cascade.py \
    --generators path/to/generators.csv \
    --static-dir data/static \
    --out out/hydro_cascade_rebuilt.csv \
    --legacy path/to/hydro_cascade.csv        # optional: enables diff/derate audits
```

Joins live fleet attributes onto the static inventory and prints its built-in
audits (completeness of constraint-binding fields, cascade order vs powerhouse
elevation, diff against the legacy table, capacity-derate reconciliation).

### Stage 3 — feed the model

Point `create_hydro_assets.main()` at the rebuilt table.

---

## Which script does what

| Script | Role | Run frequency |
|--------|------|---------------|
| `pypsa_bc/data/bootstrap_hydro_static.py` | Build the static inventory (topology + reservoirs) from the legacy file; optional GRanD storage swap | Once, then when a dam/reservoir changes |
| `build_hydro_cascade.py` | Join live CODERS generators onto the inventory → rebuilt cascade table (with audits) | Every CODERS release |
| `pypsa_bc/data/grand_reservoir_adapter.py` | Separate path: build blueprint **Table B** (`asset_id`, storage in MWh) for `hydro.load_reservoir_sites` | Only if using the Table-A/B/C blueprint pipeline |

---

## Maintaining the inventory

`hydro_reservoirs.csv` carries provenance columns so it can be trusted field by
field: `method`, `verified`, `source`, plus quality flags `flag_spill`
(placeholder/missing spill), `flag_minstor` (dead storage imputed as 0),
`flag_dup` (reservoir shared by >1 facility — correct, not an error). The
maintenance loop is: pick an unverified reservoir, confirm its storage against a
primary source (GRanD, BC Hydro, WUP), update the value, set `verified=True`.

---

## Caveat to confirm

`build_hydro_cascade.py` emits a **redesigned** schema (`facility`, `reservoir`,
`cascade_group`, `cascade_order`, `live_*_m3`, provenance columns), not the old
34-column CODERS layout. If `create_hydro_assets.main()` still expects the
legacy column names, add a thin final rename step, or update the consumer to the
new schema. Decide this before wiring Stage 3.
