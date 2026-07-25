# CODERS Schema & Column-Change Runbook

**Author:** Md Eliasinul Islam
**Schema baseline:** CODERS Data Inventory Update — 1 October 2025
**Purpose:** one place that absorbs CODERS column renames, so a CODERS data
change never means editing multiple scripts again.

---

## The problem

CODERS changes column names between releases (e.g. `summer_rating_mva` →
`ttc_summer`, `Reactance` → `Transmission Line Reactance`). Historically each
consuming script hard-coded raw column names, so every release broke several
scripts and required scattered edits.

## The fix — one control point

All CODERS column handling lives in **`config/coders.yaml`** and is applied in
**one function**, `CODERSData.normalize_columns`. Every consumer reads CODERS
through the shared client (`get_coders()`), never `pd.read_csv` on the raw file.

```
CODERS API ──▶ raw CSV cache ──▶ get_coders().load_table() ──▶ normalize_columns()
                                                                     │
                                        coders.yaml `version:` map ──┘
                                        (adds canonical alias columns)
                                                                     ▼
                                        consumer sees raw names AND canonical names
```

### Additive aliasing (why nothing breaks)

`normalize_columns` is **additive**: for each `canonical: source` mapping it
*adds* a canonical column as a copy of the source and **keeps the raw column**.
So a table ends up with a superset of names — code using either the raw name or
the canonical name works. Renaming-away is never done, so a mapping can't break
a consumer that still uses the old name.

---

## Runbook: CODERS renamed/added a column

When a release breaks something (a `KeyError` on a column, or a "mapped source
column absent" warning in `logs/pypsa_bc.log`):

1. Identify the canonical name your code uses (the one that errored) and the new
   raw name in the pull (open the cached CSV, or `get_coders().show_tables()` /
   inspect columns).
2. Add/edit one line under the table's block in `config/coders.yaml`:

   ```yaml
   coders:
     version:
       2025OCT28:
         agreement: old_to_new          # keys = canonical (your code), values = current CODERS column
         transmission_lines:
           summer_rating_mw:  ttc_summer
         generators:
           connecting_node_code: network_node_code   # <-- add the new source name here
   ```
3. Re-run. No script edits. If you want a hard guarantee, add the canonical name
   to `required:` (below) so a future disappearance fails loudly instead of
   silently.

That's the whole loop. One file, one line.

---

## Current mappings (from `coders.yaml`)

`agreement: old_to_new` → **keys are the canonical names your code uses**,
**values are the current CODERS source columns**.

**transmission_lines**

| Canonical (code) | CODERS source (Oct 2025) |
|---|---|
| `starting_node_code` | `network_node_code_starting` |
| `ending_node_code` | `network_node_code_ending` |
| `summer_rating_mw` | `ttc_summer` |
| `winter_rating_mw` | `ttc_winter` |
| `reactance_ohm` | `Transmission_Line_Reactance` |
| `segment_reactance` | `Transmission_Line_Segment_Reactance` |
| `length_km` | `line_length_km` |
| `voltage_kV` | `voltage` |

**generators**

| Canonical (code) | CODERS source (Oct 2025) |
|---|---|
| `connecting_node_code` | `network_node_code` |
| `generator_capacity` | `facility_installed_capacity` |
| `generator_fuel_type` | `gen_type` |

Raw generators columns (Oct 2025 pull), all still available post-aliasing:
`generation_unit_name, generation_unit_code, generation_facility_name,
generation_facility_code, owner, province, location, latitude, longitude,
elevation, copper_balancing_area, operating_region, network_node_name,
network_node_code, network_node_voltage, network_connection, network_contribution,
start_year, previous_renewal_year, possible_renewal_year, closure_year, gen_type,
gen_type_copper, unit_installed_capacity, facility_installed_capacity,
capacity_adjustment, unit_effective_capacity, facility_effective_capacity,
capacity_factor, unit_average_annual_energy, facility_average_annual_energy,
generation_units_represented, total_facility_generation_units, notes`.

---

## Guardrails

- **`required:` (hard fail).** Canonical columns listed under `coders.required`
  must exist after aliasing, else `load_table` raises with a clear message. Use
  it for columns a stage cannot run without (e.g. `transmission_lines:
  [starting_node_code, ending_node_code, summer_rating_mw, reactance_ohm,
  voltage_kV]`).
- **Drift warning (soft).** If a mapped source column is absent from a pull,
  `normalize_columns` logs a WARNING naming the table, the missing source, and
  the `version` key — your early signal that CODERS moved something.
- **`version:` blocks are dated.** Keep the newest block last (or set
  `active_version`); older blocks stay as a change history.

---

## The one rule for consumers

**Read CODERS through `get_coders()`, never `pd.read_csv` on the cached file.**
Only the client applies aliasing; a direct CSV read gets raw columns and
re-introduces the whack-a-mole. Pattern:

```python
from pypsa_bc.data.coders import get_coders
gens = get_coders().load_table("generators", as_gdf=False)          # aliased
gen_generic = get_coders().load_table("generation_generic", as_gdf=False)
```

---

## Provenance

Schema baseline documented from *CODERS Data Inventory Update – Release Notes –
October 1, 2025* (key changes: `Reactance` split into line/segment reactance;
interties returned to provincial inventories; costs and generation updated to
year-end 2024; no new generation attributes). Data pulls are logged with dates
in `reports/DATA_SOURCES.md`.
