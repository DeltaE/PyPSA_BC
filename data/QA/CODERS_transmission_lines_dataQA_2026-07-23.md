# CODERS `transmission_lines` — Naming Convention & Data-Gap Report

**Author:** Md Eliasinul Islam
**Check date:** 2026-07-23
**Dataset version:** CODERS 2025OCT28 pull (`province=BC`)
**File audited:** `transmission_lines.csv` — 1,239 rows × 21 columns
**API:** `http://api.sesit.ca/transmission_lines?province=BC`

---

## 1. Summary

The file is internally consistent — no missing values outside `notes` (1,189/1,239 blank, expected), unit columns cross-validate (km/mi ratio = 1.6195 ≈ 1.60934), and segment lengths never exceed full-line lengths. The dataset is usable for network build, but three structural gaps matter for PyPSA-BC and six rows need cleaning before `consistency_check()`. The single largest dependency is **absent geometry**: lines are defined only by endpoint node codes, so coordinates and regional-district bus mapping must come from a separate Nodes/substations table.

---

## 2. Node-code naming convention

Every node code (starting and ending) follows a fixed three-token grammar:

```
BC _ <XXX> _ <TYPE>
│      │       │
│      │       └── facility-type suffix (see table)
│      └────────── 3-letter location abbreviation (e.g. RUC = Ruckles, GFT = Grand Forks)
└───────────────── province prefix (always BC in this pull)
```

All 2,478 node references parse to exactly three tokens with prefix `BC`. The suffix encodes facility type, decoded from the paired `network_node_name_*` fields:

| Suffix | Count | Decoded meaning | Evidence (node names) |
|--------|-------|-----------------|-----------------------|
| `DSS` | 987 | Distribution substation (step-down / load) | Kettle Valley, Passmore, Goward |
| `JCT` | 922 | Junction / line tap (no station) | Ruckles JCT, Cascade JCT |
| `GSS` | 213 | Generating station switchyard | Kootenay Canal **GS**, South Slocan **GS** |
| `TSS` | 175 | Terminal substation (major transmission) | Grand Forks **TS**, Waneta **TS** |
| `ISS` | 154 | Industrial substation / interconnection | Houston Forest Products, BC Gas, mines |
| `SWS` | 19 | Switching station | Emerald **SS**, Brilliant **SS** |
| `INT` | 4 | International intertie (US / BPA) | INT Nelway **BPA**, INT Ingledow **BPA** |
| `IPT` | 4 | Interprovincial tie (Alberta) | IPT Natal Pocaterra, IPT Cranbrook Chapel Rock |

`INT`/`IPT` nodes and the 6 `operating_region = Alberta` rows are the cross-border ties; US connections appear at the tell-tale US nominal voltages (161, 287 kV).

### Circuit-ID convention

`transmission_circuit_id` (e.g. `60L037`, `1L365MWN`, `2L...`) follows BC Hydro line naming: the leading number is the voltage class (`60L` → 69 kV series, `1L` → 138 kV, `2L` → 230 kV, `5L` → 500 kV), optionally suffixed with a location tag (`MWN`, `GLS`). 835 distinct circuits span 1,239 rows, i.e. **many circuits are split into multiple segments** — one row per segment.

### Voltage levels present (kV)

| kV | 63 | 69 | 138 | 161 | 230 | 287 | 360 | 500 |
|----|----|----|-----|-----|-----|-----|-----|-----|
| Lines | 88 | 507 | 381 | 3 | 184 | 14 | 9 | 53 |

63/69/138/230/500 are standard BC classes; 161 and 287 are US-standard (border ties). **360 kV (9 lines) is non-standard for BC and should be verified** — likely a mislabel or a special series.

---

## 3. What is missing (modelling dependencies)

| # | Missing item | Impact on PyPSA-BC | Severity |
|---|--------------|--------------------|----------|
| 1 | **Coordinates / geometry** — only endpoint node codes given | Cannot place lines, cannot map ~900 distinct nodes to the 25 regional-district buses, cannot plot. Requires the separate **Nodes / substations** table as a hard prerequisite. | **Blocker** |
| 2 | **Resistance `r`** — only reactance `x` provided | Fine for DC/linear power flow and your transport model (uses `x` only), but no AC load flow and no ohmic loss modelling. | Medium |
| 3 | **In-service / retirement year / status** | No time-phasing of the grid by scenario year (2035/2040). Generators carry `start_year`; lines do not — all lines are implicitly always-on. | Medium |
| 4 | **True thermal rating (MVA)** — table gives `ttc_summer`/`ttc_winter` (MW) | `ttc_*` is *Total Transfer Capability*, a path limit, not a per-line MVA thermal rating. Usable as `s_nom` but semantically different (see §5). | Medium |
| 5 | **Impedance base / per-unit** — `x` in ohms, no base MVA | Need `voltage` + a base MVA to convert ohms→p.u., or feed ohms with correct `v_nom`. Manageable, not missing data per se. | Low |
| 6 | **DC lines** — `current_type` is `ac` for all 1,239 rows | The HVDC Vancouver Island link and any DC tie are **absent from this table** and must be added separately. | Medium |

---

## 4. Data-quality issues (rows to clean)

| Issue | Rows | IDs | Recommended action |
|-------|------|-----|--------------------|
| **Self-loops** (start == end) | 2 | 32350 (`BC_KI2_DSS`), 32731 (`BC_DSQ_GSS`) | Drop or collapse — zero-length self-loops break the admittance matrix. |
| **Zero TTC** (`ttc_summer = ttc_winter = 0`) | 1 | 31876 (`1L365MWN`) | Drop, or treat as out-of-service; a 0-MW `s_nom` is an infeasible flow limit. |
| **Zero reactance** | 10 | 32250, 32276, 32304, 32346, … | All are very short taps (< 0.013 km). Floor `x` at a small ε to avoid singular/near-infinite susceptance. |
| **Winter < Summer TTC** | 1 | 31607 (`1L018GLS`, 162/155) | Verify — winter rating is normally ≥ summer; likely a transcription slip. |

For reference, the healthy relationship `ttc_winter ≥ ttc_summer` holds for 1,238 of 1,239 rows (cooler ambient → higher rating), which is a good global sanity signal.

---

## 5. Notes for the base-network build

- **Segment vs line values.** `line_*` columns describe the full circuit; `line_segment_*` and `Transmission_Line_Segment_Reactance` describe *this row's* segment. To model one equivalent branch per circuit, aggregate segments (sum length, sum series reactance); to model per segment, use the segment columns. Do not mix.
- **Parallel corridors.** 165 node-pairs carry more than one line (parallel circuits or multi-segment). Your existing `aggregate_lines` step (sum `s_nom`, area-weighted impedance) is the correct handler.
- **`s_nom` mapping.** Use `ttc_summer` → `s_nom` for a summer-limited dispatch study, `ttc_winter` for winter-peak. Because these are transfer capabilities rather than MVA thermal ratings, spot-check a few known corridors (e.g. the 500 kV Peace lines) against your previous `summer_rating_mva` values before trusting the substitution wholesale.
- **Referential closure.** 381 end nodes never appear as a start (radial distribution leaves) and 42 starts never appear as an end. This is expected topology, but confirms every one of the ~900 distinct node codes must resolve in the Nodes table or the build will orphan branches.

---

## 6. Verification log

Checks run on 2026-07-23 against `transmission_lines.csv` (1,239 rows):

- Null scan: only `notes` sparse; all modelling columns 100% populated.
- Unit cross-check: `line_length_km / line_length_mi` = 1.6195 (expected 1.60934) ✓
- Monotonicity: `line_segment_length ≤ line_length` for all rows ✓
- Node-code grammar: 2,478/2,478 references parse as `BC_XXX_TYPE` ✓
- Ordering sanity: `ttc_winter ≥ ttc_summer` for 1,238/1,239 rows (1 exception logged §4)
- Duplicate scan: 0 exact duplicates on (start, end, voltage, circuit) ✓
- Anomalies logged: 2 self-loops, 1 zero-TTC, 10 zero-reactance, 1 winter<summer.
