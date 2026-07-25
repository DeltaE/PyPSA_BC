# PyPSA-BC Baseline Calibration Roadmap

**Author:** Md Eliasinul Islam  
**Purpose:** Define what a 2021 baseline model means, what data is required, what calibration checks must pass, and how the calibrated baseline feeds into scenario runs (2035, 2040+).

---

## 1. Why a 2021 Baseline?

A calibrated baseline has one purpose: **establish that the model reproduces observed system behaviour before it is used to project counterfactual futures.** Without it, you cannot distinguish genuine scenario signals from modelling artefacts.

For BC-PyPSA, 2021 is the natural baseline year because:
- BCH BalancingAuthority hourly load data is available for 2021 (`BalancingAuthorityLoad2021.xls`)
- ERA5 reanalysis weather data (for VRE and hydro inflow) has 2021 as the primary processed weather year
- BCH CEEI (Community Energy and Emissions Inventory) provides 2021 spatial load disaggregation weights
- No major generation retirements or additions occurred between the historical record and the model's generation fleet representation

---

## 2. What the 2021 Baseline Must Reproduce

The calibration target is **not** perfect replication of every hourly dispatch decision (which would require unit commitment, exact fuel prices, and BC Hydro's actual reservoir operating policy). The target is structural and operational plausibility:

| Calibration metric | Target | Data source for comparison |
|---|---|---|
| **Annual system load** | ± 2% of BCH reported 2021 provincial total (~70–75 TWh) | BCH Annual Report / CEEI 2021 |
| **Monthly load shape** | Correlation ≥ 0.95 between modelled and BCH monthly loads | BCH BalancingAuthority load data |
| **Annual hydro generation** | ± 5% of BCH reported 2021 hydro generation (~60–65 TWh) | BCH Annual Report |
| **Reservoir fill/drawdown cycle** | Qualitative agreement with USGS/WSC streamflow records | WSC hydrometric data |
| **Peak demand** | ± 5% of recorded 2021 winter peak | BCH load data |
| **Backstop dispatch** | Zero — the 2021 system was fully served | (model output check) |
| **Marginal prices** | Mean price in the $20–$80/MWh range typical of hydro systems | BC spot market / academic benchmarks |
| **VRE curtailment** | Near-zero — 2021 had minimal VRE penetration | (model output check) |

---

## 3. Data Requirements

### 3.1 Load data (already available)
- `data/downloaded_data/load/bc_hydro_load/BalancingAuthorityLoad2021.xls` — hourly provincial total
- CEEI 2021 spatial disaggregation weights — already in `disaggregate_load` module

### 3.2 Hydro inflow data (partially available)
- ERA5-derived runoff profiles for 2021 — already processed as `create_reservoir_inflows` output
- WSC (Water Survey of Canada) gauge records for major river systems:
  - Peace River at Hudson Hope (`08NE001`)
  - Columbia River at Revelstoke (`08NB005`)
  - Fraser River system gauges
- **Calibration action:** Compare the model's inflow generator `p_max_pu` profiles against WSC gauge data. If shapes differ significantly (e.g., freshet timing shifted by > 2 weeks), re-extract the ERA5 cutout with corrected basin extents.

### 3.3 Generation capacity (verify for 2021)
The current network may carry 2035/2040 capacity projections. For 2021 calibration:
- Remove CLEWs investment additions (set `capacity_choice='investment'` with 2021 CLEWs run, or manually zero out future VRE expansion sites)
- Confirm thermal generator capacities match BCH 2021 fleet (CCGT, biomass, run-of-river totals)
- Site C was under construction in 2021 — **do not include** in the 2021 baseline (commissioning: 2025)

### 3.4 Transmission (verify for 2021)
- BCH 2021 rated corridor capacities (available from BCH Open Data or IRP appendices)
- The model's current corridor ratings may reflect 2024+ planned upgrades — verify against the 2021 configuration

---

## 4. Build Steps for the 2021 Baseline

```python
# In PyPSA_builder.ipynb, Dev section:

nexus_scenario      = 'Baseline_2021'
nexus_timeslices    = 0        # no CLEWs — use actual 2021 demand directly
pypsa_simulation_year = 2021

# Key flags:
# update_load = False          — use 2021 BCH actual load directly (no CLEWs scaling)
# capacity_choice = 'committed' — only BCH 2021 actual fleet (no RESource expansion)
# tx_line_infinity = False     — use actual 2021 rated corridor capacities
```

**Step-by-step:**
1. Set snapshots to 2021 (`2021-01-01` → `2021-12-31 23:00`)
2. Load BCH 2021 actual hourly demand without CLEWs scaling — `update_load=False`, pass raw `provincial_total_load_kWh_2021` directly to `disaggregate_load.main()`
3. Remove future VRE expansion sites — do not call `add_vre_expansion_sites()` for sites with `start_year > 2021`
4. Verify reservoir `e_nom` values reflect pre-Site C capacity (omit Site C)
5. Set `e_initial` on reservoir Stores to approximate 2021 January 1 storage levels (from BCH hydro reports or WSC data) — this is the most important single calibration parameter
6. Run `n.consistency_check()` → `n.optimize()` → `n.export_to_netcdf('BCPypsa_Baseline_2021.nc')`
7. Run `build_network_doc.py --network BCPypsa_Baseline_2021.nc --scenario Baseline_2021`

---

## 5. Calibration Checks (Post-Solve)

After solving the baseline, run these checks against observed 2021 data:

```python
import pypsa, pandas as pd
n = pypsa.Network("network files/BCPypsa_Baseline_2021.nc")

# 1. Annual load check
modelled_load_TWh = n.loads_t.p_set.sum().sum() / 1e6
print(f"Modelled load: {modelled_load_TWh:.2f} TWh  | Target: ~72 TWh")

# 2. Backstop dispatch check (must be zero for a feasible baseline)
backstop = [g for g in n.generators.index if "Backstop" in g]
use = n.generators_t.p[backstop].sum().sum() / 1e6
print(f"Unserved energy: {use:.3f} TWh  | Target: 0 TWh")

# 3. Monthly load shape correlation
modelled_monthly = n.loads_t.p_set.sum(axis=1).resample("ME").sum() / 1e3  # GWh
# compare against BCH monthly actuals

# 4. Hydro generation check
disc_links = n.links.index[n.links.index.str.contains("Discharge")]
hydro_TWh = (-n.links_t.p1[disc_links]).clip(lower=0).sum().sum() / 1e6
print(f"Hydro generation: {hydro_TWh:.2f} TWh  | Target: ~62 TWh (BCH 2021 Annual Report)")

# 5. Marginal price check
avg_price = n.buses_t.marginal_price[n.buses[n.buses.carrier=='AC'].index].mean().mean()
print(f"Mean nodal marginal price: ${avg_price:.1f}/MWh  | Target: $20–$80/MWh")
```

---

## 6. Sensitivity Parameters (Calibration Levers)

If the baseline fails the checks, these are the parameters to adjust in order of impact:

| Parameter | Where | Effect | Typical adjustment |
|---|---|---|---|
| `e_initial` on reservoir Stores | `n.stores.e_initial` | Sets Jan 1 storage level; strongly affects winter dispatch | Match BCH reservoir levels at year start |
| `e_nom` on reservoir Stores | `n.stores.e_nom` | Maximum storage capacity; affects seasonal range | Verify against BCH dam specifications |
| Inflow `p_max_pu` scaling | Inflow Generator profiles | Controls total annual hydro availability | ± 10–15% to match BCH reported hydro generation |
| Load disaggregation weights | CEEI 2021 fractions | Affects spatial distribution | Use 2021 CEEI (not interpolated year) |
| `s_max_pu` on AC lines | `n.lines.s_max_pu` | N-1 security margin | Try 0.85 → 0.90 if excessive congestion in baseline |

---

## 7. Pathway to Scenario Runs

Once the baseline calibration passes all checks, the scenario workflow is:

```
Baseline_2021 (calibrated)
      │
      ├── Base_CNZ_2035    → Add CLEWs 2035 investments + scaled demand
      ├── Base_CNZ_2040    → Add CLEWs 2040 investments + scaled demand
      ├── NE_FE_2035       → NE electrification load + new resources
      └── NE_FE_2040       → NE electrification load + new resources (current focus)
```

**Key principle:** The calibrated `e_initial` values and inflow scaling factors from the 2021 baseline should be **fixed** across all scenario runs (only demand and capacity fleet change). This ensures that scenario differences in dispatch outcomes are attributable to the scenario design choices, not to uncalibrated hydro parameters.

---

## 8. Critical Thinking Prompts

Before declaring the calibration complete, interrogate:

1. **Is zero backstop in the baseline trivially true or genuinely informative?** The 2021 system was feasible, so zero backstop is correct — but it only tells you the model is not infeasible, not that the dispatch is realistic. Focus on the *quantity* of hydro spill and the seasonal reservoir trajectory.

2. **Does the model's hydro flexibility look right?** BC Hydro actively manages reservoir levels for flood control, fish flows, and export contracts — behaviours the model approximates with `e_cyclic` and `e_min` but cannot fully replicate. A systematic mismatch in winter drawdown depth is more informative than an annual energy total match.

3. **Is the spatial load distribution calibrated, not just the provincial total?** A model that gets the provincial total right but misplaces load (e.g., overloading Fraser Valley while underloading Peace River Regional District) will produce wrong congestion patterns in the scenarios.

4. **How sensitive are the scenario results to the 2021 inflow year?** ERA5 2021 may be a wet or dry year relative to the long-run mean. Consider whether your scenario conclusions hold across multiple weather years (2018, 2019, 2020, 2022) as a robustness check.
