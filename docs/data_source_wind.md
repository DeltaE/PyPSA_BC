# Data Source — Wind turbine inputs

**Author:** Md Eliasinul Islam
**Check date:** 2026-07-23

Provenance for the two wind input files read by `create_ext_wind_assets` and
`create_ext_wind_ts`.

---

## 1. `canada_turbines.xlsx` — Canadian Wind Turbine Database (CWTD)

| Field | Value |
|---|---|
| **Dataset** | Canadian Wind Turbine Database (CWTD) |
| **Provider** | Natural Resources Canada — CanmetENERGY-Ottawa, compiled with the Centre for Applied Business Research in Energy and the Environment (CABREE), University of Alberta |
| **Access** | Government of Canada Open Government Portal (`open.canada.ca`) and the Federal Geospatial Platform (`geo.ca`) |
| **Licence** | Open Government Licence – Canada (attribution required, free to reuse/redistribute) |
| **Snapshot in use** | Worksheet `Wind Turbine Apr 27` — an April-27 release; 6,698 turbine records, all 12 provinces/territories |
| **Geometry** | Point per turbine (Latitude/Longitude, WGS84) |

**Columns (15):** `OBJECTID`, `Province/Territory`, `Project name`,
`Total project capacity (MW)`, `Turbine identifier`, `Turbine number in project`,
`Turbine rated capacity (kW)`, `Rotor diameter (m)`, `Hub height (m)`,
`Manufacturer`, `Model`, `Commissioning date`, `Latitude`, `Longitude`, `Notes`.

**Use in the pipeline.** `create_ext_wind_assets` joins this to the CODERS BC
wind generators to attach turbine physical specs (hub height, rotor diameter,
rated capacity, model); `create_ext_wind_ts` uses the model to select a power
curve for the capacity-factor time series.

**Download / refresh.** Automated via `workflow/scripts/fetch_inputs.py`, driven
by the `remote:` block in `config/data.yaml`:

```bash
python -m workflow.scripts.fetch_inputs                 # fetch all (skip existing)
# or, in a notebook:  fetch_inputs.main(only=["can_turbines"], force=True)
```

Authoritative source (NRCan FTP, XLSX):
`https://ftp.cartes.canada.ca/pub/nrcan_rncan/Wind-energy_Energie-eolienne/wind_turbines_database/Wind_Turbine_Database_FGP.xlsx`

CWTD is updated periodically; re-run with `force=True` to refresh. Note the file
is saved as `canada_turbines.xlsx`; the internal sheet name (`Wind Turbine <Mon DD>`)
records which release is in use.

---

## 2. `turbine_dict.json` — CWTD Model → OEDB power-curve crosswalk

| Field | Value |
|---|---|
| **Maps** | CWTD `Model` string → OEDB entry `"<turbine_type>*<id>"` |
| **Upstream source** | **OEDB — Open Energy Database** wind turbine library, via the Open Energy Platform (OEP); consumed by `atlite.resource.get_oedb_windturbineconfig(name, id)` |
| **Example** | `"E-82/3000" → "E-82/3000"`, `"V90/3000" → "V90/3000"`; `*id` appended only when a model resolves to multiple OEDB rows |
| **Consumer** | `pypsa_bc.wind.get_config()` splits the value on `*` and feeds name + id to atlite |
| **Regenerate** | `python -m workflow.scripts.build_turbine_dict` (pulls OEDB, matches CWTD BC models, writes the JSON + a coverage report) |
| **Licence** | OEDB / OEP data under an open licence (ODbL); attribution to the Open Energy Platform |

**Not OEDI.** The *Open Energy Data Initiative* (data.openei.org) is a separate
US-DOE/NREL repository of operational turbine data, supply curves, and met
records — it does **not** provide a model→power-curve library, so it cannot
source this file. The correct source is the OEDB/OEP turbine library above.

**Caveat.** A BC `Model` with no OEDB match gets no power curve. `build_turbine_dict`
reports unmatched models (currently the BC fleet is 9 models: `3.2 M114`,
`E-141/4200`, `E-82/3000`, `GE 2.75-120`, `GE 3.2-103`, `LTW77-1500`, `V100`,
`V90`, `V90/3000`); resolve any unmatched by hand rather than defaulting silently.

---

## Sources

- [Canadian Wind Turbine Database — geo.ca record](https://app.geo.ca/en-ca/map-browser/record/79fdad93-9025-49ad-ba16-c26d718cc070)
- [Natural Resources Canada — Open Government Portal](https://open.canada.ca/data/en/organization/nrcan-rncan)
- [CanREA — "Where to find every wind turbine in Canada"](https://renewablesassociation.ca/heres-where-to-find-every-wind-turbine-in-canada/)
