# PyPSA BC Workflow Scripts Description

This document describes the workflow scripts used in the PyPSA BC modeling framework. The scripts are executed through `run_pypsa_bc.py`, which orchestrates the entire data preparation and model building process.

## Workflow Overview

The PyPSA BC workflow is organized into three main phases:

1. **Network and Asset Preparation** (optional, controlled by `update_data` flag)
2. **Load Disaggregation** (always executed)
3. **Model Building** (always executed)

---

## Phase 1: Network and Asset Preparation Scripts

### 1. `prepare_base_network.py`

**Purpose**: Creates the foundational transmission network infrastructure for PyPSA BC.

**Key Functions**:
- Loads transmission line data from CODERS database
- Processes substation information and creates buses (nodes)
- Corrects naming inconsistencies in node names
- Adds missing transmission lines (e.g., BC_WAX_GSS to BC_SEL_TSS connection)
  > Needs tracking on CODERS data updates
- Formats transmission lines with PyPSA-compatible parameters (voltage levels, bus connections, lengths)
- Creates bus dataframe with geographical coordinates (latitude/longitude)
  
**Data Sources**: 
- CODERS transmission lines and substations datasets

**Outputs**: 
- Base network structure with buses and transmission lines
- Prepared CSV files for network topology
- Currently configured at: __data/processed_data/network__
  > Files: _buses.csv, line_types.csv, lines.csv, transformer_types.csv, transformers.csv_
 

---

### 2. `create_hydro_assets.py`

**Purpose**: Prepares hydroelectric generation asset data by consolidating information from multiple sources.

**Key Functions**:
- Creates comprehensive data dictionary template for hydro assets
- Extracts features from CODERS generators database
- Integrates cascade system information (upstream/downstream relationships)
- Adds reservoir characteristics (storage capacity, water levels, discharge rates)
- Incorporates generic generation parameters (ramping rates, O&M costs, outage rates)
- Classifies hydro types (reservoir, run-of-river, pumped storage)
  
**Data Sources**: 
- CODERS generators.csv
- CODERS hydro_cascade.csv
- CODERS existing.csv
- CODERS generation_generic.csv

**Outputs**: 
- Currently configured to : _data/processed_data/hydro/existing_
  - `hydro_generation.csv` with detailed hydro generation facility parameters
  - `hydro_reservoirs.csv` with reservoir params

---

💡 **Important**: 

Once the geographic boundaries for the reservoirs are set, the cutout needs to be created using `create_cutout.py` module. __Boundaries of the basins are used as the outer extent of the cutout__. This cutout feeds the timeseries for the generation assets.

---

### 3. `create_reservoir_inflows.py`

**Purpose**: Generates hourly inflow time series for reservoir-based hydroelectric facilities.

**Key Functions**:
- Loads hydro basin data (North America and Arctic basins)
- Uses Atlite cutout for climate/weather data
- Identifies upstream basins for each reservoir
- Creates cascade-aware inflow calculations
- Applies configurable inflow calculation methods

**Outputs**: 
- Time series of reservoir inflows (m³/s or m³/hr)
- Accounts for cascade dependencies between reservoirs

**Dependencies**: 
- Hydro sites from `create_hydro_assets.py`
- Atlite cutout

**Data Sources**:
- [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins)
- ERA5

---

### 4. `create_ror_ps.py`

**Purpose**: Generates power availability time series for run-of-river (RoR) hydroelectric plants.

**Key Functions**:
- Loads RoR and RoR-water type hydro sites
- Identifies basins associated with each RoR facility
- Calculates available power based on streamflow data from Atlite cutout
- Accounts for seasonal and temporal water availability variations

**Outputs**: 
- Hourly power availability time series for RoR facilities
- Normalized capacity factors based on water flow

**Dependencies**: 
- Hydro sites from `create_hydro_assets.py`
- Atlite cutout
- HydroBasins data

---

### 5. `create_ext_wind_assets.py`

**Purpose**: Prepares existing wind generation asset data for PyPSA modeling.

**Key Functions**:
- Filters wind generators by province/region
- Merges CODERS data with Canadian Wind Turbine Database
- Maps turbine models to OEDB (Open Energy Database) configurations
- Renames and standardizes project names for data consistency
- Adds turbine specifications (hub height, rotor diameter, power curve)

**Outputs**: 
- `wind_assets.csv` with wind farm parameters
- Includes location, capacity, turbine type, and connection points

**Data Sources**:
- CODERS generators.csv
- Canadian Wind Turbine Database (.xlsx)
- OEDB turbine configurations (JSON)

**Note**: Grouse Mountain turbine is excluded as it's non-operational

---

### 6. `create_ext_wind_ts.py`

**Purpose**: Generates hourly wind power generation time series for existing wind farms.

**Key Functions**:
- Loads Atlite climate cutout with ERA5 wind speed data
- Scales ERA5 wind speeds using Global Wind Atlas data
- Matches wind farm locations to nearest ERA5 grid points
- Applies turbine power curves to calculate generation
- Overwrites cutout wind data with scaled values

**Outputs**: 
- Hourly wind generation time series for each wind farm
- Calibrated to match reported annual average generation when possible

**Dependencies**: 
- Wind assets from `create_ext_wind_assets.py`
- Atlite cutout with wind data

---

### 7. `create_ext_solar_assets.py`

**Purpose**: Prepares existing solar photovoltaic (PV) generation asset data.

**Key Functions**:
- Filters solar generators from CODERS database
- Creates unique component and asset IDs
- Adds flag column for aggregated facilities
- Renames capacity and annual average generation columns
- Sets asset_id based on generation unit codes

**Outputs**: 
- `solar_assets.csv` with solar farm parameters
- Includes location, capacity, and connection information

**Data Sources**: 
- CODERS generators.csv (solar-type generators)

---

### 8. `create_ext_solar_ts.py`

**Purpose**: Generates hourly solar PV power generation time series.

**Key Functions**:
- Loads Atlite climate cutout with solar irradiance data
- Calculates PV generation using solar_wind utility functions
- Aggregates generation for multi-unit facilities
- Applies optional calibration to match reported annual generation

**Outputs**: 
- Hourly PV generation time series for each solar facility
- Accounts for solar angle, temperature, and irradiance variations

**Dependencies**: 
- Solar assets from `create_ext_solar_assets.py`
- Atlite cutout with solar data

---

### 9. `create_ext_tpp_assets.py`

**Purpose**: Prepares existing thermal power plant (TPP) asset data.

**Key Functions**:
- Creates TPP dictionaries for PyPSA Link components
- Maps fuel types (natural gas, biomass, biogas)
- Calculates efficiency from heat rate
- Determines marginal costs (fuel + variable O&M)
- Configures ramping limits and operational constraints
- Sets up fuel bus connections (global or local gas grid)

**Outputs**: 
- TPP configuration dictionaries for model building
- Includes natural gas, biomass, and biogas generators

**Data Sources**: 
- CODERS generators.csv
- Generation parameter files

**Note**: TPPs are modeled as Links connecting fuel buses to electricity buses

---

💡 **Important**: 

The scripts with prefix `enrich_format_X.py` loads the prepared data (as mentioned in above steps) and translates them to PyPSA compatible formats.i.e. network components. 
The outputs of these scripts are configured to save at - `data/pypsa_data` as pickle files.

---


### 10. `enrich_format_hydro.py`

**Purpose**: Converts hydro asset data into PyPSA network components with proper formatting.

**Key Functions**:
- Creates multi-component representations of hydro facilities:
  - **Water Bus**: Flow control point
  - **Reservoir Bus**: Storage interface
  - **Store Link**: Water input to storage
  - **Release Link**: Water release from storage
  - **Reservoir Store**: Water storage unit
  - **Inflow Generator**: Source of natural inflows
  - **Discharge Link**: Power generation + water discharge
  - **Spill Link**: Excess water spillage
- Manages cascade relationships between reservoirs
- Configures terminal vs. intermediate stages
- Sets storage bounds, ramping limits, and efficiency parameters

**Outputs**: 
- PyPSA-compatible hydro component dictionaries
- Structured for both reservoir and RoR facilities

**Dependencies**: 
- Hydro assets
- Reservoir inflows
- RoR power series

---

### 11. `enrich_format_tpp.py`

**Purpose**: Formats thermal power plant data into PyPSA Link components.

**Key Functions**:
- Creates Link representations connecting fuel buses to electricity buses
- Calculates efficiency (inverse of heat rate)
- Determines marginal costs combining fuel and O&M
- Configures ramping constraints
- Sets minimum generation levels and extendability flags
- Handles local vs. global gas grid configurations

**Outputs**: 
- PyPSA Link component dictionaries for TPPs
- Ready for integration into network model

**Dependencies**: 
- TPP assets from `create_ext_tpp_assets.py`
- Bus mapping dictionary

---

### 12. `enrich_format_vre.py`

**Purpose**: Formats variable renewable energy (wind and solar) assets into PyPSA Generator components.

**Key Functions**:
- Creates Generator dictionaries for wind and solar assets
- Maps assets to appropriate electricity buses
- Sets nominal capacity (p_nom)
- Configures marginal costs from generic parameters
- Links generation time series to generators
- Marks as non-extendable (existing assets)

**Outputs**: 
- PyPSA Generator component dictionaries
- Separate handling for wind and solar types

**Dependencies**: 
- Wind/solar assets
- Generation time series
- Generic generation parameters

---

## Phase 2: Load Processing

### 13. `disaggregate_load.py`

**Purpose**: Distributes provincial electricity load across regional nodes based on population and economic activity.

**Key Functions**:
- Loads BC Hydro hourly load data
- Fixes and validates hourly load profiles
- Reads CEEI (Community Energy and Emissions Inventory) data
- Filters to regional districts
- Disaggregates total provincial load proportionally to regions
- Maps load to transmission network buses
- Creates load time series for each bus

**Outputs**: 
- Regional load profiles (hourly)
- Bus-level load assignments

**Data Sources**: 
- BC Hydro balancing authority load data
- CEEI regional energy data
- Network bus mapping

**Special Handling**: 
- Identifies regional districts from CEEI data
- Accounts for Metro Vancouver and other regional districts

---

## Phase 3: Model Construction

### 14. `build_model.py`

**Purpose**: Assembles all prepared data into a complete PyPSA network model and solves the optimization.

**Key Functions**:
- Initializes empty PyPSA network
- Adds all network components in sequence:
  - Buses (from base network)
  - Transmission lines
  - Hydro assets (reservoir and RoR)
  - Wind generators
  - Solar generators
  - Thermal power plants
  - Cogeneration facilities
  - Load profiles
- Creates regional bus mappings using GADM administrative boundaries
- Handles special buses not contained in standard regions
- Configures solver settings
- Runs optimal power flow optimization
- Exports results to NetCDF format

**Key Components Added**:
1. **Base network infrastructure** (buses, lines)
2. **Hydro generation** (reservoir and RoR)
3. **Variable renewables** (wind, solar)
4. **Thermal generation** (natural gas, biomass)
5. **Load** (hourly demand at each bus)
6. **Fuel infrastructure** (gas buses and stores)

**Outputs**: 
- Complete PyPSA network object
- Optimization results (generation dispatch, flows, prices)
- Exported NetCDF file for results analysis

**Special Features**:
- Can override line capacities to unlimited (testing mode)
- Handles missing GADM region assignments for special buses
- Creates busmap dictionary linking buses to administrative regions

---

## Supporting/Utility Scripts

### 15. `create_cutout.py`

**Purpose**: Creates an Atlite climate data cutout covering the study region.

**Key Functions**:
- Determines geographic bounds from hydro basins and regional polygons
- Loads GADM administrative boundaries
- Maps regions (BC, AB, SK, MB) to GADM names
- Combines hydro basin coverage with regional boundaries
- Creates ERA5 cutout with appropriate resolution and extent

**Outputs**: 
- Atlite cutout file covering study region
- Used by wind, solar, and hydro time series generation scripts

**Data Sources**: 
- GADM country boundaries (Level 1)
- HydroBasins data
- ERA5 climate reanalysis

**Note**: Currently supports ERA5; other climate sources not yet implemented

---

### 16. `plot_results_pypsa.py`

**Purpose**: Visualizes PyPSA model results through various plots and charts.

**Key Functions**:
- Extracts generator dispatch time series from solved network
- Processes data for different generator types:
  - Hydro (Discharge links)
  - RoR hydro
  - Wind
  - Solar
  - Backstop generation
  - Natural gas
  - Biomass/biogas
- Aggregates generation by type
- Creates daily peak and average resampled data
- Generates demand plots
- Creates generation stack plots
- Visualizes inter-regional transmission flows
- Plots time series data

**Output Types**:
- Generation dispatch plots
- Load duration curves
- Transmission flow diagrams
- Time series visualizations

**Dependencies**: 
- Solved PyPSA network (NetCDF)
- Plotly for interactive visualization
- BCNexus visualization utilities

---

## Workflow Execution Order

When `run_pypsa_bc.py` is executed with `update_data=True`, scripts run in this sequence:

1. `prepare_base_network` → Creates transmission network
2. `create_hydro_assets` → Prepares hydro facility data
3. `create_cutout.py` → Prepares the cutout to the extent of hydro reservoir bounds.
4. `create_reservoir_inflows` → Calculates reservoir inflows
5. `create_ror_ps` → Generates RoR power series
6. `create_ext_wind_assets` → Prepares wind farm data
7. `create_ext_wind_ts` → Generates wind time series
8. `create_ext_solar_assets` → Prepares solar facility data
9. `create_ext_solar_ts` → Generates solar time series
10. `create_ext_tpp_assets` → Prepares thermal plant data
11. `enrich_format_hydro` → Formats hydro for PyPSA
12. `enrich_format_tpp` → Formats thermal for PyPSA
13. `enrich_format_vre` → Formats wind/solar for PyPSA

Then always:

13. `disaggregate_load` → Distributes load to buses
14. `build_model` → Assembles and solves complete model

---

## Data Flow Summary

```
CODERS Database
    ↓
[Asset Creation Scripts] → Raw asset CSVs
    ↓
Atlite Cutout + HydroBasins
    ↓
[Time Series Generation] → Hourly profiles
    ↓
[Enrichment Scripts] → PyPSA component dicts
    ↓
[Build Model] → Complete PyPSA Network
    ↓
Optimization Solver
    ↓
Results (NetCDF) → [Plotting Scripts]
```

---

## Configuration

All scripts read from a central configuration managed by `AttributesParser` class, which loads:
- `config/data.yaml` for Data file paths and sources
- PyPSA-specific configuration parameters
- Regional selections (BC, AB, SK, MB)
- Output file paths
- Solver settings

---

## Notes

- **CODERS**: Canadian Open Data for Energy Regions Database
- **GADM**: Global Administrative Areas database
- **Atlite**: A lightweight Python package for renewable power system modeling
- **ERA5**: ECMWF Reanalysis v5 climate dataset
- **HydroBasins**: Global watershed boundaries database

The workflow is designed to be modular, allowing individual scripts to be run independently if needed, though dependencies exist between phases.