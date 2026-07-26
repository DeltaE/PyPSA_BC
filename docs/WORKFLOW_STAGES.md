# PyPSA-BC Workflow Stages

## Overview

The PyPSA-BC model builder uses a staged execution approach to clearly show progress and handle dependencies. Each stage loads specific data, performs transformations, and reports status.

## Notebook Execution Stages

### Stage 0: Module Reload
**Purpose**: Clear Python module cache to pick up code changes  
**Files Updated**: Any `.py` file in `workflow/scripts/` or `src/pypsa_bc/`  
**Output**: Print confirmation "✓ Modules reloaded"

```python
# Removes all pypsa_bc and workflow modules from sys.modules
# Ensures next import picks up latest code changes
```

---

### Stage 1: Configuration & Data Validation

#### 1A: Configuration Summary
**Purpose**: Display all loaded configuration parameters  
**Inputs**:
- `config/base_network.yaml` - Network electrical defaults
- `config/data.yaml` - Data paths and file structure
- `config/params.yaml` - Simulation parameters (scenario, snapshot time)
- `config/coders.yaml` - CODERS transmission line data

**Outputs**:
- Scenario name
- Time snapshot (start, end dates)
- Data root directories
- Base network parameters (frequency, transformer specs)
- Output file paths for each building stage

**Key Variables Displayed**:
```
- f_nom: Base frequency
- snapshot: Time range for optimization
- roots: Data directory structure
- output: File paths for intermediate and final data
```

#### 1B: Data Availability Check
**Purpose**: Verify all required input files exist before execution  
**Checks**:
1. **Base Network** - Network CSV files
2. **Network Assets** - Pickled asset definitions (hydro, wind, solar, thermal)
3. **Load Data** - Disaggregated load profiles (residential, CSMI, other)
4. **GADM Boundaries** - Geographic administrative regions for regionalization
5. **Profiles & Time Series** - Time-varying generation and load data

**Status Codes**:
- ✅ File exists (shows size in KB/MB)
- ❌ File missing (must be generated or obtained)

**Output Example**:
```
Base Network:
  ✅ data/processed_data/network (4.0 KB)
  → Stage ready: True

Load Data:
  ❌ data/processed_data/load/hourly_res.csv
  ❌ data/processed_data/load/hourly_csmi.csv
  → Stage ready: False
```

---

### Stage 2: Build PyPSA Model

This is the main workflow. The `build_model.main()` function builds the PyPSA network through these internal stages:

#### (0) Load Configuration
- Parse YAML config files
- Initialize AttributesParser

#### (1) Load Files
- Import base network from CSV files
- Load network topology (buses, lines, transformers)
- **Required Files**:
  - `data/processed_data/network/...` (network topology)

#### (2) Set Time-Slicing
- Configure temporal resolution
- Set snapshot period (start_time, end_time)
- From: `params.yaml`

#### (3) Start Adding Assets
- Begin adding generators and storage
- Initialize asset container structures

#### (4) Add Carriers
- Define energy carriers (electricity, hydrogen, etc.)
- Set up default carrier configurations

#### (5) Clustering
- Optional: Aggregate network nodes (if applied)
- Used for model simplification

#### (6) Get Bus-to-Region Mapping
- Fetch GADM administrative boundaries
- **File**: `config/data.yaml` → `GADM` section points to GeoJSON files
- Auto-fetches from pygadm if missing: `pygadm.Items(name="Canada", content_level=1/2)`
- Saves to: `data/downloaded_data/GADM/gadm41_CAN_*.geojson`

#### (7A) Add Regional Buses (Multi-Region Mode)
- Create new buses for each GADM region
- Map existing buses to regions using spatial intersection
- **When Active**: `copperplate=False`

#### (7B) Add Regional Buses (Single-Region Mode)
- Create single BC aggregation bus
- **When Active**: `copperplate=True`

#### (8) Add Load
- Load disaggregated demand profiles
- **Files**:
  - `data/processed_data/load/hourly_res.csv` - Residential
  - `data/processed_data/load/hourly_csmi.csv` - Commercial/Industrial
  - `data/processed_data/load/hourly_other.csv` - Other demand

#### (9) Add Generation Profiles
- Add time series for renewable generators
- Wind, solar, run-of-river hydroelectric
- **Files** (pickled):
  - `data/pypsa_data/hydro_ror_water.pickle`
  - `data/pypsa_data/wind.pickle`
  - `data/pypsa_data/solar.pickle`

#### (10) Optimize Network
- Solve linear optimal power flow (OPF) model
- **Requires**: Solver (gurobi, glpk, etc.)
- **Output**: Network solution with flows, dispatch, costs

#### (11) Generate Report
- Create markdown report: `reports/build_model_report_{year}.md`
- Contents:
  - Model statistics (buses, generators, lines, etc.)
  - Data sources
  - Optimization results (objective value, solver status)
  - Key findings

**Build Parameters**:
```python
build_model.main(
    copperplate=False,              # Regional or single-node
    capacity_choice="full_potential", # Scenario identifier
    tx_line_infinity=False,         # Transmission limits
    year=2021,                      # Analysis year
    solved_network_save_to="results/pypsa/model.nc"
)
```

---

### Stage 3: Load Disaggregation Pipeline

This stage processes raw electricity load data into disaggregated profiles. Follows pattern from `workflow/run_pypsa_bc.py`.

#### 3A: Load Raw Data
**Source**: `data/downloaded_data/load/bc_hydro_load/BalancingAuthorityLoad2021.xls`  
**Format**: Excel spreadsheet with hourly load  
**Output**: Pandas DataFrame

#### 3B: Clean & Fix Load Data
**Function**: `disaggregate_load.fix_hourly_load(load_bch_raw, year)`  
**Operations**:
- Normalize units (convert to MWh if needed)
- Fill missing values
- Align to hourly frequency
- Output: Annual total provincial load (MWh)

#### 3C: Disaggregate Load
**Function**: `disaggregate_load.main(total_load_MWh)`  
**Operations**:
- Split by sector (residential, commercial, other)
- Allocate by regional distribution (from CEEI data)
- Create hourly profiles for each sector/region
- **Output Files**:
  - `data/processed_data/load/hourly_res.csv` - Residential sector
  - `data/processed_data/load/hourly_csmi.csv` - Commercial/small-medium industrial
  - `data/processed_data/load/hourly_other.csv` - Other uses (public, industrial, etc.)

**Load Pipeline Parameters**:
```python
load_params = {
    "year": 2021,
    "raw_load_file": "data/downloaded_data/load/bc_hydro_load/BalancingAuthorityLoad2021.xls"
}
```

---

## Execution Flow

```
START
  ↓
[Stage 0] Reload Modules
  ├─ Clear sys.modules cache
  └─ Re-import code changes
  ↓
[Stage 1] Configuration & Data Validation
  ├─ Load and display configuration
  ├─ Check file availability
  └─ Report missing dependencies
  ↓
[Stage 3] Load Disaggregation (Optional)
  ├─ Load BC Hydro load data
  ├─ Clean and normalize
  ├─ Disaggregate by sector/region
  └─ Generate hourly load profiles
  ↓
[Stage 2] Build PyPSA Model
  ├─ Load base network topology
  ├─ Add renewable & thermal generators
  ├─ Create regional buses from GADM
  ├─ Load demand profiles
  ├─ Add generation time series
  ├─ Optimize network
  └─ Generate report
  ↓
COMPLETE
```

## Data Flow

```
Raw Data (Excel)
  ↓
[Load Disaggregation]
  ↓
CSV Load Profiles
  ↓
[Build Model]
  ├─ Network Topology (CSV)
  ├─ Asset Data (Pickle)
  ├─ Geographic Regions (GeoJSON)
  └─ Demand Profiles (CSV)
  ↓
PyPSA Network (.nc)
  ↓
[Optimize]
  ↓
Results & Report
```

## Configuration Structure

The workflow reads four YAML configuration files:

### `config/data.yaml`
```yaml
roots:
  downloaded: data/downloaded_data/
  processed: data/processed_data/
  
output:
  base_network: data/processed_data/network
  pypsa_dict:
    folder: data/pypsa_data/
    ror: hydro_ror.pickle
    res: hydro_reservoirs.pickle
    wind: wind.pickle
    solar: solar.pickle
    
GADM:
  country_file_L1: data/downloaded_data/GADM/gadm41_CAN_1.geojson
  country_file_L2: data/downloaded_data/GADM/gadm41_CAN_2.geojson
```

### `config/params.yaml`
```yaml
scenario: Test
snapshot:
  start: ['2021-01-01 00:00:00']
  end: ['2021-12-31 23:00:00']
```

### `config/base_network.yaml`
Network electrical defaults (frequency, transformer specs, etc.)

### `config/coders.yaml`
CODERS transmission line database configuration

## Error Handling

Each stage includes try-except blocks:

- **FileNotFoundError** → Stage dependencies not met
- **AttributeError** → Configuration key mismatch
- **ValueError** → Data format or value issue
- **Other Exceptions** → Traceback printed with full context

## Troubleshooting

### Missing Data Files
**Error**: `FileNotFoundError: data/processed_data/load/hourly_res.csv`  
**Solution**: Run Stage 3 (Load Disaggregation) or data prep scripts

### Missing Solver
**Error**: `SolverError` during optimization  
**Solution**: Install solver (e.g., `apt install glpk`) or set alternative solver

### GADM Download Fails
**Error**: `pygadm.Items()` returns empty or fails  
**Solution**: 
- Check internet connection
- Manually download from GADM website
- Place at: `data/downloaded_data/GADM/gadm41_CAN_*.geojson`

### Memory Issues
**Solution**: Reduce temporal resolution in `config/params.yaml` or use aggregation

## Extending the Workflow

To add custom stages:

1. Create new function in `workflow/scripts/`
2. Define clear inputs/outputs
3. Add markdown section to notebook
4. Create execution cell with try-except error handling
5. Document in this file

Example:
```python
# Stage N: Custom Processing
print("\n" + "=" * 80)
print("CUSTOM STAGE EXECUTION")
print("=" * 80)

try:
    print(f"\n🔄 Step 1: Processing...")
    result = custom_module.main()
    print(f"  ✅ Completed")
except Exception as e:
    print(f"  ❌ Error: {e}")
```

---

## Related Documents
- [Model Documentation](PyPSA_BC_Model_Documentation.md)
- [Setup Instructions](SETUP.md)
- [Known Issues](KNOWN_ISSUES.md)
