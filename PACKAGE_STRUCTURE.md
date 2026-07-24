# PyPSA-BC — Package Structure & Config (v2)

**Author:** Md Eliasinul Islam
**Updated:** 2026-07-23

A clean, `uv`-recognizable layout separating the reusable **package** from the
**workflow** that orchestrates it, with a config split by responsibility.

---

## Target layout

```
PyPSA/                          ← repo root = uv project root (pyproject.toml)
├── pyproject.toml              ← uv/PEP-621; src layout, deps
├── README.md
├── src/
│   └── pypsa_bc/               ← the importable library (pip/uv installable)
│       ├── __init__.py
│       ├── attributes_parser.py   ← the ONLY thing that reads config
│       ├── utils.py
│       ├── hydro.py
│       ├── data/               ← CODERS client, adapters (coders.py, …)
│       ├── network/            ← buses, lines, transformers
│       ├── reporting/          ← assumption logger
│       └── vis/
├── workflow/                   ← orchestration ONLY (imports pypsa_bc)
│   ├── run_pypsa_bc.py
│   └── scripts/                ← create_*, enrich_format_*, build_model
├── config/                     ← four single-responsibility YAMLs
│   ├── coders.yaml
│   ├── base_network.yaml
│   ├── data.yaml
│   └── params.yaml
├── data/                       ← inputs (downloaded_data, custom) + outputs
├── tests/
└── _archive/                   ← superseded code/config (recoverable)
```

## Package vs workflow — the rule

**`src/pypsa_bc/` = reusable logic, no orchestration.** Pure functions and
classes: "given these inputs, compute this." No `main()` that runs a stage, no
hard-coded file reads. Anything you'd want to unit-test or import from a
notebook lives here.

**`workflow/` = orchestration, no logic.** Each `scripts/*.py` is a thin
`main()` that: reads paths from the config (via `AttributesParser`), calls
`pypsa_bc` functions, writes an output file. If you find real algorithm code in
`workflow/`, it belongs in the package; if you find file paths hard-coded in the
package, they belong in the config.

| Concern | Home |
|---|---|
| Hydro inflow maths, cascade routing, RoR power | `pypsa_bc/hydro.py` |
| CODERS API client, schema rename, adapters | `pypsa_bc/data/` |
| Bus/line/transformer construction | `pypsa_bc/network/` |
| "Run create_hydro_assets, write the CSV" | `workflow/scripts/create_hydro_assets.py` |
| Which file goes where | `config/data.yaml` |

## Config split (no duplicated information)

| File | Responsibility | Accessor |
|---|---|---|
| `coders.yaml` | CODERS API url, tables, schema versions | `aparser.coders_cfg` |
| `base_network.yaml` | Electrical **assumptions** (pf, f_nom, reactance table, transformer) | `aparser.base_network_cfg` |
| `data.yaml` | File/dir **locations only** | `aparser.data_cfg` |
| `params.yaml` | Run **parameters** (scenario, snapshots, per-step knobs) | `aparser.params_cfg`, `get_scenario`, `get_snapshot` |

The old `pypsa_config.yaml` + `pypsa_data.yaml` (which duplicated paths and mixed
params with directories) are archived in `_archive/config_v1/`.

## Running with uv

```bash
uv sync                                   # create env from pyproject.toml
uv run python -m workflow.run_pypsa_bc    # run the pipeline
uv run pytest                             # tests
```

## Migration status & TODO

Done: `src/pypsa_bc/` scaffolded from `codebase/pypsa_bc`; `pyproject.toml`;
config split; `attributes_parser.py` rewired to the four files; CODERS API
verified live.

Still to do (deliberately left, since you're mid-restructure):

1. **Point `workflow/` at the package.** `run_pypsa_bc.py` imports
   `from pypsa_bc import ...` — with the src layout, `uv sync` makes that resolve.
   Retire the duplicate `codebase/pypsa_bc` and `scripts/pypsa_bc` trees once the
   src copy is confirmed.
2. **Migrate param access.** A few `workflow/scripts` read run-knobs that moved
   from `data.yaml` `output.*` into `params.yaml` `workflow.*`
   (inflow_method, height, flowspeed, calibration, gas_grid, UC, build_model.*).
   Switch those reads from `data_cfg["output"][step][param]` to
   `params_cfg["workflow"][step][param]`.
3. **`create_hydro_assets` regions.** It reads
   `cfg["output"]["prepare_base_network"]["regions"]`; source that from
   `base_network_cfg["regions"]` instead.
4. **Drop the network subpackage modules** (`buses.py`, `lines.py`,
   `transformers.py`) into `src/pypsa_bc/network/`.
