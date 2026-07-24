# Setup Guide

This project is packaged as `pypsa-bc` and managed with [uv](https://docs.astral.sh/uv/), a fast Python package/project manager. This guide walks through installing uv and setting up the project from a fresh clone.

## 1. Install uv

**macOS / Linux:**

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Alternatively, if you already have Python and pip/pipx:

```bash
pipx install uv
# or
pip install uv
```

Verify the install:

```bash
uv --version
```

## 2. Clone the repository

```bash
git clone https://github.com/DeltaE/PyPSA_BC.git
cd PyPSA_BC
```

## 3. Set up the environment

uv reads `pyproject.toml` and `.python-version` in the repo root. It will automatically download the pinned Python version (3.12) if it isn't already installed, create a `.venv`, and install all project dependencies plus the `pypsa-bc` package itself in editable mode:

```bash
uv sync
```

This produces:
- `.venv/` — the project's virtual environment
- `uv.lock` — a locked, reproducible set of dependency versions (commit this file)

To include development dependencies (pytest, ipykernel, etc.), `uv sync` already installs the `dev` dependency group by default. To skip them:

```bash
uv sync --no-dev
```

## 4. Run project code

Use `uv run` to execute anything inside the managed environment without manually activating it:

```bash
uv run python run_pypsa_bc.py
uv run snakemake -s workflow/snakefile
uv run jupyter lab
```

Or activate the virtual environment directly:

```bash
source .venv/bin/activate   # macOS/Linux
.venv\Scripts\activate      # Windows
```

## 5. Adding or updating dependencies

Add a new runtime dependency:

```bash
uv add <package>
```

Add a development-only dependency:

```bash
uv add --dev <package>
```

Update the lockfile after changing `pyproject.toml` by hand:

```bash
uv lock
```

Upgrade all dependencies to their latest allowed versions:

```bash
uv lock --upgrade
```

## Notes

- The project previously used a conda environment (`env/environment.yaml`) and a `setup.py`. These have been superseded by `pyproject.toml` + `uv.lock`. The conda environment file is kept for reference only.
- Python version is pinned via `.python-version` (3.12) to match the geospatial/scientific stack (`geopandas`, `rasterio`, `atlite`, etc.) most reliably.
