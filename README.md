# PyPSA_BC
Host repository for BC PyPSA work for the PICS Decarbonization project.

The project leverages BC infrastructure from CODERS.

## Setup

This project is packaged as `pypsa-bc` and managed with [uv](https://docs.astral.sh/uv/). Quick start:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh   # install uv
git clone https://github.com/DeltaE/PyPSA_BC.git
cd PyPSA_BC
uv sync                                           # creates .venv and installs pypsa-bc + dependencies
```

See [docs/SETUP.md](docs/SETUP.md) for a full walkthrough, including running scripts (`uv run ...`) and managing dependencies.

