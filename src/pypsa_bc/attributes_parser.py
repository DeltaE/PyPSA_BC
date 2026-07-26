"""
attributes_parser.py — central config access for PyPSA-BC.

Loads the four separated config files and exposes them through one object.
Nothing else in the codebase should open a YAML directly.

    config/coders.yaml        -> coders_cfg        (CODERS API pull + schema)
    config/base_network.yaml  -> base_network_cfg  (electrical assumptions)
    config/data.yaml          -> data_cfg          (file/dir locations ONLY)
    config/params.yaml        -> params_cfg         (run parameters, snapshots)

Design: each file has one responsibility, so there is no duplicated
information. Paths live only in data.yaml; assumptions only in
base_network.yaml; run knobs only in params.yaml; CODERS specifics only in
coders.yaml.
"""

import logging as log
from dataclasses import dataclass, field
from pathlib import Path

import yaml

log.basicConfig(level=log.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

CONFIG_DIR = Path("config")


@dataclass
class AttributesParser:
    """Single entry point to every config file. Child classes inherit this."""

    config_dir: str | Path = field(default=CONFIG_DIR)

    def __post_init__(self):
        d = Path(self.config_dir)
        self.coders_cfg: dict = self.load_config(d / "coders.yaml")
        self._base_network_cfg: dict = self.load_config(d / "base_network.yaml")
        self.data_cfg: dict = self.load_config(d / "data.yaml")
        self.params_cfg: dict = self.load_config(d / "params.yaml")
        self.log = log.getLogger(__name__)

    @staticmethod
    def load_config(config_file_path) -> dict:
        """Load a YAML file into a dict."""
        with open(config_file_path, "r") as file:
            return yaml.safe_load(file)

    # --- assumptions ----------------------------------------------------- #

    @property
    def base_network_cfg(self) -> dict:
        """Flat electrical assumptions: cfg['f_nom'], cfg['transformer'], ..."""
        return self._base_network_cfg

    # --- run parameters -------------------------------------------------- #

    @property
    def get_scenario(self) -> str:
        return self.params_cfg.get("scenario", "Test")

    @property
    def snapshot(self) -> tuple:
        snap = self.params_cfg["snapshot"]
        return (snap["start"][0], snap["end"][0])

    # --- reporting ------------------------------------------------------- #

    @property
    def assumption_report_save_to(self) -> Path:
        path = Path("reports/ASSUMPTIONS.md")
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def get_visual_root(self) -> Path:
        root = Path("vis/pypsa")
        root.mkdir(parents=True, exist_ok=True)
        return root

    # --- housekeeping ---------------------------------------------------- #

    @property
    def check_dirs(self):
        for directory in [
            "data/processed_data/load",
            "data/processed_data/network",
            "data/processed_data/wind",
            "data/processed_data/hydro",
            "data/processed_data/solar",
            "data/processed_data/tpp",
            "data/pypsa_data",
            "vis",
            "results/pypsa",
        ]:
            Path(directory).mkdir(parents=True, exist_ok=True)
