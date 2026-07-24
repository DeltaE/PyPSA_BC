# CODERS API loader — British Columbia only.
#
# The CODERS API trims by province server-side (?province=BC), so we fetch BC
# data directly instead of pulling all of Canada and masking in pandas.
#
#   Base URL : https://api.sesit.ca
#   Key      : ?key=<api_key>
#   Trim     : &province=BC          (codes: AB BC MB NB NL NS ON PE QC SK)
#   Extra    : &gen_type=<type>      (generators only; e.g. hydro_run, wind)
#
# Config contract (config_BC.yaml):
#
#   coders:
#     url:  https://api.sesit.ca
#     root: data/downloaded_data/CODERS      # base dir for all cached tables
#     api_config: data/downloaded_data/CODERS/coders_api.yaml
#     tables:                                # table_name: subfolder-under-root
#       generation_generic: supply
#       generators:         supply
#       storage:            supply
#       transmission_lines: network
#
# A table caches to:  <root>/<subfolder>/<table>_BC.pkl

import logging
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import requests
import yaml
from shapely.geometry import Point

from pypsa_bc.reporting.logger import get_logger, log_once

log = get_logger("coders")

PROVINCE = "BC"  # this model is BC-only; retarget the whole module here
CODERS_API_PATH = "data/downloaded_data/CODERS/coders_api.yaml"
CODERS_CONFIG_PATH = "config/coders.yaml"
NO_FILTER_TABLES = ["generation_generic"]  # tables that don't support extra filters

def load_config(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def load_api_key(file_path=CODERS_API_PATH):
    """Return (api_key, user) from a YAML file, or (None, None) if unavailable.

    Expected file structure:
        api_keys:
          your_username: your_api_key_here
        Default_user: your_username
    """
    try:
        cfg = load_config(file_path)
    except FileNotFoundError:
        print(f"API key file not found: {file_path}")
        return None, None
    if not cfg:
        print(f"API key file is empty: {file_path}")
        return None, None

    keys = cfg.get("api_keys", {})
    default_user = cfg.get("Default_user")
    if default_user and keys.get(default_user):
        return keys[default_user], default_user
    for user, key in keys.items():  # fallback: first non-empty key
        if key:
            return key, user
    return None, None


@dataclass
class CODERSData:
    """BC-only accessor for the CODERS API (https://api.sesit.ca).
    >>> coders = CODERSData("config/config_BC.yaml")
    >>> data = coders.load_tables()        # {'generators': gdf, ...}, all BC
    >>> coders.plot_inputs(data)           # spatial overview of inputs
    """

    coders_config_path: str | Path = field(default=CODERS_CONFIG_PATH)
    province: str = PROVINCE

    def __post_init__(self):
        self.cfg = load_config(self.coders_config_path).get("coders", {})
        self.url = self.cfg.get("url", "https://api.sesit.ca").rstrip("/")
        self.root = Path(self.cfg.get("root", "data/downloaded_data/CODERS"))
        self.tables = self.cfg.get("tables", {})  # {table_name: subfolder}
        self._cache: dict = {}                     # in-memory table cache (this run)

        # Column-schema management (immunity to CODERS field renames).
        self.versions = self.cfg.get("version", {})          # {date: {agreement, <table>: {..}}}
        self.active_version = self.cfg.get("active_version")  # optional; else last-defined
        self.required = self.cfg.get("required", {})          # {table: [canonical cols]}

        self.api_key, user = load_api_key(self.cfg.get("api_config", CODERS_API_PATH))
        if self.api_key is None:
            log_once("No CODERS API key found — fetching will fail until one is set.",
                     "coders", level=logging.WARNING)
        else:
            log_once(f"CODERS API key loaded (user: {user}); province = {self.province}", "coders")

    # --- API -------------------------------------------------------------- #

    def _query_withFilter(self, **filters) -> str:
        """Build '?province=BC&key=...&<filter>=<v>' for the current province."""
        params = {"province": self.province, "key": self.api_key}
        params.update({k: v for k, v in filters.items() if v is not None})
        return "&".join(f"{k}={v}" for k, v in params.items())
    
    def _query_withoutFilter(self) -> str:
        params = {"key": self.api_key}
        return "&".join(f"{k}={v}" for k, v in params.items())

    def fetch_data(self, table_name: str, **filters) -> pd.DataFrame:
        """Fetch a province-trimmed table from the CODERS API.

        Extra keyword filters are appended to the query (e.g. gen_type='wind').
        """
        if table_name in NO_FILTER_TABLES:
            r = requests.get(f"{self.url}/{table_name}?{self._query_withoutFilter()}")
        else:
            r = requests.get(f"{self.url}/{table_name}?{self._query_withFilter(**filters)}")
            
        r.raise_for_status()
        return pd.DataFrame.from_dict(r.json())

    def show_tables(self) -> list:
        """List all table names available in the database."""
        r = requests.get(f"{self.url}/tables?key={self.api_key}")
        r.raise_for_status()
        return r.json()

    # --- Loading ---------------------------------------------------------- #

    # --- Column-name normalization -------------------------------------- #

    def _rename_map(self, table_name: str) -> dict:
        """Return {current_source_column: canonical_column} for the active version.

        The config stores canonical->current (``agreement: old_to_new``); this
        inverts it so a DataFrame can be renamed *into* canonical names. Set
        ``agreement: new_to_old`` if a future block is authored the other way.
        """
        if not self.versions:
            return {}
        key = self.active_version or list(self.versions)[-1]  # newest defined last
        block = self.versions[key]
        mapping = block.get(table_name, {})
        if block.get("agreement", "old_to_new") == "old_to_new":
            return {src: canon for canon, src in mapping.items()}
        return dict(mapping)

    def normalize_columns(self, df: pd.DataFrame, table_name: str) -> pd.DataFrame:
        """Rename source columns to canonical names; warn loudly on schema drift."""
        rmap = self._rename_map(table_name)
        present = {src: canon for src, canon in rmap.items() if src in df.columns}
        missing = [src for src in rmap if src not in df.columns]
        key = self.active_version or (list(self.versions)[-1] if self.versions else "?")
        if missing:
            log.warning(f"{table_name}: mapped source column(s) absent from pull "
                        f"(CODERS schema drift vs version '{key}'?): {missing}")
        else:
            log.debug(f"{table_name}: all mapped source columns present (version '{key}').")
        df = df.rename(columns=present)

        need = self.required.get(table_name, [])
        lacking = [c for c in need if c not in df.columns]
        assert not lacking, (
            f"{table_name}: required canonical column(s) missing after rename: {lacking}. "
            f"Update coders.version['{key}']['{table_name}'] to the current field names."
        )
        return df


    @staticmethod
    def create_gdf(df: pd.DataFrame) -> gpd.GeoDataFrame:
        """Point GeoDataFrame from longitude/latitude columns (EPSG:4326)."""
        df = df.copy()
        df["geometry"] = df.apply(lambda r: Point(r["longitude"], r["latitude"]), axis=1)
        return gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")

    def table_path(self, table_name: str) -> Path:
        """<root>/<subfolder>/<table>.csv."""
        sub = self.tables.get(table_name, "")
        return self.root / sub / f"{table_name}.csv"

    def load_table(
        self,
        table_name: str,
        force_update: bool = False,
        as_gdf: bool = True,
        **filters,
    ) -> pd.DataFrame | gpd.GeoDataFrame:
        """Load one BC table, fetching from the API only if not cached.

        Args:
            table_name: key from the ``coders.tables`` mapping.
            force_update: re-fetch from the API and overwrite the cache.
            as_gdf: build a GeoDataFrame when lon/lat exist (skipped for lines).
            **filters: extra API filters, e.g. ``gen_type="hydro_run"``.
        """
        # In-memory cache: one shared instance reads each table at most once per run.
        cache_key = (table_name, as_gdf, tuple(sorted(filters.items())))
        if not force_update and cache_key in self._cache:
            return self._cache[cache_key].copy()

        path = self.table_path(table_name)
        path.parent.mkdir(parents=True, exist_ok=True)

        if path.is_file() and not force_update:
            data = pd.read_csv(path)
            log.debug(f"{table_name}: read disk cache {path}")
        else:
            data = self.fetch_data(table_name, **filters)
            data.to_csv(path, index=False)  # disk cache of RAW source columns
            log.info(f"{table_name}: fetched {len(data)} {self.province} rows -> {path}")

        # Normalize on every load so a config edit re-maps without re-fetching.
        data = self.normalize_columns(data, table_name)

        if as_gdf and "lines" not in table_name and {"longitude", "latitude"}.issubset(data.columns):
            data = self.create_gdf(data)

        self._cache[cache_key] = data
        return data.copy()  # hand callers a private copy; mutations don't poison cache

    def load_tables(
        self,
        tables: list[str] | None = None,
        force_update: bool = False,
    ) -> dict[str, pd.DataFrame | gpd.GeoDataFrame]:
        """Load every table in the ``coders.tables`` mapping (all BC)."""
        tables = tables or list(self.tables)
        out = {t: self.load_table(t, force_update=force_update) for t in tables}
        log.debug(f"Loaded {len(out)} {self.province} tables: {', '.join(out)}")
        return out

    # --- Visualization ---------------------------------------------------- #

    def plot_inputs(self, data=None, ax=None, save_path=None):
        """Spatial overview of the loaded CODERS input tables (BC).

        Point tables (generators, storage, ...) are scattered; transmission
        lines are drawn from endpoint coordinates where available.
        """
        if data is None:
            data = self.load_tables()
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 9))

        for name, tbl in data.items():  # lines underneath
            if "lines" not in name:
                continue
            for x0, y0, x1, y1 in [("lon_from", "lat_from", "lon_to", "lat_to"),
                                   ("from_lon", "from_lat", "to_lon", "to_lat")]:
                if {x0, y0, x1, y1}.issubset(tbl.columns):
                    for _, r in tbl.iterrows():
                        ax.plot([r[x0], r[x1]], [r[y0], r[y1]], color="0.6", lw=0.6, zorder=1)
                    break
            else:
                print(f">> {name}: no endpoint coords; lines not drawn.")

        for name, tbl in data.items():  # points on top
            if "lines" in name:
                continue
            if isinstance(tbl, gpd.GeoDataFrame) and tbl.geometry.notna().any():
                tbl.plot(ax=ax, markersize=14, alpha=0.75, label=name, zorder=2)
            elif {"longitude", "latitude"}.issubset(tbl.columns):
                ax.scatter(tbl["longitude"], tbl["latitude"], s=14, alpha=0.75, label=name, zorder=2)

        ax.set(title=f"CODERS inputs — {self.province}", xlabel="Longitude", ylabel="Latitude")
        ax.legend(fontsize=8, markerscale=1.5)
        if save_path:
            ax.figure.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f">> Input map saved to {save_path}")
        return ax

    @property
    def transmission_lines(self)->pd.DataFrame:
        """Return the transmission_lines table as a DataFrame."""
        return self.load_table("transmission_lines", as_gdf=False)
    
    @property
    def substations(self)->pd.DataFrame:
        """Return the substations table as a DataFrame."""
        return self.load_table("substations", as_gdf=False)
    
    @property
    def generators(self) -> gpd.GeoDataFrame:
        """Return the generators table as a GeoDataFrame."""
        return self.load_table("generators", as_gdf=True)
    
    @property
    def generation_generic(self) -> pd.DataFrame:
        """Return the generation_generic table as a DataFrame."""
        return self.load_table("generation_generic", as_gdf=False)

    @property
    def existing_hydro(self) -> pd.DataFrame:
        """Return the existing hydro generators as a DataFrame."""
        data=self.generators[self.generators.generator_fuel_type.str.contains('hydro')]
        save_to = self.table_path("existing_hydro")
        data.to_csv(save_to, index=False)  # cache RAW source columns
        log.debug(f"existing_hydro: {len(data)} {self.province} rows -> {save_to}")
        return data


# ---------------------------------------------------------------- singleton

_CODERS_SINGLETON: "CODERSData | None" = None


def get_coders(province: str = PROVINCE, config_path: str = CODERS_CONFIG_PATH) -> "CODERSData":
    """Return one shared CODERSData instance (created on first call).

    Modules should call ``get_coders()`` rather than ``CODERSData()`` so the API
    key is loaded once and every table is read/normalized at most once per run
    (cached in-memory). This removes the repeated key-load and re-parsing that
    happened when each module built its own instance.
    """
    global _CODERS_SINGLETON
    if _CODERS_SINGLETON is None:
        _CODERS_SINGLETON = CODERSData(coders_config_path=config_path, province=province)
    return _CODERS_SINGLETON
