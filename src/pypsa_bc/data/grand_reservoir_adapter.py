"""
grand_reservoir_adapter.py
==========================

Adapter: GRanD / Global Dam Watch (GDW) + HydroLAKES  ->  Blueprint Table B
(reservoir sites) consumed by ``hydro.load_reservoir_sites``.

Purpose
-------
CODERS supplies a generation-asset inventory (turbine/plant sites, capacity,
annual-average generation) but does NOT supply reservoir *storage* bounds.
For cascade dynamics and water-use-policy studies the binding inputs are the
Table B fields ``min_storage`` / ``max_storage`` (and optionally the forebay
level range). This adapter derives those fields from open, CC-BY global
reservoir data and writes a Table-B-conforming file, so the existing pipeline
runs unchanged (no edits to ``hydro.py``).

Design contract
---------------
Output columns (index = ``asset_id``), matching ``hydro.load_reservoir_sites``:

    asset_id (index) | latitude | longitude
                     | min_storage | max_storage    (MWh)
                     | min_level   | max_level       (m, forebay elevation)
                     | cascade_order
                     | geometry                       (reservoir outlet point)

Storage unit convention
------------------------
``hydro.load_reservoir_sites`` treats storage as a scalar volume/energy and the
network later maps it onto a Store ``e_nom`` (MWh). GRanD reports live/gross
storage in cubic metres (via ``CAP_MCM`` = million m^3). We convert volume to
energy with the standard hydraulic head relation:

    E [J]  = rho * g * V * H * eta
    E [MWh] = (rho * g * H * eta / 3.6e9) * V[m^3]

so a reservoir of head ``H`` (m) and round-trip/turbine efficiency ``eta``
stores ``VOLUME_TO_MWH_PER_M3(H, eta) * V`` MWh. Head must be supplied per
reservoir (dam height is a defensible proxy; see ``head_col``). If head is
unknown the row's storage is left as NaN, which the pipeline already handles
by falling back to a single ``e_nom`` default -- a documented limitation, not
a silent guess.

This module is data-only: it performs no optimisation and imports nothing from
PyPSA, so it is unaffected by the 0.28 -> 1.x upgrade.

Author: Md Eliasinul Islam
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# Physical constants (SI)
RHO_WATER = 1000.0    # kg/m^3
G = 9.81              # m/s^2
J_PER_MWH = 3.6e9     # 1 MWh = 3.6e9 J


def volume_to_mwh_per_m3(head_m: float, efficiency: float = 0.90) -> float:
    """Energy stored per cubic metre of live storage, in MWh/m^3.

    Parameters
    ----------
    head_m : float
        Effective hydraulic head (m). Dam height is an acceptable first-order
        proxy where rated head is unavailable.
    efficiency : float
        Turbine/generation efficiency (per unit). Default 0.90.
    """
    return (RHO_WATER * G * head_m * efficiency) / J_PER_MWH


def _nearest_asset(
    reservoir_pts: gpd.GeoDataFrame,
    asset_pts: gpd.GeoDataFrame,
    max_km: float = 15.0,
) -> pd.Series:
    """Match each GRanD reservoir to the nearest CODERS asset within ``max_km``.

    Returns a Series indexed like ``reservoir_pts`` giving the matched
    ``asset_id`` (or NaN if no asset falls within the tolerance). Uses a metric
    CRS (BC Albers, EPSG:3005) for a true planar nearest-join; matches beyond
    the tolerance are dropped rather than silently mis-assigned.
    """
    res_m = reservoir_pts.to_crs(3005)
    ast_m = asset_pts.to_crs(3005)

    joined = gpd.sjoin_nearest(
        res_m, ast_m[["asset_id", "geometry"]],
        how="left", max_distance=max_km * 1000.0, distance_col="_dist_m",
    )
    # sjoin_nearest can emit >1 row per reservoir on ties; keep the closest.
    joined = joined.sort_values("_dist_m").groupby(level=0).first()
    return joined["asset_id"]


def build_table_b(
    grand_path: Union[str, Path],
    coders_assets: Union[str, Path, gpd.GeoDataFrame],
    *,
    cap_mcm_col: str = "CAP_MCM",
    head_col: Optional[str] = "DAM_HGT_M",
    efficiency: float = 0.90,
    max_match_km: float = 15.0,
    region_bbox: Optional[tuple] = None,
    out_path: Optional[Union[str, Path]] = None,
) -> gpd.GeoDataFrame:
    """Produce a Table-B-conforming reservoir file from GRanD/GDW data.

    Parameters
    ----------
    grand_path : path
        GRanD / Global Dam Watch reservoir layer (shapefile or GeoPackage).
        Expected point/polygon geometry plus attribute columns ``CAP_MCM``
        (capacity, million m^3) and, ideally, ``DAM_HGT_M`` (dam height, m).
    coders_assets : path or GeoDataFrame
        CODERS reservoir/plant assets carrying an ``asset_id`` column and
        point geometry (or ``latitude``/``longitude`` columns). Used to inherit
        ``asset_id`` and ``cascade_order`` so the output keys align with the
        existing Table A the pipeline already consumes.
    cap_mcm_col, head_col : str
        Column names in the GRanD layer. Set ``head_col=None`` to leave
        storage as NaN (pipeline falls back to single ``e_nom``).
    efficiency : float
        Generation efficiency for the volume->energy conversion.
    max_match_km : float
        Spatial tolerance for reservoir<->asset matching.
    region_bbox : (minx, miny, maxx, maxy), optional
        Lon/lat clip box (e.g. BC) applied before matching to cut global data.
    out_path : path, optional
        If given, write the result as GeoPackage/GeoJSON for the pipeline.

    Returns
    -------
    GeoDataFrame indexed by ``asset_id`` with the Table B schema.
    """
    grand = gpd.read_file(grand_path)
    if grand.crs is None:
        grand = grand.set_crs(4326)
    grand = grand.to_crs(4326)

    if region_bbox is not None:
        grand = grand.cx[region_bbox[0]:region_bbox[2], region_bbox[1]:region_bbox[3]]

    # Represent each reservoir by its outlet point (polygon -> representative point).
    res_pts = grand.copy()
    res_pts["geometry"] = res_pts.geometry.representative_point()

    # --- Load CODERS assets -------------------------------------------------
    if isinstance(coders_assets, gpd.GeoDataFrame):
        assets = coders_assets.copy()
    else:
        assets = gpd.read_file(coders_assets)
    if "geometry" not in assets or assets.geometry.isna().all():
        assets = gpd.GeoDataFrame(
            assets,
            geometry=[Point(xy) for xy in zip(assets["longitude"], assets["latitude"])],
            crs=4326,
        )
    if assets.crs is None:
        assets = assets.set_crs(4326)
    assets = assets.to_crs(4326)

    # --- Match reservoirs to assets ----------------------------------------
    res_pts["asset_id"] = _nearest_asset(res_pts, assets, max_km=max_match_km).values

    matched = res_pts.dropna(subset=["asset_id"]).copy()
    n_drop = len(res_pts) - len(matched)
    if n_drop:
        print(f">> {n_drop} GRanD reservoir(s) had no CODERS asset within "
              f"{max_match_km} km and were dropped.")

    # --- Derive storage energy ---------------------------------------------
    vol_m3 = matched[cap_mcm_col].astype(float) * 1e6  # MCM -> m^3
    if head_col is not None and head_col in matched:
        head = pd.to_numeric(matched[head_col], errors="coerce")
        e_mwh = vol_m3 * head.map(lambda h: volume_to_mwh_per_m3(h, efficiency))
    else:
        print(">> No head column supplied; max_storage left NaN "
              "(pipeline falls back to single e_nom).")
        e_mwh = pd.Series(np.nan, index=matched.index)

    # --- Assemble Table B ---------------------------------------------------
    cascade_order = (
        assets.set_index("asset_id")["cascade_order"]
        if "cascade_order" in assets else pd.Series(dtype=float)
    )

    out = gpd.GeoDataFrame(
        {
            "asset_id": matched["asset_id"].values,
            "latitude": matched.geometry.y.values,
            "longitude": matched.geometry.x.values,
            "min_storage": 0.0,                 # live-storage lower bound (MWh)
            "max_storage": e_mwh.values,        # live-storage upper bound (MWh)
            "min_level": np.nan,                # forebay range: not in GRanD
            "max_level": np.nan,
            "cascade_order": matched["asset_id"].map(cascade_order).values,
        },
        geometry=matched.geometry.values,
        crs=4326,
    ).set_index("asset_id")

    # One reservoir per asset_id: keep the largest storage on collisions.
    out = out.sort_values("max_storage", ascending=False)
    out = out[~out.index.duplicated(keep="first")]

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        driver = "GeoJSON" if out_path.suffix.lower() == ".geojson" else "GPKG"
        out.reset_index().to_file(out_path, driver=driver)
        print(f">> Table B written to {out_path} ({len(out)} reservoirs).")

    return out


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="GRanD/GDW -> Blueprint Table B adapter.")
    p.add_argument("grand", help="GRanD/GDW reservoir layer (shp/gpkg).")
    p.add_argument("coders", help="CODERS assets file with asset_id + coords.")
    p.add_argument("-o", "--out", required=True, help="Output Table B file (.gpkg/.geojson).")
    p.add_argument("--head-col", default="DAM_HGT_M", help="Head/height column in GRanD.")
    p.add_argument("--efficiency", type=float, default=0.90)
    p.add_argument("--max-km", type=float, default=15.0)
    p.add_argument("--bbox", nargs=4, type=float, default=None,
                   metavar=("MINX", "MINY", "MAXX", "MAXY"),
                   help="Lon/lat clip box, e.g. BC: -139 48 -114 60")
    a = p.parse_args()

    build_table_b(
        a.grand, a.coders,
        head_col=(None if a.head_col.lower() == "none" else a.head_col),
        efficiency=a.efficiency,
        max_match_km=a.max_km,
        region_bbox=(tuple(a.bbox) if a.bbox else None),
        out_path=a.out,
    )
