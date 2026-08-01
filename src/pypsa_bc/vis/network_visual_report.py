from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import shutil

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pypsa
from shapely.geometry import LineString, Point

try:
    import contextily as ctx
except Exception:  # optional
    ctx = None

BC_ALBERS = "EPSG:3005"
CONGESTION_THRESHOLD = 0.85
# Conservative coordinate envelope used before the boundary test. It removes
# sentinel (0, 0), foreign interties, and malformed coordinates even when the
# optional BC polygon dataset has not been fetched.
BC_COORDINATE_ENVELOPE = (-139.2, -114.0, 48.2, 60.1)

# ---------------------------------------------------------------------------
# Shared presentation constants. These mirror the documentation website so the
# interactive report reads as one product: sans-serif type, forest-green
# accent, and a single, consistently placed selector for every explorer.
# ---------------------------------------------------------------------------
REPORT_FONT_FAMILY = "Arial, Helvetica, sans-serif"
REPORT_INK = "#10251f"
REPORT_MUTED = "#5c6d67"
REPORT_ACCENT = "#0b5d4b"
SELECTOR_BG = "#ffffff"
SELECTOR_BORDER = "#c9d6cf"


def _apply_selector_layout(
    fig: go.Figure,
    *,
    title: str,
    selector_label: str,
    buttons: list[dict[str, Any]],
    xaxis_title: str | None = None,
    yaxis_title: str | None = None,
    extra_layout: dict[str, Any] | None = None,
) -> go.Figure:
    """Apply the standardized explorer layout used across all dropdown figures.

    Vertical order (top to bottom): left-aligned title, small selector label,
    then the dropdown, all anchored to the left edge with a generous top margin
    so nothing overlaps regardless of figure height.
    """
    layout: dict[str, Any] = dict(
        title=dict(
            text=title,
            x=0.0,
            xanchor="left",
            xref="paper",
            y=0.97,
            yanchor="top",
            yref="container",
            font=dict(family=REPORT_FONT_FAMILY, size=17, color=REPORT_INK),
        ),
        template="plotly_white",
        font=dict(family=REPORT_FONT_FAMILY, size=13, color=REPORT_INK),
        margin=dict(l=72, r=32, t=140, b=62),
        showlegend=False,
        updatemenus=[
            dict(
                type="dropdown",
                direction="down",
                x=0.0,
                xanchor="left",
                y=1.16,
                yanchor="top",
                buttons=buttons,
                bgcolor=SELECTOR_BG,
                bordercolor=SELECTOR_BORDER,
                borderwidth=1,
                font=dict(family=REPORT_FONT_FAMILY, size=12, color=REPORT_INK),
                pad=dict(t=2, b=2, l=6, r=6),
                showactive=True,
            )
        ],
        annotations=[
            dict(
                text=f"<b>{selector_label}</b>",
                x=0.0,
                xanchor="left",
                y=1.23,
                yanchor="bottom",
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(family=REPORT_FONT_FAMILY, size=11, color=REPORT_MUTED),
            )
        ],
    )
    if xaxis_title is not None:
        layout["xaxis_title"] = xaxis_title
    if yaxis_title is not None:
        layout["yaxis_title"] = yaxis_title
    if extra_layout:
        layout.update(extra_layout)
    fig.update_layout(**layout)
    return fig


@dataclass
class PlotArtifact:
    title: str
    description: str
    relative_path: str
    kind: str  # static-image, interactive-html, or file-link
    category: str = "system"


@dataclass
class ReportResult:
    network_path: Path
    report_path: Path
    resources_dir: Path
    artifacts: list[PlotArtifact]


def _configure_plot_style() -> None:
    installed_fonts = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams["font.family"] = "Cambria" if "Cambria" in installed_fonts else "DejaVu Serif"
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["legend.frameon"] = False


def _clean_axis(ax: plt.Axes, hide_axes: bool = False) -> None:
    if hide_axes:
        ax.set_axis_off()
        return
    for side in ["top", "right"]:
        if side in ax.spines:
            ax.spines[side].set_visible(False)
    ax.grid(alpha=0.2, linewidth=0.5)


def _add_grey_basemap(ax: plt.Axes) -> None:
    """Apply a deterministic neutral map background without remote tile calls."""
    ax.set_facecolor("#e7e9e7")


def _set_bc_map_frame(
    ax: plt.Axes,
    bc_regions: gpd.GeoDataFrame | None,
    pad_fraction: float = 0.02,
) -> None:
    """Use equal projected units and a stable provincial frame."""
    ax.set_aspect("equal", adjustable="box")
    if bc_regions is None or len(bc_regions) == 0:
        return
    west, south, east, north = bc_regions.to_crs(BC_ALBERS).total_bounds
    width, height = east - west, north - south
    ax.set_xlim(west - width * pad_fraction, east + width * pad_fraction)
    ax.set_ylim(south - height * pad_fraction, north + height * pad_fraction)


def _set_geometry_map_frame(
    ax: plt.Axes,
    geometry: gpd.GeoSeries,
    pad_fraction: float = 0.07,
) -> None:
    """Frame plotted BC content without changing projected x/y scale."""
    west, south, east, north = geometry.total_bounds
    width, height = east - west, north - south
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(west - width * pad_fraction, east + width * pad_fraction)
    ax.set_ylim(south - height * pad_fraction, north + height * pad_fraction)


def _find_logo_path() -> str | None:
    logo = Path("docs/assets/delta-e.png")
    return str(logo.as_posix()) if logo.exists() else None


def _fetch_calibration_load_file(year: int, calibration_dir: Path) -> Path | None:
    target = calibration_dir / "bc_hydro" / f"BalancingAuthorityLoad{year}.xls"
    return target if target.exists() and target.stat().st_size > 0 else None


def _parse_bch_hourly_load_kwh(xls_path: Path, year: int) -> pd.Series | None:
    try:
        raw = pd.read_excel(xls_path)
    except Exception:
        return None

    numeric_col = None
    for col in raw.columns:
        s = pd.to_numeric(raw[col], errors="coerce")
        if s.notna().sum() >= 1000:
            numeric_col = s
            break
    if numeric_col is None:
        return None

    cleaned = numeric_col.dropna().reset_index(drop=True)
    idx = pd.date_range(start=f"{year}-01-01 00:00:00", end=f"{year}-12-31 23:00:00", freq="h")
    if len(cleaned) < len(idx):
        idx = idx[: len(cleaned)]
    elif len(cleaned) > len(idx):
        cleaned = cleaned.iloc[: len(idx)]

    series = pd.Series(cleaned.values * 1000.0, index=idx, name="observed_load_kWh")
    return series


def _write_load_calibration_plot(
    n: pypsa.Network,
    output_png: Path,
    calibration_dir: Path,
    year: int,
) -> bool:
    if n.loads_t.p.empty:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "No model load timeseries found.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return False

    model = _ensure_datetime_index(n.loads_t.p).sum(axis=1)
    if isinstance(model.index, pd.DatetimeIndex):
        model_monthly = model.resample("ME").mean()
    else:
        model_monthly = pd.Series(model.values)

    obs_file = _fetch_calibration_load_file(year, calibration_dir)
    obs = _parse_bch_hourly_load_kwh(obs_file, year) if obs_file else None
    has_obs = obs is not None and len(obs) > 0
    obs_monthly = obs.resample("ME").mean() if has_obs else None

    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(model_monthly.index, model_monthly.values, linewidth=1.6, label="Model load")
    if has_obs and obs_monthly is not None:
        ax.plot(obs_monthly.index, obs_monthly.values, linewidth=1.6, linestyle="--", label="Observed load")
    ax.set_title(f"Load calibration against observed balancing authority load ({year})")
    ax.set_ylabel("Load [kWh]")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)
    return bool(has_obs)


def _load_bc_gadm_regions() -> gpd.GeoDataFrame | None:
    candidate_paths = [
        Path("data/downloaded_data/GADM/gadm41_CAN_2.geojson"),
        Path("data/downloaded_data/GADM/gadm41_CAN_1.geojson"),
        Path("data/downloaded_data/GADM/gadm_can_l2.geojson"),
    ]

    for path in candidate_paths:
        if not path.exists():
            continue
        try:
            gdf = gpd.read_file(path)
        except Exception:
            continue

        if "NAME_1" in gdf.columns:
            mask = gdf["NAME_1"].astype(str).str.lower().eq("british columbia")
            subset = gdf.loc[mask].copy()
            if len(subset) > 0:
                if subset.crs is None:
                    subset = subset.set_crs("EPSG:4326")
                return subset.to_crs("EPSG:4326")
        if "NAME_0" in gdf.columns:
            subset = gdf.copy()
            if subset.crs is None:
                subset = subset.set_crs("EPSG:4326")
            return subset.to_crs("EPSG:4326")
    return None


def _get_bc_filtered_bus_coordinates(
    n: pypsa.Network,
    bc_regions: gpd.GeoDataFrame | None,
) -> tuple[pd.Series, pd.Series, pd.Index]:
    if n.buses.empty or "x" not in n.buses.columns or "y" not in n.buses.columns:
        return pd.Series(dtype=float), pd.Series(dtype=float), pd.Index([])

    x = pd.to_numeric(n.buses["x"], errors="coerce")
    y = pd.to_numeric(n.buses["y"], errors="coerce")
    west, east, south, north = BC_COORDINATE_ENVELOPE
    buses_xy = (
        x.notna()
        & y.notna()
        & x.between(west, east)
        & y.between(south, north)
    )
    x = x.loc[buses_xy]
    y = y.loc[buses_xy]

    if bc_regions is None or len(bc_regions) == 0:
        idx = x.index
        return x.loc[idx], y.loc[idx], idx

    points = gpd.GeoSeries(gpd.points_from_xy(x, y), index=x.index, crs="EPSG:4326")
    region_geom = bc_regions.to_crs("EPSG:4326").unary_union
    inside_mask = points.within(region_geom) | points.touches(region_geom)
    inside_idx = points.index[inside_mask]
    return x.loc[inside_idx], y.loc[inside_idx], inside_idx


def _iter_polygon_rings(geom):
    if geom is None or geom.is_empty:
        return
    gtype = geom.geom_type
    if gtype == "Polygon":
        yield np.asarray(geom.exterior.coords)
    elif gtype == "MultiPolygon":
        for poly in geom.geoms:
            if not poly.is_empty:
                yield np.asarray(poly.exterior.coords)


def _sanitize_filename(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in name)


def _ensure_datetime_index(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    if isinstance(df.index, pd.DatetimeIndex):
        return df
    try:
        converted = pd.to_datetime(df.index)
        out = df.copy()
        out.index = converted
        return out
    except Exception:
        return df


def _get_generator_tech_labels(n: pypsa.Network) -> pd.Series:
    labels = pd.Series("", index=n.generators.index, dtype=object)
    # PyPSA-BC currently stores its useful technology label in ``type`` while
    # ``carrier`` is blank for most generators. Prefer populated values rather
    # than merely preferring a column that exists.
    for column in ["tech", "carrier", "type"]:
        if column not in n.generators.columns:
            continue
        candidate = n.generators[column].fillna("").astype(str).str.strip()
        labels = labels.mask(labels.eq("") & candidate.ne(""), candidate)
    labels = labels.replace("", "unknown").fillna("unknown")
    return labels


def _resource_color_map() -> dict[str, str]:
    return {
        "hydro": "#1d3557",
        "ror_hydro": "#1d3557",
        "reservoir": "#264653",
        "wind": "#2a9d8f",
        "solar": "#e9c46a",
        "pv": "#e9c46a",
        "nuclear": "#6d597a",
        "tpp": "#e76f51",
        "cogen": "#f4a261",
        "backstop": "#9d0208",
        "other": "#6c757d",
        "unknown": "#6c757d",
    }


def _normalize_resource_name(name: str) -> str:
    text = str(name).lower()
    if "ror" in text:
        return "ror_hydro"
    if "hydro" in text or "reservoir" in text:
        return "hydro"
    if "wind" in text:
        return "wind"
    if "solar" in text or "pv" in text:
        return "solar"
    if "nuclear" in text:
        return "nuclear"
    if "cogen" in text:
        return "cogen"
    if "tpp" in text or "gas" in text or "coal" in text or "diesel" in text:
        return "tpp"
    if "backstop" in text:
        return "backstop"
    return "other"


def _iter_component_dicts(value):
    if isinstance(value, dict):
        if "class_name" in value:
            yield value
        else:
            for child in value.values():
                yield from _iter_component_dicts(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _iter_component_dicts(child)


def _load_prepared_generation_components() -> pd.DataFrame:
    """Read prepared PyPSA dictionaries without relying on the solved network."""
    rows: list[dict[str, Any]] = []
    for path in Path("data/pypsa_data").glob("*.pickle"):
        try:
            payload = pd.read_pickle(path)
        except Exception:
            continue
        for comp in _iter_component_dicts(payload):
            if comp.get("class_name") not in {"Generator", "Link"}:
                continue
            bus = comp.get("bus") or comp.get("bus1")
            capacity = pd.to_numeric(comp.get("p_nom", 0.0), errors="coerce")
            if not bus or pd.isna(capacity) or float(capacity) <= 0:
                continue
            resource = _normalize_resource_name(
                comp.get("type") or comp.get("carrier") or path.stem
            )
            rows.append(
                {
                    "name": comp.get("name", path.stem),
                    "bus": str(bus),
                    "capacity_mw": float(capacity),
                    "resource": resource,
                    "source": f"data/pypsa_data/{path.name}",
                }
            )
    return pd.DataFrame(rows)


def _load_processed_resource_points(buses: pd.DataFrame) -> pd.DataFrame:
    """Create a reproducible resource inventory from processed input tables."""
    inputs = [
        ("hydro", "existing", Path("data/processed_data/hydro/existing/hydro_generation.csv")),
        ("wind", "existing", Path("data/processed_data/wind/existing/bc_ext_wind_assets.csv")),
        ("solar", "existing", Path("data/processed_data/solar/existing/bc_ext_solar_assets.csv")),
        ("tpp", "existing", Path("data/processed_data/tpp/existing/bc_ext_tpp_assets.csv")),
        ("wind", "committed", Path("data/processed_data/wind/committed/BCH_CFP24_wind.csv")),
        ("solar", "committed", Path("data/processed_data/solar/committed/BCH_CFP24_solar.csv")),
        ("wind", "potential", Path("data/processed_data/wind/potential/resource_options_wind.csv")),
        ("solar", "potential", Path("data/processed_data/solar/potential/resource_options_solar.csv")),
    ]
    capacity_columns = [
        "capacity",
        "generator_capacity",
        "potential_capacity",
        "facility_installed_capacity",
        "install_capacity_in_mw",
    ]
    name_columns = [
        "project_name",
        "generation_facility_name",
        "asset_id",
        "cluster_id",
        "component_id",
    ]
    bus_columns = ["connecting_node_code", "network_node_code", "nearest_station", "node_code"]
    bus_lookup = buses.reset_index().copy()
    def station_key(value: Any) -> str | None:
        parts = str(value).upper().split("_")
        if len(parts) < 2:
            return None
        # Prepared buses use ``voltage_CODE_type`` while resource tables use
        # ``BC_CODE_type``. Their shared stable identifier is CODE.
        return parts[1] if len(parts[1]) == 3 else None

    bus_lookup["station_key"] = bus_lookup["name"].map(station_key)
    bus_lookup = bus_lookup.dropna(subset=["station_key"]).drop_duplicates("station_key")
    bus_lookup = bus_lookup.set_index("station_key")
    rows: list[dict[str, Any]] = []
    for resource, stage, path in inputs:
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        cap_col = next((c for c in capacity_columns if c in frame), None)
        name_col = next((c for c in name_columns if c in frame), None)
        bus_col = next((c for c in bus_columns if c in frame), None)
        for _, record in frame.iterrows():
            capacity = pd.to_numeric(record.get(cap_col, np.nan), errors="coerce")
            lon = pd.to_numeric(record.get("longitude", np.nan), errors="coerce")
            lat = pd.to_numeric(record.get("latitude", np.nan), errors="coerce")
            bus_name = str(record.get(bus_col, "")) if bus_col else ""
            if pd.isna(lon) or pd.isna(lat):
                key = station_key(bus_name)
                if key is not None and key in bus_lookup.index:
                    lon = pd.to_numeric(bus_lookup.at[key, "x"], errors="coerce")
                    lat = pd.to_numeric(bus_lookup.at[key, "y"], errors="coerce")
            if pd.isna(lon) or pd.isna(lat) or pd.isna(capacity) or float(capacity) <= 0:
                continue
            rows.append(
                {
                    "name": str(record.get(name_col, bus_name or path.stem)),
                    "resource": resource,
                    "stage": stage,
                    "capacity_mw": float(capacity),
                    "longitude": float(lon),
                    "latitude": float(lat),
                    "source": path.as_posix(),
                }
            )
    return pd.DataFrame(rows)


def _write_input_detailed_map(
    n: pypsa.Network,
    output_png: Path,
    *,
    full_bc_frame: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(11, 8.5))
    bc_regions = _load_bc_gadm_regions()
    prepared_buses_path = Path("data/processed_data/network/buses.csv")
    prepared_lines_path = Path("data/processed_data/network/lines.csv")
    use_prepared = prepared_buses_path.exists() and prepared_lines_path.exists()

    if use_prepared:
        buses = pd.read_csv(prepared_buses_path).set_index("name")
        lines = pd.read_csv(prepared_lines_path)
        x = pd.to_numeric(buses["x"], errors="coerce")
        y = pd.to_numeric(buses["y"], errors="coerce")
        west, east, south, north = BC_COORDINATE_ENVELOPE
        valid_mask = (
            x.notna()
            & y.notna()
            & x.between(west, east)
            & y.between(south, north)
        )
        if bc_regions is not None and len(bc_regions) > 0:
            bus_points = gpd.GeoSeries(
                gpd.points_from_xy(x, y), index=buses.index, crs="EPSG:4326"
            )
            bc_geometry = bc_regions.to_crs("EPSG:4326").geometry.union_all()
            valid_mask &= bus_points.within(bc_geometry) | bus_points.touches(bc_geometry)
        valid_buses = buses.index[valid_mask]
        x, y = x.loc[valid_buses], y.loc[valid_buses]
    else:
        lines = n.lines
        x, y, valid_buses = _get_bc_filtered_bus_coordinates(n, bc_regions)

    if len(valid_buses) == 0:
        ax.text(0.5, 0.5, "No valid BC bus coordinates available.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    if bc_regions is not None and len(bc_regions) > 0:
        base = bc_regions.to_crs(BC_ALBERS)
        base.plot(ax=ax, color="#ecebe6", edgecolor="#a7aaa7", linewidth=0.45, zorder=0)

    # Draw all retained electrical corridors in a metre-based analytical CRS.
    line_segments = []
    for _, line in lines.iterrows():
        b0, b1 = line.get("bus0"), line.get("bus1")
        if b0 in valid_buses and b1 in valid_buses:
            line_segments.append(LineString([(x.loc[b0], y.loc[b0]), (x.loc[b1], y.loc[b1])]))
    if line_segments:
        projected_lines = gpd.GeoSeries(line_segments, crs="EPSG:4326").to_crs(BC_ALBERS)
        collection = LineCollection(
            [np.asarray(geom.coords) for geom in projected_lines if not geom.is_empty],
            colors="#50575a",
            linewidths=0.7,
            alpha=0.78,
            zorder=2,
        )
        ax.add_collection(collection)
        ax.autoscale_view()

    pts = gpd.GeoDataFrame(
        {"bus": valid_buses},
        geometry=gpd.points_from_xy(x.loc[valid_buses], y.loc[valid_buses]),
        crs="EPSG:4326",
    ).to_crs(BC_ALBERS)
    if full_bc_frame:
        _set_bc_map_frame(ax, bc_regions, pad_fraction=0.045)
    else:
        _set_geometry_map_frame(ax, pts.geometry)
    _add_grey_basemap(ax)

    # Generator inventory: processed asset/resource tables are preferred because
    # they preserve technology, project stage, coordinates, and source lineage.
    resource_points = _load_processed_resource_points(buses) if use_prepared else pd.DataFrame()
    prepared_generation = _load_prepared_generation_components()
    if resource_points.empty and prepared_generation.empty:
        gen_tmp = n.generators.copy()
        gen_tmp["resource"] = _get_generator_tech_labels(n).map(_normalize_resource_name)
        gen_tmp["capacity_mw"] = pd.to_numeric(gen_tmp.get("p_nom", 0), errors="coerce").fillna(0)
        gen_tmp["bus"] = gen_tmp["bus"].astype(str)
    elif resource_points.empty:
        gen_tmp = prepared_generation
    else:
        gen_tmp = pd.DataFrame()

    cmap = _resource_color_map()
    # Substations remain visible even where no generator is connected.
    ax.scatter(
        pts.geometry.x,
        pts.geometry.y,
        s=10,
        facecolor="#f8f7f2",
        edgecolor="#4b5563",
        linewidth=0.55,
        zorder=3,
        label="substation / bus",
    )

    if not resource_points.empty:
        west, east, south, north = BC_COORDINATE_ENVELOPE
        resource_points = resource_points.loc[
            resource_points["longitude"].between(west, east)
            & resource_points["latitude"].between(south, north)
        ].copy()
        generator_pts = gpd.GeoDataFrame(
            resource_points,
            geometry=gpd.points_from_xy(
                resource_points["longitude"], resource_points["latitude"]
            ),
            crs="EPSG:4326",
        ).to_crs(BC_ALBERS)
        # Square-root scaling preserves capacity meaning without allowing large
        # potential clusters to obscure the existing transmission network.
        generator_pts["size"] = 18 + 5.5 * np.sqrt(
            generator_pts["capacity_mw"].clip(lower=0, upper=5000)
        )
        stage_markers = {"existing": "o", "committed": "D", "potential": "^"}
        for (res, stage), group in generator_pts.groupby(["resource", "stage"]):
            is_potential = stage == "potential"
            ax.scatter(
                group.geometry.x,
                group.geometry.y,
                s=group["size"],
                marker=stage_markers[stage],
                facecolors="none" if is_potential else cmap.get(str(res), "#6c757d"),
                edgecolors=cmap.get(str(res), "#6c757d") if is_potential else "white",
                linewidths=1.1 if is_potential else 0.45,
                alpha=0.88,
                label=f"{res} · {stage}",
                zorder=4,
            )
    else:
        cap = gen_tmp.groupby("bus")["capacity_mw"].sum().reindex(valid_buses).fillna(0)
        dominant_resource = pd.Series("other", index=valid_buses)
        mix = gen_tmp.groupby(["bus", "resource"])["capacity_mw"].sum().unstack(fill_value=0)
        for bus in mix.index.intersection(valid_buses):
            dominant_resource.loc[bus] = mix.loc[bus].idxmax()
        generator_pts = pts.loc[pts["bus"].map(cap).fillna(0) > 0].copy()
        generator_pts["resource"] = dominant_resource.loc[generator_pts["bus"]].values
        generator_pts["size"] = 18 + 5.5 * np.sqrt(
            generator_pts["bus"].map(cap).clip(lower=0, upper=5000)
        )
        for res, group in generator_pts.groupby("resource"):
            ax.scatter(
                group.geometry.x,
                group.geometry.y,
                s=group["size"],
                c=cmap.get(str(res), "#6c757d"),
                edgecolors="white",
                linewidths=0.4,
                label=str(res),
                zorder=4,
            )

    ax.set_title("Detailed Input Network", loc="left", fontsize=14, fontweight="bold")
    _clean_axis(ax, hide_axes=True)
    ax.legend(loc="upper left", ncol=3, fontsize=8, frameon=True, facecolor="#f7f8f3")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_ceei_heatmap_plot(output_png: Path) -> bool:
    ceei_path = Path("data/processed_data/load/CEEI_RD_ELEC_proportions.csv")
    bc_regions = _load_bc_gadm_regions()
    if not ceei_path.exists() or bc_regions is None or len(bc_regions) == 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "CEEI proportions or BC regions not found.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return False

    df = pd.read_csv(ceei_path)
    if "REGION" not in df.columns:
        df.rename(columns={df.columns[0]: "REGION"}, inplace=True)

    for c in ["PROPORTION_RES", "PROPORTION_CSMI"]:
        if c not in df.columns:
            df[c] = 0.0
    df["TOTAL_PROP"] = df["PROPORTION_RES"].fillna(0) + df["PROPORTION_CSMI"].fillna(0)

    def norm(s: str) -> str:
        return "".join(ch for ch in str(s).lower() if ch.isalnum())

    rename_map = {
        "comoxstrathcona": "strathcona",
        "greatervancouver": "metrovancouver",
        "skeenaqueencharlotte": "haida gwaii",
    }

    df["REGION_KEY"] = df["REGION"].map(norm).replace(rename_map)
    regions = bc_regions.copy().to_crs("EPSG:4326")
    if "NAME_2" not in regions.columns:
        return False
    regions["REGION_KEY"] = regions["NAME_2"].map(norm)

    merged = regions.merge(df[["REGION_KEY", "TOTAL_PROP"]], on="REGION_KEY", how="left").to_crs(BC_ALBERS)
    merged["TOTAL_PROP"] = merged["TOTAL_PROP"].fillna(0)

    fig, ax = plt.subplots(figsize=(11, 7))
    merged.plot(column="TOTAL_PROP", cmap="YlGnBu", linewidth=0.5, edgecolor="#7c8580", legend=True, ax=ax)
    _add_grey_basemap(ax)
    ax.set_title("CEEI Load Allocation by Regional District", loc="left", fontsize=14, fontweight="bold")
    _clean_axis(ax, hide_axes=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)
    return True


def _write_hydro_generation_interactive(n: pypsa.Network, output_html: Path) -> None:
    if n.generators.empty or n.generators_t.p.empty:
        go.Figure(layout=dict(title="No hydro generation data available")).write_html(output_html, include_plotlyjs="cdn")
        return

    tech = _get_generator_tech_labels(n).astype(str)
    hydro_mask = tech.str.lower().str.contains("hydro|ror|reservoir") | n.generators.index.astype(str).str.lower().str.contains("hydro|ror|reservoir")
    hydro_stations = list(n.generators.index[hydro_mask])
    if not hydro_stations:
        go.Figure(layout=dict(title="No hydro stations detected in this network")).write_html(output_html, include_plotlyjs="cdn")
        return

    ts = _ensure_datetime_index(n.generators_t.p[hydro_stations].clip(lower=0))
    fig = go.Figure()
    for i, st in enumerate(hydro_stations[:40]):
        fig.add_trace(go.Scatter(x=ts.index, y=ts[st], mode="lines", name=str(st), visible=(i == 0)))

    buttons = []
    for i, st in enumerate(hydro_stations[:40]):
        vis = [False] * min(len(hydro_stations), 40)
        vis[i] = True
        buttons.append(dict(label=str(st), method="update", args=[{"visible": vis}, {"title": f"Hydro station: {st}"}]))

    _apply_selector_layout(
        fig,
        title=f"Hydro station: {hydro_stations[0]}",
        selector_label="Select hydro station",
        buttons=buttons,
        xaxis_title="Snapshot",
        yaxis_title="Generation [MW]",
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_spill_interactive(n: pypsa.Network, output_html: Path) -> None:
    if n.links.empty or n.links_t.p0.empty:
        go.Figure(layout=dict(title="No spill-link data available")).write_html(output_html, include_plotlyjs="cdn")
        return

    spill_mask = n.links.index.astype(str).str.lower().str.contains("spill")
    if "carrier" in n.links.columns:
        spill_mask = spill_mask | n.links["carrier"].astype(str).str.lower().str.contains("spill")
    spill_links = list(n.links.index[spill_mask])
    if not spill_links:
        go.Figure(layout=dict(title="No spill links detected")).write_html(output_html, include_plotlyjs="cdn")
        return

    ts = _ensure_datetime_index(n.links_t.p0[spill_links].abs())
    fig = go.Figure()
    for i, lk in enumerate(spill_links[:40]):
        fig.add_trace(go.Scatter(x=ts.index, y=ts[lk], mode="lines", name=str(lk), visible=(i == 0)))

    buttons = []
    for i, lk in enumerate(spill_links[:40]):
        vis = [False] * min(len(spill_links), 40)
        vis[i] = True
        buttons.append(dict(label=str(lk), method="update", args=[{"visible": vis}, {"title": f"Spill link: {lk}"}]))

    _apply_selector_layout(
        fig,
        title=f"Spill link: {spill_links[0]}",
        selector_label="Select spill link",
        buttons=buttons,
        xaxis_title="Snapshot",
        yaxis_title="Spill-link flow [MW-equivalent]",
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_storage_interactive(n: pypsa.Network, output_html: Path) -> None:
    if n.stores.empty or n.stores_t.e.empty:
        go.Figure(layout=dict(title="No reservoir storage data available")).write_html(
            output_html, include_plotlyjs="cdn"
        )
        return

    capacity = pd.to_numeric(n.stores["e_nom"], errors="coerce").replace(0, np.nan)
    storage = _ensure_datetime_index(n.stores_t.e).div(capacity, axis=1)
    storage = storage.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how="all")
    stores = list(storage.columns[:50])
    if not stores:
        go.Figure(layout=dict(title="No finite reservoir state-of-charge data")).write_html(
            output_html, include_plotlyjs="cdn"
        )
        return

    fig = go.Figure()
    for i, store in enumerate(stores):
        fig.add_trace(
            go.Scatter(
                x=storage.index,
                y=storage[store],
                mode="lines",
                name=str(store),
                visible=(i == 0),
                line=dict(color="#155a78", width=1.6),
            )
        )
    buttons = []
    for i, store in enumerate(stores):
        visible = [False] * len(stores)
        visible[i] = True
        buttons.append(
            dict(
                label=str(store),
                method="update",
                args=[{"visible": visible}, {"title": f"Reservoir state of charge · {store}"}],
            )
        )
    _apply_selector_layout(
        fig,
        title=f"Reservoir state of charge · {stores[0]}",
        selector_label="Select reservoir",
        buttons=buttons,
        xaxis_title="Snapshot",
        extra_layout={
            "yaxis": dict(
                title="Stored energy / nominal capacity [p.u.]",
                range=[0, 1.05],
            )
        },
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_load_interactive_plot(n: pypsa.Network, output_html: Path) -> None:
    prepared_path = Path("data/processed_data/load/Hourly_profile_2021.csv")
    if prepared_path.exists():
        prepared = pd.read_csv(prepared_path)
        prepared["TIME"] = pd.to_datetime(prepared["TIME"], errors="coerce")
        raw = pd.to_numeric(prepared["LOAD"], errors="coerce")
        load_gw = pd.Series(raw.values / 1_000_000.0, index=prepared["TIME"], name="Load")
        source_label = "Prepared provincial load input"
    elif not n.loads_t.p.empty:
        load_gw = _ensure_datetime_index(n.loads_t.p).sum(axis=1) / 1_000.0
        source_label = "Solved-network load"
    else:
        go.Figure(layout=dict(title="No load data available")).write_html(output_html, include_plotlyjs="cdn")
        return

    load_gw = load_gw.dropna()
    duration = pd.Series(load_gw.values).sort_values(ascending=False).reset_index(drop=True)
    share = np.linspace(0, 100, len(duration))
    peak = float(load_gw.max())
    peak_time = load_gw.idxmax()
    mean = float(load_gw.mean())
    energy_twh = float(load_gw.sum() / 1_000.0)
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Hourly provincial load", "Load duration curve"),
        horizontal_spacing=0.1,
    )
    fig.add_trace(
        go.Scatter(
            x=load_gw.index,
            y=load_gw.values,
            mode="lines",
            name="Load",
            line=dict(color="#155a78", width=1.25),
            hovertemplate="%{x|%d %b %Y, %H:%M}<br>%{y:.2f} GW<extra></extra>",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=share,
            y=duration.values,
            mode="lines",
            name="Duration",
            line=dict(color="#b8860b", width=1.6),
            hovertemplate="Exceeded %{x:.1f}% of hours<br>%{y:.2f} GW<extra></extra>",
        ),
        row=1,
        col=2,
    )
    fig.add_annotation(
        x=peak_time,
        y=peak,
        text=f"Peak {peak:.2f} GW",
        showarrow=True,
        arrowcolor="#26332f",
        bgcolor="#fffdf7",
        bordercolor="#8a9691",
        row=1,
        col=1,
    )
    fig.update_xaxes(title_text="Time", row=1, col=1)
    fig.update_xaxes(title_text="Hours exceeded [%]", row=1, col=2)
    fig.update_yaxes(title_text="Load [GW]", rangemode="tozero", row=1, col=1)
    fig.update_yaxes(title_text="Load [GW]", rangemode="tozero", row=1, col=2)
    fig.update_layout(
        title=dict(
            text=(
                f"Provincial load profile · 2021"
                f"<br><sup>{source_label} · peak {peak:.2f} GW · mean {mean:.2f} GW "
                f"· annual energy {energy_twh:.1f} TWh</sup>"
            ),
            x=0.02,
        ),
        template="plotly_white",
        font=dict(family="Arial, sans-serif", size=13, color="#26332f"),
        showlegend=False,
        height=540,
        margin=dict(l=72, r=35, t=95, b=62),
        paper_bgcolor="#fbfaf6",
        plot_bgcolor="#fbfaf6",
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_price_interactive_plot(n: pypsa.Network, output_html: Path) -> None:
    if n.buses_t.marginal_price.empty:
        go.Figure(layout=dict(title="No marginal price data available")).write_html(output_html, include_plotlyjs="cdn")
        return
    prices = _ensure_datetime_index(n.buses_t.marginal_price).mean(axis=1)
    p_sorted = prices.sort_values(ascending=False).reset_index(drop=True)
    p_share = np.linspace(0, 100, len(p_sorted))

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Price duration curve", "Price over time"))
    fig.add_trace(go.Scatter(x=p_share, y=p_sorted.values, mode="lines", name="Duration"), row=1, col=1)
    fig.add_trace(go.Scatter(x=prices.index, y=prices.values, mode="lines", name="Price"), row=1, col=2)
    fig.update_layout(template="plotly_white", font=dict(family=REPORT_FONT_FAMILY, size=13, color=REPORT_INK), showlegend=False)
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_curtailment_interactive_plot(n: pypsa.Network, output_html: Path) -> None:
    if n.generators.empty or n.generators_t.p.empty or n.generators_t.p_max_pu.empty:
        go.Figure(layout=dict(title="No curtailment data available")).write_html(output_html, include_plotlyjs="cdn")
        return

    dispatch = n.generators_t.p
    p_max_pu = n.generators_t.p_max_pu.reindex(columns=dispatch.columns)
    available = p_max_pu.multiply(n.generators["p_nom"], axis=1)
    curtailed = (available - dispatch).clip(lower=0)
    tech = _get_generator_tech_labels(n)
    rel = (curtailed.sum(axis=0).groupby(tech).sum() / available.sum(axis=0).groupby(tech).sum().replace(0, np.nan) * 100).dropna()
    rel = rel.sort_values(ascending=False)

    fig = go.Figure(go.Bar(x=rel.index.astype(str), y=rel.values))
    fig.update_layout(title="Curtailment by technology", template="plotly_white", font=dict(family=REPORT_FONT_FAMILY, size=13, color=REPORT_INK))
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_static_network_map(
    n: pypsa.Network,
    output_png: Path,
    title: str,
    *,
    full_bc_frame: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 8.5))

    bc_regions = _load_bc_gadm_regions()
    x, y, valid_buses = _get_bc_filtered_bus_coordinates(n, bc_regions)
    has_coordinates = len(valid_buses) > 0

    if bc_regions is not None and len(bc_regions) > 0:
        bc_regions.to_crs(BC_ALBERS).plot(
            ax=ax,
            color="#ecebe6",
            edgecolor="#9da4a0",
            linewidth=0.5,
            alpha=0.92,
            zorder=0,
        )

    points_albers = None
    coord_map: dict[str, tuple[float, float]] = {}
    if has_coordinates:
        points_albers = gpd.GeoDataFrame(
            {"bus": valid_buses},
            geometry=gpd.points_from_xy(x.loc[valid_buses], y.loc[valid_buses]),
            crs="EPSG:4326",
        ).to_crs(BC_ALBERS).set_index("bus")
        coord_map = {
            str(bus): (float(row.geometry.x), float(row.geometry.y))
            for bus, row in points_albers.iterrows()
        }

    if has_coordinates and not n.lines.empty:
        for _, line in n.lines.iterrows():
            bus0 = line.get("bus0")
            bus1 = line.get("bus1")
            if str(bus0) in coord_map and str(bus1) in coord_map:
                p0, p1 = coord_map[str(bus0)], coord_map[str(bus1)]
                ax.plot(
                    [p0[0], p1[0]],
                    [p0[1], p1[1]],
                    color="#707070",
                    linewidth=1.0,
                    alpha=0.8,
                    zorder=1,
                )

    if has_coordinates:
        ax.scatter(
            points_albers.geometry.x,
            points_albers.geometry.y,
            s=24,
            facecolor="#f8f7f2",
            edgecolor="#4b5563",
            linewidth=0.65,
            alpha=0.9,
            zorder=2,
            label="modeled bus",
        )

        # Use the same capacity/resource visual grammar as the detailed input
        # map, but only for generators connected to geographically valid buses.
        if not n.generators.empty:
            gen = n.generators.loc[n.generators["bus"].isin(valid_buses)].copy()
            gen["resource"] = _get_generator_tech_labels(n).reindex(gen.index).map(
                _normalize_resource_name
            )
            gen["capacity_mw"] = pd.to_numeric(gen["p_nom"], errors="coerce").fillna(0)
            # Backstop is a synthetic feasibility mechanism, not a physical
            # generation asset; retain it in diagnostics but omit it here so it
            # cannot visually overwhelm the physical technology markers.
            gen = gen.loc[gen["resource"].ne("backstop")]
            mix = (
                gen.groupby(["bus", "resource"])["capacity_mw"]
                .sum()
                .reset_index()
            )
            cmap = _resource_color_map()
            for resource, group in mix.groupby("resource"):
                marker_x = [coord_map[str(bus)][0] for bus in group["bus"]]
                marker_y = [coord_map[str(bus)][1] for bus in group["bus"]]
                marker_size = 18 + 5.5 * np.sqrt(
                    group["capacity_mw"].clip(lower=0, upper=5000)
                )
                ax.scatter(
                    marker_x,
                    marker_y,
                    s=marker_size,
                    c=cmap.get(str(resource), "#6c757d"),
                    edgecolors="white",
                    linewidths=0.5,
                    alpha=0.9,
                    zorder=3,
                    label=str(resource),
                )

        for bus_name in valid_buses:
            px_, py_ = coord_map[str(bus_name)]
            ax.text(px_, py_, str(bus_name), fontsize=6, alpha=0.78)
    else:
        ax.text(
            0.5,
            0.5,
            "No bus coordinates available (x/y).",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
        ax.set_axis_off()

    if full_bc_frame:
        _set_bc_map_frame(ax, bc_regions, pad_fraction=0.045)
    elif points_albers is not None:
        _set_geometry_map_frame(ax, points_albers.geometry)
    else:
        _set_bc_map_frame(ax, bc_regions)
    _add_grey_basemap(ax)
    ax.set_title("Network Topology", loc="left", fontsize=14, fontweight="bold")
    _clean_axis(ax, hide_axes=True)
    ax.legend(loc="upper left", ncol=3, fontsize=8, frameon=True, facecolor="#f7f8f3")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_interactive_network_map(
    n: pypsa.Network,
    output_html: Path,
    title: str,
) -> None:
    bc_regions = _load_bc_gadm_regions()
    x, y, valid_buses = _get_bc_filtered_bus_coordinates(n, bc_regions)
    has_coordinates = len(valid_buses) > 0

    if not has_coordinates:
        fig = go.Figure()
        fig.add_annotation(text="No bus coordinates available for map plotting.", showarrow=False)
        fig.update_layout(title=title)
        fig.write_html(output_html, include_plotlyjs="cdn")
        return

    line_metrics = pd.DataFrame(index=n.lines.index)
    if not n.lines_t.p0.empty:
        s_nom = n.lines["s_nom"].replace(0, np.nan)
        loading = n.lines_t.p0.abs().div(s_nom, axis=1).replace([np.inf, -np.inf], np.nan)
        line_metrics["max_loading_pu"] = loading.max(axis=0)
        line_metrics["mean_loading_pu"] = loading.mean(axis=0)
        line_metrics["congestion_share"] = (loading >= 0.85).sum(axis=0) / max(len(loading.index), 1)
    else:
        line_metrics["max_loading_pu"] = np.nan
        line_metrics["mean_loading_pu"] = np.nan
        line_metrics["congestion_share"] = np.nan

    mean_price = pd.Series(dtype=float)
    if hasattr(n, "buses_t") and hasattr(n.buses_t, "marginal_price") and not n.buses_t.marginal_price.empty:
        mean_price = n.buses_t.marginal_price.mean(axis=0)

    hover_text = []
    for bus in valid_buses:
        bus_carrier = n.buses.loc[bus, "carrier"] if "carrier" in n.buses.columns else "n/a"
        msg = f"Bus: {bus}<br>Carrier: {bus_carrier}"
        if bus in mean_price.index:
            msg += f"<br>Mean marginal price: {mean_price.loc[bus]:.2f}"
        hover_text.append(msg)

    fig = go.Figure()
    if bc_regions is not None and len(bc_regions) > 0:
        bc_regions = bc_regions.to_crs("EPSG:4326")
        for _, row in bc_regions.iterrows():
            for ring in _iter_polygon_rings(row.geometry):
                fig.add_trace(
                    go.Scatter(
                        x=ring[:, 0],
                        y=ring[:, 1],
                        mode="lines",
                        fill="toself",
                        line=dict(color="#8f8f8f", width=0.7),
                        fillcolor="rgba(180,180,180,0.45)",
                        hoverinfo="skip",
                        showlegend=False,
                    )
                )

    for line_name, line in n.lines.iterrows():
        bus0 = line.get("bus0")
        bus1 = line.get("bus1")
        if bus0 not in valid_buses or bus1 not in valid_buses:
            continue

        cap = float(line.get("s_nom", np.nan))
        max_load = line_metrics.at[line_name, "max_loading_pu"] if line_name in line_metrics.index else np.nan
        mean_load = line_metrics.at[line_name, "mean_loading_pu"] if line_name in line_metrics.index else np.nan
        con_share = line_metrics.at[line_name, "congestion_share"] if line_name in line_metrics.index else np.nan

        hover = (
            f"Line: {line_name}<br>"
            f"Capacity s_nom: {cap:.1f} MW<br>"
            f"Max loading: {max_load:.2f} p.u.<br>"
            f"Mean loading: {mean_load:.2f} p.u.<br>"
            f"Congestion share (>=0.85): {100.0 * con_share:.1f}%"
        )

        color = "rgba(80,80,80,0.85)"
        if pd.notna(con_share):
            if con_share > 0.5:
                color = "rgba(179,0,0,0.95)"
            elif con_share > 0.2:
                color = "rgba(239,138,98,0.9)"

        fig.add_trace(
            go.Scatter(
                x=[float(x.loc[bus0]), float(x.loc[bus1])],
                y=[float(y.loc[bus0]), float(y.loc[bus1])],
                mode="lines",
                line=dict(width=1.8, color=color),
                hovertemplate=hover + "<extra></extra>",
                showlegend=False,
            )
        )

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(size=9, color="#1f78b4", line=dict(width=0.5, color="white")),
            text=hover_text,
            hoverinfo="text",
            name="Buses",
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=None,
        yaxis_title=None,
        template="plotly_white",
        font=dict(family=REPORT_FONT_FAMILY, size=13, color=REPORT_INK),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False, scaleanchor="x", scaleratio=1)
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_dispatch_plots(
    n: pypsa.Network,
    output_static_png: Path,
    output_interactive_html: Path,
) -> None:
    if n.generators_t.p.empty:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "No generator dispatch time series available.", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(output_static_png, dpi=200)
        plt.close(fig)
        go.Figure(layout=dict(title="No generator dispatch time series available.")).write_html(
            output_interactive_html, include_plotlyjs="cdn"
        )
        return

    dispatch = _ensure_datetime_index(n.generators_t.p.copy())
    tech = _get_generator_tech_labels(n)
    grouped = dispatch.T.groupby(tech).sum().T

    energy_totals = grouped.abs().sum(axis=0).sort_values(ascending=False)
    top_n = list(energy_totals.head(8).index)
    grouped_top = grouped[top_n].copy()
    other_cols = [c for c in grouped.columns if c not in top_n]
    if other_cols:
        grouped_top["other"] = grouped[other_cols].sum(axis=1)

    grouped_positive = grouped_top.clip(lower=0)

    if isinstance(grouped_positive.index, pd.DatetimeIndex):
        monthly = grouped_positive.resample("MS").sum() / 1e3
        x_values = monthly.index
    else:
        monthly = grouped_positive
        x_values = np.arange(len(monthly.index))

    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.stackplot(x_values, [monthly[c].values for c in monthly.columns], labels=monthly.columns)
    ax.set_title("Monthly Generation Dispatch by Technology")
    ax.set_ylabel("Energy [GWh]")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(output_static_png, dpi=220)
    plt.close(fig)

    if isinstance(grouped_positive.index, pd.DatetimeIndex):
        interactive_df = grouped_positive.copy()
        if len(interactive_df) > 24 * 31:
            interactive_df = interactive_df.iloc[: 24 * 31]
        fig_i = px.area(
            interactive_df,
            x=interactive_df.index,
            y=interactive_df.columns,
            title="Hourly Dispatch (first 31 days)",
            labels={"value": "MW", "index": "snapshot", "variable": "technology"},
            template="plotly_white",
        )
    else:
        fig_i = px.area(
            grouped_positive,
            x=grouped_positive.index,
            y=grouped_positive.columns,
            title="Dispatch",
            labels={"value": "MW", "variable": "technology"},
            template="plotly_white",
        )

    fig_i.write_html(output_interactive_html, include_plotlyjs="cdn")


def _write_load_and_price_plots(
    n: pypsa.Network,
    load_png: Path,
    price_duration_png: Path,
    price_heatmap_png: Path,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)

    if not n.loads_t.p.empty:
        load_ts = _ensure_datetime_index(n.loads_t.p).sum(axis=1)
        axes[0].plot(load_ts.index, load_ts.values, color="#023047", linewidth=1.0)
        axes[0].set_title("System Load Time Series")
        axes[0].set_ylabel("MW")

        duration = pd.Series(load_ts.values).sort_values(ascending=False).reset_index(drop=True)
        share = np.linspace(0, 100, len(duration))
        axes[1].plot(share, duration.values, color="#fb8500", linewidth=1.4)
        axes[1].set_title("Load Duration Curve")
        axes[1].set_xlabel("Share of hours [%]")
        axes[1].set_ylabel("MW")
    else:
        axes[0].text(0.5, 0.5, "No load time series available.", ha="center", va="center")
        axes[0].set_axis_off()
        axes[1].set_axis_off()

    fig.tight_layout()
    fig.savefig(load_png, dpi=220)
    plt.close(fig)

    if not n.buses_t.marginal_price.empty:
        prices = _ensure_datetime_index(n.buses_t.marginal_price).mean(axis=1)
        p_sorted = prices.sort_values(ascending=False).reset_index(drop=True)
        p_share = np.linspace(0, 100, len(p_sorted))

        fig_pdc, ax_pdc = plt.subplots(figsize=(11, 4.5))
        ax_pdc.plot(p_share, p_sorted.values, color="#7b2cbf", linewidth=1.5)
        ax_pdc.set_title("Marginal Price Duration Curve")
        ax_pdc.set_xlabel("Share of hours [%]")
        ax_pdc.set_ylabel("Price")
        fig_pdc.tight_layout()
        fig_pdc.savefig(price_duration_png, dpi=220)
        plt.close(fig_pdc)

        if isinstance(prices.index, pd.DatetimeIndex):
            price_frame = pd.DataFrame({"price": prices})
            price_frame["day"] = price_frame.index.dayofyear
            price_frame["hour"] = price_frame.index.hour
            pivot = price_frame.pivot_table(index="hour", columns="day", values="price", aggfunc="mean")

            fig_hm, ax_hm = plt.subplots(figsize=(12, 4.5))
            heat = ax_hm.imshow(pivot.values, aspect="auto", cmap="Spectral_r", origin="lower")
            ax_hm.set_title("Marginal Price Heatmap (hour x day-of-year)")
            ax_hm.set_xlabel("Day of year")
            ax_hm.set_ylabel("Hour of day")
            plt.colorbar(heat, ax=ax_hm, fraction=0.03, pad=0.02, label="Price")
            fig_hm.tight_layout()
            fig_hm.savefig(price_heatmap_png, dpi=220)
            plt.close(fig_hm)
        else:
            fig_hm, ax_hm = plt.subplots(figsize=(9, 4))
            ax_hm.text(0.5, 0.5, "Marginal price index is not datetime; skipping heatmap.", ha="center", va="center")
            ax_hm.set_axis_off()
            fig_hm.tight_layout()
            fig_hm.savefig(price_heatmap_png, dpi=220)
            plt.close(fig_hm)
    else:
        fig_pdc, ax_pdc = plt.subplots(figsize=(9, 4))
        ax_pdc.text(0.5, 0.5, "No marginal price time series available.", ha="center", va="center")
        ax_pdc.set_axis_off()
        fig_pdc.tight_layout()
        fig_pdc.savefig(price_duration_png, dpi=220)
        plt.close(fig_pdc)

        fig_hm, ax_hm = plt.subplots(figsize=(9, 4))
        ax_hm.text(0.5, 0.5, "No marginal price time series available.", ha="center", va="center")
        ax_hm.set_axis_off()
        fig_hm.tight_layout()
        fig_hm.savefig(price_heatmap_png, dpi=220)
        plt.close(fig_hm)


def _write_line_loading_plot(n: pypsa.Network, output_png: Path) -> None:
    if n.lines.empty or n.lines_t.p0.empty:
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.text(0.5, 0.5, "No line flow data available.", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    s_nom = n.lines["s_nom"].replace(0, np.nan)
    loading = n.lines_t.p0.abs().div(s_nom, axis=1).replace([np.inf, -np.inf], np.nan)
    max_loading = loading.max(axis=0).sort_values(ascending=False)
    top_lines = list(max_loading.head(8).index)

    fig, ax = plt.subplots(figsize=(11, 5))
    for line_name in top_lines:
        curve = loading[line_name].dropna().sort_values(ascending=False).reset_index(drop=True)
        if len(curve) == 0:
            continue
        share = np.linspace(0, 100, len(curve))
        ax.plot(share, curve.values, linewidth=1.2, label=str(line_name))

    ax.axhline(0.85, color="black", linestyle="--", linewidth=1, alpha=0.8)
    ax.set_title("Top Line Loading Duration Curves")
    ax.set_xlabel("Share of hours [%]")
    ax.set_ylabel("Loading [p.u.]")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_line_loading_interactive(n: pypsa.Network, output_html: Path) -> None:
    if n.lines.empty or n.lines_t.p0.empty:
        go.Figure(layout=dict(title="No line flow data available")).write_html(output_html, include_plotlyjs="cdn")
        return

    s_nom = n.lines["s_nom"].replace(0, np.nan)
    loading = n.lines_t.p0.abs().div(s_nom, axis=1).replace([np.inf, -np.inf], np.nan)
    max_loading = loading.max(axis=0).sort_values(ascending=False)
    top_lines = list(max_loading.head(10).index)

    fig = go.Figure()
    for line_name in top_lines:
        curve = loading[line_name].dropna().sort_values(ascending=False).reset_index(drop=True)
        if len(curve) == 0:
            continue
        share = np.linspace(0, 100, len(curve))
        fig.add_trace(go.Scatter(x=share, y=curve.values, mode="lines", name=str(line_name)))

    fig.add_hline(y=0.85, line_dash="dash", line_color="#303030", annotation_text="0.85 p.u.")
    fig.update_layout(
        title="Top line congestion duration curves",
        xaxis_title="Share of hours [%]",
        yaxis_title="Loading [p.u.]",
        template="plotly_white",
        font=dict(family=REPORT_FONT_FAMILY, size=13, color=REPORT_INK),
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def _write_static_capacity_congestion_map(n: pypsa.Network, output_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 8))

    bc_regions = _load_bc_gadm_regions()
    x, y, valid_buses = _get_bc_filtered_bus_coordinates(n, bc_regions)
    has_coordinates = len(valid_buses) > 0

    regions_4326 = None
    if bc_regions is not None and len(bc_regions) > 0:
        regions_4326 = bc_regions.to_crs("EPSG:4326")
        regions_4326.to_crs(BC_ALBERS).plot(
            ax=ax,
            color="#d9d9d9",
            edgecolor="#8f8f8f",
            linewidth=0.7,
            alpha=0.98,
            zorder=0,
        )

    if n.lines.empty or n.lines_t.p0.empty or not has_coordinates:
        ax.text(
            0.5,
            0.5,
            "Missing line-flow data or bus coordinates for static congestion map.",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    s_nom = n.lines["s_nom"].replace(0, np.nan)
    loading = n.lines_t.p0.abs().div(s_nom, axis=1).replace([np.inf, -np.inf], np.nan)
    congestion_share = (loading >= CONGESTION_THRESHOLD).sum(axis=0) / max(len(loading.index), 1)

    cap = s_nom.fillna(0)
    cap_max = cap.max() if len(cap) else 0
    if cap_max <= 0:
        cap_max = 1

    # Width scales with line capacity.
    widths = 0.5 + 5.5 * (cap / cap_max)

    def color_for_share(share: float) -> str:
        if share > 0.50:
            return "#c65d2e"  # heavily congested (>50% of year)
        if share > 0.20:
            return "#d9a441"  # moderately congested
        return "#315f78"  # lightly congested

    for line_name, line in n.lines.iterrows():
        bus0 = line.get("bus0")
        bus1 = line.get("bus1")
        if bus0 not in valid_buses or bus1 not in valid_buses:
            continue

        share = float(congestion_share.get(line_name, 0.0))
        color = color_for_share(share)
        width = float(widths.get(line_name, 0.8))

        line_geom = gpd.GeoSeries(
            [LineString([(x.loc[bus0], y.loc[bus0]), (x.loc[bus1], y.loc[bus1])])],
            crs="EPSG:4326",
        ).to_crs(BC_ALBERS).iloc[0]
        xy = np.asarray(line_geom.coords)
        ax.plot(
            xy[:, 0],
            xy[:, 1],
            color=color,
            linewidth=width,
            alpha=0.85,
            zorder=1,
        )

    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x, y), crs="EPSG:4326").to_crs(BC_ALBERS)
    ax.scatter(pts.geometry.x, pts.geometry.y, s=18, c="#1e495f", edgecolor="white", linewidth=0.4, zorder=2)
    _add_grey_basemap(ax)

    if regions_4326 is not None and "NAME_2" in regions_4326.columns:
        label_pts = regions_4326.to_crs(BC_ALBERS)
        for _, row in label_pts.iterrows():
            c = row.geometry.centroid
            ax.text(c.x, c.y, str(row["NAME_2"]), fontsize=6, color="#5a5a5a", alpha=0.9)

    ax.set_title("BC Line Capacity and Congestion", loc="left", fontsize=14, fontweight="bold")
    _clean_axis(ax, hide_axes=True)

    congestion_legend = [
        Line2D([0], [0], color="#c65d2e", lw=3, label="≥85% loaded for >50% of hours"),
        Line2D([0], [0], color="#d9a441", lw=3, label="≥85% loaded for 20–50% of hours"),
        Line2D([0], [0], color="#315f78", lw=3, label="≥85% loaded for <20% of hours"),
    ]

    cap_levels = np.array([0.25, 0.5, 1.0]) * cap_max
    cap_widths = 0.5 + 5.5 * (cap_levels / cap_max)
    capacity_legend = [
        Line2D([0], [0], color="#4d4d4d", lw=float(w), label=f"Capacity ~{int(c):,} MW")
        for c, w in zip(cap_levels, cap_widths)
    ]

    leg1 = ax.legend(handles=congestion_legend, loc="upper left", frameon=False)
    ax.add_artist(leg1)
    ax.legend(handles=capacity_legend, loc="upper right", frameon=False)

    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_storage_plot(n: pypsa.Network, output_png: Path) -> None:
    if n.stores.empty or n.stores_t.e.empty:
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.text(0.5, 0.5, "No store state-of-charge data available.", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    soc = _ensure_datetime_index(n.stores_t.e.copy())
    cap = n.stores["e_nom"].replace(0, np.nan)
    normalized = soc.div(cap, axis=1).replace([np.inf, -np.inf], np.nan)

    top = cap.sort_values(ascending=False).head(8).index
    normalized = normalized[top].dropna(axis=1, how="all")

    fig, ax = plt.subplots(figsize=(12, 5))
    for col in normalized.columns:
        ax.plot(normalized.index, normalized[col], linewidth=1.0, label=str(col))

    ax.set_title("Reservoir/Store State of Charge (normalized by e_nom)")
    ax.set_ylabel("State of charge [0-1]")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_curtailment_plot(n: pypsa.Network, output_png: Path) -> None:
    if n.generators.empty or n.generators_t.p.empty or n.generators_t.p_max_pu.empty:
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.text(0.5, 0.5, "No curtailment inputs available (need generators_t.p and p_max_pu).", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    dispatch = n.generators_t.p
    p_max_pu = n.generators_t.p_max_pu.reindex(columns=dispatch.columns)
    available = p_max_pu.multiply(n.generators["p_nom"], axis=1)
    available = available.where(available > 0)
    curtailed = (available - dispatch).clip(lower=0)

    tech = _get_generator_tech_labels(n)
    available_by_tech = available.sum(axis=0).groupby(tech).sum()
    curtailed_by_tech = curtailed.sum(axis=0).groupby(tech).sum()
    rel = (curtailed_by_tech / available_by_tech.replace(0, np.nan) * 100).dropna()
    rel = rel.sort_values(ascending=False).head(12)

    fig, ax = plt.subplots(figsize=(11, 4.5))
    if rel.empty:
        ax.text(0.5, 0.5, "Curtailment could not be computed from available data.", ha="center", va="center")
        ax.set_axis_off()
    else:
        ax.bar(rel.index.astype(str), rel.values, color="#ffb703")
        ax.set_title("Estimated Curtailment by Technology")
        ax.set_ylabel("Curtailment [% of available]")
        ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_resource_generation_plot(n: pypsa.Network, output_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 4.8))
    if n.generators.empty or n.generators_t.p.empty:
        ax.text(0.5, 0.5, "No generator data available.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    tech = _get_generator_tech_labels(n)
    gen = n.generators_t.p.clip(lower=0)
    energy = gen.sum(axis=0).groupby(tech).sum().sort_values(ascending=False)
    energy_gwh = energy / 1e3

    ax.bar(energy_gwh.index.astype(str), energy_gwh.values, color="#2b6cb0")
    ax.set_title("Resource-wise annual generation")
    ax.set_ylabel("Energy [GWh]")
    ax.tick_params(axis="x", rotation=35)
    _clean_axis(ax)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_station_generation_plots(
    n: pypsa.Network,
    output_static_png: Path,
    output_interactive_html: Path,
) -> None:
    if n.generators.empty or n.generators_t.p.empty:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "No station generation data available.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_static_png, dpi=220)
        plt.close(fig)
        go.Figure(layout=dict(title="No station generation data available.")).write_html(
            output_interactive_html, include_plotlyjs="cdn"
        )
        return

    gen = _ensure_datetime_index(n.generators_t.p.clip(lower=0))
    annual = gen.sum(axis=0).sort_values(ascending=False)
    top_stations = list(annual.head(15).index)

    fig, ax = plt.subplots(figsize=(11, 5))
    vals = (annual.loc[top_stations] / 1e3)
    ax.bar(vals.index.astype(str), vals.values, color="#0b5d4b")
    ax.set_title("Top stations by annual generation")
    ax.set_ylabel("Energy [GWh]")
    ax.tick_params(axis="x", rotation=60)
    _clean_axis(ax)
    fig.tight_layout()
    fig.savefig(output_static_png, dpi=220)
    plt.close(fig)

    station_ts = gen[top_stations].copy()
    if isinstance(station_ts.index, pd.DatetimeIndex) and len(station_ts) > 24 * 14:
        station_ts = station_ts.iloc[: 24 * 14]

    fig_i = go.Figure()
    for i, station in enumerate(top_stations):
        fig_i.add_trace(
            go.Scatter(
                x=station_ts.index,
                y=station_ts[station],
                mode="lines",
                name=str(station),
                visible=(i == 0),
            )
        )

    buttons = []
    for i, station in enumerate(top_stations):
        visible = [False] * len(top_stations)
        visible[i] = True
        buttons.append(
            dict(
                label=str(station),
                method="update",
                args=[{"visible": visible}, {"title": f"Station generation: {station}"}],
            )
        )

    _apply_selector_layout(
        fig_i,
        title=f"Station generation: {top_stations[0]}",
        selector_label="Select station",
        buttons=buttons,
        xaxis_title="Snapshot",
        yaxis_title="Generation [MW]",
    )
    fig_i.write_html(output_interactive_html, include_plotlyjs="cdn")


def _write_hydro_generation_plot(n: pypsa.Network, output_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    if n.generators.empty or n.generators_t.p.empty:
        ax.text(0.5, 0.5, "No generator data available.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    tech = _get_generator_tech_labels(n).astype(str).str.lower()
    hydro_mask = tech.str.contains("hydro|ror|reservoir") | n.generators.index.astype(str).str.lower().str.contains("hydro|ror|reservoir")
    hydro_stations = n.generators.index[hydro_mask]

    if len(hydro_stations) == 0:
        ax.text(0.5, 0.5, "No hydro stations matched by tech/name pattern.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    hydro = _ensure_datetime_index(n.generators_t.p[hydro_stations].clip(lower=0))
    if isinstance(hydro.index, pd.DatetimeIndex):
        monthly = hydro.sum(axis=1).resample("MS").sum() / 1e3
        ax.plot(monthly.index, monthly.values, linewidth=1.8, color="#1b8a79")
    else:
        total = hydro.sum(axis=1) / 1e3
        ax.plot(total.index, total.values, linewidth=1.8, color="#1b8a79")

    ax.set_title("Hydro generation (explicit)")
    ax.set_ylabel("Energy [GWh]")
    _clean_axis(ax)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _write_reservoir_spill_plot(n: pypsa.Network, output_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    if n.links.empty or n.links_t.p0.empty:
        ax.text(0.5, 0.5, "No link-flow data available for spill analysis.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    link_names = n.links.index.astype(str)
    spill_name_mask = link_names.str.lower().str.contains("spill")
    spill_carrier_mask = pd.Series(False, index=n.links.index)
    if "carrier" in n.links.columns:
        spill_carrier_mask = n.links["carrier"].astype(str).str.lower().str.contains("spill")

    spill_links = n.links.index[spill_name_mask | spill_carrier_mask.values]
    if len(spill_links) == 0:
        ax.text(0.5, 0.5, "No spill links found by name/carrier pattern.", ha="center", va="center")
        _clean_axis(ax, hide_axes=True)
        fig.tight_layout()
        fig.savefig(output_png, dpi=220)
        plt.close(fig)
        return

    spill_ts = _ensure_datetime_index(n.links_t.p0[spill_links].abs().sum(axis=1))
    if isinstance(spill_ts.index, pd.DatetimeIndex):
        monthly = spill_ts.resample("MS").sum() / 1e3
        ax.plot(monthly.index, monthly.values, linewidth=1.8, color="#a13e36")
    else:
        ax.plot(spill_ts.index, spill_ts.values / 1e3, linewidth=1.8, color="#a13e36")

    ax.set_title("Reservoir spill (aggregated spill-link flow)")
    ax.set_ylabel("Energy-equivalent [GWh]")
    _clean_axis(ax)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def _build_input_inventory_table() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    patterns = [
        ("Processed", "Network buses", "data/processed_data/network/buses.csv", None),
        ("Processed", "Network lines", "data/processed_data/network/lines.csv", "s_nom"),
        ("Processed", "Existing wind", "data/processed_data/wind/existing/bc_ext_wind_assets.csv", "potential_capacity"),
        ("Processed", "Committed wind", "data/processed_data/wind/committed/BCH_CFP24_wind.csv", "potential_capacity"),
        ("Processed", "Potential wind", "data/processed_data/wind/potential/resource_options_wind.csv", "potential_capacity"),
        ("Processed", "Existing solar", "data/processed_data/solar/existing/bc_ext_solar_assets.csv", "potential_capacity"),
        ("Processed", "Committed solar", "data/processed_data/solar/committed/BCH_CFP24_solar.csv", "potential_capacity"),
        ("Processed", "Potential solar", "data/processed_data/solar/potential/resource_options_solar.csv", "potential_capacity"),
        ("Processed", "Thermal assets", "data/processed_data/tpp/existing/bc_ext_tpp_assets.csv", "potential_capacity"),
        ("Processed", "Hydro generation", "data/processed_data/hydro/existing/hydro_generation.csv", "p_nom"),
    ]
    for layer, dataset, raw_path, capacity_column in patterns:
        path = Path(raw_path)
        if not path.exists():
            continue
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        capacity = ""
        if capacity_column and capacity_column in frame.columns:
            capacity = f"{pd.to_numeric(frame[capacity_column], errors='coerce').sum():,.1f}"
        rows.append(
            {
                "Layer": layer,
                "Dataset": dataset,
                "Records": f"{len(frame):,}",
                "Capacity [MW]": capacity,
            }
        )

    prepared = _load_prepared_generation_components()
    if not prepared.empty:
        for resource, group in prepared.groupby("resource"):
            rows.append(
                {
                    "Layer": "PyPSA data",
                    "Dataset": f"{resource} components",
                    "Records": f"{len(group):,}",
                    "Capacity [MW]": f"{group['capacity_mw'].sum():,.1f}",
                }
            )
    return pd.DataFrame(rows)


def _build_summary_table(n: pypsa.Network) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    def add(section: str, metric: str, value, status: str = "Info") -> None:
        rows.append({"Section": section, "Metric": metric, "Value": value, "Status": status})

    add("Scope", "Snapshots", f"{len(n.snapshots):,}")
    add("Scope", "Buses", f"{len(n.buses):,}")
    add("Scope", "Lines", f"{len(n.lines):,}")
    add("Scope", "Links", f"{len(n.links):,}")
    add("Scope", "Generators", f"{len(n.generators):,}")
    add("Scope", "Stores", f"{len(n.stores):,}")

    if {"x", "y"}.issubset(n.buses.columns):
        missing_xy = int(
            pd.to_numeric(n.buses["x"], errors="coerce").isna().sum()
            + pd.to_numeric(n.buses["y"], errors="coerce").isna().sum()
        )
        add("Input QA", "Missing bus coordinate fields", missing_xy, "Check" if missing_xy else "Pass")
    if "s_nom" in n.lines.columns:
        missing_cap = int((pd.to_numeric(n.lines["s_nom"], errors="coerce").fillna(0) <= 0).sum())
        add("Input QA", "Lines with missing/zero capacity", missing_cap, "Check" if missing_cap else "Pass")

    if not n.loads_t.p.empty:
        total_load = float(n.loads_t.p.sum().sum())
        peak_load = float(n.loads_t.p.sum(axis=1).max())
        add("Demand", "Annual load [TWh]", f"{total_load / 1e6:,.2f}")
        add("Demand", "Peak system load [MW]", f"{peak_load:,.1f}")

    if not n.generators_t.p.empty:
        generation = float(n.generators_t.p.clip(lower=0).sum().sum())
        add("Supply", "Generated energy [TWh]", f"{generation / 1e6:,.2f}")
        backstop_mask = n.generators.index.astype(str).str.lower().str.contains("backstop")
        if "tech" in n.generators.columns:
            backstop_mask |= n.generators["tech"].astype(str).str.lower().str.contains("backstop")
        backstop_cols = n.generators.index[backstop_mask].intersection(n.generators_t.p.columns)
        backstop = float(n.generators_t.p[backstop_cols].clip(lower=0).sum().sum()) if len(backstop_cols) else 0.0
        add("Adequacy", "Backstop dispatch [GWh]", f"{backstop / 1e3:,.2f}", "Check" if backstop > 1e-6 else "Pass")

    if not n.lines.empty and not n.lines_t.p0.empty and "s_nom" in n.lines.columns:
        ratings = pd.to_numeric(n.lines["s_nom"], errors="coerce").replace(0, np.nan)
        loading = n.lines_t.p0.abs().div(ratings, axis=1)
        overloaded = int((loading.max(axis=0) > 1.0).sum())
        congested = int(((loading >= CONGESTION_THRESHOLD).sum(axis=0) > 0).sum())
        add("Transmission", "Lines exceeding 1.0 p.u.", overloaded, "Check" if overloaded else "Pass")
        add("Transmission", f"Lines reaching {CONGESTION_THRESHOLD:.2f} p.u.", congested, "Check" if congested else "Pass")

    if not n.buses_t.marginal_price.empty:
        price = n.buses_t.marginal_price
        add("Prices", "Mean marginal price [$/MWh]", f"{float(price.mean().mean()):,.2f}")
        add("Prices", "Maximum marginal price [$/MWh]", f"{float(price.max().max()):,.2f}", "Check")
        negative_hours = int((price.mean(axis=1) < 0).sum())
        add("Prices", "Hours with negative mean price", negative_hours, "Check" if negative_hours else "Pass")

    return pd.DataFrame(rows)


def _render_report_html(
    network_name: str,
    artifacts: list[PlotArtifact],
    summary_table: pd.DataFrame,
    output_path: Path,
) -> None:
    logo_path = _find_logo_path()
    logo_html = f"<img src=\"../{logo_path}\" alt=\"Delta E plus\" class=\"brand-logo\">" if logo_path else ""
    summary_html = summary_table.to_html(index=False, border=0, classes="summary-table")
    input_inventory = _build_input_inventory_table()
    input_inventory_html = (
        input_inventory.to_html(index=False, border=0, classes="summary-table")
        if not input_inventory.empty
        else "<p>No prepared-input inventory was available.</p>"
    )
    validation_manifest = Path("data/validation/manifest.json")
    validation_state = (
        "Pinned validation manifest found"
        if validation_manifest.exists()
        else "Validation data not fetched — run fetch_validation_data"
    )

    category_order = ["landing", "input", "hydro", "system", "other"]
    grouped: dict[str, list[PlotArtifact]] = {k: [] for k in category_order}
    for artifact in artifacts:
        grouped.setdefault(artifact.category if artifact.category in grouped else "other", []).append(artifact)

    def render_artifact(artifact: PlotArtifact) -> str:
        if artifact.kind == "static-image":
            body = (
                f"<img src=\"{artifact.relative_path}\" alt=\"{artifact.title}\" "
                "style=\"width:100%;border-radius:10px;\" />"
            )
        elif artifact.kind == "interactive-html":
            body = (
                "<div class=\"interactive-hint\">Interactive view · use the labelled selector inside the chart, then hover, zoom, or drag to inspect.</div>"
                f"<iframe src=\"{artifact.relative_path}\" title=\"{artifact.title}\" "
                "style=\"width:100%;height:620px;border:none;\"></iframe>"
            )
        else:
            body = (
                f"<p><a href=\"{artifact.relative_path}\" target=\"_blank\" rel=\"noopener\">Open exported file</a></p>"
            )
        return f"""
        <section class=\"card\">
          <h3>{artifact.title}</h3>
          <p>{artifact.description}</p>
          {body}
        </section>
        """

    input_artifacts = grouped.get("input", [])
    detailed_map = next((a for a in input_artifacts if "Detailed input network" in a.title), None)
    topology_map = next((a for a in input_artifacts if "Network Topology" in a.title), None)
    ceei_map = next((a for a in input_artifacts if "ceei" in a.title.lower()), None)
    load_input = next((a for a in input_artifacts if "load profile" in a.title.lower()), None)

    hydro_results = [a for a in grouped.get("hydro", []) if a.kind == "interactive-html"]
    system_results = [a for a in grouped.get("system", []) if a.kind == "interactive-html"]
    landing_cards = "".join(render_artifact(a) for a in grouped.get("landing", []))

    def render_network_maps() -> str:
        if not detailed_map and not topology_map:
            return ""

        def map_card(a: PlotArtifact | None, fallback_title: str) -> str:
            if a is None:
                return f"<section class=\"card\"><h3>{fallback_title}</h3><p>Not available for this run.</p></section>"
            return (
                "<section class=\"card\">"
                f"<h3>{a.title}</h3><p>{a.description}</p>"
                f"<img src=\"{a.relative_path}\" alt=\"{a.title}\" style=\"width:100%;border-radius:10px;\"/>"
                "</section>"
            )

        return (
            "<div class=\"map-stack\">"
            + map_card(detailed_map, "Detailed Input Network")
            + map_card(topology_map, "Network Topology")
            + "</div>"
        )

    input_blocks = []
    input_blocks.append(
        "<section class=\"method-note\"><strong>Topology and demand assumptions.</strong> "
        "The first map shows the prepared substation-level input tables; the second map shows the "
        "simplified topology retained by the solved PyPSA network. Simplification is required where "
        "complete measured line, transformer, or substation data are unavailable. Provincial demand "
        "is allocated using CEEI residential and commercial/SMI regional proportions.</section>"
    )
    input_blocks.append(render_network_maps())
    input_blocks.append(
        "<section class=\"card\"><h3>Prepared input inventory</h3>"
        "<p class=\"section-intro\">Records read from <code>data/processed_data</code> and "
        "<code>data/pypsa_data</code>. Potential and committed resources remain distinct from the "
        "solved-network fleet.</p>"
        f"{input_inventory_html}</section>"
    )
    if load_input:
        input_blocks.append(render_artifact(load_input))
    if ceei_map:
        input_blocks.append(render_artifact(ceei_map))
    reserved = {id(a) for a in [detailed_map, topology_map, load_input, ceei_map] if a is not None}
    for artifact in input_artifacts:
        if id(artifact) not in reserved:
            input_blocks.append(render_artifact(artifact))

    hydro_blocks = "".join(render_artifact(a) for a in hydro_results) or "<section class=\"card\"><p>No interactive hydro artifacts available.</p></section>"
    system_blocks = "".join(render_artifact(a) for a in system_results) or "<section class=\"card\"><p>No interactive system artifacts available.</p></section>"

    html = f"""
<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width,initial-scale=1\" />
  <title>{network_name} visual report</title>
  <style>
    /* Shared design tokens, matched to the PyPSA-BC documentation website. */
    :root {{
      --ink:#10251f; --muted:#5c6d67; --paper:#f7f8f3; --card:#fff;
      --forest:#0b5d4b; --hero:#063f35; --hero2:#0b6956; --teal:#1b8a79; --lime:#cbe66b;
      --line:#d8e0d8; --amber:#b56a12; --red:#a13e36;
      --radius:14px; --shadow:0 14px 40px rgba(16,37,31,.09);
      --serif:Georgia,"Times New Roman",serif; --sans:Arial,Helvetica,sans-serif;
      --mono:Consolas,"Courier New",monospace;
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font:16px/1.65 var(--serif); color:var(--ink); background:var(--paper); }}
    h1,h2,h3,h4 {{ font-family:var(--sans); line-height:1.15; }}
    code {{ font-family:var(--mono); background:#e9eeea; padding:2px 5px; border-radius:3px; font-size:.9em; }}
    a {{ color:var(--forest); }}

    /* Sticky top bar with tab navigation. */
    .topbar {{ position:sticky; top:0; z-index:1000; background:rgba(247,248,243,.95); backdrop-filter:blur(14px); border-bottom:1px solid var(--line); }}
    .topbar-inner {{ max-width:1320px; margin:0 auto; padding:10px 24px; display:flex; justify-content:space-between; align-items:center; gap:14px; }}
    .top-tabs {{ display:flex; gap:6px; align-items:center; }}
    .top-btn {{ border:0; border-bottom:2px solid transparent; background:transparent; color:var(--muted); border-radius:0; padding:12px 12px; cursor:pointer; font:800 13px var(--sans); }}
    .top-btn:hover {{ color:var(--forest); }}
    .top-btn.active {{ color:var(--forest); border-color:var(--forest); }}
    .top-btn.icon-only {{ width:38px; padding:10px 0; text-align:center; font-size:16px; }}
    .top-right {{ display:flex; gap:12px; align-items:center; }}
    .about-link {{ border:1px solid var(--line); color:var(--ink); background:#fff; border-radius:5px; padding:8px 12px; text-decoration:none; font:800 12px var(--sans); }}
    .about-link:hover {{ border-color:var(--forest); color:var(--forest); }}
    .brand-logo {{ width:72px; max-height:34px; object-fit:contain; filter:brightness(0) saturate(100%); opacity:.72; }}

    /* Hero, mirroring the website gradient. */
    .hero {{ background:linear-gradient(125deg,var(--hero) 0%,var(--hero) 48%,var(--hero2) 150%); color:#fff; padding:34px 0 30px; }}
    .container {{ max-width:1320px; margin:0 auto; padding:24px; }}
    h1 {{ margin:0 0 .4rem; font-size:clamp(26px,3.4vw,38px); letter-spacing:-.02em; }}
    .hero-sub {{ color:#d8ebe5; margin:0; max-width:820px; }}
    .validation-state {{ display:inline-block; margin-top:12px; padding:5px 9px; background:rgba(203,230,107,.16); color:var(--lime); border:1px solid rgba(203,230,107,.4); font:800 .74rem var(--sans); letter-spacing:.03em; border-radius:4px; }}

    /* Cards, callouts and helper notes. */
    .card {{ background:var(--card); border:1px solid var(--line); border-radius:var(--radius); padding:20px; margin-bottom:18px; box-shadow:0 2px 0 rgba(16,37,31,.02); }}
    .card h3 {{ margin-top:0; font-size:20px; }}
    .btn {{ border:none; border-radius:5px; padding:9px 14px; cursor:pointer; background:var(--forest); color:#fff; font:800 13px var(--sans); }}
    .method-note {{ border-left:4px solid var(--teal); background:#eef8f4; padding:15px 18px; margin-bottom:18px; line-height:1.6; border-radius:0 6px 6px 0; }}
    .method-note strong {{ font-family:var(--sans); }}
    .section-intro {{ color:var(--muted); max-width:900px; }}
    .interactive-hint {{ background:#eef8f4; color:#24483e; border-left:3px solid var(--forest); padding:9px 12px; margin:4px 0 12px; font:700 .86rem var(--sans); border-radius:0 4px 4px 0; }}

    /* Tab panels and sub-tabs. */
    .top-panel {{ display:none; }}
    .top-panel.active {{ display:block; }}
    .subtabs {{ display:flex; gap:8px; flex-wrap:wrap; margin-bottom:16px; position:sticky; top:66px; z-index:500; background:var(--paper); padding:12px 0; border-bottom:1px solid var(--line); }}
    .subtab-btn {{ border:1px solid var(--line); background:#fff; color:var(--forest); border-radius:6px; padding:8px 14px; cursor:pointer; font:800 13px var(--sans); }}
    .subtab-btn:hover {{ border-color:var(--teal); }}
    .subtab-btn.active {{ background:var(--forest); color:#fff; border-color:var(--forest); }}
    .subtab-panel {{ display:none; }}
    .subtab-panel.active {{ display:block; }}

    /* Layout helpers for maps. */
    .two-col {{ display:grid; grid-template-columns:1fr 1fr; gap:14px; }}
    .map-stack {{ display:grid; grid-template-columns:1fr; gap:18px; }}
    .map-stack .card {{ max-width:1180px; margin-inline:auto; width:100%; }}
    .map-stack img {{ display:block; width:100%; height:auto; object-fit:contain; border-radius:10px; }}

    /* Tables, matched to the website table treatment. */
    .summary-table {{ width:100%; border-collapse:collapse; margin-top:10px; background:#fff; font-size:14px; }}
    .summary-table th {{ font:800 12px var(--sans); text-transform:uppercase; letter-spacing:.05em; background:#e7efea; }}
    .summary-table th, .summary-table td {{ border:1px solid var(--line); text-align:left; padding:10px 12px; vertical-align:top; }}
    .footer {{ color:var(--muted); font-size:.92rem; }}

    /* Summary modal. */
    .modal-bg {{ position:fixed; inset:0; background:rgba(16,37,31,.5); display:none; align-items:center; justify-content:center; z-index:2000; }}
    .modal {{ background:#fff; border-radius:12px; width:min(900px,94vw); max-height:85vh; overflow:auto; padding:20px; }}
    .close-row {{ display:flex; justify-content:space-between; align-items:center; margin-bottom:8px; }}

    @media (max-width: 900px) {{
      .two-col {{ grid-template-columns:1fr; }}
      .topbar-inner {{ flex-direction:column; align-items:flex-start; }}
      .subtabs {{ top:118px; }}
    }}
  </style>
</head>
<body>
    <div class=\"topbar\">
        <div class=\"topbar-inner\">
            <div class=\"top-tabs\">
                <button class=\"top-btn icon-only\" id=\"about-open\" title=\"About network\" aria-label=\"About network\">ⓘ</button>
                <button class=\"top-btn active\" data-top=\"input\">Input</button>
                <button class=\"top-btn\" data-top=\"results\">Results</button>
            </div>
            <div class=\"top-right\">
                <a class=\"about-link\" href=\"https://deltae.github.io/PyPSA_BC/\" target=\"_blank\" rel=\"noopener\">About project</a>
                {logo_html}
            </div>
        </div>
    </div>

    <div class=\"hero\">
        <div class=\"container\">
            <h1>PyPSA visual report: {network_name}</h1>
            <p class=\"hero-sub\">Scientific input audit and interactive diagnostics for BC system behaviour, transmission, and hydro operations.</p>
            <span class=\"validation-state\">{validation_state}</span>
            <p class=\"hero-sub\" style=\"margin-top:10px;font-size:.82rem;\">Active development: methods and publications are under preparation. Contact the model developer for full access and all features.</p>
        </div>
    </div>

    <div class=\"container\">
        <div class=\"top-panel active\" id=\"panel-input\">
            {''.join(input_blocks)}
        </div>

        <div class=\"top-panel\" id=\"panel-results\">
            {landing_cards}
            <div class=\"subtabs\">
                <button class=\"subtab-btn active\" data-subtab=\"hydro\">Hydro</button>
                <button class=\"subtab-btn\" data-subtab=\"system\">System</button>
            </div>
            <div class=\"subtab-panel active\" id=\"subtab-hydro\">{hydro_blocks}</div>
            <div class=\"subtab-panel\" id=\"subtab-system\">{system_blocks}</div>
        </div>

        <section class=\"card footer\">
            <h3>Notes</h3>
            <p>Analytical map CRS: NAD83 / BC Albers (EPSG:3005, metres). OpenStreetMap tiles are warped into that CRS when available. Congestion means line loading ≥ {CONGESTION_THRESHOLD:.2f} p.u. Hydro and spill series are inferred from carriers and component naming patterns.</p>
        </section>
  </div>

    <div class=\"modal-bg\" id=\"summary-modal\">
        <div class=\"modal\">
            <div class=\"close-row\">
                <h3 style=\"margin:0;font-family:Arial,sans-serif;\">Network Summary</h3>
                <button class=\"btn\" id=\"summary-close\">Close</button>
            </div>
            {summary_html}
        </div>
    </div>

    <script>
        const topButtons=[...document.querySelectorAll('.top-btn[data-top]')];
        const topPanels={{ input: document.getElementById('panel-input'), results: document.getElementById('panel-results') }};
        topButtons.forEach(btn=>btn.addEventListener('click',()=>{{
            topButtons.forEach(b=>b.classList.remove('active'));
            Object.values(topPanels).forEach(p=>p.classList.remove('active'));
            btn.classList.add('active');
            const key=btn.dataset.top;
            if(topPanels[key]) topPanels[key].classList.add('active');
        }}));

        const subButtons=[...document.querySelectorAll('.subtab-btn')];
        const subPanels=[...document.querySelectorAll('.subtab-panel')];
        subButtons.forEach(btn=>btn.addEventListener('click',()=>{{
            subButtons.forEach(b=>b.classList.remove('active'));
            subPanels.forEach(p=>p.classList.remove('active'));
            btn.classList.add('active');
            const panel=document.getElementById('subtab-'+btn.dataset.subtab);
            if(panel) panel.classList.add('active');
        }}));

        const modal=document.getElementById('summary-modal');
        document.getElementById('about-open').addEventListener('click',()=>modal.style.display='flex');
        document.getElementById('summary-close').addEventListener('click',()=>modal.style.display='none');
        modal.addEventListener('click',(e)=>{{ if(e.target===modal) modal.style.display='none'; }});
    </script>
</body>
</html>
"""

    output_path.write_text(html, encoding="utf-8")


def build_visual_report(
    network_path: str | Path,
    vis_dir: str | Path = "vis",
    resources_subdir: str = "vis_resources",
) -> ReportResult:
    _configure_plot_style()
    network_path = Path(network_path)
    if not network_path.exists():
        raise FileNotFoundError(f"Network file does not exist: {network_path}")

    n = pypsa.Network(str(network_path))

    stem = _sanitize_filename(network_path.stem)
    vis_dir = Path(vis_dir)
    resources_root = vis_dir / resources_subdir
    resources_dir = resources_root / f"{stem}_plots"
    resources_dir.mkdir(parents=True, exist_ok=True)
    vis_dir.mkdir(parents=True, exist_ok=True)

    artifacts: list[PlotArtifact] = []

    input_detailed_map = resources_dir / f"{stem}_input_detailed_network.svg"
    _write_input_detailed_map(n, input_detailed_map, full_bc_frame=True)
    artifacts.append(
        PlotArtifact(
            title="Detailed input network",
            description="Detailed line geometry with substation/generator markers sized by capacity and colored by resource.",
            relative_path=str(input_detailed_map.relative_to(vis_dir)),
            kind="static-image",
            category="input",
        )
    )

    map_static = resources_dir / f"{stem}_network_map_static.svg"
    _write_static_network_map(
        n,
        map_static,
        f"{stem}: network topology",
        full_bc_frame=True,
    )
    artifacts.append(
        PlotArtifact(
            title="Network Topology",
            description=(
                "Modeled PyPSA topology with physical generators sized by nominal "
                "capacity and colored by resource. Invalid/out-of-BC coordinates "
                "and synthetic backstop generators are excluded."
            ),
            relative_path=str(map_static.relative_to(vis_dir)),
            kind="static-image",
            category="input",
        )
    )

    ceei_heatmap = resources_dir / f"{stem}_ceei_disagg_heatmap.png"
    _write_ceei_heatmap_plot(ceei_heatmap)
    artifacts.append(
        PlotArtifact(
            title="CEEI load disaggregation heatmap",
            description="Regional CEEI-based load allocation proportions used for demand disaggregation.",
            relative_path=str(ceei_heatmap.relative_to(vis_dir)),
            kind="static-image",
            category="input",
        )
    )

    map_interactive = resources_dir / f"{stem}_network_map_interactive.html"
    _write_interactive_network_map(n, map_interactive, f"{stem}: interactive network map")
    artifacts.append(
        PlotArtifact(
            title="Interactive network map",
            description="Interactive map for bus inspection and high-level line topology checks.",
            relative_path=str(map_interactive.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    station_static_png = resources_dir / f"{stem}_station_generation_top.png"
    station_interactive_html = resources_dir / f"{stem}_station_generation_interactive.html"
    _write_station_generation_plots(n, station_static_png, station_interactive_html)
    artifacts.append(
        PlotArtifact(
            title="Station-wise generation explorer",
            description="Interactive station selector for hourly generation traces.",
            relative_path=str(station_interactive_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    hydro_html = resources_dir / f"{stem}_hydro_generation_interactive.html"
    _write_hydro_generation_interactive(n, hydro_html)
    artifacts.append(
        PlotArtifact(
            title="Hydro generation explorer",
            description="Interactive hydro station hourly generation explorer with dropdown selector.",
            relative_path=str(hydro_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="hydro",
        )
    )

    spill_html = resources_dir / f"{stem}_spill_interactive.html"
    _write_spill_interactive(n, spill_html)
    artifacts.append(
        PlotArtifact(
            title="Reservoir spill explorer",
            description="Interactive spill-link flow explorer with dropdown selector.",
            relative_path=str(spill_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="hydro",
        )
    )

    storage_html = resources_dir / f"{stem}_storage_interactive.html"
    _write_storage_interactive(n, storage_html)
    artifacts.append(
        PlotArtifact(
            title="Reservoir storage explorer",
            description="Interactive reservoir selector showing stored energy as a fraction of nominal capacity.",
            relative_path=str(storage_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="hydro",
        )
    )

    dispatch_static = resources_dir / f"{stem}_dispatch_monthly_static.png"
    dispatch_interactive = resources_dir / f"{stem}_dispatch_hourly_interactive.html"
    _write_dispatch_plots(n, dispatch_static, dispatch_interactive)
    artifacts.append(
        PlotArtifact(
            title="Dispatch stack (interactive)",
            description="Hourly stacked dispatch with hover/zoom (best for short-window diagnostics).",
            relative_path=str(dispatch_interactive.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    load_png = resources_dir / f"{stem}_load_profiles.png"
    price_duration_png = resources_dir / f"{stem}_price_duration_curve.png"
    price_heatmap_png = resources_dir / f"{stem}_price_heatmap.png"
    _write_load_and_price_plots(n, load_png, price_duration_png, price_heatmap_png)
    load_interactive = resources_dir / f"{stem}_load_interactive.html"
    _write_load_interactive_plot(n, load_interactive)
    artifacts.append(
        PlotArtifact(
            title="Provincial load profile",
            description=(
                "Interactive prepared 2021 provincial load and duration curve in GW, "
                "including peak, mean, annual energy, and peak timestamp."
            ),
            relative_path=str(load_interactive.relative_to(vis_dir)),
            kind="interactive-html",
            category="input",
        )
    )

    price_interactive = resources_dir / f"{stem}_price_interactive.html"
    _write_price_interactive_plot(n, price_interactive)
    artifacts.append(
        PlotArtifact(
            title="Price explorer",
            description="Interactive price duration and timeseries diagnostics.",
            relative_path=str(price_interactive.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    line_loading_png = resources_dir / f"{stem}_line_loading_duration.png"
    _write_line_loading_plot(n, line_loading_png)
    line_loading_html = resources_dir / f"{stem}_line_loading_duration_interactive.html"
    _write_line_loading_interactive(n, line_loading_html)
    artifacts.append(
        PlotArtifact(
            title="Transmission congestion duration explorer",
            description="Interactive top-line loading duration curves with 0.85 p.u. threshold.",
            relative_path=str(line_loading_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    congestion_map_png = resources_dir / f"{stem}_line_capacity_congestion_map.png"
    _write_static_capacity_congestion_map(n, congestion_map_png)
    artifacts.append(
        PlotArtifact(
            title="BC line capacity and congestion map (static)",
            description=(
                "Line width scales with capacity (s_nom), and color encodes congestion share "
                "(fraction of hours with loading >= 0.85), highlighting lines congested >50% of the year."
            ),
            relative_path=str(congestion_map_png.relative_to(vis_dir)),
            kind="static-image",
            category="landing",
        )
    )

    curtailment_html = resources_dir / f"{stem}_curtailment_interactive.html"
    _write_curtailment_interactive_plot(n, curtailment_html)
    artifacts.append(
        PlotArtifact(
            title="Curtailment explorer",
            description="Interactive curtailment share by generation technology.",
            relative_path=str(curtailment_html.relative_to(vis_dir)),
            kind="interactive-html",
            category="system",
        )
    )

    calibration_dir = Path("data/validation")
    model_year = 2021
    if isinstance(n.snapshots, pd.DatetimeIndex) and len(n.snapshots) > 0:
        model_year = int(n.snapshots[0].year)
    load_calib_png = resources_dir / f"{stem}_load_calibration_compare.png"
    has_obs = _write_load_calibration_plot(n, load_calib_png, calibration_dir, model_year)
    artifacts.append(
        PlotArtifact(
            title="Model vs observed load",
            description=(
                f"Monthly load comparison against fetched BC Hydro balancing authority data ({model_year})"
                if has_obs
                else "Observed calibration data unavailable; model-only monthly load shown."
            ),
            relative_path=str(load_calib_png.relative_to(vis_dir)),
            kind="static-image",
            category="input",
        )
    )

    # Build and register the reproducible standalone input suite. Importing here
    # avoids a module-level cycle because the script reuses shared map helpers.
    from workflow.scripts.build_input_visuals import build_all as build_input_visuals

    input_visual_dir = vis_dir / "input_visuals"
    build_input_visuals(input_visual_dir)
    extra_inputs = [
        (
            "Provincial Load Disaggregation",
            "CEEI-based residential and commercial/small-industry allocation, with a BC regional map and ranked shares.",
            input_visual_dir / "provincial_load_disaggregation.svg",
            "static-image",
        ),
        (
            "CODERS Generation Parameter Explorer",
            "Interactive comparison of generic technology size, cost, efficiency, emissions, flexibility, and outage assumptions.",
            input_visual_dir / "coders_generation_parameters.html",
            "interactive-html",
        ),
        (
            "Hydrological Input Basins",
            "HydroBASINS level-12 catchments intersecting BC, with modeled hydro assets sized by capacity.",
            input_visual_dir / "hydrobasins_bc_inputs.svg",
            "static-image",
        ),
        (
            "Modeled Hydro Cascade Structure",
            "Major modeled cascades, downstream sequence, plant type, installed capacity, and total cascade capacity.",
            input_visual_dir / "hydro_cascade_schematic.svg",
            "static-image",
        ),
        (
            "Hydro Cascade Asset Flow",
            (
                "Interactive reservoir-to-generation-to-downstream-reservoir connectivity. "
                "Direction represents modeled water routing; link width is an installed-capacity proxy."
            ),
            input_visual_dir / "hydro_cascade_asset_flow.html",
            "interactive-html",
        ),
        (
            "Demand Input Calibration Explorer",
            "Interactive monthly magnitude and shape comparison with pinned 2021 BC Hydro balancing-authority observations.",
            input_visual_dir / "demand_calibration_2021.html",
            "interactive-html",
        ),
    ]
    for title, description, path, kind in extra_inputs:
        artifacts.append(
            PlotArtifact(
                title=title,
                description=description,
                relative_path=str(path.relative_to(vis_dir)),
                kind=kind,
                category="input",
            )
        )

    summary = _build_summary_table(n)
    summary_csv = resources_dir / f"{stem}_summary_metrics.csv"
    summary.to_csv(summary_csv, index=False)
    artifacts.append(
        PlotArtifact(
            title="Summary metrics table (CSV)",
            description="Tabular summary exported for downstream analytics.",
            relative_path=str(summary_csv.relative_to(vis_dir)),
            kind="file-link",
            category="other",
        )
    )

    report_path = vis_dir / f"{stem}_visual_report.html"
    _render_report_html(stem, artifacts, summary, report_path)

    return ReportResult(
        network_path=network_path,
        report_path=report_path,
        resources_dir=resources_dir,
        artifacts=artifacts,
    )
