"""Build publication-oriented PyPSA-BC input visualizations.

Examples
--------
python -m workflow.scripts.build_input_visuals
python -m workflow.scripts.build_input_visuals --output-dir vis/input_visuals

All inputs are local project data. The script does not download or modify source
datasets, making the visual outputs reproducible from a checked data snapshot.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from pypsa_bc.vis.network_visual_report import (
    BC_ALBERS,
    _add_grey_basemap,
    _load_bc_gadm_regions,
    _parse_bch_hourly_load_kwh,
    _set_bc_map_frame,
    _write_input_detailed_map,
    _write_static_network_map,
)

OFF_WHITE = "#fbfaf6"
INK = "#26332f"
BLUE = "#155a78"
GOLD = "#b8860b"
TEAL = "#2a9d8f"
ORANGE = "#e76f51"


def _save_map(fig: plt.Figure, output_stem: Path) -> list[Path]:
    paths = [output_stem.with_suffix(".png"), output_stem.with_suffix(".svg")]
    fig.savefig(paths[0], dpi=320, bbox_inches="tight", facecolor=OFF_WHITE)
    fig.savefig(paths[1], bbox_inches="tight", facecolor=OFF_WHITE)
    plt.close(fig)
    return paths


def build_network_map_previews(output_dir: Path) -> list[Path]:
    """Save uncropped full-BC previews without registering them in the report."""
    import pypsa

    network_path = Path("results/PyPSA_BC_network_s_2021.nc")
    if not network_path.exists():
        raise FileNotFoundError(f"Solved network unavailable: {network_path}")
    network = pypsa.Network(network_path)
    outputs = [
        output_dir / "detailed_input_network_full_bc.png",
        output_dir / "detailed_input_network_full_bc.svg",
        output_dir / "network_topology_full_bc.png",
        output_dir / "network_topology_full_bc.svg",
    ]
    _write_input_detailed_map(network, outputs[0], full_bc_frame=True)
    _write_input_detailed_map(network, outputs[1], full_bc_frame=True)
    _write_static_network_map(
        network,
        outputs[2],
        "Network Topology",
        full_bc_frame=True,
    )
    _write_static_network_map(
        network,
        outputs[3],
        "Network Topology",
        full_bc_frame=True,
    )
    return outputs


def build_regional_load_map(output_dir: Path) -> list[Path]:
    """Map CEEI-based provincial load allocation by regional district."""
    c_path = Path("data/processed_data/load/CEEI_RD_ELEC_proportions.csv")
    regions = _load_bc_gadm_regions()
    if regions is None or not c_path.exists():
        raise FileNotFoundError("BC GADM boundary or CEEI proportion table is unavailable")
    shares = pd.read_csv(c_path)

    def key(value: object) -> str:
        return "".join(ch for ch in str(value).lower() if ch.isalnum())

    aliases = {
        "greatervancouver": "metrovancouver",
        "comoxstrathcona": "strathcona",
        "skeenaqueencharlotte": "haidagwaii",
    }
    shares["key"] = shares["REGION"].map(key).replace(aliases)
    shares["total_share"] = shares["PROPORTION_RES"] + shares["PROPORTION_CSMI"]
    regions = regions.copy()
    regions["key"] = regions["NAME_2"].map(key)
    mapped = regions.merge(
        shares[["REGION", "key", "PROPORTION_RES", "PROPORTION_CSMI", "total_share"]],
        on="key",
        how="left",
    ).to_crs(BC_ALBERS)

    fig = plt.figure(figsize=(15, 9), facecolor=OFF_WHITE)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.08)
    ax = fig.add_subplot(grid[0, 0])
    rank_ax = fig.add_subplot(grid[0, 1])
    mapped.plot(
        column="total_share",
        cmap="Blues",
        edgecolor="#707b76",
        linewidth=0.45,
        legend=False,
        missing_kwds={"color": "#e5e5e1", "label": "No CEEI match"},
        ax=ax,
        zorder=2,
    )
    _set_bc_map_frame(ax, regions)
    top_labels = mapped.nlargest(7, "total_share")
    for _, row in top_labels.iterrows():
        point = row.geometry.representative_point()
        ax.text(
            point.x,
            point.y,
            f"{row['total_share']:.1%}",
            ha="center",
            va="center",
            fontsize=8,
            color="white" if row["total_share"] > 0.04 else INK,
            weight="bold",
            zorder=3,
        )
    ax.set_title("Regional allocation map", loc="left", fontsize=13, weight="bold")
    ax.set_axis_off()

    ranking = shares.nlargest(12, "total_share").sort_values("total_share")
    rank_ax.barh(
        ranking["REGION"],
        ranking["PROPORTION_RES"],
        color=BLUE,
        label="Residential",
    )
    rank_ax.barh(
        ranking["REGION"],
        ranking["PROPORTION_CSMI"],
        left=ranking["PROPORTION_RES"],
        color=GOLD,
        label="Commercial / SMI",
    )
    for y_pos, total in enumerate(ranking["total_share"]):
        rank_ax.text(total + 0.002, y_pos, f"{total:.1%}", va="center", fontsize=8)
    rank_ax.set_xlabel("Share of provincial electricity allocation")
    rank_ax.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    rank_ax.set_title("Largest regional allocations", loc="left", fontsize=13, weight="bold")
    rank_ax.legend(loc="lower right", frameon=False)
    rank_ax.grid(axis="x", color="#d8ddd9", linewidth=0.6)
    rank_ax.spines[["top", "right", "left"]].set_visible(False)
    fig.suptitle(
        "Provincial Load Disaggregation",
        x=0.025,
        y=0.98,
        ha="left",
        fontsize=19,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.025,
        0.93,
        "CEEI 2021 residential and commercial/small-industry electricity shares allocate the provincial hourly profile to regional districts.",
        fontsize=10,
        color="#53615b",
    )
    return _save_map(fig, output_dir / "provincial_load_disaggregation")


def build_coders_parameter_explorer(output_dir: Path) -> Path:
    """Interactive comparison of generic CODERS technology assumptions."""
    path = Path("data/downloaded_data/CODERS/data-pull/supply/generation_generic.csv")
    data = pd.read_csv(path)
    metrics = [
        ("capital_cost_CAD_per_kW", "Capital cost [CAD/kW]"),
        ("annualized_capital_cost_CAD_per_kwyear", "Annualized capital cost [CAD/kW-year]"),
        ("variable_om_costs", "Variable O&M [CAD/MWh]"),
        ("efficiency", "Efficiency [p.u.]"),
        ("carbon_emissions", "Carbon emissions [tCO₂e/MWh]"),
        ("forced_outage_rate", "Forced outage rate [p.u.]"),
        ("ramp_rate_percent_per_min", "Ramp rate [%/min]"),
    ]
    palette = px.colors.qualitative.Safe
    classes = sorted(data["gen_type_copper"].fillna("other").astype(str).unique())
    colors = {name: palette[i % len(palette)] for i, name in enumerate(classes)}
    fig = go.Figure()
    for metric_index, (metric, label) in enumerate(metrics):
        for technology_class in classes:
            group = data.loc[data["gen_type_copper"].fillna("other").astype(str) == technology_class]
            fig.add_trace(
                go.Scatter(
                    x=group["generic_generation_size_MW"],
                    y=pd.to_numeric(group[metric], errors="coerce"),
                    mode="markers",
                    name=technology_class,
                    legendgroup=technology_class,
                    showlegend=metric_index == 0,
                    visible=metric_index == 0,
                    marker=dict(
                        size=11,
                        color=colors[technology_class],
                        line=dict(color="white", width=0.8),
                    ),
                    customdata=np.column_stack(
                        [
                            group["generation_type_description"].fillna(group["gen_type"]),
                            group["economic_life"],
                            group["min_capacity_factor"],
                            group["max_capacity_factor"],
                        ]
                    ),
                    hovertemplate=(
                        "<b>%{customdata[0]}</b><br>Generic size: %{x:.1f} MW"
                        f"<br>{label}: %{{y:.3g}}"
                        "<br>Economic life: %{customdata[1]} years"
                        "<br>Capacity factor range: %{customdata[2]}–%{customdata[3]}"
                        "<extra></extra>"
                    ),
                )
            )
    buttons = []
    traces_per_metric = len(classes)
    for metric_index, (_, label) in enumerate(metrics):
        visible = [False] * (len(metrics) * traces_per_metric)
        start = metric_index * traces_per_metric
        visible[start : start + traces_per_metric] = [True] * traces_per_metric
        buttons.append(
            dict(
                label=label,
                method="update",
                args=[{"visible": visible}, {"yaxis": {"title": label, "rangemode": "tozero"}}],
            )
        )
    fig.update_layout(
        title=(
            "CODERS Generic Generation Parameters"
            "<br><sup>Technology assumptions used to enrich PyPSA-BC generation assets</sup>"
        ),
        xaxis_title="Generic generation size [MW]",
        yaxis_title=metrics[0][1],
        template="plotly_white",
        paper_bgcolor=OFF_WHITE,
        plot_bgcolor=OFF_WHITE,
        font=dict(family="Arial, sans-serif", color=INK),
        height=620,
        margin=dict(l=75, r=35, t=135, b=65),
        legend=dict(orientation="h", y=-0.2),
        updatemenus=[
            dict(
                buttons=buttons,
                x=0,
                y=1.13,
                xanchor="left",
                yanchor="top",
                bgcolor="#fffdf7",
                bordercolor="#7b8b84",
            )
        ],
        annotations=[
            dict(
                text="<b>Select parameter</b>",
                x=0,
                y=1.22,
                xref="paper",
                yref="paper",
                showarrow=False,
                xanchor="left",
            )
        ],
    )
    output = output_dir / "coders_generation_parameters.html"
    fig.write_html(output, include_plotlyjs="cdn")
    return output


def _read_bc_basins() -> gpd.GeoDataFrame:
    bbox = (-139.2, 48.2, -114.0, 60.1)
    paths = [
        Path("data/downloaded_data/HydroBASINS/hybas_na_lev12_v1c.shp"),
        Path("data/downloaded_data/HydroBASINS_arctic/hybas_ar_lev12_v1c.shp"),
    ]
    frames = [gpd.read_file(path, bbox=bbox) for path in paths if path.exists()]
    if not frames:
        raise FileNotFoundError("HydroBASINS North America and Arctic shapefiles are unavailable")
    basins = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=frames[0].crs)
    return basins.drop_duplicates("HYBAS_ID")


def build_basin_map(output_dir: Path) -> list[Path]:
    """Show level-12 HydroBASINS units and hydro assets used by the model."""
    regions = _load_bc_gadm_regions()
    if regions is None:
        raise FileNotFoundError("BC boundary unavailable")
    basins = _read_bc_basins().to_crs(BC_ALBERS)
    bc_union = regions.to_crs(BC_ALBERS).geometry.union_all()
    basins = basins.loc[basins.intersects(bc_union)].copy()
    hydro = pd.read_csv("data/processed_data/hydro/existing/hydro_generation.csv")
    hydro = gpd.GeoDataFrame(
        hydro,
        geometry=gpd.points_from_xy(hydro["longitude"], hydro["latitude"]),
        crs="EPSG:4326",
    ).to_crs(BC_ALBERS)
    hydro = hydro.loc[hydro.within(bc_union)]

    fig, ax = plt.subplots(figsize=(9.2, 10.2), facecolor=OFF_WHITE)
    basins.plot(
        column=np.log10(pd.to_numeric(basins["UP_AREA"], errors="coerce").clip(lower=1)),
        cmap="Greys",
        linewidth=0.08,
        edgecolor="#a4aaa7",
        alpha=0.72,
        legend=True,
        legend_kwds={
            "label": "log₁₀ upstream area [km²]",
            "orientation": "horizontal",
            "shrink": 0.55,
            "pad": 0.025,
        },
        ax=ax,
        zorder=1,
    )
    regions.to_crs(BC_ALBERS).boundary.plot(ax=ax, color="#47544e", linewidth=0.7, zorder=2)
    ax.scatter(
        hydro.geometry.x,
        hydro.geometry.y,
        s=12 + 2.6 * np.sqrt(pd.to_numeric(hydro["capacity"], errors="coerce").clip(0, 5000)),
        c=BLUE,
        edgecolors="white",
        linewidths=0.45,
        alpha=0.88,
        label="Hydro generation asset",
        zorder=3,
    )
    _set_bc_map_frame(ax, regions)
    ax.set_title("Hydrological Input Basins", loc="left", fontsize=16, weight="bold")
    ax.text(
        0,
        1.005,
        f"HydroBASINS level 12 units intersecting BC (n={len(basins):,}); marker area scales with hydro capacity",
        transform=ax.transAxes,
        fontsize=9,
        color="#53615b",
    )
    ax.legend(loc="lower left", frameon=True, facecolor=OFF_WHITE)
    ax.set_axis_off()
    return _save_map(fig, output_dir / "hydrobasins_bc_inputs")


def build_cascade_schematic(output_dir: Path) -> list[Path]:
    """Create small-multiple schematics of the major modeled hydro cascades."""
    hydro = pd.read_csv("data/processed_data/hydro/existing/hydro_generation.csv")
    hydro["capacity"] = pd.to_numeric(hydro["capacity"], errors="coerce").fillna(0)
    hydro["cascade_order"] = pd.to_numeric(hydro["cascade_order"], errors="coerce")
    grouped = (
        hydro.loc[hydro["cascade_group"].ne("DEFAULT")]
        .groupby("cascade_group")["capacity"]
        .sum()
        .sort_values(ascending=False)
    )
    cascades = grouped.head(6).index.tolist()
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), facecolor=OFF_WHITE)
    for ax, cascade in zip(axes.flat, cascades):
        group = hydro.loc[hydro["cascade_group"].eq(cascade)].copy()
        group = group.sort_values(["cascade_order", "longitude"])
        group["plot_order"] = group["cascade_order"].fillna(
            pd.Series(np.arange(1, len(group) + 1), index=group.index)
        )
        x = np.linspace(0.2, 0.8, len(group)) if len(group) > 1 else np.array([0.5])
        y = -group["plot_order"].to_numpy(dtype=float)
        for i in range(len(group) - 1):
            ax.annotate(
                "",
                xy=(x[i + 1], y[i + 1] + 0.08),
                xytext=(x[i], y[i] - 0.08),
                arrowprops=dict(arrowstyle="->", color="#6096ba", lw=1.5),
            )
        for i, (_, row) in enumerate(group.iterrows()):
            marker = "s" if "reservoir" in str(row["hydro_type"]).lower() else "o"
            ax.scatter(
                x[i],
                y[i],
                s=90 + 2.8 * np.sqrt(min(row["capacity"], 5000)) * 10,
                marker=marker,
                facecolor="#dbeaf2" if marker == "s" else OFF_WHITE,
                edgecolor=BLUE,
                linewidth=1.5,
                zorder=3,
            )
            ax.text(
                x[i],
                y[i] - 0.18,
                f"{row['asset_id'].replace('_GSS', '')}\n{row['capacity']:,.0f} MW",
                ha="center",
                va="top",
                fontsize=8,
            )
        ax.set_title(f"{cascade} · {group['capacity'].sum():,.0f} MW", color=BLUE, weight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(min(y) - 0.75, max(y) + 0.55)
        ax.set_axis_off()
    fig.suptitle(
        "Modeled Hydro Cascade Structure",
        x=0.04,
        ha="left",
        fontsize=18,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.04,
        0.94,
        "Arrows indicate modeled downstream sequence; squares are reservoir plants and circles are run-of-river/water-linked plants.",
        fontsize=10,
        color="#53615b",
    )
    fig.tight_layout(rect=[0.03, 0.03, 0.98, 0.91])
    return _save_map(fig, output_dir / "hydro_cascade_schematic")


def build_cascade_asset_flow(output_dir: Path) -> Path:
    """Interactive reservoir → generation asset → downstream reservoir flow."""
    hydro = pd.read_csv("data/processed_data/hydro/existing/hydro_generation.csv")
    hydro["capacity"] = pd.to_numeric(hydro["capacity"], errors="coerce").fillna(0)
    hydro = hydro.loc[
        hydro["cascade_group"].notna()
        & hydro["cascade_group"].ne("DEFAULT")
        & hydro["upper_reservoir_id"].notna()
        & hydro["upper_reservoir_id"].ne("DEFAULT")
    ].copy()
    cascade_capacity = (
        hydro.groupby("cascade_group")["capacity"].sum().sort_values(ascending=False)
    )
    cascades = cascade_capacity.index.tolist()
    reservoir_color = "#b9d7e8"
    plant_colors = {
        "reservoir": "#155a78",
        "reservoir-impute": "#477c98",
        "ror-water": "#b8860b",
        "ror": "#b8860b",
    }

    fig = go.Figure()
    trace_meta: list[tuple[str, float, int, int]] = []

    # Combined view: identical reservoir IDs are deliberately reused as one
    # node, so connections shared across cascade labels remain visible.
    all_group = hydro.sort_values(["cascade_group", "cascade_order", "asset_id"])
    all_keys: list[str] = []
    all_labels: dict[str, str] = {}
    all_colors: dict[str, str] = {}

    def add_all_node(key: str, label: str, color: str) -> None:
        if key not in all_labels:
            all_keys.append(key)
            all_labels[key] = label
            all_colors[key] = color

    for _, row in all_group.iterrows():
        cascade = str(row["cascade_group"])
        upper = str(row["upper_reservoir_id"])
        plant = f"PLANT::{row['component_id']}"
        lower_raw = str(row.get("lower_reservoir_id", "DEFAULT"))
        lower = (
            f"OUTLET::{cascade}"
            if lower_raw in {"DEFAULT", "nan", "None", ""}
            else lower_raw
        )
        add_all_node(
            upper,
            f"Reservoir - {upper.replace('BC_', '').replace('_RES', '')}",
            reservoir_color,
        )
        add_all_node(
            plant,
            (
                f"Plant - {str(row['asset_id']).replace('BC_', '').replace('_GSS', '')}"
                f" [{cascade}]<br>{row['capacity']:,.0f} MW"
            ),
            plant_colors.get(str(row["hydro_type"]).lower(), "#155a78"),
        )
        add_all_node(
            lower,
            (
                f"Outlet - {cascade}"
                if lower.startswith("OUTLET::")
                else f"Reservoir - {lower.replace('BC_', '').replace('_RES', '')}"
            ),
            "#d9ddd9" if lower.startswith("OUTLET::") else reservoir_color,
        )

    all_index = {key: i for i, key in enumerate(all_keys)}
    all_source: list[int] = []
    all_target: list[int] = []
    all_value: list[float] = []
    all_link_color: list[str] = []
    all_link_label: list[str] = []
    for _, row in all_group.iterrows():
        cascade = str(row["cascade_group"])
        upper = str(row["upper_reservoir_id"])
        plant = f"PLANT::{row['component_id']}"
        lower_raw = str(row.get("lower_reservoir_id", "DEFAULT"))
        lower = (
            f"OUTLET::{cascade}"
            if lower_raw in {"DEFAULT", "nan", "None", ""}
            else lower_raw
        )
        capacity = max(float(row["capacity"]), 0.1)
        all_source.extend([all_index[upper], all_index[plant]])
        all_target.extend([all_index[plant], all_index[lower]])
        all_value.extend([capacity, capacity])
        all_link_color.extend(["rgba(21,90,120,0.26)", "rgba(184,134,11,0.26)"])
        all_link_label.extend(
            [
                f"{cascade}: reservoir intake to generating asset",
                f"{cascade}: generating asset to downstream reservoir/outlet",
            ]
        )

    fig.add_trace(
        go.Sankey(
            visible=True,
            arrangement="snap",
            orientation="h",
            node=dict(
                pad=17,
                thickness=18,
                line=dict(color="#52615a", width=0.7),
                label=[all_labels[key] for key in all_keys],
                color=[all_colors[key] for key in all_keys],
                hovertemplate="%{label}<br>Connected capacity: %{value:,.0f} MW<extra></extra>",
            ),
            link=dict(
                source=all_source,
                target=all_target,
                value=all_value,
                color=all_link_color,
                label=all_link_label,
                hovertemplate="%{label}<br>Capacity proxy: %{value:,.0f} MW<extra></extra>",
            ),
            name="All cascades",
        )
    )
    trace_meta.append(
        (
            "All cascades",
            float(all_group["capacity"].sum()),
            len(all_group),
            len(all_keys),
        )
    )

    for trace_index, cascade in enumerate(cascades):
        group = hydro.loc[hydro["cascade_group"].eq(cascade)].sort_values(
            ["cascade_order", "asset_id"]
        )
        node_keys: list[str] = []
        node_labels: dict[str, str] = {}
        node_colors: dict[str, str] = {}

        def add_node(key: str, label: str, color: str) -> None:
            if key not in node_labels:
                node_keys.append(key)
                node_labels[key] = label
                node_colors[key] = color

        for _, row in group.iterrows():
            upper = str(row["upper_reservoir_id"])
            plant = f"PLANT::{row['component_id']}"
            lower_raw = str(row.get("lower_reservoir_id", "DEFAULT"))
            lower = (
                f"OUTLET::{cascade}"
                if lower_raw in {"DEFAULT", "nan", "None", ""}
                else lower_raw
            )
            add_node(
                upper,
                f"Reservoir · {upper.replace('BC_', '').replace('_RES', '')}",
                reservoir_color,
            )
            add_node(
                plant,
                (
                    f"Plant · {str(row['asset_id']).replace('BC_', '').replace('_GSS', '')}"
                    f"<br>{row['capacity']:,.0f} MW"
                ),
                plant_colors.get(str(row["hydro_type"]).lower(), "#155a78"),
            )
            add_node(
                lower,
                (
                    "Cascade outlet"
                    if lower.startswith("OUTLET::")
                    else f"Reservoir · {lower.replace('BC_', '').replace('_RES', '')}"
                ),
                "#d9ddd9" if lower.startswith("OUTLET::") else reservoir_color,
            )

        index = {key: i for i, key in enumerate(node_keys)}
        source: list[int] = []
        target: list[int] = []
        value: list[float] = []
        link_color: list[str] = []
        link_label: list[str] = []
        for _, row in group.iterrows():
            upper = str(row["upper_reservoir_id"])
            plant = f"PLANT::{row['component_id']}"
            lower_raw = str(row.get("lower_reservoir_id", "DEFAULT"))
            lower = (
                f"OUTLET::{cascade}"
                if lower_raw in {"DEFAULT", "nan", "None", ""}
                else lower_raw
            )
            capacity = max(float(row["capacity"]), 0.1)
            source.extend([index[upper], index[plant]])
            target.extend([index[plant], index[lower]])
            value.extend([capacity, capacity])
            link_color.extend(["rgba(21,90,120,0.30)", "rgba(184,134,11,0.30)"])
            link_label.extend(
                [
                    "Reservoir intake → generating asset",
                    "Generating asset → downstream reservoir/outlet",
                ]
            )

        fig.add_trace(
            go.Sankey(
                visible=False,
                arrangement="snap",
                orientation="h",
                node=dict(
                    pad=24,
                    thickness=22,
                    line=dict(color="#52615a", width=0.8),
                    label=[node_labels[key] for key in node_keys],
                    color=[node_colors[key] for key in node_keys],
                    hovertemplate="%{label}<br>Connected capacity: %{value:,.0f} MW<extra></extra>",
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_color,
                    label=link_label,
                    hovertemplate="%{label}<br>Capacity proxy: %{value:,.0f} MW<extra></extra>",
                ),
                name=cascade,
            )
        )
        trace_meta.append(
            (cascade, float(group["capacity"].sum()), len(group), len(node_keys))
        )

    buttons = []
    for i, (cascade, capacity, plants, nodes) in enumerate(trace_meta):
        visible = [False] * len(trace_meta)
        visible[i] = True
        buttons.append(
            dict(
                label=f"{cascade} · {capacity:,.0f} MW",
                method="update",
                args=[
                    {"visible": visible},
                    {
                        "title": (
                            f"Hydro Cascade Asset Flow · {cascade}"
                            f"<br><sup>{plants} generating assets · {capacity:,.0f} MW installed "
                            f"· {nodes} connected reservoir/plant/outlet nodes</sup>"
                        )
                    },
                ],
            )
        )

    first_name, first_capacity, first_plants, first_nodes = trace_meta[0]
    fig.update_layout(
        title=(
            f"Hydro Cascade Asset Flow · {first_name}"
            f"<br><sup>{first_plants} generating assets · {first_capacity:,.0f} MW installed "
            f"· {first_nodes} connected reservoir/plant/outlet nodes</sup>"
        ),
        template="plotly_white",
        paper_bgcolor=OFF_WHITE,
        font=dict(family="Arial, sans-serif", color=INK, size=13),
        height=820,
        margin=dict(l=35, r=35, t=155, b=95),
        updatemenus=[
            dict(
                buttons=buttons,
                x=0,
                y=1.14,
                xanchor="left",
                yanchor="top",
                bgcolor="#fffdf7",
                bordercolor="#7b8b84",
                borderwidth=1,
            )
        ],
        annotations=[
            dict(
                text="<b>Select cascade</b>",
                x=0,
                y=1.25,
                xref="paper",
                yref="paper",
                xanchor="left",
                showarrow=False,
            ),
            dict(
                text=(
                    "Direction represents modeled water routing. Link width represents installed "
                    "generation capacity—not water volume or simulated dispatch. In the combined "
                    "view, repeated reservoir IDs are merged into shared nodes."
                ),
                x=0,
                y=-0.11,
                xref="paper",
                yref="paper",
                xanchor="left",
                showarrow=False,
                font=dict(size=11, color="#59655f"),
            ),
        ],
    )
    output = output_dir / "hydro_cascade_asset_flow.html"
    fig.write_html(output, include_plotlyjs="cdn")
    return output


def build_calibration_explorer(output_dir: Path) -> Path:
    """Interactive monthly load magnitude and shape calibration."""
    model = pd.read_csv("data/processed_data/load/Hourly_profile_2021.csv")
    model["TIME"] = pd.to_datetime(model["TIME"])
    model_series = pd.Series(
        pd.to_numeric(model["LOAD"], errors="coerce").values / 1_000_000.0,
        index=model["TIME"],
    )
    obs_path = Path("data/validation/bc_hydro/BalancingAuthorityLoad2021.xls")
    observed_kwh = _parse_bch_hourly_load_kwh(obs_path, 2021)
    if observed_kwh is None:
        raise ValueError("BC Hydro validation load could not be parsed")
    observed_gw = observed_kwh / 1_000_000.0
    monthly_model = model_series.resample("MS").sum()  # GWh because GW × hour
    monthly_obs = observed_gw.resample("MS").sum()
    aligned = pd.concat([monthly_model.rename("Model"), monthly_obs.rename("Observed")], axis=1).dropna()
    correlation = aligned.corr().iloc[0, 1]
    magnitude_bias = 100 * (aligned["Model"].sum() / aligned["Observed"].sum() - 1)

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Monthly electricity", "Monthly shape calibration"),
        horizontal_spacing=0.12,
    )
    fig.add_trace(
        go.Bar(x=aligned.index, y=aligned["Model"], name="Model", marker_color=BLUE),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=aligned.index,
            y=aligned["Observed"],
            name="BC Hydro BA",
            mode="lines+markers",
            line=dict(color=ORANGE, width=2),
        ),
        row=1,
        col=1,
    )
    low = float(aligned.min().min())
    high = float(aligned.max().max())
    fig.add_trace(
        go.Scatter(
            x=[low, high],
            y=[low, high],
            mode="lines",
            name="1:1",
            line=dict(color="#777", dash="dash"),
            hoverinfo="skip",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=aligned["Observed"],
            y=aligned["Model"],
            mode="markers+text",
            text=aligned.index.strftime("%b"),
            textposition="top center",
            marker=dict(size=10, color=TEAL),
            name="Month",
            hovertemplate="Observed %{x:.0f} GWh<br>Model %{y:.0f} GWh<extra></extra>",
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="Month", row=1, col=1)
    fig.update_yaxes(title_text="Monthly electricity [GWh]", rangemode="tozero", row=1, col=1)
    fig.update_xaxes(title_text="Observed BC Hydro [GWh]", row=1, col=2)
    fig.update_yaxes(title_text="Modeled [GWh]", scaleanchor="x2", scaleratio=1, row=1, col=2)
    fig.update_layout(
        title=(
            "Demand Input Calibration · 2021"
            f"<br><sup>Pearson r = {correlation:.3f} · annual magnitude bias = {magnitude_bias:+.1f}% "
            "· model scope differs from all-sector BA load</sup>"
        ),
        template="plotly_white",
        paper_bgcolor=OFF_WHITE,
        plot_bgcolor=OFF_WHITE,
        font=dict(family="Arial, sans-serif", color=INK),
        height=570,
        margin=dict(l=75, r=35, t=105, b=65),
        legend=dict(orientation="h", y=-0.18),
    )
    output = output_dir / "demand_calibration_2021.html"
    fig.write_html(output, include_plotlyjs="cdn")
    return output


def build_all(output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    outputs.extend(build_network_map_previews(output_dir))
    outputs.extend(build_regional_load_map(output_dir))
    outputs.append(build_coders_parameter_explorer(output_dir))
    outputs.extend(build_basin_map(output_dir))
    outputs.extend(build_cascade_schematic(output_dir))
    outputs.append(build_cascade_asset_flow(output_dir))
    outputs.append(build_calibration_explorer(output_dir))
    manifest = pd.DataFrame(
        [{"artifact": path.name, "path": path.as_posix()} for path in outputs]
    )
    manifest.to_csv(output_dir / "manifest.csv", index=False)
    outputs.append(output_dir / "manifest.csv")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("vis/input_visuals"))
    args = parser.parse_args()
    outputs = build_all(args.output_dir)
    print(f"Generated {len(outputs) - 1} input visualization files")
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
