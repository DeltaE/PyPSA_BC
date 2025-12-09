import pandas as pd
from pathlib import Path
from typing import Optional
import plotly.express as px
from pypsa_bc import utils
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import contextily as ctx
import geopandas as gpd
import matplotlib as mpl
import pypsa
from pypsa_bc import calc

"""
def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years, fill_value=0)
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df

""" 

def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years) #, fill_value=0
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df

def visualize_timeseries(data:pd.DataFrame,
                         year:int,
                         type='area',
                         plot_title:str=None,
                         yaxis_title:str=None,
                         xaxis_title:str=None,
                         legend_title_text:str=None,
                         save_to:Optional[Path]=None,
                         show:Optional[bool]=False):
    data.index = data.index.map(lambda x: x.replace(year=year))
    if type=='area':
        fig = px.area(data, title="Title" if plot_title is None else plot_title)
    if type=='line':
        fig = px.line(data, title="Title" if plot_title is None else plot_title)
    
    fig.update_layout(xaxis_title=None if xaxis_title is None else xaxis_title , 
                      yaxis_title='Variable' if yaxis_title is None else yaxis_title,
                      legend_title_text="" if legend_title_text is None else legend_title_text)
  
    if save_to is not None:
        save_to = Path(save_to)
        save_to.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_to)
        utils.print_update(level=2,message=f'Plot save to :{save_to}')
        
    if show:
        fig.show()
        
    return fig


def plot_inter_region_link_usage(
    year: int,
    pypsa_network: pypsa.Network,
    regional_boundaries_GADM_L2: gpd.GeoDataFrame,
    usage_statistics: str = "usage_avg",
    plot_dpi: int = 300,
    line_cmap: str = "Reds",
    plot_title: str | None = None,
    show: bool = False,
    plot_save_to: str | Path | None = None,
):
    """Plot inter-regional link usage with line width ∝ s_nom and color ∈ [0, 1]."""

    usage_type = "Average" if "avg" in usage_statistics.lower() else "Max"
    title = (
        f"Inter Region Link {usage_type} Usage [Simulated for {year}]"
        if plot_title is None
        else plot_title
    )

    # 1) Get inter-region line usage (GeoDataFrame)
    inter_region_line_with_usage = calc.get_line_usage(
        pypsa_network, regional_boundaries_GADM_L2
    )

    # 2) Prepare region boundary frame and ensure a 'Region' column exists
    boundary = regional_boundaries_GADM_L2.copy()
    if "Region" not in boundary.columns:
        boundary = boundary.reset_index().rename(columns={"index": "Region"})

    boundary["Region_Number"] = range(1, len(boundary) + 1)

    # 3) Project to Web Mercator
    boundary = boundary.to_crs(epsg=3857)
    inter_region_lines = inter_region_line_with_usage.to_crs(epsg=3857)

    # 4) Line widths from s_nom (guard against constant s_nom)
    min_width, max_width = 2, 10
    s_nom = inter_region_lines["s_nom"]

    if s_nom.max() == s_nom.min():
        line_width = (min_width + max_width) / 2.0
    else:
        scaled_width = (s_nom - s_nom.min()) / (s_nom.max() - s_nom.min())
        line_width = min_width + scaled_width * (max_width - min_width)

    inter_region_lines["line_width"] = line_width

    # 5) Fixed colorbar ∈ [0, 1]
    if usage_statistics not in inter_region_lines.columns:
        raise KeyError(
            f"Column '{usage_statistics}' not found in inter_region_lines. "
            f"Available: {list(inter_region_lines.columns)}"
        )

    vmin, vmax = 0.0, 1.0
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # 6) Plot
    fig, ax = plt.subplots(figsize=(12, 8))

    boundary.plot(ax=ax, color="none", edgecolor="k", linewidth=0.3, alpha=0.5)

    inter_region_lines.plot(
        column=usage_statistics,
        cmap=line_cmap,
        linewidth=inter_region_lines["line_width"],
        alpha=1,
        legend=False,
        ax=ax,
        norm=norm,
    )

    # Region labels
    for _, row in boundary.iterrows():
        centroid = row.geometry.centroid
        ax.annotate(
            f"{row['Region_Number']}",
            xy=(centroid.x, centroid.y),
            ha="center",
            fontsize=8,
            color="black",
            bbox=dict(facecolor="none", edgecolor="none"),
        )

    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)
    ax.axis("off")

    # Colorbar 0–1
    sm = plt.cm.ScalarMappable(cmap=line_cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(
        sm,
        ax=ax,
        orientation="horizontal",
        shrink=0.8,
        pad=0.001,
        anchor=(0.5, 0.8),
    )
    cbar.set_label(f"{usage_type} Usage (0–1)", rotation=0, labelpad=5, fontsize=12)
    cbar.outline.set_visible(False)
    cbar.ax.patch.set_alpha(0.7)

    # Optional legend (number → Region name)
    handles = [
        Line2D(
            [0],
            [0],
            linestyle="",
            marker="",
            label=f"{row['Region_Number']}: {row['Region']}",
        )
        for _, row in boundary.iterrows()
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1, 0.98),
        frameon=False,
    )

    ax.set_title(title)
    plt.tight_layout()

    if plot_save_to is not None:
        plot_save_to = Path(plot_save_to)
        plot_save_to.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_save_to, dpi=plot_dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig
def plot_inter_region_link_usage_with_extendables(
    year: int,
    pypsa_network: pypsa.Network,
    regional_boundaries_GADM_L2: gpd.GeoDataFrame,
    usage_statistics: str = "usage_avg",
    plot_dpi: int = 300,
    line_cmap: str = "Reds",
    plot_title: str | None = None,
    show: bool = False,
    plot_save_to: str | Path | None = None,
):
    """
    Left: usage (assumed to be ∈ [0, 1]), width ∝ existing capacity (s_nom).
    Right: width ∝ capacity increase (%_change), colour = SAME usage_statistics (0–1).
    """

    usage_type = "Average" if "avg" in usage_statistics.lower() else "Max"
    title = (
        f"Inter Region Link {usage_type} Usage & Extendables [Simulated for {year}]"
        if plot_title is None
        else plot_title
    )

    # 1) Get inter-region line usage
    inter_region_line_with_usage = calc.get_line_usage(
        pypsa_network, regional_boundaries_GADM_L2
    )

    # 2) % change in s_nom  (from network lines)
    comparison = pypsa_network.lines.loc[:, ["s_nom", "s_nom_opt"]].copy()
    comparison["%_change"] = (
        (comparison["s_nom_opt"] - comparison["s_nom"]) / comparison["s_nom"] * 100
    )
    comparison = comparison.reset_index().rename(columns={"index": "Line"})

    inter_region_line_with_usage = inter_region_line_with_usage.merge(
        comparison[["Line", "%_change"]],
        left_index=True,
        right_on="Line",
        how="left",
    )

    # 3) Regions
    boundary = regional_boundaries_GADM_L2.copy()
    if "Region" not in boundary.columns:
        boundary = boundary.reset_index().rename(columns={"index": "Region"})
    boundary["Region_Number"] = range(1, len(boundary) + 1)

    # 4) Project
    boundary = boundary.to_crs(epsg=3857)
    inter_region_line_with_usage = inter_region_line_with_usage.to_crs(epsg=3857)

    # 5) Line widths (subplot 1) – based on existing capacity s_nom
    min_width, max_width = 2, 10
    s_nom = inter_region_line_with_usage["s_nom"]

    if s_nom.max() == s_nom.min():
        line_width = (min_width + max_width) / 2.0
    else:
        scaled_width = (s_nom - s_nom.min()) / (s_nom.max() - s_nom.min())
        line_width = min_width + scaled_width * (max_width - min_width)

    inter_region_line_with_usage["line_width"] = line_width

    # 6) Usage color norm: fixed 0–1
    if usage_statistics not in inter_region_line_with_usage.columns:
        raise KeyError(
            f"Column '{usage_statistics}' not found in inter_region_line_with_usage. "
            f"Available: {list(inter_region_line_with_usage.columns)}"
        )
    usage_norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

    # 7) Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(20, 8))

    # --- Subplot 1: Usage ---
    ax = axs[0]
    boundary.plot(ax=ax, color="none", edgecolor="k", linewidth=0.3, alpha=0.2)

    inter_region_line_with_usage.plot(
        column=usage_statistics,
        linewidth=inter_region_line_with_usage["line_width"],
        cmap=line_cmap,
        alpha=1,
        legend=False,
        ax=ax,
        norm=usage_norm,
    )

    # Region labels
    for _, row in boundary.iterrows():
        centroid = row.geometry.centroid
        ax.annotate(
            f"{row['Region_Number']}",
            xy=(centroid.x, centroid.y),
            ha="center",
            fontsize=8,
            color="black",
            bbox=dict(facecolor="none", edgecolor="none"),
        )

    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)
    ax.axis("off")
    ax.set_title(f"{usage_type} Line Usage")

    usage_sm = plt.cm.ScalarMappable(cmap=line_cmap, norm=usage_norm)
    usage_sm.set_array([])
    cbar = fig.colorbar(
        usage_sm,
        ax=ax,
        orientation="horizontal",
        shrink=0.8,
        pad=0.01,
    )
    cbar.set_label(f"{usage_type} Usage (0–1)", rotation=0, labelpad=5, fontsize=12)
    cbar.outline.set_visible(False)

    # --- Subplot 2: capacity extensions (width) + SAME usage metric (colour) ---
    ax2 = axs[1]

    change_df = inter_region_line_with_usage.dropna(subset=["%_change"]).copy()
    # positive change only for width
    change_df["%_change_pos"] = change_df["%_change"].clip(lower=0)

    # width ∝ %_change_pos
    if change_df["%_change_pos"].max() > 0:
        w_scaled = change_df["%_change_pos"] / change_df["%_change_pos"].max()
        line_width2 = min_width + w_scaled * (max_width - min_width)
    else:
        line_width2 = (min_width + max_width) / 2.0

    change_df["line_width_ext"] = line_width2

    # colour uses SAME usage_statistics, with fixed 0–1 norm
    if usage_statistics not in change_df.columns:
        raise KeyError(
            f"Column '{usage_statistics}' not found in change_df. "
            f"Available: {list(change_df.columns)}"
        )
    usage_norm_ext = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

    change_df.plot(
        column=usage_statistics,
        cmap=line_cmap,
        linewidth=change_df["line_width_ext"],
        alpha=0.9,
        legend=False,
        ax=ax2,
        norm=usage_norm_ext,
    )

    boundary.plot(ax=ax2, color="none", edgecolor="black", linewidth=0.3, alpha=0.2)
    ctx.add_basemap(ax2, source=ctx.providers.CartoDB.Positron)

    ax2.axis("off")
    ax2.set_title(f"{usage_type} Usage (0–1) & Capacity Extensions (line width)")

    change_sm = plt.cm.ScalarMappable(cmap=line_cmap, norm=usage_norm_ext)
    change_sm.set_array([])
    cbar2 = fig.colorbar(
        change_sm,
        ax=ax2,
        orientation="horizontal",
        shrink=0.8,
        pad=0.01,
    )
    cbar2.set_label(f"{usage_type} Usage (0–1)", rotation=0, labelpad=5, fontsize=12)
    cbar2.outline.set_visible(False)

    # Overall title
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()

    if plot_save_to is not None:
        plot_save_to = Path(plot_save_to)
        plot_save_to.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_save_to, dpi=plot_dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig
