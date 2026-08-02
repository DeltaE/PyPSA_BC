from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import rasterio
from rasterio.transform import from_origin
from workflow.scripts.validate_model_network import _validate_dynamic_inputs

from pypsa_bc.wind import sample_gwa_wind_speeds


def test_gwa_sampling_uses_raster_georeferencing(tmp_path: Path) -> None:
    raster_path = tmp_path / "wind.tif"
    values = np.array([[1.0, 2.0], [3.0, 4.0]], dtype="float32")
    with rasterio.open(
        raster_path,
        "w",
        driver="GTiff",
        height=2,
        width=2,
        count=1,
        dtype=values.dtype,
        crs="EPSG:4326",
        transform=from_origin(-130.0, 52.0, 1.0, 1.0),
    ) as dataset:
        dataset.write(values, 1)

    assets = pd.DataFrame(
        {
            "asset_id": ["northwest", "southeast"],
            "longitude": [-129.5, -128.5],
            "latitude": [51.5, 50.5],
        }
    )
    speeds = sample_gwa_wind_speeds(assets, raster_path)

    assert speeds.tolist() == [1.0, 4.0]


def test_dynamic_validation_rejects_nan_generator_availability() -> None:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=2, freq="h")
    network.set_snapshots(snapshots)
    network.add("Bus", "bus")
    network.add(
        "Generator",
        "wind",
        bus="bus",
        p_nom=10.0,
        p_max_pu=pd.Series([0.5, np.nan], index=snapshots),
    )

    errors = _validate_dynamic_inputs(network)

    assert any("generators_t.p_max_pu" in error and "wind" in error for error in errors)


def test_dynamic_validation_rejects_fixed_dispatch_outside_bounds() -> None:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=2, freq="h")
    network.set_snapshots(snapshots)
    network.add("Bus", "external")
    network.add(
        "Generator",
        "trade",
        bus="external",
        p_nom=100.0,
        p_min_pu=0.0,
        p_max_pu=1.0,
        p_set=pd.Series([-20.0, 10.0], index=snapshots),
    )

    errors = _validate_dynamic_inputs(network)

    assert any("p_set lies outside dispatch bounds" in error and "trade" in error for error in errors)
