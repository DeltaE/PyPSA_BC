from __future__ import annotations

import pandas as pd
import pypsa

from pypsa_bc.studies.station_validation import (
    STATION_ENERGY_META_KEY,
    add_named_station_energy_constraints,
    apply_named_station_capacities,
    attach_named_station_energy_calibration,
    compare_named_station_generation,
)


def test_station_capacity_alignment_and_comparison() -> None:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    network.add("Bus", "water")
    network.add("Bus", "electricity")
    network.add(
        "Link",
        "station Discharge Link",
        bus0="water",
        bus1="electricity",
        p_nom=20,
        efficiency=0.5,
    )
    network.add("Link", "station Spill Link", bus0="water", bus1="water", p_nom=20)
    evidence = pd.DataFrame(
        {
            "facility": ["Station"],
            "model_component_type": ["link"],
            "model_component_name": ["station Discharge Link"],
            "capacity_mw": [8.0],
            "energy_gwh": [0.01],
            "source_printed_page": [134],
        }
    )
    audit = apply_named_station_capacities(network, evidence)
    assert audit.iloc[0]["capacity_before_mw"] == 10
    assert network.links.at["station Discharge Link", "p_nom"] == 16
    network.links_t.p0["station Discharge Link"] = [8.0, 8.0]
    network.links_t.p1["station Discharge Link"] = [-4.0, -4.0]
    network.links_t.p0["station Spill Link"] = [16.0, 16.0]
    comparison = compare_named_station_generation(network, evidence)
    assert comparison.iloc[0]["capacity_check"] == "PASS"
    assert comparison.iloc[0]["model_calendar_2021_energy_gwh"] == 0.008
    assert comparison.iloc[0]["turbine_capture_fraction"] == 1 / 3
    assert comparison.iloc[0]["water_use_diagnostic"] == "generation_route_used"


def test_cross_period_station_energy_band_is_explicit_and_solvable() -> None:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    network.add("Bus", "water")
    network.add("Bus", "electricity")
    network.add("Generator", "water supply", bus="water", p_nom=20, marginal_cost=1)
    network.add(
        "Link",
        "station Discharge Link",
        bus0="water",
        bus1="electricity",
        p_nom=10,
        efficiency=1,
    )
    network.add("Load", "load", bus="electricity", p_set=5)
    evidence = pd.DataFrame(
        {
            "facility": ["Station"],
            "model_component_type": ["link"],
            "model_component_name": ["station Discharge Link"],
            "capacity_mw": [10],
            "energy_gwh": [0.01],
            "source_printed_page": [134],
        }
    )
    audit = attach_named_station_energy_calibration(
        network,
        evidence,
        "bc_hydro_fiscal_2021_screening_band",
        relative_tolerance=0.25,
    )

    assert len(audit) == 1
    assert STATION_ENERGY_META_KEY in network.meta
    assert audit.iloc[0]["lower_mwh"] == 7.5
    status, termination = network.optimize(
        solver_name="highs",
        extra_functionality=add_named_station_energy_constraints,
    )
    assert (status, termination) == ("ok", "optimal")
    assert network.links_t.p0["station Discharge Link"].sum() == 10
