import pandas as pd
import pypsa

from pypsa_bc.studies.rolling_diagnostics import (
    CommittedWindow,
    cyclic_store_closure_audit,
    dispatch_timeseries,
    external_interchange_audit,
    line_utilization_audit,
    minimum_transfer_utilization_audit,
    release_compliance_audit,
    reservoir_state_audit,
)


def _diagnostic_network() -> pypsa.Network:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=2, freq="h")
    network.set_snapshots(snapshots)
    network.add("Bus", "ac", carrier="AC")
    network.add("Bus", "water", carrier="Water")
    network.add("Load", "load", bus="ac", p_set=[10.0, 10.0])
    network.add(
        "Generator",
        "US External Market",
        bus="ac",
        carrier="external_trade",
        p_nom=10.0,
        p_min_pu=-1.0,
    )
    network.add("Generator", "Backstop Region", bus="ac", p_nom=10.0)
    network.add("Generator", "Plant RoR Generator", bus="ac", p_nom=10.0)
    network.add("Link", "Turbine", bus0="water", bus1="ac", p_nom=10.0, carrier="Water")
    network.links.loc["Turbine", "hydro_release_role"] = "turbine"
    network.links.loc["Turbine", "hydro_station"] = "Station"
    network.links.loc["Turbine", "water_output_factor"] = 1.0
    network.links.loc["Turbine", "minimum_release_m3_per_hour"] = 4.0
    network.links.loc["Turbine", "minimum_release_schedule_id"] = ""
    network.links.loc["Turbine", "minimum_release_evidence_class"] = ""
    network.add("Store", "Reservoir", bus="water", e_nom=20.0, e_cyclic=False)

    network.generators_t.p = pd.DataFrame(
        {
            "US External Market": [-2.0, 3.0],
            "Backstop Region": [1.0, 0.0],
            "Plant RoR Generator": [4.0, 4.0],
        },
        index=snapshots,
    )
    network.generators_t.p_set = pd.DataFrame(
        {"US External Market": [-2.0, 3.0]}, index=snapshots
    )
    network.loads_t.p = pd.DataFrame({"load": [10.0, 10.0]}, index=snapshots)
    network.links_t.p0 = pd.DataFrame({"Turbine": [4.0, 5.0]}, index=snapshots)
    network.links_t.p1 = pd.DataFrame({"Turbine": [-4.0, -5.0]}, index=snapshots)
    network.stores_t.e = pd.DataFrame({"Reservoir": [8.0, 9.0]}, index=snapshots)
    network.meta["water_volume_unit_m3"] = 1.0
    return network


def test_dispatch_and_external_boundary_use_committed_sign_convention() -> None:
    network = _diagnostic_network()
    window = CommittedWindow(1, network, pd.DatetimeIndex(network.snapshots))

    dispatch = dispatch_timeseries([window])
    external = external_interchange_audit([window])

    assert dispatch["exports_mw"].tolist() == [2.0, 0.0]
    assert dispatch["imports_mw"].tolist() == [0.0, 3.0]
    assert dispatch["backstop_mw"].sum() == 1.0
    assert dispatch["reservoir_hydro_mw"].tolist() == [4.0, 5.0]
    assert external["maximum_absolute_error_mw"].max() == 0.0


def test_release_and_physical_reservoir_audits_preserve_units() -> None:
    network = _diagnostic_network()
    window = CommittedWindow(1, network, pd.DatetimeIndex(network.snapshots))
    source = _diagnostic_network()
    source.stores.loc["Reservoir", "e_cyclic"] = True

    releases = release_compliance_audit([window])
    trajectory, bounds = reservoir_state_audit([window], source)

    assert releases.loc[0, "minimum_slack_m3_per_hour"] == 0.0
    assert releases.loc[0, "violation_hours"] == 0
    assert trajectory["storage_m3"].tolist() == [8.0, 9.0]
    assert bounds.loc[0, "minimum_storage_fraction"] == 0.4
    assert not bounds.loc[0, "lower_bound_violation"]
    assert not bounds.loc[0, "upper_bound_violation"]


def test_cyclic_closure_uses_first_hour_dispatch_not_equal_end_states() -> None:
    network = _diagnostic_network()
    network.stores.loc["Reservoir", "e_cyclic"] = True
    network.stores_t.p = pd.DataFrame(
        {"Reservoir": [1.0, -1.0]}, index=network.snapshots
    )

    closure = cyclic_store_closure_audit(network)

    assert closure.loc[0, "first_storage_m3"] == 8.0
    assert closure.loc[0, "last_storage_m3"] == 9.0
    assert closure.loc[0, "closure_residual_m3"] == 0.0
    assert closure.loc[0, "normalized_closure_residual"] == 0.0


def test_directional_transport_utilization_uses_net_corridor_flow() -> None:
    network = _diagnostic_network()
    network.add("Bus", "ac2", carrier="AC")
    for direction, bus0, bus1 in (
        ("forward", "ac", "ac2"),
        ("reverse", "ac2", "ac"),
    ):
        name = f"Transport corridor {direction}"
        network.add(
            "Link",
            name,
            bus0=bus0,
            bus1=bus1,
            carrier="electricity_transport",
            p_nom=10.0,
        )
        network.links.at[name, "source_line_name"] = "corridor"
        network.links.at[name, "transport_direction"] = direction
    network.links_t.p0.loc[:, "Transport corridor forward"] = [6.0, 4.0]
    network.links_t.p0.loc[:, "Transport corridor reverse"] = [1.0, 0.0]
    window = CommittedWindow(1, network, pd.DatetimeIndex(network.snapshots))

    utilization = line_utilization_audit([window])

    assert utilization.loc[0, "line"] == "corridor"
    assert utilization.loc[0, "maximum_utilization_pu"] == 0.5


def test_minimum_transfer_postprocessing_removes_zero_cost_cycle() -> None:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=1, freq="h"))
    for bus in ("a", "b", "c"):
        network.add("Bus", bus, carrier="AC")
    for label, bus0, bus1 in (
        ("a-b", "a", "b"),
        ("b-c", "b", "c"),
        ("c-a", "c", "a"),
    ):
        name = f"Transport {label}"
        network.add(
            "Link",
            name,
            bus0=bus0,
            bus1=bus1,
            carrier="electricity_transport",
            p_nom=10.0,
            p_min_pu=-1.0,
        )
        network.links.at[name, "source_line_name"] = label
    network.links_t.p0 = pd.DataFrame(
        [[6.0, 5.0, 5.0]],
        index=network.snapshots,
        columns=["Transport a-b", "Transport b-c", "Transport c-a"],
    )

    corridors, hourly = minimum_transfer_utilization_audit(network)

    assert hourly.loc[0, "raw_absolute_transfer_mw"] == 16.0
    assert hourly.loc[0, "identified_absolute_transfer_mw"] == 1.0
    assert hourly.loc[0, "circulating_transfer_removed_mw"] == 15.0
    assert hourly.loc[0, "maximum_nodal_divergence_residual_mw"] < 1e-9
    assert corridors["maximum_utilization_pu"].max() == 0.1
