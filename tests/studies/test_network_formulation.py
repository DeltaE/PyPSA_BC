import pandas as pd
import pypsa

from pypsa_bc.studies.network_formulation import apply_network_flow_formulation


def _network() -> pypsa.Network:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    network.add("Bus", "a", carrier="AC")
    network.add("Bus", "b", carrier="AC")
    network.add("Line", "a-b", bus0="a", bus1="b", s_nom=125.0, x=0.2, r=0.01)
    return network


def test_transport_formulation_preserves_corridor_and_capacity() -> None:
    network = _network()
    audit = apply_network_flow_formulation(network, "transport")

    assert network.lines.empty
    assert "Transport a-b" in network.links.index
    assert network.links.at["Transport a-b", "p_nom"] == 125.0
    assert network.links.at["Transport a-b", "p_min_pu"] == -1.0
    assert network.links.at["Transport a-b", "efficiency"] == 1.0
    assert audit.loc[0, "source_reactance"] == 0.2
    assert not audit.loc[0, "kvl_enforced"]
    assert network.meta["network_flow_formulation"] == "transport"
    assert "minimum-transfer" in network.meta["transport_flow_reporting_contract"]


def test_dc_formulation_retains_lines_and_records_kvl() -> None:
    network = _network()
    audit = apply_network_flow_formulation(network, "dc_load_flow")

    assert "a-b" in network.lines.index
    assert network.links.empty
    assert audit.loc[0, "modeled_branch"] == "a-b"
    assert audit.loc[0, "kvl_enforced"]
