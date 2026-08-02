import pandas as pd
import pypsa
import pytest

from pypsa_bc.studies.hydro_dispatch_policy import (
    apply_hydro_dispatch_cost_policy,
)


def _network() -> pypsa.Network:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    network.add("Bus", "water", carrier="Water")
    network.add("Bus", "electricity")
    network.add(
        "Link",
        "turbine",
        bus0="water",
        bus1="electricity",
        efficiency=4,
        p_nom=5,
        marginal_cost=27.52,
        hydro_release_role="turbine",
    )
    return network


def test_uniform_policy_replaces_electric_output_cost() -> None:
    network = _network()
    audit = apply_hydro_dispatch_cost_policy(
        network, "uniform_reservoir_vom", uniform_cost_cad_per_mwh=1.97
    )

    assert network.links.at["turbine", "marginal_cost"] == pytest.approx(7.88)
    assert audit.iloc[0]["source_cost_cad_per_mwh"] == pytest.approx(6.88)
    assert audit.iloc[0]["applied_cost_cad_per_mwh"] == pytest.approx(1.97)
    assert network.meta["hydro_dispatch_cost_policy"] == "uniform_reservoir_vom"


def test_source_policy_preserves_cost() -> None:
    network = _network()
    apply_hydro_dispatch_cost_policy(network, "source_generic_vom")

    assert network.links.at["turbine", "marginal_cost"] == pytest.approx(27.52)
