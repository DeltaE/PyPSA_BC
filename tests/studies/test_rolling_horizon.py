import pandas as pd
import pypsa
from workflow.scripts.solve_rolling_horizon import _maximum_nodal_residual

from pypsa_bc.studies.rolling_horizon import (
    build_rolling_windows,
    configure_window_network,
    initial_store_state,
)


def _network() -> pypsa.Network:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=10, freq="h")
    network.set_snapshots(snapshots)
    network.add("Bus", "water", carrier="Water")
    network.add("Store", "reservoir", bus="water", e_nom=100.0, e_cyclic=True)
    network.add("Store", "sink", bus="water", e_nom=1_000.0, e_cyclic=False)
    return network


def test_windows_commit_every_snapshot_once_and_add_lookahead() -> None:
    snapshots = pd.date_range("2021-01-01", periods=10, freq="h")
    windows = build_rolling_windows(snapshots, window_hours=4, lookahead_hours=2)

    committed = windows[0].commit_snapshots.append(
        [window.commit_snapshots for window in windows[1:]]
    )
    assert committed.equals(snapshots)
    assert [len(window.solve_snapshots) for window in windows] == [6, 6, 2]
    assert [len(window.commit_snapshots) for window in windows] == [4, 4, 2]


def test_final_closure_window_is_not_left_with_a_short_remainder() -> None:
    snapshots = pd.date_range("2021-01-01", periods=25, freq="h")
    windows = build_rolling_windows(
        snapshots,
        window_hours=6,
        lookahead_hours=2,
        final_window_hours=12,
    )

    committed = windows[0].commit_snapshots.append(
        [window.commit_snapshots for window in windows[1:]]
    )
    assert committed.equals(snapshots)
    assert [len(window.commit_snapshots) for window in windows] == [6, 6, 13]
    assert windows[-1].solve_end == snapshots[-1]


def test_window_configuration_uses_scenario_initial_state_and_final_target() -> None:
    network = _network()
    initial, cyclic = initial_store_state(network, cyclic_fraction=0.5)
    window = build_rolling_windows(
        network.snapshots, window_hours=4, lookahead_hours=2
    )[0]

    configure_window_network(
        network,
        window,
        initial,
        annual_closure_target=initial.loc[cyclic],
    )

    assert not network.stores.e_cyclic.any()
    assert network.stores.at["reservoir", "e_initial"] == 50.0
    assert network.stores.at["sink", "e_initial"] == 0.0
    assert network.stores_t.e_set.at[window.solve_end, "reservoir"] == 50.0
    assert network.stores_t.e_set.iloc[:-1].isna().all().all()


def test_nodal_residual_includes_passive_branch_terminal_flows() -> None:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=2, freq="h")
    network.set_snapshots(snapshots)
    network.add("Bus", "supply")
    network.add("Bus", "demand")
    network.add("Line", "line", bus0="supply", bus1="demand", x=0.1, s_nom=100.0)
    network.add("Generator", "generator", bus="supply", p_nom=20.0, marginal_cost=1.0)
    network.add("Load", "load", bus="demand", p_set=10.0)
    status, termination = network.optimize(solver_name="highs")

    assert status == "ok"
    assert termination == "optimal"
    assert _maximum_nodal_residual(network) < 1e-8
