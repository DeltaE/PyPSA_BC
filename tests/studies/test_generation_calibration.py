from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pypsa

from pypsa_bc.studies.generation_calibration import (
    META_KEY,
    add_monthly_generation_calibration_constraints,
    apply_ror_dispatch_policy,
    attach_monthly_generation_calibration,
)


def _network() -> pypsa.Network:
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", "2021-12-31 23:00", freq="h")
    network.set_snapshots(snapshots)
    for bus in ("fuel", "electricity"):
        network.add("Bus", bus)
    network.add("Carrier", "Biomass")
    network.add(
        "Link",
        "bio",
        bus0="fuel",
        bus1="electricity",
        carrier="Biomass",
        p_nom=20,
        efficiency=0.5,
    )
    network.add("Generator", "fuel", bus="fuel", p_nom=100, marginal_cost=1)
    network.add(
        "Generator", "backstop", bus="electricity", p_nom=100, marginal_cost=100
    )
    network.add("Load", "load", bus="electricity", p_set=5)
    return network


def test_ror_fixed_profile_sets_matching_lower_bound() -> None:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    network.add("Bus", "electricity")
    network.add(
        "Generator",
        "ror",
        bus="electricity",
        type="RoR",
        p_nom=10,
        p_max_pu=[0.2, 0.4],
    )
    audit = apply_ror_dispatch_policy(network, "fixed_source_profile")
    pd.testing.assert_series_equal(
        network.generators_t.p_min_pu["ror"],
        network.generators_t.p_max_pu["ror"],
        check_names=False,
    )
    assert audit.iloc[0]["annual_profile_energy_mwh"] == 6


def test_monthly_constraint_matches_target(tmp_path: Path) -> None:
    network = _network()
    rows = []
    for month in range(1, 13):
        for category in ("biomass", "nonrenewable_combustible"):
            rows.append(
                {"month": f"2021-{month:02d}", "category": category, "observed_mwh": 120}
            )
    observed = tmp_path / "observed.csv"
    pd.DataFrame(rows).to_csv(observed, index=False)
    network.add("Carrier", "NG")
    network.add(
        "Link", "gas", bus0="fuel", bus1="electricity", carrier="NG", p_nom=20,
        efficiency=0.5,
    )
    audit = attach_monthly_generation_calibration(
        network, observed, "statcan_monthly_nonhydro"
    )
    assert len(audit) == 24
    assert len(json.loads(network.meta[META_KEY])) == 24
    status, termination = network.optimize(
        solver_name="highs",
        extra_functionality=add_monthly_generation_calibration_constraints,
    )
    assert (status, termination) == ("ok", "optimal")
    biomass_output = float((network.links_t.p0["bio"] * 0.5).sum())
    assert abs(biomass_output - 1440) < 1e-7


def test_archived_statcan_wide_table_is_accepted(tmp_path: Path) -> None:
    network = _network()
    network.add("Carrier", "NG")
    network.add(
        "Link", "gas", bus0="fuel", bus1="electricity", carrier="NG", p_nom=20,
        efficiency=0.5,
    )
    observed = tmp_path / "observed-wide.csv"
    pd.DataFrame(
        {
            "month": [f"2021-{month:02d}" for month in range(1, 13)],
            "biomass_mwh": [120] * 12,
            "nonrenewable_combustible_mwh": [80] * 12,
            "hydro_mwh": [1000] * 12,
        }
    ).to_csv(observed, index=False)

    audit = attach_monthly_generation_calibration(
        network, observed, "statcan_monthly_nonhydro"
    )

    assert len(audit) == 24
    assert set(audit["category"]) == {"biomass", "nonrenewable_combustible"}
