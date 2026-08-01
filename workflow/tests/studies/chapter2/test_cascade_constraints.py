from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pypsa

from pypsa_bc.studies.cascade.constraints import (
    add_cascade_constraints,
    attach_route_release_schedules,
)


def _network(snapshots: int) -> pypsa.Network:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=snapshots, freq="h"))
    network.add("Carrier", "water")
    network.add("Carrier", "electricity")
    return network


class CascadeDelayTests(unittest.TestCase):
    def test_electricity_is_immediate_and_bus2_water_is_delayed(self):
        network = _network(4)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "electric", carrier="electricity")
        network.add("Bus", "downstream", carrier="water")
        pulse = pd.Series([1.0, 0.0, 0.0, 0.0], index=network.snapshots)
        network.add(
            "Generator",
            "inflow",
            bus="upstream",
            p_nom=10.0,
            p_min_pu=pulse,
            p_max_pu=pulse,
            marginal_cost=0.01,
        )
        network.add(
            "Link",
            "turbine",
            bus0="upstream",
            bus1="electric",
            bus2="downstream",
            p_nom=10.0,
            efficiency=0.5,
            efficiency2=1.0,
            delay2=2,
            cyclic_delay2=False,
            marginal_cost=0.01,
        )
        network.add(
            "Load",
            "electric demand",
            bus="electric",
            p_set=pd.Series([5.0, 0.0, 0.0, 0.0], index=network.snapshots),
        )
        network.add(
            "Load",
            "downstream water use",
            bus="downstream",
            p_set=pd.Series([0.0, 0.0, 10.0, 0.0], index=network.snapshots),
        )
        status, condition = network.optimize(
            solver_name="highs",
            log_to_console=False,
            include_objective_constant=False,
        )
        self.assertEqual((status, condition), ("ok", "optimal"))
        self.assertAlmostEqual(network.links_t.p0.loc[network.snapshots[0], "turbine"], 10.0)
        self.assertAlmostEqual(network.links_t.p1.loc[network.snapshots[0], "turbine"], -5.0)
        self.assertAlmostEqual(network.links_t.p2.loc[network.snapshots[0], "turbine"], 0.0)
        self.assertAlmostEqual(network.links_t.p2.loc[network.snapshots[2], "turbine"], -10.0)


class MinimumReleaseTests(unittest.TestCase):
    def test_turbine_plus_spill_meets_verified_release(self):
        network = _network(2)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "electric", carrier="electricity")
        network.add("Bus", "downstream", carrier="water")
        network.add(
            "Generator",
            "inflow",
            bus="upstream",
            p_nom=10.0,
            marginal_cost=0.01,
        )
        network.add(
            "Store",
            "downstream storage",
            bus="downstream",
            e_nom=100.0,
            e_cyclic=False,
        )
        common = {
            "p_nom": 10.0,
            "hydro_station": "TEST_GSS",
            "minimum_release_m3_per_hour": 5.0,
            "water_output_factor": 1.0,
        }
        network.add(
            "Link",
            "turbine",
            bus0="upstream",
            bus1="electric",
            bus2="downstream",
            efficiency=0.5,
            efficiency2=1.0,
            marginal_cost=1.0,
            hydro_release_role="turbine",
            **common,
        )
        network.add(
            "Link",
            "spill",
            bus0="upstream",
            bus1="downstream",
            efficiency=1.0,
            marginal_cost=0.1,
            hydro_release_role="spill",
            **common,
        )
        status, condition = network.optimize(
            solver_name="highs",
            extra_functionality=add_cascade_constraints,
            log_to_console=False,
            include_objective_constant=False,
        )
        self.assertEqual((status, condition), ("ok", "optimal"))
        combined = network.links_t.p0[["turbine", "spill"]].sum(axis=1)
        self.assertTrue((combined.sub(5.0).abs() < 1e-7).all())
        self.assertTrue((network.links_t.p0["turbine"].abs() < 1e-7).all())
        self.assertTrue((network.links_t.p0["spill"].sub(5.0).abs() < 1e-7).all())


class ScheduledRouteReleaseTests(unittest.TestCase):
    def test_schedule_applies_only_to_explicit_release_route(self):
        network = _network(4)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "electric", carrier="electricity")
        network.add("Bus", "downstream", carrier="water")
        network.add(
            "Generator",
            "inflow",
            bus="upstream",
            p_nom=10.0,
            marginal_cost=0.01,
        )
        network.add(
            "Store",
            "downstream storage",
            bus="downstream",
            e_nom=100.0,
            e_cyclic=False,
        )
        network.add(
            "Link",
            "turbine",
            bus0="upstream",
            bus1="electric",
            bus2="downstream",
            p_nom=10.0,
            efficiency=0.5,
            efficiency2=1.0,
            marginal_cost=1.0,
            hydro_station="TEST_GSS",
            hydro_release_role="turbine",
            water_output_factor=1.0,
            minimum_release_m3_per_hour=0.0,
        )
        network.add(
            "Link",
            "river outlet",
            bus0="upstream",
            bus1="downstream",
            p_nom=10.0,
            efficiency=1.0,
            marginal_cost=0.1,
            hydro_station="TEST_GSS",
            hydro_release_role="environmental_release",
            water_output_factor=1.0,
            minimum_release_m3_per_hour=0.0,
        )
        expected = [1.0, 3.0, 0.0, 2.0]
        schedule = pd.DataFrame(
            {
                "link": "river outlet",
                "snapshot": network.snapshots,
                "minimum_release_m3_per_hour": expected,
                "schedule_id": "test_seasonal_route",
                "source": "analytic_fixture",
                "evidence_class": "observed",
            }
        )
        self.assertEqual(attach_route_release_schedules(network, schedule), 1)
        status, condition = network.optimize(
            solver_name="highs",
            extra_functionality=add_cascade_constraints,
            log_to_console=False,
            include_objective_constant=False,
        )
        self.assertEqual((status, condition), ("ok", "optimal"))
        self.assertTrue(
            (
                network.links_t.p0["river outlet"]
                .sub(pd.Series(expected, index=network.snapshots))
                .abs()
                < 1e-7
            ).all()
        )
        self.assertTrue((network.links_t.p0["turbine"].abs() < 1e-7).all())

    def test_schedule_requires_complete_snapshot_coverage(self):
        network = _network(2)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "downstream", carrier="water")
        network.add(
            "Link",
            "river outlet",
            bus0="upstream",
            bus1="downstream",
            hydro_release_role="environmental_release",
            water_output_factor=1.0,
        )
        incomplete = pd.DataFrame(
            {
                "link": ["river outlet"],
                "snapshot": [network.snapshots[0]],
                "minimum_release_m3_per_hour": [1.0],
                "schedule_id": ["incomplete"],
                "source": ["analytic_fixture"],
                "evidence_class": ["observed"],
            }
        )
        with self.assertRaisesRegex(ValueError, "does not exactly cover"):
            attach_route_release_schedules(network, incomplete)

    def test_later_route_schedule_preserves_existing_routes(self):
        network = _network(2)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "downstream", carrier="water")
        for link in ("observed route", "scenario route"):
            network.add(
                "Link",
                link,
                bus0="upstream",
                bus1="downstream",
                hydro_release_role="environmental_release",
                water_output_factor=1.0,
                minimum_release_m3_per_hour=0.0,
            )
        for link, schedule_id, evidence_class, values in (
            ("observed route", "observed", "observed", [1.0, 2.0]),
            ("scenario route", "scenario", "scenario", [3.0, 4.0]),
        ):
            attach_route_release_schedules(
                network,
                pd.DataFrame(
                    {
                        "link": link,
                        "snapshot": network.snapshots,
                        "minimum_release_m3_per_hour": values,
                        "schedule_id": schedule_id,
                        "source": "analytic_fixture",
                        "evidence_class": evidence_class,
                    }
                ),
            )
        self.assertEqual(
            network.links_t.minimum_release_m3_per_hour.columns.tolist(),
            ["observed route", "scenario route"],
        )
        self.assertEqual(
            network.links_t.minimum_release_m3_per_hour["observed route"].tolist(),
            [1.0, 2.0],
        )

    def test_schedule_survives_netcdf_round_trip(self):
        network = _network(2)
        network.add("Bus", "upstream", carrier="water")
        network.add("Bus", "downstream", carrier="water")
        network.add(
            "Link",
            "river outlet",
            bus0="upstream",
            bus1="downstream",
            hydro_release_role="environmental_release",
            water_output_factor=1.0,
            minimum_release_m3_per_hour=0.0,
        )
        schedule = pd.DataFrame(
            {
                "link": "river outlet",
                "snapshot": network.snapshots.astype(str),
                "minimum_release_m3_per_hour": [1.0, 2.0],
                "schedule_id": "round_trip",
                "source": "analytic_fixture",
                "evidence_class": "observed",
            }
        )
        attach_route_release_schedules(network, schedule)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scheduled-release.nc"
            network.export_to_netcdf(path)
            restored = pypsa.Network(path)
        self.assertEqual(
            restored.links.at["river outlet", "minimum_release_schedule_id"],
            "round_trip",
        )
        self.assertEqual(
            restored.links_t.minimum_release_m3_per_hour["river outlet"].tolist(),
            [1.0, 2.0],
        )


if __name__ == "__main__":
    unittest.main()
