import unittest

import pandas as pd
import pypsa

from pypsa_bc.studies.cascade.uncertainty import (
    apply_release_uncertainty_case,
    load_release_uncertainty_protocol,
)


def _protocol(stations):
    return {
        "protocol_version": "test",
        "evidence_class": "scenario",
        "stations": stations,
    }


def _definition(scope, lower, upper, link=None):
    result = {
        "upper_reservoir_id": "R",
        "constraint_scope": scope,
        "source": "test_protocol",
        "locator": "unit test",
        "rationale": "Test-only structural sensitivity.",
        "lower_bound": lower,
        "upper_bound": upper,
    }
    if link:
        result["link"] = link
    return result


def _network():
    network = pypsa.Network()
    snapshots = pd.date_range("2021-01-01", periods=3, freq="h")
    network.set_snapshots(snapshots)
    for bus in ("water", "power", "downstream"):
        network.add("Bus", bus)
    network.add("Generator", "R Inflow Generator", bus="water", p_nom=1000)
    network.generators_t.p_set["R Inflow Generator"] = [100.0, 200.0, 300.0]
    network.add("Link", "Turbine", bus0="water", bus1="power", p_nom=1000)
    network.add("Link", "Spill", bus0="water", bus1="downstream", p_nom=1000)
    network.links["hydro_station"] = "Station"
    network.links["hydro_release_role"] = ["turbine", "spill"]
    network.links["water_output_factor"] = 1.0
    network.links["minimum_release_m3_per_hour"] = 0.0
    return network


class ReleaseUncertaintyTests(unittest.TestCase):
    def test_absolute_combined_release_is_converted_once(self):
        network = _network()
        protocol = _protocol(
            {
                "Station": _definition(
                    "combined_station_release",
                    {"method": "absolute_m3_per_s", "value": 5.7},
                    {"method": "absolute_m3_per_s", "value": 18.4},
                )
            }
        )

        audit = apply_release_uncertainty_case(network, protocol, "upper_bound")

        expected = 18.4 * 3600
        self.assertTrue((network.links["minimum_release_m3_per_hour"] == expected).all())
        self.assertEqual(audit.iloc[0]["minimum_release_m3_per_hour"], expected)
        self.assertEqual(network.meta["release_uncertainty_case"], "upper_bound")

    def test_fraction_uses_mean_natural_inflow(self):
        network = _network()
        protocol = _protocol(
            {
                "Station": _definition(
                    "combined_station_release",
                    {"method": "fraction_of_mean_natural_inflow", "value": 0.0},
                    {"method": "fraction_of_mean_natural_inflow", "value": 0.1},
                )
            }
        )

        apply_release_uncertainty_case(network, protocol, "upper_bound")

        self.assertTrue((network.links["minimum_release_m3_per_hour"] == 20.0).all())

    def test_named_route_schedule_does_not_constrain_turbine(self):
        network = _network()
        network.links.loc[:, "minimum_release_m3_per_hour"] = 9.0
        protocol = _protocol(
            {
                "Station": _definition(
                    "named_route",
                    {"method": "fraction_of_mean_natural_inflow", "value": 0.0},
                    {"method": "fraction_of_mean_natural_inflow", "value": 0.1},
                    link="Spill",
                )
            }
        )

        apply_release_uncertainty_case(network, protocol, "upper_bound")

        self.assertEqual(network.links.at["Turbine", "minimum_release_m3_per_hour"], 0.0)
        self.assertEqual(network.links.at["Spill", "minimum_release_m3_per_hour"], 0.0)
        self.assertEqual(network.links.at["Spill", "minimum_release_evidence_class"], "scenario")
        self.assertNotEqual(network.links.at["Spill", "minimum_release_schedule_id"], "")
        self.assertEqual(
            network.links_t.minimum_release_m3_per_hour["Spill"].tolist(),
            [20.0, 20.0, 20.0],
        )
        self.assertNotIn("Turbine", network.links_t.minimum_release_m3_per_hour.columns)

    def test_reversed_bounds_are_rejected(self):
        protocol = _protocol(
            {
                "Station": _definition(
                    "combined_station_release",
                    {"method": "absolute_m3_per_s", "value": 2.0},
                    {"method": "absolute_m3_per_s", "value": 1.0},
                )
            }
        )
        with self.assertRaisesRegex(ValueError, "reversed"):
            load_release_uncertainty_protocol(protocol)


if __name__ == "__main__":
    unittest.main()
