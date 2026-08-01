import unittest

import pandas as pd
import pypsa

from pypsa_bc.studies.cascade.storage_uncertainty import (
    apply_storage_uncertainty_case,
    load_storage_uncertainty_protocol,
)


def _protocol(lower=0.0, upper=24.0):
    return {
        "protocol_version": "test",
        "evidence_class": "scenario",
        "cases": {
            "zero_pondage": {
                "hours_of_mean_natural_inflow": lower,
                "description": "test lower",
            },
            "upper_pondage": {
                "hours_of_mean_natural_inflow": upper,
                "description": "test upper",
            },
        },
        "stations": {
            "Station": {
                "upper_reservoir_id": "R",
                "source": ["test"],
                "locator": "unit test",
                "rationale": "test-only structural sensitivity",
            }
        },
    }


def _network():
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=3, freq="h"))
    network.add("Bus", "R Water Bus", carrier="Water")
    network.add("Generator", "R Inflow Generator", bus="R Water Bus", p_nom=1000)
    network.generators_t.p_set["R Inflow Generator"] = [100.0, 200.0, 300.0]
    return network


class StorageUncertaintyTests(unittest.TestCase):
    def test_zero_case_retains_no_store(self):
        network = _network()
        audit = apply_storage_uncertainty_case(network, _protocol(), "zero_pondage")
        self.assertTrue(network.stores.empty)
        self.assertEqual(audit.iloc[0]["controllable_storage_m3"], 0.0)

    def test_upper_case_scales_mean_inflow_once(self):
        network = _network()
        audit = apply_storage_uncertainty_case(network, _protocol(), "upper_pondage")
        expected = 24 * 200.0
        self.assertEqual(network.stores.iloc[0]["e_nom"], expected)
        self.assertEqual(audit.iloc[0]["controllable_storage_m3"], expected)
        self.assertEqual(network.meta["storage_uncertainty_case"], "upper_pondage")

    def test_existing_store_is_rejected(self):
        network = _network()
        network.add("Store", "existing", bus="R Water Bus", e_nom=1.0)
        with self.assertRaisesRegex(ValueError, "already has a Store"):
            apply_storage_uncertainty_case(network, _protocol(), "upper_pondage")

    def test_reversed_bounds_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "reversed"):
            load_storage_uncertainty_protocol(_protocol(lower=24.0, upper=0.0))


if __name__ == "__main__":
    unittest.main()
