from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pypsa

from pypsa_bc.studies.external_boundary import apply_external_boundary


class ExternalBoundaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.snapshots = pd.date_range("2021-01-01", periods=3, freq="h")
        self._write_inputs()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _write_inputs(self) -> None:
        dates = ["01/01/2021"] * 3
        hours = [1, 1, 3]  # Deliberate duplicate/missing label; row order is valid.
        ttc = {
            "ABBCTTC2021.xlsx": [100.0, 90.0, 80.0],
            "BCABTTC2021.xlsx": [70.0, 60.0, 50.0],
            "USBCTTC2021.xlsx": [300.0, 290.0, 280.0],
            "BCUSTTC2021.xlsx": [250.0, 240.0, 230.0],
        }
        for filename, values in ttc.items():
            pd.DataFrame(
                {"Date": dates, "HE": hours, "Path": "test", "TTC": values}
            ).to_excel(self.root / filename, index=False)
        actual = pd.DataFrame(
            {
                "Date": dates,
                "HE": hours,
                "US Tielines": [25.0, -30.0, 10.0],
                "AB Tielines": [-20.0, 15.0, 5.0],
            }
        )
        with pd.ExcelWriter(self.root / "HourlyTielineData2021.xlsx") as writer:
            actual.to_excel(writer, index=False, startrow=1)

    def _network(self) -> pypsa.Network:
        network = pypsa.Network()
        network.set_snapshots(self.snapshots)
        network.add("Carrier", "AC")
        for bus in [
            "BC West",
            "BC East",
            "138_ABBC01_IPT",
            "500_BCUS01_INT",
        ]:
            network.add("Bus", bus, carrier="AC")
        network.add(
            "Line",
            "AB intertie",
            bus0="BC East",
            bus1="138_ABBC01_IPT",
            x=0.1,
            s_nom=40.0,
        )
        network.add(
            "Line",
            "US intertie",
            bus0="BC West",
            bus1="500_BCUS01_INT",
            x=0.1,
            s_nom=40.0,
        )
        return network

    def test_observed_interchange_uses_publisher_sign_convention(self) -> None:
        network = self._network()
        audit = apply_external_boundary(
            network, self.root, 2021, "observed_interchange"
        )

        self.assertNotIn("138_ABBC01_IPT", network.buses.index)
        self.assertNotIn("500_BCUS01_INT", network.buses.index)
        self.assertEqual(network.lines.at["AB intertie", "bus1"], "AB Trade")
        self.assertEqual(network.lines.at["US intertie", "bus1"], "US Trade")
        pd.testing.assert_series_equal(
            network.generators_t.p_set["AB External Market"],
            pd.Series([20.0, -15.0, -5.0], index=self.snapshots),
            check_names=False,
        )
        self.assertEqual(network.generators.at["AB External Market", "p_min_pu"], -1.0)
        self.assertEqual(network.generators.at["AB External Market", "p_max_pu"], 1.0)
        self.assertEqual(audit["duplicate_hour_labels"].tolist(), [1, 1])
        self.assertTrue(audit["maximum_import_ttc_mw"].isna().all())
        self.assertTrue(audit["observed_hours_outside_ttc"].isna().all())
        self.assertTrue((audit["evidence_class"] == "observed").all())

    def test_hourly_ttc_sets_directional_market_bounds(self) -> None:
        network = self._network()
        apply_external_boundary(network, self.root, 2021, "hourly_ttc")

        p_nom = network.generators.at["AB External Market", "p_nom"]
        self.assertEqual(p_nom, 100.0)
        pd.testing.assert_series_equal(
            network.generators_t.p_max_pu["AB External Market"],
            pd.Series([1.0, 0.9, 0.8], index=self.snapshots),
            check_names=False,
        )
        pd.testing.assert_series_equal(
            network.generators_t.p_min_pu["AB External Market"],
            pd.Series([-0.7, -0.6, -0.5], index=self.snapshots),
            check_names=False,
        )

    def test_closed_policy_does_not_mutate_network(self) -> None:
        network = self._network()
        before_buses = network.buses.index.tolist()
        audit = apply_external_boundary(network, self.root, 2021, "closed")
        self.assertEqual(before_buses, network.buses.index.tolist())
        self.assertEqual(audit.iloc[0]["evidence_class"], "counterfactual")


if __name__ == "__main__":
    unittest.main()
