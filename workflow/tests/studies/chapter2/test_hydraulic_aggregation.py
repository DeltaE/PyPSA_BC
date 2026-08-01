from __future__ import annotations

import unittest

import pandas as pd
import pypsa

from pypsa_bc.studies.cascade.representation import aggregate_cascade_representation


def _topology() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cascade_group": ["Test", "Test", "Test"],
            "cascade_order": [1, 2, 3],
            "station_asset_id": ["S1", "S2", "S3"],
            "upper_reservoir_id": ["R1", "R2", "R3"],
        }
    )


def _network() -> pypsa.Network:
    network = pypsa.Network()
    network.set_snapshots(pd.date_range("2021-01-01", periods=2, freq="h"))
    for bus in [
        "R1 Water Bus",
        "R2 Water Bus",
        "R3 Water Bus",
        "S1 Water Bus",
        "S2 Water Bus",
        "S3 Water Bus",
        "Test Water Bus",
        "External Creek Water Bus",
    ]:
        network.add("Bus", bus, carrier="Water")
    for bus in ["E1", "E2", "E3"]:
        network.add("Bus", bus, carrier="AC")

    for number, storage in enumerate([100.0, 200.0, 300.0], start=1):
        network.add(
            "Store",
            f"R{number} Reservoir Store",
            bus=f"R{number} Water Bus",
            e_nom=storage,
            e_initial=storage / 2,
            e_cyclic=True,
        )
        network.add(
            "Generator",
            f"R{number} Inflow Generator",
            bus=f"R{number} Water Bus",
            p_nom=number * 10,
            p_set=[number, number + 1],
        )
        network.add(
            "Link",
            f"S{number} Store Link",
            bus0=f"S{number} Water Bus",
            bus1=f"R{number} Water Bus",
            p_nom=storage,
        )
        network.add(
            "Link",
            f"S{number} Release Link",
            bus0=f"R{number} Water Bus",
            bus1=f"S{number} Water Bus",
            p_nom=storage,
        )

    downstream = ["R2 Water Bus", "R3 Water Bus", "Test Water Bus"]
    for number, destination in enumerate(downstream, start=1):
        network.add(
            "Link",
            f"S{number} Discharge Link",
            bus0=f"S{number} Water Bus",
            bus1=f"E{number}",
            bus2=destination,
            efficiency=0.5,
            efficiency2=1.0,
            p_nom=number * 100.0,
            hydro_station=f"S{number}",
            hydro_release_role="turbine",
            water_output_factor=1.0,
        )
        spill_destination = (
            "External Creek Water Bus" if number == 2 else destination
        )
        network.add(
            "Link",
            f"S{number} Spill Link",
            bus0=f"S{number} Water Bus",
            bus1=spill_destination,
            p_nom=number * 1000.0,
            hydro_station=f"S{number}",
            hydro_release_role="spill",
            water_output_factor=1.0,
        )
    return network


class HydraulicAggregationTests(unittest.TestCase):
    def test_full_cascade_is_unchanged(self):
        network = _network()
        stores = network.stores.copy()
        links = network.links.copy()
        result = aggregate_cascade_representation(network, _topology(), "A_full_cascade")
        self.assertTrue(result.summary.empty)
        pd.testing.assert_frame_equal(network.stores, stores)
        pd.testing.assert_frame_equal(network.links, links)

    def test_two_store_conserves_budgets_and_removes_internal_self_loops(self):
        network = _network()
        inflows_before = network.generators_t.p_set.copy()
        electric_buses_before = network.links.loc[
            ["S1 Discharge Link", "S2 Discharge Link", "S3 Discharge Link"], "bus1"
        ].copy()
        result = aggregate_cascade_representation(network, _topology(), "B_two_store")

        aggregate_stores = network.stores[
            network.stores.index.str.contains("Aggregate Reservoir Store")
        ]
        self.assertEqual(len(aggregate_stores), 2)
        self.assertEqual(aggregate_stores["e_nom"].sum(), 600.0)
        self.assertEqual(aggregate_stores["e_initial"].sum(), 300.0)
        self.assertEqual(set(result.summary["bucket"]), {"Upper", "Lower"})
        self.assertEqual(
            network.links.at["S1 Discharge Link", "bus2"],
            "Test Lower Aggregate Water Bus",
        )
        self.assertEqual(
            network.links.at["S2 Discharge Link", "bus2"],
            "Test Lower Aggregate Water Bus",
        )
        self.assertEqual(network.links.at["S3 Discharge Link", "bus2"], "Test Water Bus")
        self.assertEqual(
            network.links.at["S2 Spill Link", "bus1"], "External Creek Water Bus"
        )
        self.assertFalse(network.links.index.str.endswith(("Store Link", "Release Link")).any())
        self.assertEqual(network.links.loc["S1 Discharge Link", "p_nom"], 100.0)
        self.assertEqual(network.links.loc["S2 Discharge Link", "p_nom"], 200.0)
        self.assertEqual(network.links.loc["S3 Discharge Link", "p_nom"], 300.0)
        pd.testing.assert_series_equal(
            network.links.loc[electric_buses_before.index, "bus1"], electric_buses_before
        )
        pd.testing.assert_frame_equal(network.generators_t.p_set, inflows_before)
        self.assertEqual(network.meta["reservoir_representation"], "B_two_store")

    def test_single_bucket_conserves_storage_and_external_routes(self):
        network = _network()
        aggregate_cascade_representation(network, _topology(), "C_single_bucket")
        aggregate_stores = network.stores[
            network.stores.index.str.contains("Aggregate Reservoir Store")
        ]
        self.assertEqual(len(aggregate_stores), 1)
        self.assertEqual(aggregate_stores.iloc[0]["e_nom"], 600.0)
        for number in (1, 2, 3):
            link = f"S{number} Discharge Link"
            self.assertEqual(network.links.at[link, "bus0"], "Test Single Aggregate Water Bus")
            self.assertEqual(network.links.at[link, "bus2"], "Test Water Bus")
        self.assertEqual(
            network.links.at["S2 Spill Link", "bus1"], "External Creek Water Bus"
        )
        self.assertEqual(network.meta["reservoir_representation"], "C_single_bucket")


if __name__ == "__main__":
    unittest.main()
