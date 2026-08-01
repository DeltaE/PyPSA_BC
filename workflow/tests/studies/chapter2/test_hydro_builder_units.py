from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.modules.setdefault("atlite", types.ModuleType("atlite"))

from workflow.scripts.create_hydro_assets import (
    apply_hydraulic_route_overrides,
    custom_bridge_agg,
)
from workflow.scripts.enrich_format_hydro import (
    get_reservoir_dict,
    get_ror_water_dict,
)


class ReservoirBuilderUnitTests(unittest.TestCase):
    def test_inflow_nominal_capacity_is_not_converted_twice(self):
        source = Path("workflow/scripts/enrich_format_hydro.py").read_text(encoding="utf-8")
        compact = source.replace(" ", "")
        self.assertIn("max_inflow=max(inflow)", compact)
        self.assertNotIn("max_inflow=max(inflow)*3600", compact)

    @patch("workflow.scripts.enrich_format_hydro.utils.get_gen_bus", return_value="ELC")
    def test_terminal_ror_water_station_creates_sink_and_delays_only_water(self, _):
        site = pd.Series(
            {
                "asset_id": "BC_WAN_GSS",
                "upper_reservoir_id": "BC_WAN_RES",
                "lower_reservoir_id": "DEFAULT",
                "cascade_group": "Seven Mile",
                "connecting_node_code": "NODE",
                "capacity": 480.0,
                "max_water_discharge": 838.0,
                "max_spill": "DEFAULT",
                "variable_om_costs": 0.0,
                "travel_time_h": 2,
                "minimum_release_m3_per_s": pd.NA,
            }
        )
        result = get_ror_water_dict(
            site,
            pd.Series([100.0, 200.0]),
            {},
            ["BC_WAN_RES", "DEFAULT"],
        )
        self.assertEqual(result["terminal bus"]["name"], "Seven Mile Water Bus")
        self.assertFalse(result["terminal store"]["e_cyclic"])
        self.assertEqual(result["discharge link"]["delay2"], 2)
        self.assertNotIn("delay", result["discharge link"])
        self.assertEqual(result["discharge link"]["efficiency2"], 1.0)
        self.assertEqual(result["discharge link"]["water_output_factor"], 1.0)
        self.assertAlmostEqual(
            result["discharge link"]["efficiency"],
            480.0 / (838.0 * 3600.0),
        )
        self.assertEqual(result["discharge link"]["p_nom"], 838.0 * 3600.0)
        self.assertEqual(result["spill link"]["delay"], 2)
        self.assertEqual(result["spill link"]["efficiency"], 1.0)
        self.assertGreater(result["spill link"]["p_nom"], 0)

    @patch("workflow.scripts.enrich_format_hydro.utils.get_gen_bus", return_value="ELC")
    def test_ror_water_headpond_retains_documented_live_storage(self, _):
        site = pd.Series(
            {
                "asset_id": "BC_WHN_GSS",
                "upper_reservoir_id": "BC_WHN_RES",
                "lower_reservoir_id": "BC_ALH_RES",
                "cascade_group": "Mica/Columbia",
                "connecting_node_code": "NODE",
                "capacity": 8.0,
                "max_water_discharge": 4.33,
                "max_spill": 11.0,
                "variable_om_costs": 0.0,
                "travel_time_h": 0,
                "minimum_release_m3_per_s": pd.NA,
            }
        )
        reservoir = pd.Series({"min_storage": 0.0, "max_storage": 459350.0})
        result = get_ror_water_dict(
            site,
            pd.Series([100.0, 200.0]),
            {},
            ["BC_WHN_RES", "BC_ALH_RES"],
            reservoir=reservoir,
        )
        self.assertEqual(result["reservoir store"]["e_nom"], 459350.0)
        self.assertTrue(result["reservoir store"]["e_cyclic"])
        self.assertEqual(result["reservoir store"]["bus"], "BC_WHN_RES Water Bus")

    @patch("workflow.scripts.enrich_format_hydro.utils.get_gen_bus", return_value="ELC")
    def test_turbine_and_dam_release_can_have_different_downstream_routes(self, _):
        site = pd.Series(
            {
                "asset_id": "BC_BR1_GSS",
                "upper_reservoir_id": "BC_BR1_RES",
                "lower_reservoir_id": "BC_SON_RES",
                "spill_lower_reservoir_id": "BC_LOWER_BRIDGE_RIVER_SINK",
                "cascade_group": "Bridge",
                "connecting_node_code": "NODE",
                "capacity": 508.0,
                "max_water_discharge": 160.0,
                "max_spill": 1171.2,
                "variable_om_costs": 0.0,
                "travel_time_h": 0,
                "minimum_release_m3_per_s": pd.NA,
            }
        )
        reservoir = pd.Series(
            {
                "max_storage": 1_012_500_000.0,
                "min_storage": 0.0,
            }
        )
        result = get_reservoir_dict(
            site,
            reservoir,
            pd.Series([100.0, 200.0]),
            ["BC_BR1_RES", "BC_SON_RES"],
            {},
        )
        self.assertEqual(
            result["discharge link"]["bus2"],
            "BC_SON_RES Water Bus",
        )
        self.assertEqual(
            result["spill link"]["bus1"],
            "BC_LOWER_BRIDGE_RIVER_SINK Water Bus",
        )
        self.assertFalse(result["spill terminal store"]["e_cyclic"])


class HydroTopologyPreparationTests(unittest.TestCase):
    def test_route_override_restores_lajoie_and_splits_bridge_routes(self):
        hydro = pd.DataFrame(
            {
                "asset_id": ["BC_LAJ_GSS", "BC_BR1_GSS"],
                "lower_reservoir_id": ["DEFAULT", "BC_SON_RES"],
                "spill_lower_reservoir_id": ["DEFAULT", "DEFAULT"],
                "hydraulic_route_source": ["", ""],
                "hydraulic_route_locator": ["", ""],
                "hydraulic_route_note": ["", ""],
            }
        )
        overrides = pd.DataFrame(
            {
                "asset_id": ["BC_LAJ_GSS", "BC_BR1_GSS"],
                "turbine_lower_reservoir_id": ["BC_BR1_RES", "BC_SON_RES"],
                "spill_lower_reservoir_id": [
                    "BC_BR1_RES",
                    "BC_LOWER_BRIDGE_RIVER_SINK",
                ],
                "source": ["bridge_wup_2011", "bridge_wup_2011"],
                "locator": ["plan page 6", "plan page 7"],
                "note": ["Middle Bridge River", "Terzaghi route"],
            }
        )
        result = apply_hydraulic_route_overrides(hydro, overrides).set_index("asset_id")
        self.assertEqual(
            result.at["BC_LAJ_GSS", "lower_reservoir_id"],
            "BC_BR1_RES",
        )
        self.assertEqual(
            result.at["BC_BR1_GSS", "spill_lower_reservoir_id"],
            "BC_LOWER_BRIDGE_RIVER_SINK",
        )
        self.assertEqual(
            result.at["BC_BR1_GSS", "hydraulic_route_source"],
            "bridge_wup_2011",
        )

    def test_parallel_bridge_turbines_share_one_dam_release_capacity(self):
        frame = pd.DataFrame(
            {
                "capacity": [208.0, 300.0],
                "annual_avg_energy": [597.6, 794.88],
                "max_water_discharge": [65.0, 95.0],
                "max_spill": [1171.2, 1171.2],
                "num_of_units": [4, 4],
                "ramp_up": [10.0, 12.0],
                "ramp_down": [10.0, 12.0],
            },
            index=["BC_BR1_GSS", "BC_BR2_GSS"],
        )
        result = custom_bridge_agg(frame)
        self.assertEqual(result.index.tolist(), ["BC_BR1_GSS"])
        self.assertEqual(result.at["BC_BR1_GSS", "capacity"], 508.0)
        self.assertEqual(
            result.at["BC_BR1_GSS", "max_water_discharge"],
            160.0,
        )
        self.assertEqual(result.at["BC_BR1_GSS", "max_spill"], 1171.2)


if __name__ == "__main__":
    unittest.main()
