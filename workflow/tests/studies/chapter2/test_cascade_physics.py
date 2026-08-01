from __future__ import annotations

import unittest

import pandas as pd

from pypsa_bc.studies.cascade.physics import (
    coerce_delay_hours,
    coerce_minimum_release_m3_per_hour,
    merge_cascade_evidence,
)


class HydraulicParameterTests(unittest.TestCase):
    def test_integer_delay_is_required(self):
        self.assertEqual(coerce_delay_hours("3"), 3)
        self.assertEqual(coerce_delay_hours(pd.NA), 0)
        with self.assertRaises(ValueError):
            coerce_delay_hours(1.5)
        with self.assertRaises(ValueError):
            coerce_delay_hours(-1)

    def test_minimum_release_is_converted_once(self):
        self.assertEqual(coerce_minimum_release_m3_per_hour(2.5), 9000.0)
        self.assertEqual(coerce_minimum_release_m3_per_hour("DEFAULT"), 0.0)

    def test_evidence_merge_updates_only_verified_release(self):
        generation = pd.DataFrame(
            {
                "asset_id": ["A", "B"],
                "min_water_discharge": pd.Series(["1.0", "2.0"], dtype="str"),
            }
        )
        register = pd.DataFrame(
            {
                "station_asset_id": ["A", "B"],
                "travel_time_h": [2, pd.NA],
                "travel_time_evidence_class": ["derived", "unknown"],
                "travel_time_source": ["distance/velocity", ""],
                "minimum_release_m3_per_s": [3.0, pd.NA],
                "minimum_release_evidence_class": ["observed", "unknown"],
                "minimum_release_source": ["WUP", ""],
                "resolution_status": ["resolved", "evidence_required"],
            }
        )
        result = merge_cascade_evidence(generation, register).set_index("asset_id")
        self.assertEqual(result.loc["A", "min_water_discharge"], 3.0)
        self.assertEqual(result.loc["B", "min_water_discharge"], "2.0")
        self.assertEqual(result.loc["A", "travel_time_h"], 2)


if __name__ == "__main__":
    unittest.main()
