import unittest

import pandas as pd

from pypsa_bc.studies.cascade.schedules import build_release_schedule


class ReleaseScheduleTests(unittest.TestCase):
    def test_seasonal_minimum_overrides_default_and_converts_units(self):
        snapshots = pd.date_range("2021-04-14 23:00", "2021-06-15 00:00", freq="h")
        rules = {
            "route_release_rules": [
                {
                    "link": "spill",
                    "schedule_id": "spring",
                    "source": "wup",
                    "evidence_class": "observed",
                    "default_m3_per_s": 1.5,
                    "seasonal_minima": [
                        {"start": "04-15", "end": "06-14", "minimum_m3_per_s": 3.0}
                    ],
                }
            ]
        }
        result = build_release_schedule(snapshots, rules).set_index("snapshot")
        self.assertEqual(result.loc[pd.Timestamp("2021-04-14 23:00"), "minimum_release_m3_per_hour"], 5400)
        self.assertEqual(result.loc[pd.Timestamp("2021-04-15 00:00"), "minimum_release_m3_per_hour"], 10800)
        self.assertEqual(result.loc[pd.Timestamp("2021-06-14 23:00"), "minimum_release_m3_per_hour"], 10800)
        self.assertEqual(result.loc[pd.Timestamp("2021-06-15 00:00"), "minimum_release_m3_per_hour"], 5400)

    def test_duplicate_route_is_rejected(self):
        snapshots = pd.date_range("2021-01-01", periods=2, freq="h")
        rule = {
            "link": "spill",
            "schedule_id": "one",
            "source": "wup",
            "evidence_class": "observed",
            "default_m3_per_s": 1.0,
        }
        with self.assertRaisesRegex(ValueError, "unique"):
            build_release_schedule(snapshots, {"route_release_rules": [rule, rule]})


if __name__ == "__main__":
    unittest.main()
