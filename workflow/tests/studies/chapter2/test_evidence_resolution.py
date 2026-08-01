from __future__ import annotations

import unittest

import pandas as pd
from workflow.scripts.studies.chapter2.resolve_cascade_evidence import resolve


class EvidenceResolutionTests(unittest.TestCase):
    def test_scenario_delay_is_not_mislabelled_as_observed(self):
        register = pd.DataFrame(
            {
                "station_asset_id": ["A"],
                "minimum_release_evidence_class": ["unknown"],
                "resolution_status": ["evidence_required"],
            }
        )
        decisions = {
            "travel_time_default": {
                "travel_time_h": 0,
                "travel_time_evidence_class": "scenario",
                "travel_time_source": "method_precedent",
                "travel_time_low_h": 0,
                "travel_time_high_h": 12,
                "travel_time_method": "structural sensitivity",
                "travel_time_sensitivity_cases": "0|1|3|6|12",
            },
            "minimum_release_decisions": {},
        }
        result = resolve(register, decisions).iloc[0]
        self.assertEqual(result["travel_time_h"], 0)
        self.assertEqual(result["travel_time_evidence_class"], "scenario")
        self.assertEqual(result["resolution_status"], "minimum_release_evidence_required")

    def test_schedule_rule_remains_open_until_implemented(self):
        register = pd.DataFrame(
            {
                "station_asset_id": ["A"],
                "minimum_release_evidence_class": ["unknown"],
                "resolution_status": ["evidence_required"],
            }
        )
        decisions = {
            "travel_time_default": {
                "travel_time_h": 0,
                "travel_time_evidence_class": "scenario",
                "travel_time_source": "method_precedent",
            },
            "minimum_release_decisions": {
                "A": {
                    "minimum_release_evidence_class": "observed",
                    "minimum_release_source": "accepted_plan",
                    "minimum_release_rule_type": "seasonal_route_specific",
                    "minimum_release_schedule_id": "schedule_a",
                }
            },
        }
        result = resolve(register, decisions).iloc[0]
        self.assertEqual(result["minimum_release_schedule_id"], "schedule_a")
        self.assertEqual(result["resolution_status"], "schedule_implementation_required")


if __name__ == "__main__":
    unittest.main()
