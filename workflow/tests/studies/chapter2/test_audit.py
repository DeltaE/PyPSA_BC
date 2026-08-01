from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from pypsa_bc.studies.cascade.audit import (
    _inflow_checks,
    _topology,
    _verified_release_schedule_map,
    find_cycle,
    is_missing,
    monthly_flow_reconciliation,
)


class MissingValueTests(unittest.TestCase):
    def test_protocol_tokens_are_missing(self):
        for value in ("", "  ", "DEFAULT", None, pd.NA, float("nan")):
            with self.subTest(value=value):
                self.assertTrue(is_missing(value))

    def test_zero_and_named_reservoir_are_not_missing(self):
        self.assertFalse(is_missing(0))
        self.assertFalse(is_missing("BC_GMS_RES"))


class ReleaseScheduleAuditTests(unittest.TestCase):
    def test_complete_source_labelled_schedule_is_verified(self):
        snapshots = pd.date_range("2021-01-01 00:00", periods=3, freq="h")
        frame = pd.DataFrame(
            {
                "link": "station Spill Link",
                "snapshot": snapshots,
                "minimum_release_m3_per_hour": [360.0, 360.0, 720.0],
                "schedule_id": "route_rule",
                "source": "accepted_wup",
                "evidence_class": "observed",
            }
        )
        issues = []
        result = _verified_release_schedule_map(
            frame,
            {"start": str(snapshots[0]), "end": str(snapshots[-1])},
            issues,
        )
        self.assertEqual(issues, [])
        self.assertEqual(result["route_rule"]["minimum_release_m3_per_hour_min"], 360.0)
        self.assertEqual(result["route_rule"]["minimum_release_m3_per_hour_max"], 720.0)

    def test_incomplete_schedule_is_not_counted(self):
        frame = pd.DataFrame(
            {
                "link": ["station Spill Link"],
                "snapshot": ["2021-01-01 00:00"],
                "minimum_release_m3_per_hour": [360.0],
                "schedule_id": ["route_rule"],
                "source": ["accepted_wup"],
                "evidence_class": ["observed"],
            }
        )
        issues = []
        result = _verified_release_schedule_map(
            frame,
            {"start": "2021-01-01 00:00", "end": "2021-01-01 01:00"},
            issues,
        )
        self.assertEqual(result, {})
        self.assertEqual(issues[0]["check"], "release_schedule_coverage")


class TopologyTests(unittest.TestCase):
    def test_acyclic_chain(self):
        self.assertEqual(find_cycle([("A", "B"), ("B", "C")]), [])

    def test_cycle_is_reported_with_closed_path(self):
        cycle = find_cycle([("A", "B"), ("B", "C"), ("C", "A")])
        self.assertEqual(cycle[0], cycle[-1])
        self.assertEqual(set(cycle[:-1]), {"A", "B", "C"})

    def test_evidence_backed_non_generating_reservoir_is_not_an_orphan(self):
        generation = pd.DataFrame(
            {
                "asset_id": ["GEN"],
                "capacity": [1.0],
                "cascade_group": ["Example"],
                "cascade_order": [1],
                "hydro_type": ["reservoir"],
                "upper_reservoir_id": ["GENERATING_RES"],
                "lower_reservoir_id": ["DEFAULT"],
            }
        )
        reservoirs = pd.DataFrame(
            {
                "asset_id": ["GENERATING_RES", "CONTROL_RES"],
                "cascade_group": ["Example", "Example"],
            }
        )
        decisions = pd.DataFrame(
            {
                "reservoir_id": ["CONTROL_RES"],
                "role": ["non_generating_control_reservoir"],
                "generation_expected": [False],
                "source": ["official_source"],
                "locator": ["facility entry"],
                "note": ["No generation at the dam."],
            }
        )
        issues = []
        _, _, roles = _topology(generation, reservoirs, decisions, {}, issues)
        self.assertFalse(
            any(issue["check"] == "cascade_reservoir_station_coverage" for issue in issues)
        )
        control = roles.set_index("reservoir_id").loc["CONTROL_RES"]
        self.assertEqual(control["status"], "resolved_non_generating")


class UnitReconciliationTests(unittest.TestCase):
    def test_hourly_m3_per_hour_matches_monthly_m3_per_second(self):
        index = pd.date_range("2021-01-01", "2021-12-31 23:00", freq="h")
        source = pd.Series({month: float(month) for month in range(1, 13)})
        hourly = pd.Series([timestamp.month * 3600.0 for timestamp in index], index=index)
        result = monthly_flow_reconciliation(hourly, source)
        self.assertAlmostEqual(float(result["relative_difference"].max()), 0.0)

    def test_wrong_unit_is_detected(self):
        index = pd.date_range("2021-01-01", "2021-01-31 23:00", freq="h")
        hourly = pd.Series(10.0, index=index)
        source = pd.Series({1: 10.0})
        result = monthly_flow_reconciliation(hourly, source)
        self.assertGreater(
            float(result.loc[result["month"] == 1, "relative_difference"].iloc[0]), 0.9
        )

    def test_evidence_backed_non_generating_inflow_table_is_not_warned(self):
        index = pd.date_range("2021-01-01", "2021-01-01 01:00", freq="h")
        inflows = pd.DataFrame({"snapshot": index, "ACTIVE_RES": [3600.0, 3600.0]})
        with TemporaryDirectory() as directory:
            stats_dir = Path(directory)
            pd.DataFrame(
                {
                    "Month": ["January"],
                    "Mean Monthly Inflow": [1.0],
                }
            ).to_csv(stats_dir / "CONTROL_RES.csv", index=False)
            issues = []
            checks = _inflow_checks(
                inflows,
                stats_dir,
                {"start": str(index[0]), "end": str(index[-1])},
                {
                    "unit_scale_max_relative_difference": 0.02,
                    "monthly_calibration_max_relative_difference": 0.01,
                },
                issues,
                {"CONTROL_RES"},
            )
        control = checks.loc[checks["subject"] == "CONTROL_RES"].iloc[0]
        self.assertEqual(control["status"], "pass")
        self.assertFalse(
            any(issue["check"] == "unused_source_inflow_statistics" for issue in issues)
        )


if __name__ == "__main__":
    unittest.main()
