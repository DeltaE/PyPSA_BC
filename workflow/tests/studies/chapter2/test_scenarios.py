import unittest

from pypsa_bc.studies.cascade.scenarios import expand_scenarios, runnable_case_ids


class ScenarioExpansionTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "scenario_axes": {
                "reservoir_representations": {
                    "A": {"implemented": True},
                    "B": {"implemented": False},
                },
                "water_use_policies": {
                    "evidence": {"implemented": True, "constraint_mode": "evidence_based"},
                    "off": {"implemented": True, "constraint_mode": "no_minimum_release"},
                },
            }
        }

    def test_cartesian_product_is_explicit(self):
        cases = expand_scenarios(self.config)
        self.assertEqual(len(cases), 4)
        self.assertEqual(
            {case["case_id"] for case in cases},
            {"A__evidence", "A__off", "B__evidence", "B__off"},
        )

    def test_unimplemented_representation_is_not_runnable(self):
        self.assertEqual(runnable_case_ids(self.config), ["A__evidence", "A__off"])


if __name__ == "__main__":
    unittest.main()
