from __future__ import annotations

import unittest

import pandas as pd
from workflow.scripts.studies.chapter2.normalize_cascade_representation import normalize


class CascadeRepresentationTests(unittest.TestCase):
    def test_only_cascade_ror_is_converted_to_water_passing(self):
        frame = pd.DataFrame(
            {
                "asset_id": ["CASCADE", "STANDALONE", "RESERVOIR"],
                "cascade_group": ["Example", "DEFAULT", "Example"],
                "hydro_type": ["ror", "ror", "reservoir"],
            }
        )
        result, changed = normalize(frame)
        self.assertEqual(changed, ["CASCADE"])
        self.assertEqual(result.loc[0, "hydro_type"], "ror-water")
        self.assertEqual(result.loc[1, "hydro_type"], "ror")
        self.assertEqual(result.loc[2, "hydro_type"], "reservoir")


if __name__ == "__main__":
    unittest.main()
