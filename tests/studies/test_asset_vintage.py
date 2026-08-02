from __future__ import annotations

import pandas as pd
import pytest

from pypsa_bc.studies.asset_vintage import (
    active_inventory_codes,
    filter_component_dictionary,
)


def test_active_inventory_codes_excludes_future_and_closed_assets() -> None:
    inventory = pd.DataFrame(
        {
            "code": ["active", "future", "closed", "closes_this_year"],
            "start_year": [2010, 2025, 2000, 2000],
            "closure_year": [None, None, 2020, 2021],
        }
    )
    assert active_inventory_codes(inventory, 2021, code_column="code") == {
        "active",
        "closes_this_year",
    }


def test_filter_component_dictionary_records_vintage_decisions() -> None:
    components = {
        "active": {"name": "active Solar Generator"},
        "future": {"name": "future Solar Generator"},
    }
    retained, audit = filter_component_dictionary(
        components, {"active"}, match="key", resource="solar"
    )
    assert list(retained) == ["active"]
    assert audit.set_index("inventory_code").at["future", "decision"] == (
        "outside_model_year"
    )


def test_name_prefix_mapping_must_be_unambiguous() -> None:
    with pytest.raises(ValueError, match="uniquely map"):
        filter_component_dictionary(
            {0: {"name": "unknown biomass Link"}},
            {"BC_A_GEN"},
            match="name_prefix",
            resource="thermal",
        )
