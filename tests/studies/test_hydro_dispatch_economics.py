import pandas as pd
from workflow.scripts.studies.chapter2.audit_hydro_dispatch_economics import (
    build_economic_table,
)


def test_economic_table_compares_output_vom_with_local_price() -> None:
    links = pd.DataFrame(
        {
            "bus1": ["Electric"],
            "efficiency": [2.0],
            "marginal_cost": [12.0],
        },
        index=["Turbine"],
    )
    output = pd.DataFrame({"Turbine": [-2.0, -4.0]})
    prices = pd.DataFrame({"Electric": [5.0, 7.0]})
    mapping = pd.DataFrame(
        {
            "facility": ["Station"],
            "model_component_type": ["link"],
            "model_component_name": ["Turbine"],
        }
    )

    result = build_economic_table(links, output, prices, mapping).iloc[0]

    assert result["effective_turbine_vom_cad_per_mwh"] == 6.0
    assert result["hours_price_strictly_above_vom"] == 1
    assert result["generation_gwh"] == 0.006
    assert result["economic_ordering"] == "price_exceeds_vom_in_some_hours"
