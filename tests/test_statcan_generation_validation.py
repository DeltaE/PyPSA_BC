import pandas as pd
from workflow.scripts.fetch_statcan_generation_validation import (
    TOTAL_CLASS,
    TYPE_MAP,
    filter_bc_generation,
)


def test_filter_bc_generation_requires_complete_monthly_components() -> None:
    rows = []
    for month in range(1, 13):
        values = {
            "Total all types of electricity generation": 21.0,
            "Hydraulic turbine": 10.0,
            "Total electricity production from non-renewable combustible fuels": 2.0,
            "Total electricity production from biomass": 3.0,
            "Wind power turbine": 4.0,
            "Solar": 1.0,
            "Other types of electricity generation": 1.0,
        }
        for category, value in values.items():
            rows.append(
                {
                    "REF_DATE": f"2021-{month:02d}",
                    "GEO": "British Columbia",
                    "Class of electricity producer": TOTAL_CLASS,
                    "Type of electricity generation": category,
                    "UOM": "Megawatt hours",
                    "SCALAR_FACTOR": "units",
                    "VALUE": value,
                }
            )
    result = filter_bc_generation(pd.DataFrame(rows), 2021)

    assert len(result) == 12
    assert set(TYPE_MAP.values()).issubset(result.columns)
    assert result["component_reconciliation_error_mwh"].eq(0.0).all()
