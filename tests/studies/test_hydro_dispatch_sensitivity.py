import pandas as pd
from workflow.scripts.studies.chapter2.compare_hydro_dispatch_sensitivities import (
    build_comparison,
)


def test_station_sensitivity_comparison_preserves_evidence_columns() -> None:
    reference = pd.DataFrame(
        {
            "facility": ["Station"],
            "model_calendar_2021_energy_gwh": [1.0],
            "observed_fiscal_2021_energy_gwh": [4.0],
            "turbine_capture_fraction": [0.1],
            "energy_screening_flag": ["REVIEW"],
        }
    )
    sensitivity = reference.copy()
    sensitivity["model_calendar_2021_energy_gwh"] = 3.0
    sensitivity["turbine_capture_fraction"] = 0.8

    result = build_comparison(reference, sensitivity)

    assert result.iloc[0]["energy_change_gwh"] == 2.0
    assert result.iloc[0]["energy_change_percent_of_observed"] == 0.5
