import pandas as pd

from pypsa_bc.studies.baseline_calibration import (
    compare_generation,
    compare_load,
    evaluate_calibration_gate,
)


def test_load_comparison_distinguishes_uniform_scaling_from_chronology() -> None:
    index = pd.date_range("2021-01-01", periods=4, freq="h")
    observed = pd.Series([10.0, 20.0, 30.0, 20.0], index=index)
    model = observed * 0.99

    _, metrics = compare_load(model, observed)

    assert abs(metrics["annual_relative_error"] + 0.01) < 1e-12
    assert metrics["hourly_pearson_r"] == 1.0
    assert metrics["peak_position_error_hours"] == 0
    assert metrics["evidence_role"] == "input_fidelity_not_independent_validation"


def test_generation_gate_rejects_hydro_substitution_for_material_thermal() -> None:
    months = [f"2021-{month:02d}" for month in range(1, 13)]
    observed = pd.DataFrame(
        {
            "month": months,
            "hydro_mwh": [80.0] * 12,
            "nonrenewable_combustible_mwh": [5.0] * 12,
            "biomass_mwh": [5.0] * 12,
            "wind_mwh": [10.0] * 12,
            "solar_mwh": [0.0] * 12,
            "total_generation_mwh": [100.0] * 12,
        }
    )
    model = observed.copy()
    model["hydro_mwh"] = 90.0
    model["nonrenewable_combustible_mwh"] = 0.0
    model["biomass_mwh"] = 0.0
    _, annual, metrics = compare_generation(model, observed)
    load = {
        "annual_relative_error": 0.0,
        "peak_relative_error": 0.0,
        "hourly_pearson_r": 1.0,
        "peak_position_error_hours": 0,
    }

    checks = evaluate_calibration_gate(load, annual, metrics)

    assert checks["total_generation_magnitude"] == "PASS"
    assert checks["hydro_generation_magnitude"] == "FAIL"
    assert checks["material_nonhydro_dispatch"] == "FAIL"
    assert checks["generation_mix"] == "FAIL"
