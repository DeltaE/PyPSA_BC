import pandas as pd
import pytest

from pypsa_bc.hydro import align_inflow_index


def test_align_inflow_index_precedes_monthly_calibration():
    utc_cutout = pd.DataFrame(
        {"BC_TEST_RES": range(8760)},
        index=pd.date_range("2023-01-01 07:00", periods=8760, freq="h"),
    )

    aligned = align_inflow_index(
        utc_cutout,
        "2021-01-01 00:00",
        "2021-12-31 23:00",
    )

    assert aligned.index[0] == pd.Timestamp("2021-01-01 00:00")
    assert aligned.index[-1] == pd.Timestamp("2021-12-31 23:00")
    assert aligned.loc["2021-02"].shape[0] == 28 * 24
    assert aligned.iloc[:, 0].tolist() == utc_cutout.iloc[:, 0].tolist()


def test_align_inflow_index_rejects_length_mismatch():
    cutout = pd.DataFrame(
        {"BC_TEST_RES": [1.0, 2.0]},
        index=pd.date_range("2023-01-01", periods=2, freq="h"),
    )

    with pytest.raises(ValueError, match="cutout has 2 samples"):
        align_inflow_index(
            cutout,
            "2021-01-01 00:00",
            "2021-01-01 02:00",
        )
