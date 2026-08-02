from workflow.scripts.studies.chapter2.audit_publication_readiness import (
    evaluate_readiness,
)


def test_calibration_pass_does_not_equal_independent_validation() -> None:
    summary = evaluate_readiness(
        {"computational_integrity_gate_status": "PASS"},
        {
            "calibration_fit_status": "PASS",
            "independent_validation_status": "NOT_ESTABLISHED",
        },
        {
            "capacity_gate_status": "PASS",
            "energy_screening_status": "REVIEW",
            "energy_screening_facilities": ["Seven Mile"],
            "predominantly_spilled_facilities": ["Seven Mile"],
        },
    )

    assert summary["matched_representation_authorization"] == "HOLD"
    assert summary["publication_claim_readiness"] == "HOLD"
    assert summary["representation_blockers"] == ["named_station_energy_screen"]
    assert "independent_validation" in summary["publication_blockers"]


def test_all_independent_layers_can_release_gate() -> None:
    summary = evaluate_readiness(
        {"computational_integrity_gate_status": "PASS"},
        {"calibration_fit_status": "PASS", "independent_validation_status": "PASS"},
        {
            "capacity_gate_status": "PASS",
            "energy_screening_status": "WITHIN_SCREENING_BAND",
        },
    )

    assert summary["matched_representation_authorization"] == "GO"
    assert summary["publication_claim_readiness"] == "READY"
