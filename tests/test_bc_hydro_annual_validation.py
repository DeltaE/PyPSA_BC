from pathlib import Path

from workflow.scripts.fetch_bc_hydro_annual_validation import _sha256


def test_sha256_is_stable(tmp_path: Path) -> None:
    source = tmp_path / "evidence.pdf"
    source.write_bytes(b"%PDF-test")
    assert _sha256(source) == (
        "3c87d37f1dbea6909f917ce437c390fb8e655a774387d9e69301c0b2283d5b63"
    )
