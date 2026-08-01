"""Fetch external validation datasets used by the visual report.

Examples
--------
uv run python -m workflow.scripts.fetch_validation_data --year 2021
uv run python -m workflow.scripts.fetch_validation_data --year 2021 --force

Files are written under ``data/validation`` with a machine-readable manifest.
The report never downloads data implicitly; it reads only these pinned files.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import requests

VALIDATION_ROOT = Path("data/validation")
BC_HYDRO_URLS = {
    2021: (
        "https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/"
        "corporate/suppliers/transmission-system/balancing_authority_load_data/"
        "Historical%20Transmission%20Data/BalancingAuthorityLoad2021.xls"
    ),
    2022: (
        "https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/"
        "corporate/suppliers/transmission-system/balancing_authority_load_data/"
        "Historical%20Transmission%20Data/BalancingAuthorityLoad2022.xls"
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(year: int, force: bool = False) -> Path:
    if year not in BC_HYDRO_URLS:
        raise ValueError(f"No registered BC Hydro validation source for {year}")

    output_dir = VALIDATION_ROOT / "bc_hydro"
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / f"BalancingAuthorityLoad{year}.xls"
    url = BC_HYDRO_URLS[year]

    if force or not target.exists() or target.stat().st_size == 0:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        target.write_bytes(response.content)

    manifest_path = VALIDATION_ROOT / "manifest.json"
    manifest = {"datasets": []}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    record = {
        "dataset": "BC Hydro balancing authority hourly load",
        "year": year,
        "source_url": url,
        "local_path": target.as_posix(),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "bytes": target.stat().st_size,
        "sha256": _sha256(target),
        "use": "Independent comparison of model system-load shape and magnitude",
    }
    records = [
        item
        for item in manifest.get("datasets", [])
        if not (item.get("dataset") == record["dataset"] and item.get("year") == year)
    ]
    records.append(record)
    manifest_path.write_text(
        json.dumps({"datasets": records}, indent=2) + "\n",
        encoding="utf-8",
    )
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch pinned report-validation data.")
    parser.add_argument("--year", type=int, default=2021)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    path = fetch(args.year, force=args.force)
    print(f"Validation data: {path}")
    print(f"Manifest: {VALIDATION_ROOT / 'manifest.json'}")


if __name__ == "__main__":
    main()
