"""Fetch external validation datasets used by the visual report.

Examples
--------
uv run python -m workflow.scripts.fetch_validation_data --year 2021
uv run python -m workflow.scripts.fetch_validation_data --year 2021 --dataset all
uv run python -m workflow.scripts.fetch_validation_data --year 2021 --force

Files are written under ``data/validation`` with a machine-readable manifest.
The report never downloads data implicitly; it reads only these pinned files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import requests

VALIDATION_ROOT = Path("data/validation")
TRANSMISSION_DATA_ROOT = (
    "https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/"
    "corporate/suppliers/transmission-system/balancing_authority_load_data/"
    "Historical%20Transmission%20Data"
)
ACTUAL_FLOW_ROOT = (
    "https://www.bchydro.com/content/dam/BCHydro/customer-portal/documents/"
    "corporate/suppliers/transmission-system/actual_flow_data/historical_data"
)

DATASET_SPECS = {
    "load": {
        "filename": "BalancingAuthorityLoad{year}.xls",
        "url": TRANSMISSION_DATA_ROOT + "/BalancingAuthorityLoad{year}.xls",
        "dataset": "BC Hydro balancing authority hourly load",
        "use": "Independent comparison of model system-load shape and magnitude",
        "claim_boundary": "Gross telemetered provincial load is not regional demand.",
    },
    "bc_ab_ttc": {
        "filename": "BCABTTC{year}.xls",
        "url": TRANSMISSION_DATA_ROOT + "/BCABTTC{year}.xls",
        "dataset": "BC Hydro hourly TTC: British Columbia to Alberta",
        "use": "Observed study-year upper bound for the BC-to-Alberta interface",
        "claim_boundary": "Directional TTC is a path limit, not realized flow or an internal line rating.",
    },
    "ab_bc_ttc": {
        "filename": "ABBCTTC{year}.xls",
        "url": TRANSMISSION_DATA_ROOT + "/ABBCTTC{year}.xls",
        "dataset": "BC Hydro hourly TTC: Alberta to British Columbia",
        "use": "Observed study-year upper bound for the Alberta-to-BC interface",
        "claim_boundary": "Directional TTC is a path limit, not realized flow or an internal line rating.",
    },
    "bc_us_ttc": {
        "filename": "BCUSTTC{year}.xls",
        "url": TRANSMISSION_DATA_ROOT + "/BCUSTTC{year}.xls",
        "dataset": "BC Hydro hourly TTC: British Columbia to United States",
        "use": "Observed study-year upper bound for the BC-to-US interface",
        "claim_boundary": "Directional TTC is a path limit, not realized flow or an internal line rating.",
    },
    "us_bc_ttc": {
        "filename": "USBCTTC{year}.xls",
        "url": TRANSMISSION_DATA_ROOT + "/USBCTTC{year}.xls",
        "dataset": "BC Hydro hourly TTC: United States to British Columbia",
        "use": "Observed study-year upper bound for the US-to-BC interface",
        "claim_boundary": "Directional TTC is a path limit, not realized flow or an internal line rating.",
    },
    "net_actual_flow": {
        "filename": "HourlyTielineData{year}.xls",
        "url": ACTUAL_FLOW_ROOT + "/HourlyTielineData{year}.xls",
        "dataset": "BC Hydro hourly net actual tie-line flow",
        "use": "Independent validation of model import-export chronology and magnitude",
        "claim_boundary": "Net actual flow requires an inspected and documented sign convention.",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(year: int, force: bool = False, dataset: str = "load") -> Path:
    if dataset not in DATASET_SPECS:
        raise ValueError(f"Unknown BC Hydro validation dataset: {dataset}")
    if year < 2007 and dataset != "load":
        raise ValueError(f"{dataset} is not registered before 2007")
    if year < 2001:
        raise ValueError("BC Hydro hourly load is not registered before 2001")

    output_dir = VALIDATION_ROOT / "bc_hydro"
    output_dir.mkdir(parents=True, exist_ok=True)
    spec = DATASET_SPECS[dataset]
    target = output_dir / str(spec["filename"]).format(year=year)
    url = str(spec["url"]).format(year=year)

    if force or not target.exists() or target.stat().st_size == 0:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        target.write_bytes(response.content)

    manifest_path = VALIDATION_ROOT / "manifest.json"
    manifest = {"datasets": []}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    record = {
        "dataset_id": dataset,
        "dataset": spec["dataset"],
        "year": year,
        "source_url": url,
        "local_path": target.as_posix(),
        "retrieved_utc": datetime.now(UTC).isoformat(),
        "bytes": target.stat().st_size,
        "sha256": _sha256(target),
        "use": spec["use"],
        "evidence_class": "observed",
        "claim_boundary": spec["claim_boundary"],
    }
    records = [
        item
        for item in manifest.get("datasets", [])
        if not (
            item.get("dataset") == record["dataset"]
            and item.get("year") == year
        )
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
    parser.add_argument(
        "--dataset",
        choices=[*DATASET_SPECS, "all"],
        default="load",
        help="BC Hydro dataset to pin; 'all' includes load, four TTC directions, and net actual flow.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    dataset_ids = list(DATASET_SPECS) if args.dataset == "all" else [args.dataset]
    paths = [fetch(args.year, force=args.force, dataset=item) for item in dataset_ids]
    for path in paths:
        print(f"Validation data: {path}")
    print(f"Manifest: {VALIDATION_ROOT / 'manifest.json'}")


if __name__ == "__main__":
    main()
