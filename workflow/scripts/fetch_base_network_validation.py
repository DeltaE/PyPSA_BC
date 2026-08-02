"""Fetch and fingerprint public evidence used by Rule 1 manual corrections."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import requests

from pypsa_bc.attributes_parser import AttributesParser


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(force: bool = False) -> tuple[Path, Path]:
    """Download the configured BC Hydro map when needed and write metadata."""
    spec = AttributesParser().data_cfg["remote"]["bc_hydro_transmission_map"]
    destination = Path(spec["dest"])
    metadata_path = Path(spec["metadata"])
    destination.parent.mkdir(parents=True, exist_ok=True)

    downloaded = False
    if force or not destination.exists():
        response = requests.get(spec["url"], timeout=120)
        response.raise_for_status()
        destination.write_bytes(response.content)
        downloaded = True

    metadata = {
        "source_name": "BC Hydro Provincial Transmission System Map",
        "source_url": spec["url"],
        "source_date": "2025-04",
        "local_path": str(destination),
        "sha256": sha256(destination),
        "bytes": destination.stat().st_size,
        "metadata_generated_at_utc": datetime.now(UTC).isoformat(),
        "downloaded_this_run": downloaded,
        "use_limitation": (
            "Schematic evidence only; manually interpreted points are representative "
            "and must not be used for engineering distance or asset-location claims."
        ),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return destination, metadata_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    destination, metadata = fetch(force=args.force)
    print(f"Evidence: {destination}")
    print(f"Metadata: {metadata}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
