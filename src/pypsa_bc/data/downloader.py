"""
downloader.py — fetch public raw inputs into their local paths.

Config-driven: `config/data.yaml` holds a `remote:` map of
`name -> {url, dest}`. Downloads stream to a `.part` file and are moved into
place only on success (so a broken download never leaves a half file that looks
valid). Existing files are skipped unless `force=True`.
"""

from __future__ import annotations

from pathlib import Path

import requests

from pypsa_bc.reporting.logger import get_logger
from pypsa_bc.reporting.provenance import record_source

log = get_logger("downloader")


def download_file(url: str, dest: str | Path, *, force: bool = False,
                  chunk: int = 1 << 16, timeout: int = 60) -> Path:
    """Stream-download `url` to `dest`; skip if present unless `force`."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not force:
        log.info(f"{dest.name}: already present ({dest.stat().st_size/1e6:.1f} MB), skipping")
        return dest

    log.info(f"{dest.name}: downloading from {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for block in r.iter_content(chunk):
                if block:
                    f.write(block)
    tmp.replace(dest)  # atomic move on success
    size_mb = dest.stat().st_size / 1e6
    log.info(f"{dest.name}: saved {size_mb:.1f} MB -> {dest}")
    record_source(dest.stem, url, dest, detail=f"{size_mb:.1f} MB")
    return dest


def download_all(remote_cfg: dict, *, only: list[str] | None = None,
                 force: bool = False) -> dict[str, Path]:
    """Download every entry in a `remote:` config (or the subset in `only`)."""
    out: dict[str, Path] = {}
    for name, spec in remote_cfg.items():
        if only and name not in only:
            continue
        out[name] = download_file(spec["url"], spec["dest"], force=force)
    return out
