"""
provenance.py — dated log of what external data was sourced, and when.

Maintains one Markdown table (reports/DATA_SOURCES.md), upserted by dataset, so
it always reflects the current provenance: source URL, last-retrieved date,
local path, and a short detail (rows / size). Called automatically when the
downloader saves a file or CODERSData fetches a table from the API.

    from pypsa_bc.reporting.provenance import record_source
    record_source("CODERS:generators", url, path, detail="268 BC rows")
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

DEFAULT_REPORT = "reports/DATA_SOURCES.md"
COLUMNS = ["Dataset", "Source", "Retrieved", "Path", "Detail"]


def _esc(x) -> str:
    return str(x).replace("|", "\\|").replace("\n", " ").strip()


def _parse(report: Path) -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    if not report.exists():
        return rows
    for line in report.read_text(encoding="utf-8").splitlines():
        if line.startswith("|") and "---" not in line:
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) == len(COLUMNS) and cells != COLUMNS:
                rows[cells[0]] = cells
    return rows


def record_source(dataset: str, source: str, path, detail: str = "",
                  report: str | Path = DEFAULT_REPORT) -> None:
    """Upsert one dataset's provenance row (keyed by dataset; date = today)."""
    report = Path(report)
    report.parent.mkdir(parents=True, exist_ok=True)

    rows = _parse(report)
    key = _esc(dataset)
    rows[key] = [key, _esc(source), datetime.now().strftime("%Y-%m-%d"),
                 _esc(path), _esc(detail)]

    header = [
        "# Data Sources",
        "",
        "Provenance of external inputs (auto-maintained by pypsa_bc). "
        "One row per dataset; `Retrieved` is the last fetch date.",
        "",
        "| " + " | ".join(COLUMNS) + " |",
        "|" + "|".join([" --- "] * len(COLUMNS)) + "|",
    ]
    body = ["| " + " | ".join(r) + " |" for r in sorted(rows.values(), key=lambda r: r[0])]
    report.write_text("\n".join(header + body) + "\n", encoding="utf-8")
