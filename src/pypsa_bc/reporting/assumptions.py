"""
assumptions.py — append-only assumption logger for PyPSA-BC.

Keeps a single Markdown table (ASSUMPTION_REPORT.md) as a running audit trail.
Call ``log_assumption(...)`` wherever a modelling assumption is plugged in; a
timestamped row is appended. The file and header are created on first use, and
exact duplicates are skipped so re-running a notebook does not spam the log.

Quick use
---------
>>> from pypsa_bc.assumptions import log_assumption
>>> log_assumption("line s_nom source", "ttc_summer", unit="MW",
...                rationale="CODERS TTC used as thermal proxy (pf=1)",
...                source="transmission_lines.csv 2025OCT28",
...                module="base_network", scenario="Base_CNZ_2035")

Or hold one logger bound to a scenario:

>>> from pypsa_bc.assumptions import AssumptionLog
>>> log = AssumptionLog(scenario="NE_FE_2040")
>>> log.add("reservoir e_initial", 0, unit="MWh", rationale="Cold-start; WSC Jan-1 storage not yet integrated")

Author: Md Eliasinul Islam
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

DEFAULT_REPORT = "docs/ASSUMPTION_REPORT.md"
COLUMNS = ["Date", "Scenario", "Module", "Parameter", "Value", "Unit", "Rationale", "Source"]


class AssumptionLog:
    """A single Markdown assumption table, appended to one row at a time."""

    def __init__(self, path: str | Path = DEFAULT_REPORT, scenario: str = "-"):
        self.path = Path(path)
        self.scenario = scenario
        self._ensure_header()

    # --- internals ------------------------------------------------------- #

    def _ensure_header(self) -> None:
        if self.path.exists():
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        head = (
            "# Assumption Report\n\n"
            "Running log of modelling assumptions for PyPSA-BC. One row per\n"
            "assumption, appended automatically when it is plugged into the model.\n"
            "Do not edit rows by hand — append via `pypsa_bc.assumptions.log_assumption`.\n\n"
            "| " + " | ".join(COLUMNS) + " |\n"
            "|" + "|".join([" --- "] * len(COLUMNS)) + "|\n"
        )
        self.path.write_text(head, encoding="utf-8")

    @staticmethod
    def _esc(x) -> str:
        return str(x).replace("|", "\\|").replace("\n", " ").strip()

    def _data_rows(self) -> list[list[str]]:
        rows = []
        for line in self.path.read_text(encoding="utf-8").splitlines():
            if line.startswith("|") and "---" not in line:
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if cells != COLUMNS:  # skip the header row
                    rows.append(cells)
        return rows

    # --- public API ------------------------------------------------------ #

    def add(
        self,
        parameter: str,
        value,
        unit: str = "",
        rationale: str = "",
        source: str = "",
        module: str = "",
        scenario: str | None = None,
        date: str | None = None,
        dedupe: bool = True,
    ) -> None:
        """Append one assumption row.

        Args:
            parameter: what is being assumed (e.g. "line s_nom source").
            value: the assumed value (e.g. "ttc_summer", 0.10, "pf=1").
            unit, rationale, source, module: context columns.
            scenario: overrides the logger's scenario for this row.
            date: ISO date; defaults to today.
            dedupe: skip if an identical (scenario, module, parameter, value)
                row already exists — keeps re-runs clean while still recording
                a genuinely changed value as a new row.
        """
        scenario = scenario if scenario is not None else self.scenario
        date = date or datetime.now().strftime("%Y-%m-%d")
        cells = [date, scenario, module, parameter, value, unit, rationale, source]
        cells = [self._esc(c) for c in cells]

        if dedupe:
            key = lambda r: (r[1], r[2], r[3], r[4])  # scenario, module, param, value
            if any(key(r) == key(cells) for r in self._data_rows()):
                return

        with self.path.open("a", encoding="utf-8") as f:
            f.write("| " + " | ".join(cells) + " |\n")


# --- module-level convenience (default report path) --------------------- #

_default_log: AssumptionLog | None = None


def log_assumption(parameter: str, value, *, path: str | Path = DEFAULT_REPORT,
                   **kwargs) -> None:
    """Append an assumption to the default report (``docs/ASSUMPTION_REPORT.md``).

    Accepts the same keyword arguments as :meth:`AssumptionLog.add`
    (unit, rationale, source, module, scenario, date, dedupe).
    """
    global _default_log
    if _default_log is None or str(_default_log.path) != str(path):
        _default_log = AssumptionLog(path)
    _default_log.add(parameter, value, **kwargs)
