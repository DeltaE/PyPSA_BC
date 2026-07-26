"""
logger.py — staged, checklist-style progress logger for PyPSA-BC.

Console view per run:

    Base network — 3 stages: lines, buses, transformers
    ✓ lines  (0.2s)
        ⚠ [buses] no information to create bus for: BC_GST_JCT
        ⚠ [buses] substation code WEY is missing from the lines dataset
    ✓ buses  (13.1s)
        ⚠ [transformers] skipped 1 bus row(s) with non-string name: [812]
    ✓ transformers  (3.0s)
    → (live bar for whatever is currently running)

    Base network preparation completed. Delivered:
    ┌────────────────┬───────────────────────────────────────┬──────┐
    │ File           │ Path                                  │ Rows │
    └────────────────┴───────────────────────────────────────┴──────┘

Warnings/errors are printed above the bar (via tqdm.write) and tagged with the
stage they came from, so they link to the right step. Everything (INFO/DEBUG up)
is written in full to logs/pypsa_bc.log.
"""

from __future__ import annotations

import logging
import time
import csv
import re
import pickle
from contextlib import contextmanager
from pathlib import Path

try:
    from tqdm import tqdm  # plain text bar; works in terminal + notebook, no ipywidgets
except Exception:
    tqdm = None

_RESET = "\033[0m"
_COLOR = {logging.WARNING: "\033[33m", logging.ERROR: "\033[1;31m"}
_PREFIX = {logging.INFO: "›", logging.WARNING: "⚠", logging.ERROR: "✖", logging.DEBUG: "·"}
_ROOT = "pypsa_bc"


def _write(msg: str) -> None:
    """Print above an active tqdm bar without tearing it."""
    if tqdm is not None:
        tqdm.write(msg)
    else:
        print(msg, flush=True)


class _ConsoleFormatter(logging.Formatter):
    def format(self, record):
        color = _COLOR.get(record.levelno, "")
        prefix = _PREFIX.get(record.levelno, "")
        tag = record.name.split(".")[-1]
        tag = "" if tag == _ROOT else f"[{tag}] "
        return f"    {color}{prefix} {tag}{record.getMessage()}{_RESET}"  # indented under stage


class _TqdmHandler(logging.StreamHandler):
    def emit(self, record):
        try:
            _write(self.format(record))
        except Exception:
            self.handleError(record)


class _DedupFilter(logging.Filter):
    def __init__(self):
        super().__init__()
        self._seen: set = set()

    def filter(self, record) -> bool:
        if getattr(record, "once", False):
            key = (record.name, record.getMessage())
            if key in self._seen:
                return False
            self._seen.add(key)
        return True


def get_logger(name: str = _ROOT, *, logfile: str | None = "logs/pypsa_bc.log",
               console_level: int = logging.WARNING) -> logging.Logger:
    root = logging.getLogger(_ROOT)
    if not root.handlers:
        root.setLevel(logging.DEBUG)
        root.propagate = False
        dedup = _DedupFilter()

        console = _TqdmHandler()
        console.setLevel(console_level)
        console.setFormatter(_ConsoleFormatter())
        console.addFilter(dedup)
        root.addHandler(console)

        if logfile:
            Path(logfile).parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(logfile)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter(
                "%(asctime)s %(levelname)-7s %(name)-26s | %(message)s"))
            fh.addFilter(dedup)
            root.addHandler(fh)

    return root if name == _ROOT else logging.getLogger(f"{_ROOT}.{name}")


def log_once(msg: str, name: str = _ROOT, level: int = logging.INFO) -> None:
    get_logger(name).log(level, msg, extra={"once": True})


def set_console_level(level: int) -> None:
    get_logger()
    for h in logging.getLogger(_ROOT).handlers:
        if isinstance(h, _TqdmHandler):
            h.setLevel(level)


@contextmanager
def stage(title: str, name: str = _ROOT):
    log = get_logger(name)
    log.info(f"[{title}] start")
    t0 = time.perf_counter()
    try:
        yield log
    finally:
        log.info(f"[{title}] done in {time.perf_counter() - t0:.1f}s")


def _render_table(headers: list[str], rows: list[tuple]) -> str:
    cols = list(zip(*([tuple(headers)] + rows))) if rows else [(h,) for h in headers]
    w = [max(len(str(c)) for c in col) for col in cols]
    line = lambda l, m, r: l + m.join("─" * (x + 2) for x in w) + r
    row = lambda r: "│" + "│".join(f" {str(v):<{x}} " for v, x in zip(r, w)) + "│"
    return "\n".join([line("┌", "┬", "┐"), row(headers), line("├", "┼", "┤"),
                      *[row(r) for r in rows], line("└", "┴", "┘")])


def _format_size(num_bytes: int) -> str:
    if num_bytes < 1024:
        return f"{num_bytes} B"
    if num_bytes < 1024**2:
        return f"{num_bytes / 1024:.2f} KB"
    if num_bytes < 1024**3:
        return f"{num_bytes / (1024**2):.2f} MB"
    return f"{num_bytes / (1024**3):.2f} GB"


def _shape_text(rows: int, cols: int) -> str:
    return f"Rows {rows}, Columns {cols}"


def _csv_shape(path: Path) -> str | None:
    try:
        with path.open("r", newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                return _shape_text(0, 0)
            ncols = len(header)
            nrows = sum(1 for _ in reader)
            return _shape_text(nrows, ncols)
    except Exception:
        return None


def _shape_from_object(obj) -> str | None:
    """Best-effort rows/cols for common serialized objects."""
    try:
        # pandas DataFrame / Series
        if hasattr(obj, "shape") and isinstance(getattr(obj, "shape"), tuple):
            shp = obj.shape
            if len(shp) == 2:
                return _shape_text(int(shp[0]), int(shp[1]))
            if len(shp) == 1:
                return _shape_text(int(shp[0]), 1)

        if isinstance(obj, list):
            nrows = len(obj)
            if nrows == 0:
                return _shape_text(0, 0)
            first = obj[0]
            if isinstance(first, dict):
                return _shape_text(nrows, len(first))
            if isinstance(first, (list, tuple)):
                return _shape_text(nrows, len(first))
            return _shape_text(nrows, 1)

        if isinstance(obj, dict):
            nrows = len(obj)
            if nrows == 0:
                return _shape_text(0, 0)
            first_val = next(iter(obj.values()))
            if isinstance(first_val, dict):
                return _shape_text(nrows, len(first_val))
            if isinstance(first_val, (list, tuple)):
                return _shape_text(nrows, len(first_val))
            return _shape_text(nrows, 1)
    except Exception:
        return None
    return None


def _pickle_shape(path: Path) -> str | None:
    try:
        with path.open("rb") as f:
            obj = pickle.load(f)
        return _shape_from_object(obj)
    except Exception:
        return None


def _parse_legacy_size_note(note: str) -> tuple[str | None, str | None]:
    """Parse strings like '1.48 MB (exists locally)' into (size, status)."""
    m = re.match(
        r"^\s*([0-9]+(?:\.[0-9]+)?\s*(?:B|KB|MB|GB))\s*(?:\(([^)]+)\))?\s*$",
        note,
        flags=re.IGNORECASE,
    )
    if not m:
        return None, None
    size = m.group(1)
    status = m.group(2)
    return size, status


class Pipeline:
    """A staged run: checklist of ticks, a live bar for the current stage, and a
    delivered-files table at the end."""

    def __init__(self, title: str, stages: list[str], name: str = _ROOT):
        self.title = title
        self.stages = stages
        self.n = len(stages)
        self.log = get_logger(name)
        self._bar = None
        self._done = 0
        self._delivered: list[tuple] = []

    def __enter__(self):
        self.log.info(f"=== {self.title}: {self.n} stages ===")
        _write(f"\033[1m{self.title}\033[0m — {self.n} stages: {', '.join(self.stages)}")
        if tqdm is not None:
            self._bar = tqdm(total=self.n, desc=f"→ {self.stages[0]}", leave=True,
                             bar_format="{desc} |{bar}| {elapsed}{postfix}")
        return self

    @contextmanager
    def stage(self, label: str):
        current = self._done + 1
        if self._bar is not None:
            self._bar.set_description_str(f"→ {label} ({current}/{self.n})")
            self._bar.set_postfix_str("")
            self._bar.refresh()
        self.log.info(f"[{label}] start ({current}/{self.n})")
        t0 = time.perf_counter()
        try:
            yield self.log
        except Exception as exc:
            _write(f"\033[1;31m✗ {label} FAILED: {exc}\033[0m")
            self.log.error(f"[{label}] FAILED: {exc}")
            raise
        else:
            dt = time.perf_counter() - t0
            _write(f"\033[32m✓ {label}\033[0m  ({dt:.1f}s)")
            self.log.info(f"[{label}] done in {dt:.1f}s")
        finally:
            self._done += 1
            if self._bar is not None:
                self._bar.update(1)

    def deliver(self, filename: str, path, rows: int | str | None = None) -> None:
        """Record an output file to list in the completion table."""
        p = Path(path)
        shape = ""
        size = ""

        if p.exists() and p.is_file():
            size = _format_size(p.stat().st_size)
            if p.suffix.lower() == ".csv":
                shape = _csv_shape(p) or ""
            elif p.suffix.lower() in {".pickle", ".pkl"}:
                shape = _pickle_shape(p) or ""
        elif not p.exists():
            shape = "missing"

        if isinstance(rows, int):
            if not shape:
                shape = f"({rows}, ?)"
        elif isinstance(rows, str):
            note = rows.strip()
            lower = note.lower()
            legacy_size, legacy_status = _parse_legacy_size_note(note)
            if legacy_size is not None:
                size = legacy_size
                if legacy_status:
                    shape = f"{shape} [{legacy_status}]" if shape else legacy_status
            elif "skip" in lower or "missing" in lower or "exists locally" in lower:
                shape = f"{shape} [{note}]" if shape else note
            elif not size:
                size = note
            elif not shape:
                shape = note

        self._delivered.append((filename, str(path), shape or "-", size or "-"))
        self.log.debug(f"delivered {filename} (shape={shape}, size={size}) -> {path}")

    def status(self, msg: str) -> None:
        if self._bar is not None:
            self._bar.set_postfix_str(msg)
        self.log.debug(msg)

    def __exit__(self, *exc):
        if self._bar is not None:
            self._bar.set_description_str(f"{self.title}: complete")
            self._bar.close()
        self.log.info(f"=== {self.title}: complete ===")
        if self._delivered:
            table = _render_table(["File", "Path", "Data Shape", "Size"], self._delivered)
            _write(f"\n\033[1m{self.title} preparation completed. Delivered:\033[0m\n{table}")
        return False
