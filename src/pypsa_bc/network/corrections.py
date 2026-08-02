"""Load, apply, and audit explicit Rule 1 correction registers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pypsa_bc.reporting.logger import get_logger

log = get_logger("network.corrections")

NODE_REQUIRED = {
    "correction_id",
    "node_code",
    "action",
    "longitude",
    "latitude",
    "method",
    "source",
    "rationale",
    "confidence",
    "manual_treatment",
    "active",
}
LINE_REQUIRED = {
    "correction_id",
    "transmission_line_id",
    "action",
    "field",
    "original_value",
    "replacement_value",
    "method",
    "source",
    "rationale",
    "confidence",
    "manual_treatment",
    "active",
}


def _active_rows(path: str | Path, required: set[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"transmission_line_id": "string"})
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Correction register {path} is missing columns: {missing}")
    duplicate_ids = frame["correction_id"].duplicated(keep=False)
    if duplicate_ids.any():
        repeated = sorted(frame.loc[duplicate_ids, "correction_id"].astype(str).unique())
        raise ValueError(f"Correction register {path} has duplicate IDs: {repeated}")
    return frame[frame["active"].astype(str).str.lower().eq("true")].copy()


def load_node_corrections(path: str | Path) -> pd.DataFrame:
    """Return active node corrections indexed by source node code."""
    frame = _active_rows(path, NODE_REQUIRED)
    duplicate_nodes = frame["node_code"].duplicated(keep=False)
    if duplicate_nodes.any():
        repeated = sorted(frame.loc[duplicate_nodes, "node_code"].astype(str).unique())
        raise ValueError(f"Active node corrections are not unique: {repeated}")
    frame["longitude"] = pd.to_numeric(frame["longitude"], errors="raise")
    frame["latitude"] = pd.to_numeric(frame["latitude"], errors="raise")
    return frame.set_index("node_code", drop=False)


def load_line_corrections(path: str | Path) -> pd.DataFrame:
    """Return active line corrections in declared register order."""
    return _active_rows(path, LINE_REQUIRED)


def _line_id(value: object) -> str:
    try:
        numeric = float(value)
        if numeric.is_integer():
            return str(int(numeric))
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def apply_line_corrections(lines: pd.DataFrame, register: pd.DataFrame) -> pd.DataFrame:
    """Apply active line patches while checking their declared source values."""
    corrected = lines.copy()
    corrected_ids = corrected["transmission_line_id"].map(_line_id)

    for _, correction in register.iterrows():
        correction_id = str(correction["correction_id"])
        line_id = _line_id(correction["transmission_line_id"])
        action = str(correction["action"])
        matches = corrected_ids.eq(line_id)
        if matches.sum() != 1:
            raise ValueError(
                f"{correction_id} expected exactly one line {line_id}; found {matches.sum()}"
            )

        if action == "exclude_line":
            corrected = corrected.loc[~matches].copy()
            corrected_ids = corrected["transmission_line_id"].map(_line_id)
            log.warning("%s excluded line %s through the correction register", correction_id, line_id)
            continue

        if action == "replace_value":
            field = str(correction["field"])
            if field not in corrected.columns:
                raise ValueError(f"{correction_id} references missing field {field}")
            observed = corrected.loc[matches, field].iloc[0]
            expected = correction["original_value"]
            if pd.notna(expected):
                numeric_match = False
                try:
                    numeric_match = bool(np.isclose(float(observed), float(expected)))
                except (TypeError, ValueError):
                    pass
                if not numeric_match and str(observed).strip() != str(expected).strip():
                    raise ValueError(
                        f"{correction_id} source drift: expected {field}={expected}, "
                        f"found {observed}"
                    )
            corrected.loc[matches, field] = correction["replacement_value"]
            log.warning(
                "%s replaced line %s %s: %s -> %s",
                correction_id,
                line_id,
                field,
                observed,
                correction["replacement_value"],
            )
            continue

        if action not in {"retain_boundary_intertie", "retain_after_node_patch"}:
            raise ValueError(f"{correction_id} has unsupported action {action}")

    return corrected.reset_index(drop=True)


def write_correction_audit(
    node_register: pd.DataFrame,
    line_register: pd.DataFrame,
    buses: pd.DataFrame,
    lines: pd.DataFrame,
    output_path: str | Path,
) -> pd.DataFrame:
    """Verify registered treatments are reflected in prepared tables."""
    rows: list[dict[str, object]] = []
    bus_names = buses["name"].dropna().astype(str)
    prepared_line_ids = lines["transmission_line_id"].map(_line_id)

    for _, correction in node_register.iterrows():
        suffix = "_" + "_".join(str(correction["node_code"]).split("_")[1:])
        observed = sorted(bus_names[bus_names.str.endswith(suffix)].tolist())
        rows.append(
            {
                "correction_id": correction["correction_id"],
                "record_type": "node",
                "record_id": correction["node_code"],
                "action": correction["action"],
                "status": "PASS" if observed else "FAIL",
                "observed": ";".join(observed),
                "source": correction["source"],
                "rationale": correction["rationale"],
            }
        )

    for _, correction in line_register.iterrows():
        line_id = _line_id(correction["transmission_line_id"])
        selected = lines.loc[prepared_line_ids.eq(line_id)]
        action = str(correction["action"])
        status = "PASS"
        observed = "excluded" if selected.empty else "retained"
        if action == "exclude_line":
            status = "PASS" if selected.empty else "FAIL"
        elif selected.empty:
            status = "FAIL"
        elif action == "replace_value":
            field = str(correction["field"])
            observed_value = selected[field].iloc[0]
            observed = f"{field}={observed_value}"
            try:
                matches = np.isclose(
                    float(observed_value), float(correction["replacement_value"])
                )
            except (TypeError, ValueError):
                matches = str(observed_value) == str(correction["replacement_value"])
            status = "PASS" if matches else "FAIL"
        rows.append(
            {
                "correction_id": correction["correction_id"],
                "record_type": "line",
                "record_id": line_id,
                "action": action,
                "status": status,
                "observed": observed,
                "source": correction["source"],
                "rationale": correction["rationale"],
            }
        )

    audit = pd.DataFrame(rows)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, index=False)
    if audit["status"].ne("PASS").any():
        failed = audit.loc[audit["status"].ne("PASS"), "correction_id"].tolist()
        raise ValueError(f"Prepared network failed correction audit: {failed}")
    return audit
