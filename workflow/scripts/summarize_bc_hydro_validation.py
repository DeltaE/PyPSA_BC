"""Summarize pinned BC Hydro study-year validation workbooks.

The output is descriptive evidence only. Directional TTC is not realized flow,
and the sign of published actual interchange is deliberately left uninterpreted
until it is documented in the validation protocol.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _series_summary(series: pd.Series) -> dict[str, float | int]:
    values = pd.to_numeric(series, errors="coerce")
    valid = values.dropna()
    return {
        "rows": len(values),
        "valid_rows": len(valid),
        "missing_rows": int(values.isna().sum()),
        "minimum_mw": float(valid.min()),
        "p05_mw": float(valid.quantile(0.05)),
        "median_mw": float(valid.median()),
        "mean_mw": float(valid.mean()),
        "p95_mw": float(valid.quantile(0.95)),
        "maximum_mw": float(valid.max()),
    }


def summarize(year: int, root: Path) -> dict[str, object]:
    load_path = root / f"BalancingAuthorityLoad{year}.xls"
    flow_path = root / f"HourlyTielineData{year}.xls"
    ttc_paths = {
        "bc_to_ab": root / f"BCABTTC{year}.xls",
        "ab_to_bc": root / f"ABBCTTC{year}.xls",
        "bc_to_us": root / f"BCUSTTC{year}.xls",
        "us_to_bc": root / f"USBCTTC{year}.xls",
    }
    required = [load_path, flow_path, *ttc_paths.values()]
    missing = [path.as_posix() for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing pinned validation files: " + ", ".join(missing))

    load_frame = pd.read_excel(load_path, header=1)
    load = pd.to_numeric(load_frame.iloc[:, 2], errors="coerce")
    peak_index = load.idxmax()
    peak_row = load_frame.loc[peak_index]

    ttc = {}
    for direction, path in ttc_paths.items():
        frame = pd.read_excel(path)
        summary = _series_summary(frame["TTC"])
        summary.update({"path": path.as_posix(), "sha256": _sha256(path)})
        ttc[direction] = summary

    flow_frame = pd.read_excel(flow_path, header=1)
    actual_flow = {}
    for column in flow_frame.columns[2:]:
        summary = _series_summary(flow_frame[column])
        summary["signed_sum_gwh"] = float(
            pd.to_numeric(flow_frame[column], errors="coerce").sum() / 1_000
        )
        actual_flow[str(column)] = summary

    return {
        "year": year,
        "claim_boundaries": {
            "load": "Gross telemetered provincial load is not regional demand.",
            "ttc": "Directional TTC is a path limit, not realized flow or an internal line rating.",
            "actual_flow": "Signed interchange is reported without assigning import/export meaning until the publisher sign convention is documented.",
        },
        "load": {
            **_series_summary(load),
            "annual_energy_twh": float(load.sum() / 1_000_000),
            "peak_date": str(peak_row.iloc[0]),
            "peak_hour_ending": int(peak_row.iloc[1]),
            "path": load_path.as_posix(),
            "sha256": _sha256(load_path),
        },
        "directional_ttc": ttc,
        "actual_interchange": {
            "path": flow_path.as_posix(),
            "sha256": _sha256(flow_path),
            "series": actual_flow,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2021)
    parser.add_argument("--root", type=Path, default=Path("data/validation/bc_hydro"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.root / f"summary_{args.year}.json"
    payload = summarize(args.year, args.root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"BC Hydro validation summary: {output}")


if __name__ == "__main__":
    main()
