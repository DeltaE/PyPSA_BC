"""
prepare_datahandovertopypsa - stage 4 of the PyPSA-BC workflow.

Creates PyPSA component dictionaries from prepared assets and profiles by
orchestrating the enrich_format_* scripts in staged order:

    1) VRE dicts   (wind + solar)
    2) TPP dicts   (thermal + fossil-fuel infrastructure)
    3) Hydro dicts (reservoir, ror, ror-water)

Behavior mirrors other prepare_* orchestrators:
    - preflight checks for required input artifacts
    - per-stage enable/disable switches
    - skip execution when outputs already exist (unless force_update=True)
    - staged progress bar + delivered-files table
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.logger import Pipeline


def _rows(path: Path) -> int | str:
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ""


def _artifact_note(path: Path) -> int | str:
    if not path.exists():
        return "missing"
    if path.suffix.lower() == ".csv":
        return _rows(path)
    return f"{path.stat().st_size / 1e6:.2f} MB"


def _dict_path(out_cfg: dict, key: str) -> Path:
    return Path(out_cfg["pypsa_dict"]["folder"]) / out_cfg["pypsa_dict"][key]


def _all_exist(paths: list[Path]) -> bool:
    return all(p.exists() for p in paths)


def main(
    run_vre_dict: bool = True,
    run_tpp_dict: bool = True,
    run_hydro_dict: bool = True,
    force_update: bool = False,
):
    from workflow.scripts import enrich_format_hydro, enrich_format_tpp, enrich_format_vre

    out = AttributesParser().data_cfg["output"]

    stage_cfg = {
        "vre dict": {
            "enabled": run_vre_dict,
            "runner": enrich_format_vre.main,
            "outputs": [
                ("wind.pickle", _dict_path(out, "wind")),
                ("solar.pickle", _dict_path(out, "solar")),
            ],
            "needs": {
                "buses.csv": Path(out["base_network"]) / "buses.csv",
                "wind assets": Path(out["create_ext_wind_assets"]["fname"]),
                "wind ts": Path(out["create_ext_wind_ts"]["fname"]),
                "solar assets": Path(out["create_ext_solar_assets"]["fname"]),
                "solar ts": Path(out["create_ext_solar_ts"]["fname"]),
            },
        },
        "tpp dict": {
            "enabled": run_tpp_dict,
            "runner": enrich_format_tpp.main,
            "outputs": [
                ("tpp.pickle", _dict_path(out, "tpp")),
                ("ff_infra.pickle", _dict_path(out, "ff_infrastructure")),
            ],
            "needs": {
                "buses.csv": Path(out["base_network"]) / "buses.csv",
                "tpp assets": Path(out["create_ext_tpp_assets"]["fname"]),
            },
        },
        "hydro dict": {
            "enabled": run_hydro_dict,
            "runner": enrich_format_hydro.main,
            "outputs": [
                ("hydro_reservoirs.pickle", _dict_path(out, "res")),
                ("hydro_ror.pickle", _dict_path(out, "ror")),
                ("hydro_ror_water.pickle", _dict_path(out, "ror_water")),
            ],
            "needs": {
                "buses.csv": Path(out["base_network"]) / "buses.csv",
                "hydro generation": Path(out["create_hydro_assets"]["hydro_generation"]),
                "hydro reservoirs": Path(out["create_hydro_assets"]["hydro_reservoir"]),
                "reservoir inflows": Path(out["reservoir_inflows"]["fname"]),
                "ror power": Path(out["ror_ps"]["fname"]),
            },
        },
    }

    missing: dict[str, str] = {}
    for label, cfg in stage_cfg.items():
        if not cfg["enabled"]:
            continue

        out_paths = [p for _, p in cfg["outputs"]]
        if _all_exist(out_paths) and not force_update:
            continue

        for name, path in cfg["needs"].items():
            if not path.exists():
                missing[f"{label}:{name}"] = str(path)

    if missing:
        raise FileNotFoundError(
            "prepare_datahandovertopypsa: missing required inputs for selected stages. "
            f"Missing: {missing}"
        )

    stages = ["vre dict", "tpp dict", "hydro dict"]
    with Pipeline("Data handover to PyPSA", stages) as pipe:
        for label in stages:
            cfg = stage_cfg[label]
            outputs = cfg["outputs"]
            out_paths = [p for _, p in outputs]

            with pipe.stage(label):
                if not cfg["enabled"]:
                    for artifact, path in outputs:
                        pipe.deliver(artifact, path, "skipped (disabled)")
                    continue

                if _all_exist(out_paths) and not force_update:
                    for artifact, path in outputs:
                        pipe.deliver(artifact, path, f"{_artifact_note(path)} (exists locally)")
                    continue

                cfg["runner"]()
                for artifact, path in outputs:
                    pipe.deliver(artifact, path, _artifact_note(path))


if __name__ == "__main__":
    main()
