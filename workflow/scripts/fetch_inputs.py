"""
fetch_inputs — download public raw inputs listed under `remote:` in data.yaml.

    python -m workflow.scripts.fetch_inputs            # fetch all (skip existing)
    fetch_inputs.main(only=["can_turbines"])           # a subset
    fetch_inputs.main(force=True)                       # re-download everything

Shows the staged checklist + a delivered-files table; detail in logs/pypsa_bc.log.
"""

from pathlib import Path

from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.downloader import download_file
from pypsa_bc.reporting.logger import Pipeline


def main(only: list[str] | None = None, force: bool = False):
    remote = AttributesParser().data_cfg.get("remote", {})
    names = list(remote) if only is None else [n for n in only if n in remote]
    if not names:
        print("fetch_inputs: nothing to download (no matching 'remote:' entries).")
        return {}

    got: dict[str, Path] = {}
    with Pipeline("Fetch inputs", names) as pipe:
        for name in names:
            with pipe.stage(name):
                spec = remote[name]
                p = download_file(spec["url"], spec["dest"], force=force)
                size = f"{p.stat().st_size/1e6:.1f} MB" if p.exists() else "—"
                pipe.deliver(p.name, p, size)
                got[name] = p
    return got


if __name__ == "__main__":
    main()
