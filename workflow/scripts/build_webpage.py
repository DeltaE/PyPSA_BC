"""Rebuild the static PyPSA-BC documentation website.

The site is authored as a single template at ``docs/site-src/index.template.html``
and published to ``docs/index.html``. Brand assets are copied from a logos
directory (overridable with the ``PYPSA_BC_LOGOS`` environment variable) into
``docs/assets``. Missing logos are reported as warnings rather than fatal errors
so the HTML can still be regenerated on machines without the brand assets.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

# repo root is two levels up from workflow/scripts/
ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
ASSETS = DOCS / "assets"
TEMPLATE = DOCS / "site-src" / "index.template.html"
OUTPUT = DOCS / "index.html"

# Default brand-asset location; override with PYPSA_BC_LOGOS if needed.
LOGOS = Path(os.environ.get("PYPSA_BC_LOGOS", r"E:\CoWork\PROJECTS\logos"))

# source filename -> published asset filename
LOGO_MAP = {
    "pypsa_bc_logo.svg": "pypsa-bc.svg",
    "DeltaE_white.png": "delta-e.png",
}


def _publish_template() -> None:
    if not TEMPLATE.exists():
        raise FileNotFoundError(f"Site template not found: {TEMPLATE}")
    OUTPUT.write_text(TEMPLATE.read_text(encoding="utf-8"), encoding="utf-8")


def _copy_logos() -> list[str]:
    """Copy available brand assets; return warnings for any that are missing."""
    warnings: list[str] = []
    for source, target in LOGO_MAP.items():
        src = LOGOS / source
        if src.exists():
            shutil.copy2(src, ASSETS / target)
        else:
            warnings.append(f"  - missing logo (kept existing if any): {src}")
    return warnings


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    _publish_template()
    warnings = _copy_logos()

    print(f"Built {OUTPUT}")
    if warnings:
        print("Warning: some brand assets were not copied:")
        print("\n".join(warnings))


if __name__ == "__main__":
    main()
