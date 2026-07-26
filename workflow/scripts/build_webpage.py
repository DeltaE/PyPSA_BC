"""Rebuild the static PyPSA-BC documentation website."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
ASSETS = DOCS / "assets"
LOGOS = Path(r"E:\CoWork\PROJECTS\logos")


def main() -> None:
    ASSETS.mkdir(exist_ok=True)
    template = DOCS / "site-src" / "index.template.html"
    (DOCS / "index.html").write_text(
        template.read_text(encoding="utf-8"), encoding="utf-8"
    )
    for source, target in {
        "pypsa_bc_logo.svg": "pypsa-bc.svg",
        "DeltaE_white.png": "delta-e.png",
    }.items():
        src = LOGOS / source
        if not src.exists():
            raise FileNotFoundError(f"Required logo not found: {src}")
        shutil.copy2(src, ASSETS / target)
    print(f"Built {DOCS / 'index.html'}")


if __name__ == "__main__":
    main()
