"""Pinned TLEs from tests/reference/tle/ (fetched once by scripts/fetch_tle.py)."""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from skyfield.api import EarthSatellite, load

REF_DIR = Path(__file__).resolve().parents[3] / "tests" / "reference"


@lru_cache(maxsize=1)
def timescale():
    """Offline timescale (bundled leap-second and Delta-T data)."""
    return load.timescale(builtin=True)


@dataclass(frozen=True)
class PinnedTle:
    path: Path
    name: str
    norad: int
    sat: EarthSatellite


def pinned_path(which):
    """'starlink' -> the single tle/starlink_<norad>.txt; 'iss' -> tle/iss.txt."""
    if which == "iss":
        return REF_DIR / "tle" / "iss.txt"
    found = sorted((REF_DIR / "tle").glob(f"{which}_*.txt"))
    if len(found) != 1:
        raise FileNotFoundError(f"expected exactly one pinned {which}_*.txt, found {len(found)}")
    return found[0]


def load_tle(path):
    name, l1, l2 = Path(path).read_text(encoding="ascii").splitlines()[:3]
    sat = EarthSatellite(l1, l2, name.strip(), timescale())
    return PinnedTle(Path(path), name.strip(), int(l1[2:7]), sat)


def load_pinned(which):
    return load_tle(pinned_path(which))
