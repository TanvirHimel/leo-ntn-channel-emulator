"""Pin three passes of the main (Starlink-class) TLE to tests/reference/passes.json.

The passes are scenario inputs (window epochs for tests, demos and M5), not expected values.
Search window: TLE epoch +- 3 days, mask 10 deg, UE = rounded constants. In each
max-elevation band (low 20-35, mid 45-60, high 75-90 deg) the pass whose TCA is closest to
the TLE epoch is taken. Stops if a band is empty.

    .\\.venv\\Scripts\\python.exe scripts\\pin_passes.py
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from ntn_emu.constants import UE_ALT_M, UE_LAT_DEG, UE_LON_DEG  # noqa: E402
from ntn_emu.orbit import geometry, passes, tle  # noqa: E402

MASK_DEG = 10.0
SEARCH_DAYS = 3.0
BANDS = {"low": (20.0, 35.0), "mid": (45.0, 60.0), "high": (75.0, 90.0)}


def main():
    main_tle = tle.load_pinned("starlink")
    sat, ue = main_tle.sat, geometry.ue_position()
    epoch = sat.epoch
    t0 = geometry.times_from_offsets(epoch, -SEARCH_DAYS * 86_400.0)
    t1 = geometry.times_from_offsets(epoch, SEARCH_DAYS * 86_400.0)
    found = passes.find_passes(sat, ue, t0, t1, mask_deg=MASK_DEG)
    print(f"{len(found)} complete passes above {MASK_DEG} deg within +-{SEARCH_DAYS} d of epoch")

    chosen = []
    for label, (lo, hi) in BANDS.items():
        band = [p for p in found if lo <= p.max_el_deg <= hi]
        if not band:
            sys.exit(f"no pass with max elevation in [{lo}, {hi}] deg; stopping")
        p = min(band, key=lambda q: abs(q.tca - epoch))
        chosen.append((label, p))
        print(f"  {label:4s}: max el {p.max_el_deg:6.2f} deg, TCA {geometry.time_to_iso(p.tca)}, "
              f"{(p.los - p.aos) * 86_400:.0f} s above mask, {len(band)} in band")

    doc = {
        "description": "Pinned passes of the main TLE: scenario epochs (inputs), not expected values.",
        "generated_by": "scripts/pin_passes.py",
        "tle_file": main_tle.path.relative_to(tle.REF_DIR).as_posix(),
        "tle_name": main_tle.name,
        "tle_epoch_utc": geometry.time_to_iso(epoch),
        "ue_lat_deg": UE_LAT_DEG, "ue_lon_deg": UE_LON_DEG, "ue_alt_m": UE_ALT_M,
        "search": {"start_utc": geometry.time_to_iso(t0), "stop_utc": geometry.time_to_iso(t1),
                   "rule": "per band, pass with TCA closest to the TLE epoch"},
        "passes": [{
            "label": label,
            "mask_deg": p.mask_deg,
            "max_el_deg": round(p.max_el_deg, 6),
            "aos_utc": geometry.time_to_iso(p.aos),
            "t_max_el_utc": geometry.time_to_iso(p.t_max_el),
            "tca_utc": geometry.time_to_iso(p.tca),
            "los_utc": geometry.time_to_iso(p.los),
        } for label, p in chosen],
    }
    out = tle.REF_DIR / "passes.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
