"""Pass finder: AOS/LOS at an elevation mask, max-elevation epoch and geometric TCA.

Skyfield find_events (geometric, no refraction) gives first guesses; its timing is too
coarse (0.02 deg at the mask was seen), so AOS/LOS are refined as roots of elevation - mask,
the max-elevation epoch by a bounded scalar search, and TCA (minimum geometric range) as a
root of rho_dot, all to ~1 us. load_pinned_passes() reads tests/reference/passes.json.
"""

import json
from dataclasses import dataclass

from scipy.optimize import brentq, minimize_scalar

from ntn_emu.orbit.geometry import range_geometric, time_from_iso, times_from_offsets
from ntn_emu.orbit.tle import REF_DIR


@dataclass
class Pass:
    aos: object       # Skyfield Time
    t_max_el: object
    tca: object       # minimum geometric range
    los: object
    max_el_deg: float
    mask_deg: float


def elevation_deg(sat, ue, t):
    alt, _, _ = (sat - ue).at(t).altaz()
    return alt.degrees


def _refine_max_el(sat, ue, t_guess, half_s=20.0):
    res = minimize_scalar(lambda s: -elevation_deg(sat, ue, times_from_offsets(t_guess, s)),
                          bounds=(-half_s, half_s), method="bounded", options={"xatol": 1e-6})
    t = times_from_offsets(t_guess, res.x)
    return t, float(elevation_deg(sat, ue, t))


def _refine_crossing(sat, ue, t_guess, mask_deg, half_s=5.0):
    def f(s):
        return float(elevation_deg(sat, ue, times_from_offsets(t_guess, s))) - mask_deg
    return times_from_offsets(t_guess, brentq(f, -half_s, half_s, xtol=1e-7))


def _tca(sat, ue, t_guess, half_s=120.0):
    def rate(s):
        return float(range_geometric(sat, ue, times_from_offsets(t_guess, s)).rho_dot)
    return times_from_offsets(t_guess, brentq(rate, -half_s, half_s, xtol=1e-7))


def find_passes(sat, ue, t0, t1, mask_deg=10.0):
    """Complete passes (rise, culminate, set) above mask_deg between t0 and t1."""
    times, events = sat.find_events(ue, t0, t1, altitude_degrees=mask_deg)
    out = []
    aos = culm = None
    for t, ev in zip(times, events):
        if ev == 0:
            aos, culm = t, None
        elif ev == 1 and aos is not None:
            culm = t
        elif ev == 2 and aos is not None and culm is not None:
            t_max, el_max = _refine_max_el(sat, ue, culm)
            out.append(Pass(_refine_crossing(sat, ue, aos, mask_deg), t_max,
                            _tca(sat, ue, t_max), _refine_crossing(sat, ue, t, mask_deg),
                            el_max, mask_deg))
            aos = culm = None
    return out


def load_pinned_passes(path=REF_DIR / "passes.json"):
    """{label: entry} with Skyfield times for aos/t_max_el/tca/los, plus '_doc' (raw JSON)."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    out = {"_doc": doc}
    for entry in doc["passes"]:
        e = dict(entry)
        for key in ("aos", "t_max_el", "tca", "los"):
            e[key] = time_from_iso(entry[f"{key}_utc"])
        out[entry["label"]] = e
    return out
