"""Satellite-UE range in GCRS: geometric and light-time (SPEC §2 "Geometry").

Light-time range at receive time t_rx:
    rho(t_rx) = |r_sat(t_rx - tau) - r_ue(t_rx)|,   tau = rho / c,
solved by fixed-point iteration (contraction ~|v_sat|/c ~ 2.5e-5) until |d rho| < 1 um.
The geometric mode uses t_tx = t_rx.

rho_dot is the 4th-order central difference of that same range function (light-time included),
in both modes. SGP4 velocity is never used for rho_dot: it is inconsistent with the
derivative of SGP4 position by ~1 cm/s (measured in scripts/measure_m1a.py), which would put
range and Doppler on different tracks. The stencil samples the orbit up to 2h beyond the
requested times, so knot slopes at window edges are valid.

Times are Skyfield Time objects; offsets are applied to the two-part TT Julian date
(whole + fraction), giving ~5e-12 s time resolution.
"""

from dataclasses import dataclass
from datetime import datetime

import numpy as np
from skyfield.api import wgs84

from ntn_emu.constants import C_LIGHT, UE_ALT_M, UE_LAT_DEG, UE_LON_DEG
from ntn_emu.orbit.tle import timescale

LIGHTTIME_TOL_M = 1e-6
LIGHTTIME_MAX_ITER = 10
# Stencil step for rho_dot. The closed-form circular-orbit test gives 1.3e-7 m/s at h = 0.2 s.
# On the real track, h = 0.1 s vs 0.2 s differ by ~5e-7 m/s (median) and ~2e-5 m/s (max).
# The max is set by SGP4 Kepler-solver steps (several um of range, see measure_m1a.py), and
# time-quantization noise of the two-part Julian date scales as 1/h. The larger step,
# h = 0.2 s, keeps both small.
RATE_STEP_S = 0.2


def ue_position(lat_deg=UE_LAT_DEG, lon_deg=UE_LON_DEG, alt_m=UE_ALT_M):
    return wgs84.latlon(lat_deg, lon_deg, elevation_m=alt_m)


def times_from_offsets(t0, offsets_s):
    """t0 + offsets (seconds, scalar or array), on the two-part TT Julian date."""
    return timescale().tt_jd(t0.whole, t0.tt_fraction + np.asarray(offsets_s) / 86_400.0)


def time_from_iso(text):
    """Parse an ISO-8601 UTC timestamp ('...Z', microseconds) into a Skyfield Time."""
    return timescale().from_datetime(datetime.fromisoformat(text.replace("Z", "+00:00")))


def time_to_iso(t):
    return t.utc_datetime().strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def central_diff4(fn, x, h):
    """4th-order central difference of fn at x with step h."""
    return (fn(x - 2 * h) - 8.0 * fn(x - h) + 8.0 * fn(x + h) - fn(x + 2 * h)) / (12.0 * h)


@dataclass
class RangeState:
    rho: np.ndarray      # m
    rho_dot: np.ndarray  # m/s, derivative of the rho track (never SGP4 velocity)
    u_hat: np.ndarray    # 3 x N, UE -> satellite (at transmit position)
    v_sat: np.ndarray    # 3 x N, SGP4 GCRS velocity at transmit time. Informational only.


def _range(sat, ue, t_rx, lighttime):
    """rho at receive times t_rx, with unit vector and SGP4 velocity at the transmit time."""
    r_ue = ue.at(t_rx).position.m
    s = sat.at(t_rx)
    d = s.position.m - r_ue
    rho = np.linalg.norm(d, axis=0)
    if lighttime:
        for _ in range(LIGHTTIME_MAX_ITER):
            s = sat.at(times_from_offsets(t_rx, -rho / C_LIGHT))
            d = s.position.m - r_ue
            new = np.linalg.norm(d, axis=0)
            converged = np.max(np.abs(new - rho)) < LIGHTTIME_TOL_M
            rho = new
            if converged:
                break
        else:
            raise RuntimeError("light-time iteration did not converge")
    return rho, d / rho, s.velocity.m_per_s


def _state(sat, ue, t, lighttime, rate_step_s):
    rho, u_hat, v_sat = _range(sat, ue, t, lighttime)

    def rho_at(offset_s):
        return _range(sat, ue, times_from_offsets(t, offset_s), lighttime)[0]

    rate = central_diff4(rho_at, 0.0, rate_step_s)
    return RangeState(rho, rate, u_hat, v_sat)


def range_geometric(sat, ue, t, rate_step_s=RATE_STEP_S):
    """Instantaneous range |r_sat(t) - r_ue(t)| (no light-time) and its derivative."""
    return _state(sat, ue, t, False, rate_step_s)


def range_lighttime(sat, ue, t_rx, rate_step_s=RATE_STEP_S):
    """Light-time range at receive times t_rx and its derivative."""
    return _state(sat, ue, t_rx, True, rate_step_s)
