"""Orbit geometry: acceptance test 6 (light-time) plus closed-form and definitional checks.

- Pinned TLEs carry valid mod-10 checksums.
- The derivative routine reproduces the analytic rho_dot of a synthetic circular orbit.
- rho_dot is insensitive to the stencil step (h = 0.1 s vs 0.2 s) on a real pass, within a
  limit derived from the test-5 budget.
- The light-time solution satisfies its own definition: |r_sat(t_rx - rho/c) - r_ue(t_rx)| = rho.
- Pinned passes: elevation at AOS/LOS equals the mask, the max-elevation and TCA epochs are
  extrema, and the three passes are within 3 days of the TLE epoch with distinct max elevations.

Measurement functions return the metric; tests assert on it, and scripts/measure_m1a.py prints it.
"""

import numpy as np
import pytest

from conftest import PASS_LABELS
from ntn_emu.constants import C_LIGHT
from ntn_emu.orbit import geometry, passes, tle

MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6_378_137.0      # m


def mod10(line):
    return sum(int(ch) if ch.isdigit() else (1 if ch == "-" else 0) for ch in line[:68]) % 10


def pass_offsets(p, step_s):
    duration = np.floor((p["los"] - p["aos"]) * 86_400.0)
    return np.arange(0.0, duration + step_s / 2, step_s)


def circular_diff_error(h=geometry.RATE_STEP_S):
    """max |central_diff4(rho) - analytic rho_dot| for a circular 550 km orbit over a fixed UE.

    rho(t) = sqrt(R^2 + Re^2 - 2 R Re cos(w t)),  rho_dot = R Re w sin(w t) / rho.
    """
    r_orb = R_EARTH + 550e3
    w = np.sqrt(MU_EARTH / r_orb**3)

    def rho(t):
        return np.sqrt(r_orb**2 + R_EARTH**2 - 2.0 * r_orb * R_EARTH * np.cos(w * t))

    t = np.linspace(-300.0, 300.0, 601)
    exact = r_orb * R_EARTH * w * np.sin(w * t) / rho(t)
    return np.max(np.abs(geometry.central_diff4(rho, t, h) - exact))


def rate_step_diffs(sat, ue, p, mode, h1=0.1, h2=0.2):
    """|rho_dot(h1) - rho_dot(h2)| over the whole pass at 0.25 s steps (no exclusions)."""
    fn = geometry.range_geometric if mode == "geometric" else geometry.range_lighttime
    t = geometry.times_from_offsets(p["aos"], pass_offsets(p, 0.25))
    return np.abs(fn(sat, ue, t, rate_step_s=h1).rho_dot - fn(sat, ue, t, rate_step_s=h2).rho_dot)


def lighttime_first_order_error(sat, ue, p):
    """Test 6: max |(rho_lt - rho_geo) - (-(v_sat . u_hat) tau)| over the pass at 1 s steps.

    The first-order predictor uses SGP4 velocity; its ~1 cm/s inconsistency moves the
    prediction by ~1e-2 m/s * tau (< 7 ms) < 1e-4 m, far below the 1 cm limit.
    """
    t = geometry.times_from_offsets(p["aos"], pass_offsets(p, 1.0))
    geo = geometry.range_geometric(sat, ue, t)
    lt = geometry.range_lighttime(sat, ue, t)
    predicted = -np.sum(geo.v_sat * geo.u_hat, axis=0) * geo.rho / C_LIGHT
    return np.max(np.abs((lt.rho - geo.rho) - predicted))


def lighttime_definition_error(sat, ue, p):
    """max | |r_sat(t_rx - rho/c) - r_ue(t_rx)| - rho | recomputed from scratch (7 s steps)."""
    t = geometry.times_from_offsets(p["aos"], pass_offsets(p, 7.0))
    lt = geometry.range_lighttime(sat, ue, t)
    t_tx = geometry.times_from_offsets(t, -lt.rho / C_LIGHT)
    r_sat = sat.at(t_tx).position.m
    r_ue = ue.at(t).position.m
    return np.max(np.abs(np.linalg.norm(r_sat - r_ue, axis=0) - lt.rho))


@pytest.mark.parametrize("which", ["starlink", "iss"])
def test_pinned_tle_checksums(which):
    text = tle.pinned_path(which).read_text(encoding="ascii").splitlines()
    assert len(text) == 3
    for line, tag in ((text[1], "1"), (text[2], "2")):
        assert len(line) == 69 and line[0] == tag
        assert mod10(line) == int(line[68])


def test_central_diff4_circular_orbit():
    """Closed form: derivative routine at the production step, error < 1 um/s."""
    assert circular_diff_error() < 1e-6


# Budget-derived limit. Test 5 allows 1 mm/s of rho_dot error, which is
# f_c/c * 1e-3 m/s = 6.7e-3 Hz of Doppler at 2 GHz. The stencil-step sensitivity may use 10% of
# that tolerance: 1e-4 m/s (6.7e-4 Hz). The real track is not smooth to 1 um/s: SGP4
# Kepler-solver steps of several um and Julian-date time quantization (both measured in
# scripts/measure_m1a.py) dominate, so the derivative routine is checked separately in
# closed form (test_central_diff4_circular_orbit, < 1e-6 m/s).
RATE_STEP_LIMIT_M_S = 1e-4


@pytest.mark.parametrize("mode", ["geometric", "lighttime"])
def test_rate_step_consistency(mode, starlink, ue, pinned_passes):
    """h = 0.1 s vs 0.2 s on the high pass, all epochs: max difference < 1e-4 m/s."""
    diffs = rate_step_diffs(starlink.sat, ue, pinned_passes["high"], mode)
    assert np.max(diffs) < RATE_STEP_LIMIT_M_S


@pytest.mark.parametrize("label", PASS_LABELS)
def test_lighttime_first_order(label, starlink, ue, pinned_passes):
    """Test 6: rho_lt - rho_geo = -(v_sat,GCRS . u_hat) * tau within 1 cm over the pass."""
    assert lighttime_first_order_error(starlink.sat, ue, pinned_passes[label]) < 0.01


@pytest.mark.parametrize("label", PASS_LABELS)
def test_lighttime_solution_definition(label, starlink, ue, pinned_passes):
    """The light-time range satisfies its defining equation within 1 um."""
    assert lighttime_definition_error(starlink.sat, ue, pinned_passes[label]) < 1e-6


def test_pinned_pass_set(starlink, pinned_passes):
    """Three passes of the main TLE, within 3 days of its epoch, low/mid/high max elevation."""
    doc = pinned_passes["_doc"]
    assert doc["tle_file"] == tle.pinned_path("starlink").relative_to(tle.REF_DIR).as_posix()
    bands = {"low": (20.0, 35.0), "mid": (45.0, 60.0), "high": (75.0, 90.0)}
    for label, (lo, hi) in bands.items():
        p = pinned_passes[label]
        assert lo <= p["max_el_deg"] <= hi
        assert abs(p["tca"] - starlink.sat.epoch) < 3.0
        assert p["aos"] < p["t_max_el"] < p["los"] and p["aos"] < p["tca"] < p["los"]


@pytest.mark.parametrize("label", PASS_LABELS)
def test_pass_epochs_are_definitional(label, starlink, ue, pinned_passes):
    p = pinned_passes[label]
    mask = p["mask_deg"]
    el = lambda t: passes.elevation_deg(starlink.sat, ue, t)  # noqa: E731
    assert abs(el(p["aos"]) - mask) < 0.01
    assert abs(el(p["los"]) - mask) < 0.01
    around = geometry.times_from_offsets(p["t_max_el"], np.array([-1.0, 0.0, 1.0]))
    e = el(around)
    assert e[0] < e[1] and e[2] < e[1]
    assert abs(e[1] - p["max_el_deg"]) < 1e-6
    around = geometry.times_from_offsets(p["tca"], np.array([-1.0, 0.0, 1.0]))
    g = geometry.range_geometric(starlink.sat, ue, around)
    assert g.rho[0] > g.rho[1] and g.rho[2] > g.rho[1]
    assert abs(g.rho_dot[1]) < 1e-3
