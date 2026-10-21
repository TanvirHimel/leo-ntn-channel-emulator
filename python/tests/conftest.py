"""Shared fixtures: pinned Starlink TLE, rounded UE, pinned passes (tests/reference/)."""

import pytest

from ntn_emu.orbit import geometry, passes, tle

PASS_LABELS = ["low", "mid", "high"]


@pytest.fixture(scope="session")
def starlink():
    return tle.load_pinned("starlink")


@pytest.fixture(scope="session")
def ue():
    return geometry.ue_position()


@pytest.fixture(scope="session")
def pinned_passes():
    return passes.load_pinned_passes()
