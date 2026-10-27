"""Import-time CPU check of ntn_emu._core (SPEC §6, M2 step 5).

Module init calls ntn::require_avx2() before anything else and raises ImportError naming the
missing features. The bindings are compiled without /arch:AVX2, so this check runs before any
AVX2 instruction. On a capable machine the test-only switch NTN_TEST_HIDE_AVX2 (cpu_check.hpp)
hides one feature; each case runs `import ntn_emu._core` in a fresh interpreter.
"""

import os
import subprocess
import sys

import pytest

from ntn_emu import _core

REQUIRED = "required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV)"
HIDE_CASES = [
    ("1", "AVX2"),
    ("AVX2", "AVX2"),
    ("FMA3", "FMA3"),
    ("BMI1", "BMI1"),
    ("BMI2", "BMI2"),
    ("OSAVX", "OS AVX state (XGETBV)"),
]

pytestmark = pytest.mark.skipif(not _core.cpu_has_avx2(), reason="CPU without the AVX2 feature set")


def import_core(hide=None):
    env = dict(os.environ)
    env.pop("NTN_TEST_HIDE_AVX2", None)
    if hide is not None:
        env["NTN_TEST_HIDE_AVX2"] = hide
    return subprocess.run([sys.executable, "-c", "import ntn_emu._core"], env=env, capture_output=True,
                          text=True, timeout=60)


@pytest.mark.parametrize("hide,name", HIDE_CASES)
def test_import_fails_with_feature_list(hide, name):
    r = import_core(hide)
    assert r.returncode != 0
    assert "ImportError" in r.stderr, r.stderr
    assert f"missing: {name};" in r.stderr, r.stderr
    assert REQUIRED in r.stderr, r.stderr


@pytest.mark.parametrize("value", ["0", "fma3", "FMA3 ", "1,FMA3", "AVX512"])
def test_import_ignores_other_switch_values(value):
    r = import_core(value)
    assert r.returncode == 0, r.stderr


def test_import_without_switch():
    r = import_core()
    assert r.returncode == 0, r.stderr
