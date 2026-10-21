"""M1a acceptance test 7: SigMF write -> read round trip and metadata placement."""

import json
import struct
from datetime import datetime, timezone

import numpy as np
import pytest

from ntn_emu import sigmf_io

FS = 7.68e6
FC = 2.0e9
T0 = datetime(2026, 1, 2, 3, 4, 5, 678901, tzinfo=timezone.utc)


def awkward_samples():
    """Random complex64 plus +-0, subnormals, extremes, inf and NaN (compared as raw bytes)."""
    rng = np.random.default_rng(7)
    x = (rng.standard_normal(10_000) + 1j * rng.standard_normal(10_000)).astype(np.complex64)
    f32 = np.finfo(np.float32)
    specials = np.array([0.0, -0.0, f32.smallest_subnormal, -f32.smallest_subnormal,
                         f32.tiny, f32.max, -f32.max, np.inf, -np.inf, np.nan], dtype=np.float32)
    x[:len(specials)].real = specials
    x[:len(specials)].imag = specials[::-1]
    return x


def write_example(base, x, **kw):
    sigmf_io.write(base, x, fs=FS, fc=FC, datetime_utc=T0,
                   ntn={"tau0_s": 2e-3, "tau0_samples": 15_360, "tle_name": "TEST"},
                   annotations=[{"core:sample_start": 100, "core:sample_count": 512,
                                 "core:label": "PSS"}], **kw)


def test_round_trip_bit_identical(tmp_path):
    """Test 7: write -> read is bit-identical, little-endian cf32 on disk."""
    x = awkward_samples()
    base = tmp_path / "rt"
    write_example(base, x)
    y, _ = sigmf_io.read(base)
    assert y.dtype == np.complex64
    assert y.tobytes() == x.tobytes()
    raw = (tmp_path / "rt.sigmf-data").read_bytes()
    assert len(raw) == 8 * len(x)
    assert raw[8 * 20:8 * 21] == struct.pack("<ff", x[20].real, x[20].imag)


def test_metadata_placement(tmp_path):
    """Test 7: core:datatype == cf32_le; core:frequency in captures, absent from global."""
    base = tmp_path / "meta"
    write_example(base, awkward_samples())
    meta = json.loads((tmp_path / "meta.sigmf-meta").read_text())
    g, cap = meta["global"], meta["captures"]
    assert g["core:datatype"] == "cf32_le"
    assert g["core:sample_rate"] == FS
    assert "core:frequency" not in g
    assert cap[0]["core:frequency"] == FC
    assert cap[0]["core:sample_start"] == 0
    assert cap[0]["core:datetime"] == "2026-01-02T03:04:05.678901Z"
    for ext in g["core:extensions"]:
        assert isinstance(ext["name"], str)
        assert isinstance(ext["version"], str)
        assert isinstance(ext["optional"], bool)
    assert any(e["name"] == "ntn" for e in g["core:extensions"])
    assert g["ntn:tau0_samples"] == 15_360 and g["ntn:tle_name"] == "TEST"
    assert "ntn:git_sha" in g
    assert meta["annotations"][0]["core:label"] == "PSS"


def test_complex128_is_rounded_to_cf32(tmp_path):
    x = np.exp(1j * np.linspace(0, 10, 1000))
    write_example(tmp_path / "c128", x)
    y, _ = sigmf_io.read(tmp_path / "c128")
    assert y.tobytes() == x.astype(np.complex64).tobytes()


def test_read_rejects_corrupt_data(tmp_path):
    base = tmp_path / "bad"
    write_example(base, awkward_samples())
    data = tmp_path / "bad.sigmf-data"
    raw = bytearray(data.read_bytes())
    raw[0] ^= 1
    data.write_bytes(bytes(raw))
    with pytest.raises(ValueError, match="sha512"):
        sigmf_io.read(base)


def test_read_rejects_other_datatype(tmp_path):
    base = tmp_path / "dt"
    write_example(base, awkward_samples())
    p = tmp_path / "dt.sigmf-meta"
    meta = json.loads(p.read_text())
    meta["global"]["core:datatype"] = "ci16_le"
    p.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="datatype"):
        sigmf_io.read(base)
