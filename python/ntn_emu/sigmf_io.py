"""Minimal SigMF writer/reader for cf32_le recordings (SPEC §5.3).

- Data: explicit little-endian complex64 ('<c8'). complex128 input is rounded to complex64.
- global: core:datatype "cf32_le", core:sample_rate, core:version, core:sha512,
  core:extensions declaring "ntn", and the ntn:* keys (git_sha/git_dirty are added).
- captures[0]: core:sample_start 0, core:frequency (= f_c, never in global), core:datetime.
- annotations: passed through, sorted by core:sample_start.
"""

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SIGMF_VERSION = "1.2.0"
NTN_EXT_VERSION = "0.1.0"
DATATYPE = "cf32_le"
_DTYPE = np.dtype("<c8")
_REPO = Path(__file__).resolve().parents[2]


def _paths(base):
    base = Path(base)
    return base.with_name(base.name + ".sigmf-data"), base.with_name(base.name + ".sigmf-meta")


def _git_state():
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO, capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO, capture_output=True,
                               text=True, check=True).stdout.strip() != ""
        return sha, dirty
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def _iso_utc(dt):
    if dt.tzinfo is None:
        raise ValueError("datetime_utc must be timezone-aware")
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def write(base, samples, *, fs, fc, datetime_utc: datetime, ntn=None, annotations=None):
    """Write <base>.sigmf-data and <base>.sigmf-meta. ntn keys are given without the prefix."""
    data_path, meta_path = _paths(base)
    data = np.ascontiguousarray(samples).astype(_DTYPE, copy=False).tobytes()

    sha, dirty = _git_state()
    glob = {
        "core:datatype": DATATYPE,
        "core:sample_rate": float(fs),
        "core:version": SIGMF_VERSION,
        "core:sha512": hashlib.sha512(data).hexdigest(),
        "core:recorder": "ntn_emu",
        "core:extensions": [{"name": "ntn", "version": NTN_EXT_VERSION, "optional": True}],
        "ntn:git_sha": sha,
        "ntn:git_dirty": dirty,
    }
    for key, value in (ntn or {}).items():
        if ":" in key:
            raise ValueError(f"ntn key {key!r} must be given without a namespace prefix")
        glob[f"ntn:{key}"] = value

    meta = {
        "global": glob,
        "captures": [{
            "core:sample_start": 0,
            "core:frequency": float(fc),
            "core:datetime": _iso_utc(datetime_utc),
        }],
        "annotations": sorted(annotations or [], key=lambda a: a["core:sample_start"]),
    }
    data_path.parent.mkdir(parents=True, exist_ok=True)
    data_path.write_bytes(data)
    meta_path.write_text(json.dumps(meta, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read(base):
    """Read a cf32_le recording; returns (complex64 samples, metadata dict)."""
    data_path, meta_path = _paths(base)
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    glob = meta["global"]
    if glob.get("core:datatype") != DATATYPE:
        raise ValueError(f"unsupported core:datatype {glob.get('core:datatype')!r}")
    data = data_path.read_bytes()
    if len(data) % _DTYPE.itemsize:
        raise ValueError("data file size is not a multiple of 8 bytes")
    if "core:sha512" in glob and hashlib.sha512(data).hexdigest() != glob["core:sha512"]:
        raise ValueError("core:sha512 mismatch: data file does not match metadata")
    return np.frombuffer(data, dtype=_DTYPE).astype(np.complex64, copy=True), meta
