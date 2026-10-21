"""Fetch and pin the reference TLEs (maintenance script; tests never fetch).

One request per source: the CelesTrak Starlink group (one satellite is selected from it) and
the ISS by NORAD id. Raw downloads are cached in --download-dir with a .fetch.json sidecar
(URL, UTC fetch time). A rerun reuses the cache and never downloads again.

Selection of the Starlink-class satellite (deterministic):
  inclination 53.0-53.3 deg, mean-motion altitude 530-560 km, eccentricity < 0.001,
  epoch within 2 days of the fetch; closest to 550 km, ties broken by lowest NORAD id.
  (The window was 540-560 km until the 2026-10-21 fetch showed no 53-deg shell above
  539.8 km; the 53-deg shells were lowered, most to ~463 km. The nominal 540 km shell is
  the closest to SPEC's ~550 km.)

Writes tests/reference/tle/{starlink_<norad>,iss}.txt and tests/reference/SOURCES.md.

    .\\.venv\\Scripts\\python.exe scripts\\fetch_tle.py --download-dir <empty scratch folder>
"""

import argparse
import hashlib
import json
import math
import sys
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "tests" / "reference"
MU_KM3_S2 = 398_600.4418
R_EARTH_KM = 6_378.137

SOURCES = {
    "starlink": "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle",
    "iss": "https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=tle",
}
ISS_NORAD = 25544


def fetch_once(name, url, dl_dir):
    raw_path = dl_dir / f"{name}.tle.txt"
    side_path = dl_dir / f"{name}.fetch.json"
    if raw_path.exists() and side_path.exists():
        print(f"[{name}] reusing cached download {raw_path}")
        return raw_path.read_bytes(), json.loads(side_path.read_text())
    fetched = datetime.now(timezone.utc)
    req = urllib.request.Request(url, headers={"User-Agent": "ntn-emulator-fetch_tle/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            status = resp.status
            body = resp.read()
    except OSError as exc:
        sys.exit(f"[{name}] FETCH FAILED for {url}: {exc}")
    if status != 200:
        sys.exit(f"[{name}] FETCH FAILED for {url}: HTTP {status}")
    side = {"url": url, "fetched_utc": fetched.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "http_status": status, "bytes": len(body),
            "sha256": hashlib.sha256(body).hexdigest()}
    raw_path.write_bytes(body)
    side_path.write_text(json.dumps(side, indent=2))
    print(f"[{name}] fetched {len(body)} bytes from {url}")
    return body, side


def checksum_ok(line):
    total = sum(int(ch) if ch.isdigit() else (1 if ch == "-" else 0) for ch in line[:68])
    return len(line) == 69 and line[68].isdigit() and total % 10 == int(line[68])


def parse_tles(body):
    """Yield (name, line1, line2) for every entry with valid structure and checksums."""
    text = body.decode("ascii")  # TLEs are ASCII; anything else is rejected
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    good, bad = [], 0
    for i in range(0, len(lines) - 2, 3):
        name, l1, l2 = lines[i].strip(), lines[i + 1], lines[i + 2]
        if (l1.startswith("1 ") and l2.startswith("2 ") and checksum_ok(l1) and checksum_ok(l2)
                and l1[2:7] == l2[2:7]):
            good.append((name, l1, l2))
        else:
            bad += 1
    if len(lines) % 3:
        bad += 1
    return good, bad


def elements(l1, l2):
    yy = int(l1[18:20])
    year = 2000 + yy if yy < 57 else 1900 + yy
    epoch = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=float(l1[20:32]) - 1.0)
    incl = float(l2[8:16])
    ecc = float("0." + l2[26:33])
    n_rad_s = float(l2[52:63]) * 2.0 * math.pi / 86_400.0
    a_km = (MU_KM3_S2 / n_rad_s**2) ** (1.0 / 3.0)
    return {"norad": int(l1[2:7]), "epoch": epoch, "incl": incl, "ecc": ecc,
            "alt_km": a_km - R_EARTH_KM}


def write_tle(path, name, l1, l2):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{name}\n{l1}\n{l2}\n", encoding="ascii", newline="\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--download-dir", type=Path, required=True)
    args = ap.parse_args()
    args.download_dir.mkdir(parents=True, exist_ok=True)

    pinned = []
    # Starlink-class: select one from the group.
    body, side = fetch_once("starlink", SOURCES["starlink"], args.download_dir)
    tles, bad = parse_tles(body)
    print(f"[starlink] {len(tles)} valid TLEs, {bad} rejected")
    fetched = datetime.strptime(side["fetched_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    cands = []
    for name, l1, l2 in tles:
        e = elements(l1, l2)
        if (53.0 <= e["incl"] <= 53.3 and 530.0 <= e["alt_km"] <= 560.0 and e["ecc"] < 1e-3
                and abs((fetched - e["epoch"]).total_seconds()) < 2 * 86_400):
            cands.append((abs(e["alt_km"] - 550.0), e["norad"], name, l1, l2, e))
    if not cands:
        sys.exit("[starlink] no TLE matches the selection criteria; stopping")
    cands.sort(key=lambda c: (c[0], c[1]))
    _, norad, name, l1, l2, e = cands[0]
    print(f"[starlink] {len(cands)} candidates; selected {name} ({norad}), "
          f"alt {e['alt_km']:.1f} km, incl {e['incl']:.4f} deg")
    write_tle(REF / "tle" / f"starlink_{norad}.txt", name, l1, l2)
    pinned.append((f"tle/starlink_{norad}.txt", name, e, side,
                   "Starlink-class main case. Selected from the group file: inclination "
                   "53.0-53.3 deg, mean-motion altitude 530-560 km, e < 0.001, epoch within "
                   f"2 days of fetch; closest to 550 km ({len(cands)} candidates). "
                   f"Mean-motion altitude {e['alt_km']:.1f} km, inclination {e['incl']:.4f} deg. "
                   "No 53-deg Starlink was above 539.8 km in this fetch (the 53-deg shells "
                   "were lowered, most to ~463 km), so the nominal 540 km shell stands in "
                   "for SPEC's ~550 km. The pre-approved plan used 540-560 km; the widening "
                   "to 530-560 km was decided after this fetch returned no candidates. The "
                   "user did not pre-approve it and accepts the result as a scenario input."))

    # ISS: secondary case.
    body, side = fetch_once("iss", SOURCES["iss"], args.download_dir)
    tles, bad = parse_tles(body)
    iss = [t for t in tles if int(t[1][2:7]) == ISS_NORAD]
    if len(iss) != 1:
        sys.exit(f"[iss] expected exactly one valid TLE for NORAD {ISS_NORAD}, got {len(iss)}")
    name, l1, l2 = iss[0]
    e = elements(l1, l2)
    write_tle(REF / "tle" / "iss.txt", name, l1, l2)
    pinned.append(("tle/iss.txt", name, e, side, "ISS, secondary case."))

    out = ["# Reference sources", "",
           "Pinned inputs for the tests and demos. Protected from tag `golden-v1` onward.",
           "Fetched once with `scripts/fetch_tle.py`; tests never fetch.", "",
           "## TLEs", "",
           "| File | Object | NORAD | TLE epoch (UTC) | Source URL | Fetched (UTC) | Notes |",
           "|---|---|---|---|---|---|---|"]
    for path, name, e, side, note in pinned:
        out.append(f"| `{path}` | {name} | {e['norad']} | "
                   f"{e['epoch'].strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]} | {side['url']} | "
                   f"{side['fetched_utc']} | {note} Download sha256 `{side['sha256'][:16]}...` |")
    out += ["", "Checksums of both TLE lines were validated (mod-10 rule) before pinning.", ""]
    (REF / "SOURCES.md").write_text("\n".join(out), encoding="utf-8", newline="\n")
    print(f"wrote {REF / 'SOURCES.md'}")


if __name__ == "__main__":
    main()
