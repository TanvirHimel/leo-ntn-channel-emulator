# Reference sources

Pinned inputs for the tests and demos. Protected from tag `golden-v1` onward.
Fetched once with `scripts/fetch_tle.py`; tests never fetch.

## TLEs

| File | Object | NORAD | TLE epoch (UTC) | Source URL | Fetched (UTC) | Notes |
|---|---|---|---|---|---|---|
| `tle/starlink_53934.txt` | STARLINK-4775 | 53934 | 2026-10-20 04:06:17.748 | https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle | 2026-10-21T00:00:39Z | Starlink-class main case. Selected from the group file: inclination 53.0-53.3 deg, mean-motion altitude 530-560 km, e < 0.001, epoch within 2 days of fetch; closest to 550 km (31 candidates). Mean-motion altitude 539.8 km, inclination 53.1597 deg. No 53-deg Starlink was above 539.8 km in this fetch (the 53-deg shells were lowered, most to ~463 km), so the nominal 540 km shell stands in for SPEC's ~550 km. The pre-approved plan used 540-560 km; the widening to 530-560 km was decided after this fetch returned no candidates. The user did not pre-approve it and accepts the result as a scenario input. Download sha256 `ccf3aaa4f3ccfe36...` |
| `tle/iss.txt` | ISS (ZARYA) | 25544 | 2026-10-20 12:45:45.460 | https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=tle | 2026-10-21T00:28:25Z | ISS, secondary case. Download sha256 `0cc254aa792e9342...` |

Checksums of both TLE lines were validated (mod-10 rule) before pinning.
