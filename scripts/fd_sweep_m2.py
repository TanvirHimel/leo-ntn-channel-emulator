"""M2 fractional-delay design sweep: run the gtest, save its table, plot it.

Run from the repo root after building (build/ntn_tests.exe):
    .\\.venv\\Scripts\\python.exe scripts\\fd_sweep_m2.py
Outputs:
    out/fd_design_sweep_m2.txt        raw [ measure ] lines of FracDelayDesign.Sweep + header
    docs/figures/m2_fd_design_sweep.png   worst-case error vs tap count, 128 / 256 phases
The numbers come only from the test's output; nothing is typed in here.
"""

from __future__ import annotations

import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / "ntn_tests.exe"
OUT_TXT = ROOT / "out" / "fd_design_sweep_m2.txt"
OUT_PNG = ROOT / "docs" / "figures" / "m2_fd_design_sweep.png"

ROW_RE = re.compile(r"^\[ measure  \] \| (\d+) \| (\d+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.e+]+) \| ([-\d.]+) \|$")
TARGET_DB = -60.0  # SPEC §2 target over |f| < 0.4 fs
PRODUCTION = (32, 256)


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], check=True, capture_output=True, text=True).stdout.strip()


def main() -> int:
    run = subprocess.run([str(EXE), "--gtest_filter=FracDelayDesign.Sweep"], capture_output=True, text=True)
    if run.returncode != 0 or "[  PASSED  ] 1 test" not in run.stdout:
        print(run.stdout, run.stderr, sep="\n")
        return 1
    measure = [ln for ln in run.stdout.splitlines() if ln.startswith("[ measure  ]")]
    rows = [m.groups() for ln in measure if (m := ROW_RE.match(ln))]
    if len(rows) != 6:
        print(f"expected 6 design rows, parsed {len(rows)}", file=sys.stderr)
        return 1

    dirty = " + uncommitted changes" if git("status", "--porcelain") else ""
    header = [
        f"M2 fractional-delay design sweep ({datetime.now().astimezone():%Y-%m-%d %H:%M %Z})",
        f"git: {git('rev-parse', '--short', 'HEAD')}{dirty}",
        f"command: build\\ntn_tests.exe --gtest_filter=FracDelayDesign.Sweep",
        "",
    ]
    OUT_TXT.parent.mkdir(exist_ok=True)
    OUT_TXT.write_text("\n".join(header + measure) + "\n", encoding="utf-8")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    series = {128: ("#2a78d6", "o"), 256: ("#eb6834", "s")}  # validated default palette, slots 1-2
    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=150)
    fig.patch.set_facecolor(surface)
    ax.set_facecolor(surface)
    for phases, (color, marker) in series.items():
        pts = sorted((int(r[0]), float(r[4])) for r in rows if int(r[1]) == phases)
        xs, ys = zip(*pts)
        ax.plot(xs, ys, color=color, lw=2, marker=marker, ms=8, mec=surface, mew=2, label=f"{phases} phases", zorder=3)
        ax.annotate(f"{phases} phases", (xs[-1], ys[-1]), xytext=(10, -4 if phases == 256 else 6),
                    textcoords="offset points", color=ink2, fontsize=9, va="center")
    prod = next(r for r in rows if (int(r[0]), int(r[1])) == PRODUCTION)
    ax.annotate(f"production 32 x 256: {float(prod[4]):.1f} dB", (32, float(prod[4])), xytext=(-14, 0),
                textcoords="offset points", ha="right", va="center", color=ink, fontsize=9)
    ax.axhline(TARGET_DB, color=ink2, lw=1, ls=(0, (4, 3)), zorder=2)
    ax.text(24.5, TARGET_DB + 1.5, "target -60 dB (SPEC)", color=ink2, fontsize=9)
    ax.set_xticks([16, 24, 32])
    ax.set_xlim(14.5, 36.5)
    ax.set_xlabel("taps", color=ink2)
    ax.set_ylabel("worst error over 200 delays (dB)", color=ink2)
    ax.set_title("Fractional delay: worst-case error vs design (|f| < 0.4 fs)", color=ink, fontsize=11, loc="left")
    ax.grid(axis="y", color=grid, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=ink2)
    ax.set_ylim(-106, -50)
    ax.legend(frameon=False, loc="lower left", labelcolor=ink2)
    fig.tight_layout()
    fig.savefig(OUT_PNG, facecolor=surface)
    print(f"wrote {OUT_TXT.relative_to(ROOT)} and {OUT_PNG.relative_to(ROOT)}")
    print("\n".join(measure))
    return 0


if __name__ == "__main__":
    sys.exit(main())
