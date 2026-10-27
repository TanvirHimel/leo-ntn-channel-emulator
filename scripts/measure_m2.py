"""Re-run the M2 measurements and store their raw outputs in out/ (input to write_results_m2.py).

Run inside the MSVC x64 environment (ctest, dumpbin and cl must be on PATH), from the repo root,
after building build/ (Release) and installing ntn_emu._core:
    .\\.venv\\Scripts\\python.exe scripts\\measure_m2.py [--skip-slow]

Steps (each output overwritten):
  toolchain + pinned dependency tags   -> out/toolchain_m2.json
  ctest (fast) / ctest -L slow          -> out/ctest_fast_m2.txt, out/ctest_slow_m2.txt
  build/ntn_tests.exe, NTN_SLOW=1       -> out/gtest_measure_m2.txt (every [ measure ] line, with its test)
  pytest fast / -m slow (with -s)       -> out/pytest_fast_m2.txt, out/pytest_slow_m2.txt
  scripts/check_isa.py                  -> out/isa_check_m2.txt
  scripts/check_fma.py                  -> out/fma_contraction_m2.txt
  scripts/fd_sweep_m2.py                -> out/fd_design_sweep_m2.txt, docs/figures/m2_fd_design_sweep.png
  scripts/demo_m2.py                    -> out/demo_m2.txt, out/m2_results.json
The benchmarks are separate (scripts/bench_m2.py), so they can run on a quiet, plugged-in machine.
Nothing here is typed in: every number in the M2 section comes from these files.
"""

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out"
PY = sys.executable


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def tool_version(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()[0].strip()
    except (OSError, IndexError):
        return "unknown"


def msvc_version():
    for f in sorted((ROOT / "build" / "CMakeFiles").glob("*/CMakeCXXCompiler.cmake")):
        m = re.search(r'set\(CMAKE_CXX_COMPILER_VERSION "([^"]+)"\)', f.read_text())
        if m:
            return m.group(1)
    return "unknown"


def dep_tags():
    text = (ROOT / "cmake" / "Deps.cmake").read_text(encoding="utf-8")
    tags = {}
    for m in re.finditer(r"FetchContent_Declare\(\s*(\w+)\s+GIT_REPOSITORY\s+(\S+)\s+GIT_TAG\s+(\S+)", text):
        tags[m.group(1)] = {"repository": m.group(2), "tag": m.group(3)}
    return tags


def step(name, cmd, out_file, env=None, cwd=ROOT):
    t = time.perf_counter()
    print(f"[measure_m2] {name} ...", flush=True)
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env, errors="replace")
    text = r.stdout + (("\n[stderr]\n" + r.stderr) if r.stderr.strip() else "")
    if out_file is not None:
        (OUT / out_file).write_text(text, encoding="utf-8")
    dt = time.perf_counter() - t
    print(f"[measure_m2] {name}: exit {r.returncode}, {dt:.1f} s" + (f" -> out/{out_file}" if out_file else ""),
          flush=True)
    return r.returncode, dt


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--skip-slow", action="store_true", help="skip the slow suites (not for the milestone record)")
    args = ap.parse_args()
    for tool in ("ctest", "dumpbin", "cl"):
        if shutil.which(tool) is None:
            sys.exit(f"measure_m2: {tool} not on PATH; run inside the VsDevCmd environment (CLAUDE.md)")
    OUT.mkdir(exist_ok=True)
    t0 = time.perf_counter()
    dirty = bool(git("status", "--porcelain"))
    sys.path.insert(0, str(ROOT / "python"))
    from ntn_emu import _core

    toolchain = {
        "time": datetime.now().astimezone().isoformat(timespec="seconds"),
        "commit": git("rev-parse", "HEAD"), "commit_short": git("rev-parse", "--short", "HEAD"), "dirty": dirty,
        "msvc": msvc_version(), "cmake": tool_version(["cmake", "--version"]),
        "ninja": tool_version(["ninja", "--version"]), "python": platform.python_version(),
        "packages": {p: metadata.version(p) for p in ("pybind11", "scikit_build_core", "numpy", "scipy", "skyfield")},
        "deps": dep_tags(), "core_module": Path(_core.__file__).name, "cpu_has_avx2": bool(_core.cpu_has_avx2()),
    }
    (OUT / "toolchain_m2.json").write_text(json.dumps(toolchain, indent=2), encoding="utf-8")

    env_fast = {k: v for k, v in os.environ.items() if k not in ("NTN_SLOW", "NTN_TEST_HIDE_AVX2")}
    env_slow = dict(env_fast, NTN_SLOW="1")
    results = {}
    results["ctest_fast"] = step("ctest fast", ["ctest", "--test-dir", "build", "--output-on-failure"],
                                 "ctest_fast_m2.txt", env_fast)
    if not args.skip_slow:
        results["ctest_slow"] = step("ctest slow", ["ctest", "--test-dir", "build", "-L", "slow", "--output-on-failure"],
                                     "ctest_slow_m2.txt", env_slow)
    results["gtest"] = step("gtest measure lines", [str(ROOT / "build" / "ntn_tests.exe")], "gtest_measure_m2.txt",
                            env_slow if not args.skip_slow else env_fast)
    results["pytest_fast"] = step("pytest fast", [PY, "-m", "pytest", "python/tests", "-q", "-s"],
                                  "pytest_fast_m2.txt", env_fast)
    if not args.skip_slow:
        results["pytest_slow"] = step("pytest slow", [PY, "-m", "pytest", "python/tests", "-q", "-s", "-m", "slow"],
                                      "pytest_slow_m2.txt", env_fast)
    results["check_isa"] = step("check_isa", [PY, "scripts/check_isa.py"], None, env_fast)
    results["check_fma"] = step("check_fma", [PY, "scripts/check_fma.py"], None, env_fast)
    results["fd_sweep"] = step("design sweep", [PY, "scripts/fd_sweep_m2.py"], None, env_fast)
    results["demo"] = step("demo_m2", [PY, "scripts/demo_m2.py"], "demo_m2.txt", env_fast)

    total = time.perf_counter() - t0
    (OUT / "measure_m2_run.json").write_text(json.dumps(
        {"toolchain": toolchain, "steps": {k: {"exit": v[0], "seconds": v[1]} for k, v in results.items()},
         "slow": not args.skip_slow, "total_s": total}, indent=2), encoding="utf-8")
    failed = [k for k, v in results.items() if v[0] != 0]
    print(f"[measure_m2] done in {total:.1f} s; " + (f"FAILED: {', '.join(failed)}" if failed else "all steps exit 0"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
