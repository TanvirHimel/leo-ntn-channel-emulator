"""M2 throughput benchmarks: run build/ntn_bench.exe and summarise (SPEC §8 M2, "Bench").

Run from the repo root (Release build with -DNTN_BUILD_BENCH=ON):
    .\\.venv\\Scripts\\python.exe scripts\\bench_m2.py                 # timed run
    .\\.venv\\Scripts\\python.exe scripts\\bench_m2.py --preflight-only  # machine/power state only

Records CPU, clock (before and after), Windows power plan (powercfg) and power mode (the
Windows 11 slider, PowerGetEffectiveOverlayScheme), AC/battery state, commit and dirty state,
cpu_has_avx2(), the default kernel, and the benchmark block sizes. Prints a prominent warning if
the laptop is on battery or not set to best performance.

The benchmarks run single-threaded with --benchmark_repetitions=5 and
--benchmark_enable_random_interleaving=true, so scalar and AVX2 repetitions are interleaved and
thermal drift does not favour one side. Rates use wall-clock time (UseRealTime in
bench/bench_m2.cpp): process CPU time on Windows has 15.625 ms resolution. For each benchmark: median Msps over the 5 repetitions
and spread = (max - min) / median, flagged above 5%. Real-time factor = Msps / fs (fs in Msps).
Throughput is reported, never asserted.

Outputs: out/bench_m2.json (Google Benchmark JSON), out/bench_m2_console.txt,
out/bench_m2.txt (readable report), out/bench_m2_summary.json (input to write_results_m2.py).
"""

from __future__ import annotations

import argparse
import ctypes
import json
import re
import statistics
import subprocess
import sys
import uuid
from ctypes import wintypes
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
OUT = ROOT / "out"
EXE = ROOT / "build" / "ntn_bench.exe"
REPETITIONS = 5
SPREAD_FLAG = 0.05

# Windows 10/11 power-mode overlay schemes (the "Power mode" slider).
OVERLAY_NAMES = {
    "00000000-0000-0000-0000-000000000000": "Balanced",
    "ded574b5-45a0-4f42-8737-46345c09c238": "Best performance",
    "961cc777-2547-4f9d-8174-7d86181b8a7a": "Best power efficiency",
    "3af9b8d9-7c97-431d-ad78-34a8bfea439f": "Better performance",
}
PERFORMANCE_PLANS = ("High performance", "Ultimate Performance")


def ps(cmd: str) -> str:
    r = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", cmd],
                       capture_output=True, text=True)
    return r.stdout.strip()


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True).stdout.strip()


def cpu_info() -> dict:
    raw = ps("Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors,"
             "MaxClockSpeed,L2CacheSize,L3CacheSize | ConvertTo-Json")
    d = json.loads(raw)
    d = d[0] if isinstance(d, list) else d
    return {"name": d["Name"].strip(), "cores": d["NumberOfCores"], "threads": d["NumberOfLogicalProcessors"],
            "max_clock_mhz": d["MaxClockSpeed"], "l2_total_kib": d["L2CacheSize"], "l3_kib": d["L3CacheSize"]}


def clock_now() -> dict:
    """Effective clock = Processor Frequency x % Processor Performance (3 samples, 1 s apart)."""
    raw = ps(r"$c = Get-Counter -Counter '\Processor Information(_Total)\Processor Frequency',"
             r"'\Processor Information(_Total)\% Processor Performance' -SampleInterval 1 -MaxSamples 3;"
             r" $c.CounterSamples | Select-Object Path,CookedValue | ConvertTo-Json")
    try:
        rows = json.loads(raw)
    except json.JSONDecodeError:
        return {"available": False, "note": "performance counters unavailable"}
    freq = [r["CookedValue"] for r in rows if r["Path"].endswith("processor frequency")]
    perf = [r["CookedValue"] for r in rows if r["Path"].endswith("% processor performance")]
    if not freq or not perf:
        return {"available": False, "note": "performance counters unavailable"}
    f, p = statistics.mean(freq), statistics.mean(perf)
    return {"available": True, "nominal_mhz": f, "performance_pct": p, "effective_mhz": f * p / 100.0,
            "time": datetime.now().astimezone().isoformat(timespec="seconds")}


def power_plan() -> dict:
    out = subprocess.run(["powercfg", "/getactivescheme"], capture_output=True, text=True).stdout.strip()
    m = re.search(r"([0-9a-fA-F-]{36})\s+\((.+)\)", out)
    return {"raw": out, "guid": m.group(1) if m else None, "name": m.group(2) if m else None}


def power_mode() -> dict:
    """Effective power-mode overlay (Windows 10 1709+ / 11): the 'Power mode' slider."""
    try:
        g = (ctypes.c_byte * 16)()
        rc = ctypes.WinDLL("powrprof").PowerGetEffectiveOverlayScheme(ctypes.byref(g))
        if rc != 0:
            return {"guid": None, "name": None, "note": f"PowerGetEffectiveOverlayScheme returned {rc}"}
        guid = str(uuid.UUID(bytes_le=bytes(bytearray(g))))
        return {"guid": guid, "name": OVERLAY_NAMES.get(guid, "unknown overlay")}
    except (OSError, AttributeError) as e:
        return {"guid": None, "name": None, "note": f"unavailable: {e}"}


class SYSTEM_POWER_STATUS(ctypes.Structure):
    _fields_ = [("ACLineStatus", wintypes.BYTE), ("BatteryFlag", wintypes.BYTE),
                ("BatteryLifePercent", wintypes.BYTE), ("SystemStatusFlag", wintypes.BYTE),
                ("BatteryLifeTime", wintypes.DWORD), ("BatteryFullLifeTime", wintypes.DWORD)]


def ac_power() -> dict:
    s = SYSTEM_POWER_STATUS()
    if not ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(s)):
        return {"ac_online": None, "note": "GetSystemPowerStatus failed"}
    ac = {0: False, 1: True}.get(s.ACLineStatus & 0xFF)
    return {"ac_online": ac, "battery_percent": None if s.BatteryLifePercent & 0xFF == 255 else s.BatteryLifePercent & 0xFF}


def build_type() -> str:
    cache = ROOT / "build" / "CMakeCache.txt"
    m = re.search(r"^CMAKE_BUILD_TYPE:\w+=(\S*)$", cache.read_text(encoding="utf-8"), re.M) if cache.exists() else None
    return m.group(1) if m else "unknown"


def engine_state() -> dict:
    from ntn_emu import _core  # import-time CPU check runs here
    import numpy as np
    t = _core.RangeTrack(7.68e6, np.array([0], np.int64), np.array([1000], np.int64), np.array([[1e6, 0, 0, 0]]))
    default = _core.Pipeline(t, _core.GeoConfig(fc_hz=2e9, fs_hz=7.68e6)).kernel
    return {"cpu_has_avx2": bool(_core.cpu_has_avx2()), "default_kernel": default.name}


def preflight() -> dict:
    plan, mode, ac = power_plan(), power_mode(), ac_power()
    warnings = []
    if ac.get("ac_online") is not True:
        warnings.append("NOT ON AC POWER (battery or unknown): throughput will be limited and noisy")
    best = mode.get("name") == "Best performance" or plan.get("name") in PERFORMANCE_PLANS
    if not best:
        warnings.append(f"POWER MODE IS NOT 'Best performance' (mode: {mode.get('name')}, plan: {plan.get('name')})")
    bt = build_type()
    if bt != "Release":
        warnings.append(f"BUILD TYPE IS {bt}, NOT Release")
    dirty = bool(git("status", "--porcelain"))
    return {"time": datetime.now().astimezone().isoformat(timespec="seconds"), "cpu": cpu_info(),
            "clock_before": clock_now(), "power_plan": plan, "power_mode": mode, "power": ac,
            "commit": git("rev-parse", "HEAD"), "commit_short": git("rev-parse", "--short", "HEAD"),
            "dirty": dirty, "build_type": bt, "engine": engine_state(), "warnings": warnings}


def warning_box(warnings: list[str]) -> list[str]:
    if not warnings:
        return []
    bar = "#" * 96
    return [bar, "#  WARNING: benchmark conditions are not as required", *[f"#  - {w}" for w in warnings],
            "#  Numbers from this run are not representative; fix the above and rerun.", bar]


def parse_name(name: str) -> dict:
    d = {"bench": name.split("/")[0]}
    for part in name.split("/")[1:]:
        if ":" in part:
            k, v = part.split(":", 1)
            d[k] = int(v)
    return d


def summarise(bench_json: dict, cpu: dict) -> list[dict]:
    runs: dict[str, list[dict]] = {}
    for b in bench_json["benchmarks"]:
        if b.get("run_type") == "iteration":
            runs.setdefault(b["run_name"], []).append(b)
    l2_core = cpu["l2_total_kib"] / cpu["cores"] if cpu.get("cores") else None
    rows = []
    for name, reps in runs.items():
        msps = sorted(r["items_per_second"] / 1e6 for r in reps)
        med = statistics.median(msps)
        spread = (msps[-1] - msps[0]) / med
        p = parse_name(name)
        ws = reps[0].get("ws_KiB")
        fits = ("L2" if l2_core and ws <= l2_core else "L3" if cpu.get("l3_kib") and ws <= cpu["l3_kib"] else "beyond L3")
        row = {"name": name, **p, "kernel_name": reps[0].get("label"), "repetitions": len(msps),
               "msps": msps, "median_msps": med, "spread": spread, "spread_flag": spread > SPREAD_FLAG,
               "ws_kib": ws, "fits": fits, "block": int(reps[0].get("block", p.get("block", 0)))}
        if "fs_Msps" in reps[0]:
            row["fs_msps"] = reps[0]["fs_Msps"]
            row["rtf"] = med / reps[0]["fs_Msps"]
        rows.append(row)
    return sorted(rows, key=lambda r: (r["bench"], r.get("fs_kHz", 0), r["block"], r.get("kernel", 0)))


def report(pre: dict, post_clock: dict, rows: list[dict], flags: list[str], ctx: dict) -> list[str]:
    L = [*warning_box(pre["warnings"])]
    c = pre["cpu"]
    L += [f"M2 benchmarks ({pre['time']})",
          f"commit: {pre['commit_short']}{' + uncommitted changes' if pre['dirty'] else ' (clean)'}",
          f"CPU: {c['name']}, {c['cores']} cores / {c['threads']} threads, max clock {c['max_clock_mhz']} MHz, "
          f"L2 {c['l2_total_kib']} KiB total ({c['l2_total_kib'] / c['cores']:.0f} KiB/core), L3 {c['l3_kib']} KiB",
          f"power plan (powercfg): {pre['power_plan']['raw']}",
          f"power mode (overlay): {pre['power_mode'].get('name')} ({pre['power_mode'].get('guid')})",
          f"AC power: {pre['power'].get('ac_online')} (battery {pre['power'].get('battery_percent')}%)"]
    for label, ck in (("clock before", pre["clock_before"]), ("clock after", post_clock)):
        L.append(f"{label}: " + (f"effective {ck['effective_mhz']:.0f} MHz (nominal {ck['nominal_mhz']:.0f} MHz x "
                                 f"{ck['performance_pct']:.1f}% performance, {ck['time']})"
                                 if ck.get("available") else ck.get("note", "unavailable")))
    e = pre["engine"]
    L += [f"cpu_has_avx2(): {e['cpu_has_avx2']}; default kernel: {e['default_kernel']}; build: {pre['build_type']}",
          f"benchmark context: {json.dumps(ctx, sort_keys=True)}",
          f"flags: {' '.join(flags)}",
          "timing: wall clock (UseRealTime; Windows process CPU time has 15.625 ms resolution)",
          "single-threaded; throughput is reported, not asserted", ""]
    if not rows:
        return L
    L.append("Per stage (Msps = 1e6 complex samples/s; median of repetitions; spread = (max - min) / median)")
    L.append(f"{'benchmark':48s} {'kernel':7s} {'block':>6s} {'ws KiB':>8s} {'fits':>9s} {'median':>9s} "
             f"{'spread':>7s} {'RTF':>7s}")
    for r in rows:
        flag = " (>5%)" if r["spread_flag"] else ""
        rtf = f"{r['rtf']:7.2f}" if "rtf" in r else "      -"
        L.append(f"{r['name']:48s} {r['kernel_name']:7s} {r['block']:6d} {r['ws_kib']:8.0f} {r['fits']:>9s} "
                 f"{r['median_msps']:9.2f} {100 * r['spread']:6.1f}%{rtf}{flag}")
    flagged = [r["name"] for r in rows if r["spread_flag"]]
    L.append("")
    L.append(f"spread > 5%: {', '.join(flagged) if flagged else 'none'}")
    L.append("'full chain' (BM_Pipeline) = geometric stage only in M2 (TDL arrives in M3, AWGN in M4)")
    return L + ["", *warning_box(pre["warnings"])]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preflight-only", action="store_true", help="record machine/power state; run no benchmark")
    ap.add_argument("--min-time", default="1s", help="--benchmark_min_time per repetition (default 1s)")
    args = ap.parse_args()

    pre = preflight()
    for ln in warning_box(pre["warnings"]):
        print(ln)
    if args.preflight_only:
        print(json.dumps(pre, indent=2))
        return 0
    if not EXE.exists():
        print(f"{EXE} not found: configure with -DNTN_BUILD_BENCH=ON (Release) and build", file=sys.stderr)
        return 2

    OUT.mkdir(exist_ok=True)
    flags = [f"--benchmark_repetitions={REPETITIONS}", "--benchmark_enable_random_interleaving=true",
             f"--benchmark_min_time={args.min_time}", f"--benchmark_out={OUT / 'bench_m2.json'}",
             "--benchmark_out_format=json"]
    run = subprocess.run([str(EXE), *flags], capture_output=True, text=True)
    (OUT / "bench_m2_console.txt").write_text(run.stdout + run.stderr, encoding="utf-8")
    if run.returncode != 0:
        print(run.stdout, run.stderr, sep="\n")
        return run.returncode
    post_clock = clock_now()
    bj = json.loads((OUT / "bench_m2.json").read_text(encoding="utf-8"))
    rows = summarise(bj, pre["cpu"])
    ctx = {k: v for k, v in bj["context"].items()
           if k in ("cpu_has_avx2", "default_kernel", "ntn_build_type", "frac_delay_design", "full_chain",
                    "library_build_type", "num_cpus", "mhz_per_cpu", "cpu_scaling_enabled")}
    lines = report(pre, post_clock, rows, flags, ctx)
    (OUT / "bench_m2.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "bench_m2_summary.json").write_text(json.dumps(
        {"preflight": pre, "clock_after": post_clock, "flags": flags, "context": ctx, "rows": rows,
         "spread_flag_threshold": SPREAD_FLAG}, indent=2), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
