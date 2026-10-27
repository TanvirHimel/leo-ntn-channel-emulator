"""Write the M2 section of docs/results.md from the files in out/ (no typed numbers).

    .\\.venv\\Scripts\\python.exe scripts\\write_results_m2.py [--dry-run]

Inputs (all produced by scripts/measure_m2.py and scripts/bench_m2.py):
  out/toolchain_m2.json, out/measure_m2_run.json, out/ctest_{fast,slow}_m2.txt,
  out/pytest_{fast,slow}_m2.txt, out/gtest_measure_m2.txt, out/fd_design_sweep_m2.txt,
  out/isa_check_m2.txt, out/fma_contraction_m2.txt, out/bench_m2_summary.json, out/m2_results.json
Only the text between "## M2 — Accuracy, AVX2, benchmarks" and "## M3" is replaced; the script
checks that everything before the M2 header (Machine, M1a, M1b) and from "## M3" on is unchanged.
The section is also written to out/results_m2_section.md. --dry-run writes only that file.
Measured values are copied or computed from the inputs; the only literals are limits and labels.
"""

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out"
RESULTS_MD = ROOT / "docs" / "results.md"
M2_HEADER = "## M2 — Accuracy, AVX2, benchmarks"
M3_HEADER = "## M3 — "
INPUTS = ["toolchain_m2.json", "measure_m2_run.json", "ctest_fast_m2.txt", "ctest_slow_m2.txt", "pytest_fast_m2.txt",
          "pytest_slow_m2.txt", "gtest_measure_m2.txt", "fd_design_sweep_m2.txt", "isa_check_m2.txt",
          "fma_contraction_m2.txt", "bench_m2_summary.json", "m2_results.json"]


def read(name):
    return (OUT / name).read_text(encoding="utf-8", errors="replace")


def last_match(name, pattern):
    hits = [ln.strip() for ln in read(name).splitlines() if re.search(pattern, ln)]
    return hits[-1] if hits else f"(no match in {name})"


def md_cell(s):
    return str(s).replace("|", "\\|")


def gtest_by_test():
    """{test name: [measure lines]} from the gtest run, plus its summary lines."""
    cur, by, summary = None, {}, []
    for ln in read("gtest_measure_m2.txt").splitlines():
        if ln.startswith("[ RUN      ] "):
            cur = ln[len("[ RUN      ] "):].strip()
        elif ln.startswith("[ measure  ] ") and cur:
            by.setdefault(cur, []).append(ln[len("[ measure  ] "):].rstrip())
        elif ln.startswith(("[  PASSED  ]", "[  FAILED  ]", "[  SKIPPED ]")) and "test" in ln:
            summary.append(ln.strip())
    return by, summary


def code_block(by, prefixes):
    lines = []
    for test, ms in by.items():
        if any(test.startswith(p) for p in prefixes):
            lines += [f"{test}: {m}" for m in ms]
    return ["```", *lines, "```"] if lines else ["_(no measure lines)_"]


def golden_rows():
    rows = []
    rx = re.compile(r"\[ measure \] core vs golden, (\w+), ([\d.]+) s: (-?[\d.]+) dB \(limit (-?[\d.]+) dB\)")
    for f in ("pytest_fast_m2.txt", "pytest_slow_m2.txt"):
        for m in rx.finditer(read(f)):
            rows.append((m.group(1), m.group(2), m.group(3), m.group(4), f))
    return rows


def sweep_rows():
    rx = re.compile(r"\| (\d+) \| (\d+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.e+]+) \| ([-\d.]+) \|$")
    return [m.groups() for ln in read("fd_design_sweep_m2.txt").splitlines() if (m := rx.search(ln))]


def isa_excerpt():
    lines = read("isa_check_m2.txt").splitlines()
    keep, section = [], None
    for ln in lines:
        if ln.startswith("== "):
            section = ln[3:4]
        if ln.startswith(("ISA check", "git:", "scope:", "       ")) or ln.startswith("RESULT"):
            keep.append(ln)
        elif section == "1" and ln.startswith("[build]"):
            keep.append(ln.rstrip())
        elif section in ("2", "4", "5", "6") and (ln.startswith(("==", "build ", "[build]", "(", " ")) or "build dir" in ln):
            keep.append(ln.rstrip())
        elif ln.startswith(("INFO list", "  - [build]")):
            keep.append(ln.rstrip())
    return keep


def fma_excerpt():
    lines = read("fma_contraction_m2.txt").splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("== FMAs per"))
    return lines[:2] + [""] + lines[start:]


def fmt_ms(r):
    flag = " ⚠ spread > 5%" if r["spread_flag"] else ""
    return f"{r['median_msps']:.2f} ({100 * r['spread']:.1f}%{flag})"


def bench_section(b):
    pre, rows = b["preflight"], b["rows"]
    c, e = pre["cpu"], pre["engine"]
    L = []
    if pre["warnings"]:
        L += ["> **WARNING: the benchmark conditions were not as required:** " + "; ".join(pre["warnings"]) +
              ". These throughput numbers are not representative.", ""]
    ck = lambda k: (f"{k['effective_mhz']:.0f} MHz effective ({k['nominal_mhz']:.0f} MHz nominal x "  # noqa: E731
                    f"{k['performance_pct']:.1f}% processor performance, {k['time']})" if k.get("available")
                    else k.get("note", "unavailable"))
    L += ["| Item | Value |", "|---|---|",
          f"| Commit | `{pre['commit_short']}`{' + uncommitted changes' if pre['dirty'] else ' (clean)'} |",
          f"| Time | {pre['time']} |",
          f"| CPU | {c['name']}, {c['cores']} cores / {c['threads']} threads, max clock {c['max_clock_mhz']} MHz |",
          f"| Caches | L2 {c['l2_total_kib']} KiB total ({c['l2_total_kib'] / c['cores']:.0f} KiB per core), L3 {c['l3_kib']} KiB |",
          f"| Clock before / after | {ck(pre['clock_before'])} / {ck(b['clock_after'])} |",
          f"| Power plan (powercfg) | {md_cell(pre['power_plan']['name'])} (`{pre['power_plan']['guid']}`) |",
          f"| Power mode (Windows overlay) | {md_cell(pre['power_mode'].get('name'))} (`{pre['power_mode'].get('guid')}`) |",
          f"| AC power | {pre['power'].get('ac_online')} (battery {pre['power'].get('battery_percent')}%) |",
          f"| Build | {pre['build_type']}; Google Benchmark library build: {b['context'].get('library_build_type')} |",
          f"| cpu_has_avx2() / default kernel | {e['cpu_has_avx2']} / {e['default_kernel']} |",
          f"| Flags | `{' '.join(f for f in b['flags'] if not f.startswith('--benchmark_out'))}` |",
          "", "Msps = 1e6 complex samples per second; median of the repetitions, spread = (max - min) / median "
          "in parentheses (⚠ marks spread > 5%). Single-threaded. Throughput is reported, not pass/fail.", ""]

    def pick(bench, kernel, block, fs=None):
        for r in rows:
            if r["bench"] == bench and r.get("kernel") == kernel and r["block"] == block and (fs is None or r.get("fs_kHz") == fs):
                return r
        return None

    L += ["#### Per stage, block 4096 samples", "",
          "| Stage | Scalar Msps | AVX2 Msps | AVX2 / scalar | Working set (est.), scalar / AVX2 |", "|---|---|---|---|---|"]
    for bench, label in (("BM_FracDelay", "Fractional delay (kernel)"), ("BM_Rotator", "Rotator (kernel)"),
                         ("BM_GeoChannel", "Geometric channel (fs 7.68 Msps)")):
        s, a = pick(bench, 0, 4096), pick(bench, 1, 4096)
        if s and a:
            L.append(f"| {label} | {fmt_ms(s)} | {fmt_ms(a)} | {a['median_msps'] / s['median_msps']:.2f}x | "
                     f"{s['ws_kib']:.0f} KiB ({s['fits']}) / {a['ws_kib']:.0f} KiB ({a['fits']}) |")
    L += ["", "#### Full chain (M2: geometric stage only; TDL arrives in M3, AWGN in M4)", "",
          "| fs (Msps) | Block | Scalar Msps | Scalar RTF | AVX2 Msps | AVX2 RTF | AVX2 / scalar | Working set (est.), AVX2 |",
          "|---|---|---|---|---|---|---|---|"]
    for fs in sorted({r["fs_kHz"] for r in rows if r["bench"] == "BM_Pipeline"}):
        for block in sorted({r["block"] for r in rows if r["bench"] == "BM_Pipeline"}):
            s, a = pick("BM_Pipeline", 0, block, fs), pick("BM_Pipeline", 1, block, fs)
            if s and a:
                L.append(f"| {s['fs_msps']:.2f} | {block} | {fmt_ms(s)} | {s['rtf']:.2f} | {fmt_ms(a)} | {a['rtf']:.2f} | "
                         f"{a['median_msps'] / s['median_msps']:.2f}x | {a['ws_kib']:.0f} KiB ({a['fits']}) |")
    flagged = [r["name"] for r in rows if r["spread_flag"]]
    L += ["", "Notes:",
          "- Real-time factor (RTF) = Msps / fs; below 1 means slower than real time.",
          "- \"Full chain\" is `Pipeline::process`, which in M2 contains only the geometric stage (range evaluation, "
          "fractional delay, carrier rotation, input history). TDL fading (M3) and AWGN (M4) will lower these numbers.",
          "- Each benchmark processes the same block repeatedly, so the data stay in cache as far as the working set "
          "allows. The working set is an estimate (input, output and history buffers, per-4096-sample scratch, filter "
          "table), compared with the per-core L2 and the shared L3; the 65536-sample block shows the effect of a "
          "working set beyond L2.",
          "- Scalar and AVX2 repetitions were interleaved in random order (`--benchmark_enable_random_interleaving`).",
          "- Rates use wall-clock time (`UseRealTime`). A first run used Google Benchmark's default basis, process CPU "
          "time, which on Windows advances in 15.625 ms ticks; its medians and spreads were snapped to that grid "
          "(kept as `out/bench_m2_cputime_quantized.*`, not used here).",
          f"- Spread > 5%: {', '.join(flagged) if flagged else 'none'}.",
          "- The scalar kernels are a separate target (`ntn_scalar`) compiled without `/arch:AVX2` (verified by "
          "`scripts/check_isa.py`), so the comparison is not inflated by compiler auto-vectorisation of the scalar path."]
    return L


def build_section():
    tc = json.loads(read("toolchain_m2.json"))
    run = json.loads(read("measure_m2_run.json"))
    bench = json.loads(read("bench_m2_summary.json"))
    demo = json.loads(read("m2_results.json"))
    by, gsummary = gtest_by_test()
    S = [M2_HEADER,
         f"Code commit: `{tc['commit_short']}`{' + uncommitted changes' if tc['dirty'] else ' (clean)'} "
         f"(full hash `{tc['commit']}`), measurements started {tc['time']}; benchmarks at commit "
         f"`{bench['preflight']['commit_short']}`{' + uncommitted changes' if bench['preflight']['dirty'] else ' (clean)'}.",
         "",
         "Commands (repo root, PowerShell; build and C++ commands inside the VsDevCmd environment, see CLAUDE.md):",
         "- `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DNTN_BUILD_BENCH=ON`, `cmake --build build -j 2`",
         "- `.\\.venv\\Scripts\\python.exe scripts\\bench_m2.py` -> `out/bench_m2.txt`, `out/bench_m2.json`, `out/bench_m2_summary.json`",
         f"- `.\\.venv\\Scripts\\python.exe scripts\\measure_m2.py` -> suites, gtest measure lines, ISA/FMA checks, design "
         f"sweep, demo (`out/*_m2*`; {run['total_s']:.0f} s)",
         "- `.\\.venv\\Scripts\\python.exe scripts\\write_results_m2.py` -> this section", "",
         "Notes:",
         "- Scalar kernels (`frac_delay`, `rotator`) are a separate CMake target (`ntn_scalar`) compiled without "
         "`/arch:AVX2`; the AVX2 kernels (xsimd) are in `ntn_avx2`, compiled with `/arch:AVX2`. Proof: compile flags and "
         "disassembly below (`out/isa_check_m2.txt`).",
         "- Departure from SPEC §6 (recorded in SPEC.md): the bindings (`_core`) and `ntn_scalar` (which holds the CPU "
         "check) are built without `/arch:AVX2`, and the `.pyd` without LTCG, so the import-time check runs before any "
         "AVX2 instruction.",
         "- CPU requirement: AVX2, FMA3, BMI1, BMI2 and OS support for AVX state (XGETBV). `ntn_core` is compiled with "
         "`/arch:AVX2` (MSVC emits BMI1/BMI2 there; xsimd uses FMA3), so there is no runtime fallback: `import "
         "ntn_emu._core` raises ImportError and `GeometricChannel`/`Pipeline` throw for either kernel when a feature is "
         "missing. The scalar kernels are a test oracle and benchmark baseline.",
         "- Fractional delay: 32 taps x 256 phases kept (decision from the design sweep below); Kaiser beta from the "
         "window design formula (not tuned).",
         "- Derived limit (new check, not a loosened tolerance): the vectorised float32 sincos alone must stay below "
         "1e-6 per sample (|(c + js) - e^{ja}| and phase), half of the 2e-6 rad tone-phase budget of M1a test 1 on the "
         "C++ path; the float32 angle rounding and the complex multiply take the rest.",
         "- \"Full chain\" in the benchmarks is the geometric stage only (TDL arrives in M3, AWGN in M4).",
         "- Throughput is not pass/fail (SPEC §3).", ""]

    S += ["### Toolchain and dependencies", "| Item | Value |", "|---|---|",
          f"| MSVC (cl) | {tc['msvc']} |", f"| CMake | {tc['cmake']} |", f"| Ninja | {tc['ninja']} |",
          f"| Python | {tc['python']} |"]
    S += [f"| {k} | {v} |" for k, v in tc["packages"].items()]
    S += [f"| {name} (FetchContent) | `{d['tag']}` ({d['repository']}) |" for name, d in tc["deps"].items()]
    S += [f"| ntn_emu._core | {tc['core_module']}; cpu_has_avx2() = {tc['cpu_has_avx2']} |", ""]

    S += ["### Test suites (from out/ logs)", "| Suite | Result |", "|---|---|",
          f"| pytest fast | `{last_match('pytest_fast_m2.txt', r'passed|failed')}` |",
          f"| pytest slow | `{last_match('pytest_slow_m2.txt', r'passed|failed')}` |",
          f"| ctest fast | `{last_match('ctest_fast_m2.txt', r'tests passed|tests failed')}`, "
          f"`{last_match('ctest_fast_m2.txt', r'Total Test time')}` |",
          f"| ctest slow (NTN_SLOW=1, -L slow) | `{last_match('ctest_slow_m2.txt', r'tests passed|tests failed')}`, "
          f"`{last_match('ctest_slow_m2.txt', r'Total Test time')}` |",
          f"| gtest, all tests with NTN_SLOW=1 (source of the measure lines below) | `{'; '.join(gsummary)}` |", ""]

    sw = sweep_rows()
    S += ["### Fractional-delay design sweep (`out/fd_design_sweep_m2.txt`)",
          "200 delays (100 grid fractions, 100 random in [-8, 8)), band-limited periodic noise |f| <= 0.4 fs. "
          "Worst/median: error against the FFT-domain reference. Peak: max |H(f) - exp(-j2 pi f d)| over |f| <= 0.4 fs "
          "and all delays, measured through the kernel.", "",
          "| Taps | Phases | Kaiser beta | Design A (dB) | Worst (dB) | Median (dB) | Peak \\|H - ideal\\| | Peak (dB) |",
          "|---|---|---|---|---|---|---|---|"]
    for r in sw:
        prod = " (production)" if (r[0], r[1]) == ("32", "256") else ""
        S.append(f"| {r[0]}{prod} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]} | {r[6]} | {r[7]} |")
    S += ["", "Figure: `docs/figures/m2_fd_design_sweep.png` (worst-case error vs taps, 128 and 256 phases, -60 dB target).",
          "", "![Design sweep](figures/m2_fd_design_sweep.png)", ""]

    S += ["### Acceptance tests 1-2: fractional-delay accuracy and slip continuity (gtest, scalar kernel)",
          "Limits: -60 dB for each of 50 random delays; -60 dB overall and within +-64 samples of every slip. "
          "The negative control must exceed -60 dB.", *code_block(by, ("FracDelayM2.", "SlipContinuity.")), ""]

    S += ["### Acceptance test 3: AVX2 vs scalar (limit -100 dB)",
          "Lengths 0, 1, 7, 8, 9, 31, 33, 4095, 4097; buffers aligned and misaligned by 8 bytes; integer and fractional "
          "delays; NaN sentinels next to the allowed read range; output guard cells.",
          *code_block(by, ("FracDelayAvx2.", "RotatorAvx2.IsNot", "RotatorAvx2.Matches")), ""]

    S += ["### Vectorised sincos accuracy (derived limit 1e-6, see notes)",
          "Reference: float64 sin/cos of the same float32 argument; range [-pi_f, pi_f] (the rotator's wrapped "
          "argument). Fast: 2^22 + 1 grid points plus special values; slow: every float32 in the range.",
          *code_block(by, ("SincosAvx2.",)), ""]

    S += ["### M1a closed-form tests on both kernels (gtest, channel level; same limits as M1b)",
          "Phase < 2e-6 rad with no drift (last 1e5 samples), FFT peak within 0.1 Hz, integer shift bit-exact, "
          "quadratic phase within 0.01 Hz, time-varying delay < -90 dB. Kernel in the test name (`/Scalar`, `/Avx2`); "
          "`RotatorAvx2.TonePhase*` is the same tone test on the AVX2 rotator alone.",
          *code_block(by, ("Kernels/", "RotatorAvx2.TonePhase")), ""]

    S += ["### Acceptance test 4: AVX2 vs golden, block invariance",
          "Golden: float64 model on the M1b band-limited OFDM stimulus, real pass, window crossing a segment knot.", "",
          "| Kernel | Window | Error power | Limit | Log |", "|---|---|---|---|---|"]
    S += [f"| {k} | {w} s | {e} dB | {lim} dB | `{f}` |" for k, w, e, lim, f in golden_rows()]
    S += ["", "Block sizes 1, 7 and 4096 (and a mixed pattern); limits: scalar bit-identical, AVX2 < -120 dB.",
          *code_block(by, ("PipelineBlocks.",)), ""]

    S += ["### CPU feature check", *code_block(by, ("CpuCheck.",)),
          "Import-time check: `python/tests/test_cpu_guard.py` (ImportError naming each hidden feature, in a fresh "
          "interpreter; part of the fast pytest suite above).", ""]

    S += ["### Instruction-set checks (`out/isa_check_m2.txt`, `out/fma_contraction_m2.txt`)",
          "ISA check (object files of our sources; the full report also covers the `.pyd` build tree):",
          "```", *isa_excerpt(), "```", "",
          "FMA contraction under `/arch:AVX2 /fp:precise` (each FMA attributed via a `/FAs` listing to an xsimd header "
          "or to our own code):", "```", *fma_excerpt(), "```", ""]

    S += ["### Benchmarks (`out/bench_m2.txt`, `out/bench_m2_summary.json`)", *bench_section(bench), ""]

    S += ["### End-to-end demo (`out/m2_results.json`, `out/demo_m2.txt`)",
          f"Window: pass '{demo['pass_label']}', {demo['window_start_utc']} + {demo['window_s']} s (same as M1a/M1b); "
          f"kernel {demo['kernel']}; cpu_has_avx2() = {demo['cpu_has_avx2']}; latency L = {demo['latency_samples']}, "
          f"tau0 = {demo['tau0_samples']} samples.", "",
          "| Sub | Centre (s) | rho_dot (m/s) | Measured tone offset (Hz) | Predicted -(fc+f_tone)rho_dot/c (Hz) | Error (Hz) |",
          "|---|---|---|---|---|---|"]
    S += [f"| {r['sub']} | {r['centre_s']:.3f} | {r['rho_dot_mps']:.4f} | {r['measured_hz']:.4f} | "
          f"{r['predicted_hz']:.4f} | {r['error_hz']:+.4f} |" for r in demo["tone_rows"]]
    S += ["", f"Result: max \\|measured - predicted\\| = {demo['doppler_max_err_hz']:.4f} Hz over "
          f"{len(demo['tone_rows'])} sub-windows (limit {demo['doppler_limit_hz']:g} Hz) -> "
          f"{'PASS' if demo['doppler_max_err_hz'] < demo['doppler_limit_hz'] else 'FAIL'}.",
          f"AVX2 vs scalar on the same 1 s input: {demo['avx2_vs_scalar_db']:.1f} dB, "
          f"{demo['avx2_vs_scalar_differing_samples']} of {demo['n_samples']} samples differ in bits. "
          f"Wall time AVX2 / scalar (informational, not a benchmark): {demo['wall_s']['avx2']:.2f} s / "
          f"{demo['wall_s']['scalar']:.2f} s.", ""]
    return "\n".join(S).rstrip("\n") + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true", help="write out/results_m2_section.md only")
    args = ap.parse_args()
    missing = [f for f in INPUTS if not (OUT / f).exists()]
    if missing:
        sys.exit("write_results_m2: missing inputs in out/: " + ", ".join(missing))
    section = build_section()
    (OUT / "results_m2_section.md").write_text(section, encoding="utf-8")
    if args.dry_run:
        print(f"wrote out/results_m2_section.md ({len(section.splitlines())} lines)")
        return
    text = RESULTS_MD.read_text(encoding="utf-8")
    a, b = text.index(M2_HEADER), text.index("\n" + M3_HEADER)
    new = text[:a] + section + text[b:]
    assert new[:a] == text[:a] and new.endswith(text[b:]), "only the M2 section may change"
    RESULTS_MD.write_text(new, encoding="utf-8", newline="\n")
    print(f"docs/results.md: M2 section written ({len(section.splitlines())} lines); M1a/M1b unchanged")


if __name__ == "__main__":
    main()
