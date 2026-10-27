"""Does MSVC contract a*b + c to FMA under /arch:AVX2 /fp:precise? (M2 step 5)

Run inside the MSVC x64 environment (needs cl and dumpbin), from the repo root, after building:
    .\\.venv\\Scripts\\python.exe scripts\\check_fma.py  ->  out/fma_contraction_m2.txt

Method: every source compiled with /arch:AVX2 (targets ntn_avx2, ntn_core, ntn_isa_probe) is
recompiled with its exact command from build/compile_commands.json plus /FAs (assembly listing
annotated with source file and line), into a scratch directory; the build tree is not touched.
Each FMA instruction (vfmadd*/vfmsub*/vfnmadd*/vfnmsub*) is attributed to the source file of the
listing's most recent "; File" marker:
  - an xsimd header: an explicit xsimd::fma or an xsimd math routine (e.g. the sincos polynomial)
    -- deliberate FMA, not a contraction;
  - one of our own files: the compiler fused a separate multiply and add -- a contraction.
Per-function counts for the probe (tests/cpp/isa_probe_avx2.cpp) answer the question directly:
plain C++ a*b + c (float, double, loop), xsimd operators a*b + c, and the explicit xsimd::fma
positive control. Consistency check: the listing's FMA count per object must equal the count in
the real build object (dumpbin /disasm); otherwise the script fails.
"""

from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "fma_contraction_m2.txt"
AVX2_TARGETS = ("ntn_avx2", "ntn_core", "ntn_isa_probe")
FMA_RE = re.compile(r"^\s+(vfn?m(?:add|sub)(?:sub)?\d{3}[ps][sd])\b")
FILE_RE = re.compile(r"^; File (.+)$")
PROC_RE = re.compile(r"^(\S+)\s+PROC\b(.*)$")
ENDP_RE = re.compile(r"^\S+\s+ENDP\b")
XSIMD_DIR = "xsimd-src"


def target_of(output: str) -> str:
    m = re.search(r"CMakeFiles[\\/]([^\\/]+)\.dir[\\/]", output)
    return m.group(1) if m else "?"


def run(cmd, cwd=None) -> str:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, errors="replace", shell=isinstance(cmd, str))
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}\n{r.stdout}\n{r.stderr}")
    return r.stdout


def git(*args: str) -> str:
    return run(["git", "-C", str(ROOT), *args]).strip()


def listing_fmas(asm: Path):
    """Yields (function demangled name, source file, mnemonic) for each FMA in the listing."""
    cur_file, cur_func = "?", "?"
    for line in asm.read_text(encoding="utf-8", errors="replace").splitlines():
        if m := FILE_RE.match(line):
            cur_file = m.group(1).strip()
        elif m := PROC_RE.match(line):
            # "name PROC ; demangled name, COMDAT" (C++) or "name PROC ; COMDAT" (extern "C").
            demangled = m.group(2).strip().lstrip(";").split(",")[0].strip()
            cur_func = m.group(1) if demangled in ("", "COMDAT") else demangled
        elif ENDP_RE.match(line):
            cur_func = "?"
        elif m := FMA_RE.match(line):
            yield cur_func, cur_file, m.group(1)


def dumpbin_fma_count(obj: Path) -> int:
    n = 0
    for line in run(["dumpbin", "/nologo", "/disasm", str(obj)]).splitlines():
        parts = line.split()
        if any(re.fullmatch(r"vfn?m(add|sub)(sub)?\d{3}[ps][sd]", p) for p in parts[1:12]):
            n += 1
    return n


def main() -> int:
    db = json.loads((ROOT / "build" / "compile_commands.json").read_text(encoding="utf-8"))
    entries = [e for e in db if target_of(e.get("output", "")) in AVX2_TARGETS]
    if not entries:
        print("no /arch:AVX2 sources in build/compile_commands.json", file=sys.stderr)
        return 2

    lines: list[str] = []
    failures: list[str] = []
    say = lines.append
    dirty = " + uncommitted changes" if git("status", "--porcelain") else ""
    say(f"FMA contraction check, M2 ({datetime.now().astimezone():%Y-%m-%d %H:%M %Z})")
    say(f"git: {git('rev-parse', '--short', 'HEAD')}{dirty}")
    say("Question: does MSVC contract a*b + c to FMA in our own code under /arch:AVX2 /fp:precise?")
    say("Method: recompile each /arch:AVX2 source with its build command + /FAs; attribute every FMA")
    say("to the source file of the preceding '; File' marker (xsimd header = deliberate; our file = contraction).")
    say("")

    per_obj = []
    probe_funcs: dict[str, Counter] = defaultdict(Counter)
    contracted_sites: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        for e in sorted(entries, key=lambda x: x["file"]):
            src = Path(e["file"])
            rel = src.resolve().relative_to(ROOT).as_posix()
            tgt = target_of(e["output"])
            stem = f"{tgt}_{src.stem}"
            obj_arg = f"/Fo{tmpd / (stem + '.obj')}"
            cmd = re.sub(r"/Fo\S+", lambda _: obj_arg, e["command"])
            cmd += f" /FAs /Fa{tmpd / (stem + '.asm')}"
            if "/arch:AVX2" not in shlex.split(cmd, posix=False):
                failures.append(f"{rel}: expected /arch:AVX2 in its compile command")
            run(cmd, cwd=e["directory"])
            counts = Counter()
            for func, f, mnem in listing_fmas(tmpd / (stem + ".asm")):
                origin = "xsimd" if XSIMD_DIR in f else "own"
                counts[origin] += 1
                if tgt == "ntn_isa_probe":
                    probe_funcs[func][origin] += 1
                if origin == "own":
                    contracted_sites.append(f"{rel}: {mnem} in {func} (file {f})")
            real = dumpbin_fma_count((Path(e["directory"]) / e["output"]).resolve())
            total = counts["xsimd"] + counts["own"]
            if real != total:
                failures.append(f"{rel}: listing has {total} FMAs, build object has {real}")
            per_obj.append((tgt, rel, total, counts["xsimd"], counts["own"], real))
            flags = " ".join(t for t in e["command"].split()[1:]
                             if t.startswith(("/O", "-O", "/fp", "/arch", "-std", "/std")))
            say(f"{rel}: flags {flags}")

    say("")
    say("== FMAs per /arch:AVX2 object ==")
    say(f"{'target':14s} {'source':34s} {'total':>6s} {'xsimd':>6s} {'own (contracted)':>17s} {'build obj':>10s}")
    for tgt, rel, total, xs, own, real in per_obj:
        say(f"{tgt:14s} {rel:34s} {total:6d} {xs:6d} {own:17d} {real:10d}")
    say("")
    say("== Probe functions (tests/cpp/isa_probe_avx2.cpp) ==")
    expect = {
        "ntn_isa_probe_contract_f32": "plain C++ float a*b + c",
        "ntn_isa_probe_contract_f64": "plain C++ double a*b + c",
        "ntn_isa_probe_contract_loop": "plain C++ loop acc = a*b + acc",
        "ntn_isa_probe_xsimd_mul_add": "xsimd operators a*b + c",
        "ntn_isa_probe_xsimd_fma": "explicit xsimd::fma (positive control)",
        "ntn_isa_probe_xsimd_sincos": "xsimd::sincos",
    }
    for fn, what in expect.items():
        c = probe_funcs.get(fn, Counter())
        say(f"{fn:30s} {what:40s} FMAs: {c['xsimd'] + c['own']:3d} (xsimd {c['xsimd']}, own {c['own']})")
    if sum(probe_funcs.get("ntn_isa_probe_xsimd_fma", Counter()).values()) == 0:
        failures.append("positive control: explicit xsimd::fma produced no FMA in the listing")
    say("")
    contracting = [fn for fn in ("ntn_isa_probe_contract_f32", "ntn_isa_probe_contract_f64",
                                 "ntn_isa_probe_contract_loop", "ntn_isa_probe_xsimd_mul_add")
                   if sum(probe_funcs.get(fn, Counter()).values())]
    n_own = sum(o[4] for o in per_obj if o[0] != "ntn_isa_probe")
    say("== Conclusion ==")
    say("MSVC contracts a*b + c to FMA under /arch:AVX2 /fp:precise: "
        + ("YES (" + ", ".join(contracting) + ")" if contracting else "NO (all four probe forms stay mul + add)"))
    say(f"FMAs attributed to our own code in ntn_avx2/ntn_core (contractions): {n_own}")
    for site in contracted_sites:
        say("  - " + site)
    say("All FMAs in ntn_avx2/ntn_core come from xsimd (explicit xsimd::fma, xsimd math routines)."
        if n_own == 0 else "Some FMAs come from our own code: see the list above.")
    say("")
    say("RESULT: " + ("PASS (consistency checks)" if not failures else "FAIL"))
    for f in failures:
        say("  - " + f)

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
