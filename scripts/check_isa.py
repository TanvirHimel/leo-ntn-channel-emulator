"""Instruction-set check for the scalar / AVX2 target split (SPEC §6, M2 step 1).

Run inside the MSVC x64 environment (needs dumpbin), from the repo root:
    .\\.venv\\Scripts\\python.exe scripts\\check_isa.py  ->  out/isa_check_m2.txt

Checks, each PASS/FAIL:
 1. Compile flags (compile_commands.json): /arch:AVX2 exactly on ntn_avx2, ntn_core and
    ntn_isa_probe; absent from ntn_scalar and the bindings (_core).
 2. Objects of the non-AVX2 targets contain no VEX/EVEX-encoded instruction and no ymm/zmm
    register (dumpbin /disasm). /arch:AVX2 VEX-encodes even scalar float code, so "no VEX" is
    stricter than "no ymm".
 3. Positive control: the AVX2 probe object does contain ymm instructions, so check 2 would
    have seen them.
 4. COMDAT hazard: an inline/template function defined (External, COMDAT) in both a non-AVX2
    object and an AVX2 object is folded by the linker to one copy, chosen arbitrarily. If the
    AVX2 copy contains VEX code, the scalar path (or the code that runs before the import-time
    AVX2 check) may execute it. Every such function shared with ntn_scalar or the bindings and
    referenced by a relocation in the link unit must have a VEX-free AVX2 copy.
 5. No dynamic initialisers (.CRT$XCU) in AVX2 objects, so loading the .pyd runs no AVX2 code
    before the import-time check.
 6. Every instruction family beyond AVX/AVX2 in the AVX2 objects (FMA3, BMI1/2, LZCNT, ...) is
    one that cpu_check tests; the AVX2 kernel objects contain ymm code.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "isa_check_m2.txt"

AVX2_TARGETS = {"ntn_avx2", "ntn_core", "ntn_isa_probe"}
NO_AVX2_TARGETS = {"ntn_scalar", "_core"}  # _core: only bindings.cpp is compiled in it
PROBE_TARGET = "ntn_isa_probe"
# Implemented AVX2 kernels: their objects must contain ymm code (really vectorised).
MUST_HAVE_YMM = {"src/ntn/frac_delay_avx2.cpp", "src/ntn/rotator_avx2.cpp"}
# Instruction families beyond AVX/AVX2 that /arch:AVX2 code may contain. Every family that
# appears in an AVX2 object must be one cpu_check tests (src/ntn/cpu_check.hpp); the VEX/AVX/AVX2
# encodings themselves are covered by its AVX2 and OS-AVX-state checks.
CPU_CHECK_TESTED = {"FMA3", "BMI1", "BMI2"}
EXT_FAMILIES = {
    "FMA3": re.compile(r"^vfn?m(add|sub|addsub|subadd)"),
    "BMI1": re.compile(r"^(andn|bextr|blsi|blsmsk|blsr|tzcnt)$"),
    "BMI2": re.compile(r"^(bzhi|mulx|pdep|pext|rorx|sarx|shlx|shrx)$"),
    "LZCNT": re.compile(r"^lzcnt$"),
    "POPCNT": re.compile(r"^popcnt$"),
    "MOVBE": re.compile(r"^movbe$"),
    "F16C": re.compile(r"^vcvt(ph2ps|ps2ph)$"),
    "ADX": re.compile(r"^(adcx|adox)$"),
    "AES": re.compile(r"^v?aes"),
    "PCLMUL": re.compile(r"^v?pclmul"),
    "SHA": re.compile(r"^sha"),
    "RDRAND/RDSEED": re.compile(r"^rd(rand|seed)$"),
    "AVX-VNNI": re.compile(r"^vpdp(bus|ws)d"),
}

# VMX instructions start with "v" but are not VEX-encoded.
VMX = {"vmcall", "vmclear", "vmfunc", "vmlaunch", "vmptrld", "vmptrst", "vmread", "vmresume",
       "vmwrite", "vmxoff", "vmxon", "verr", "verw"}

# dumpbin /disasm instruction line: "  0000000000000010: C5 FC 10 04 0A     vmovups  ymm0,..."
INSN_RE = re.compile(r"^\s+[0-9A-F]{8,16}:\s+(?:[0-9A-F]{2}\s)+\s*([a-z][a-z0-9]*)\s*(.*)$")
YMM_RE = re.compile(r"\bymm\d+\b")
ZMM_RE = re.compile(r"\bzmm\d+\b")
LABEL_RE = re.compile(r"^(\S.*):$")
# dumpbin /symbols line: "008 00000000 SECT3  notype ()    External     | ?f@@YAXXZ (void __cdecl f(void))"
SYM_RE = re.compile(r"^[0-9A-F]{3,}\s+[0-9A-F]{8}\s+(SECT[0-9A-F]+|UNDEF|ABS|DEBUG)\s+\w+\s+(\(\))?\s+(\w+)\s+\|\s+(\S+)")


def target_of(output: str) -> str:
    m = re.search(r"CMakeFiles[\\/]([^\\/]+)\.dir[\\/]", output)
    return m.group(1) if m else "?"


def short_flags(command: str) -> str:
    keep = []
    for tok in command.split()[1:]:
        if tok.startswith(("-I", "/I", "-external:I", "/external:I", "/Fo", "/Fd", "/FS", "-D", "/D")):
            continue
        if tok.endswith((".cpp", ".cc")) or tok == "-c" or tok == "/c":
            continue
        keep.append(tok)
    return " ".join(keep)


def run(cmd: list[str]) -> str:
    return subprocess.run(cmd, check=True, capture_output=True, text=True, errors="replace").stdout


def disasm_by_function(obj: Path) -> dict[str, list[tuple[str, str]]]:
    funcs: dict[str, list[tuple[str, str]]] = defaultdict(list)
    cur = "<none>"
    for line in run(["dumpbin", "/nologo", "/disasm",str(obj)]).splitlines():
        m = LABEL_RE.match(line)
        if m and not line.startswith(("Dump of", "File Type")):
            cur = m.group(1).strip()
            continue
        m = INSN_RE.match(line)
        if m:
            funcs[cur].append((m.group(1), m.group(2)))
    return funcs


def is_vex(mnem: str) -> bool:
    return mnem.startswith("v") and mnem not in VMX


def isa_counts(insns: list[tuple[str, str]]) -> dict[str, int]:
    c = {"insns": len(insns), "vex": 0, "ymm": 0, "zmm": 0}
    for mnem, ops in insns:
        c["vex"] += is_vex(mnem)
        c["ymm"] += bool(YMM_RE.search(ops))
        c["zmm"] += bool(ZMM_RE.search(ops))
    return c


def external_functions(obj: Path) -> set[str]:
    out = set()
    for line in run(["dumpbin", "/nologo", "/symbols", str(obj)]).splitlines():
        m = SYM_RE.match(line.strip())
        if m and m.group(1).startswith("SECT") and m.group(2) and m.group(3) == "External":
            out.add(m.group(4))
    return out


def referenced_symbols(obj: Path) -> set[str]:
    """Symbols that any relocation in obj points at (calls, jumps, address-of)."""
    out = set()
    for line in run(["dumpbin", "/nologo", "/relocations", str(obj)]).splitlines():
        parts = line.split()
        if len(parts) >= 5 and re.fullmatch(r"[0-9A-F]{8}", parts[0]) and parts[1].startswith(("REL", "ADDR", "SECREL", "SECTION")):
            out.add(parts[4])
    return out


def dynamic_initialisers(obj: Path) -> int:
    text = run(["dumpbin", "/nologo", "/headers", str(obj)])
    return len(re.findall(r"^\s*\.CRT\$XC[A-Z]", text, flags=re.M))


def main() -> int:
    if shutil.which("dumpbin") is None:
        print("check_isa: dumpbin not on PATH; run inside the VsDevCmd environment", file=sys.stderr)
        return 2

    dbs = [ROOT / "build" / "compile_commands.json"] + sorted((ROOT / "build").glob("py-*/compile_commands.json"))
    entries = []  # (target, source, obj, flags, has_avx2, build_dir)
    for db in dbs:
        if not db.exists():
            continue
        for e in json.loads(db.read_text(encoding="utf-8")):
            src = Path(e["file"])
            try:
                rel = src.resolve().relative_to(ROOT)
            except ValueError:
                continue  # dependency sources
            if rel.parts[0] == "build":
                continue
            out = e.get("output") or re.search(r"/Fo(\S+)", e["command"]).group(1)
            obj = (Path(e["directory"]) / out).resolve()
            tgt = target_of(out)
            entries.append((tgt, rel.as_posix(), obj, short_flags(e["command"]),
                            "/arch:AVX2" in e["command"].split(), db.parent.name))

    lines: list[str] = []
    failures: list[str] = []
    info: list[str] = []  # unreferenced (or test-only) shared functions with VEX code
    say = lines.append
    say(f"ISA check, M2 ({datetime.now().astimezone():%Y-%m-%d %H:%M %Z})")
    say(f"git: {run(['git', '-C', str(ROOT), 'rev-parse', '--short', 'HEAD']).strip()}"
        f"{' + uncommitted changes' if run(['git', '-C', str(ROOT), 'status', '--porcelain']).strip() else ''}")
    say(f"compile databases: {', '.join(str(d.relative_to(ROOT)) for d in dbs if d.exists())}")
    say("scope: inspects the object files (.obj) of our sources, NOT the final linked _core.pyd or")
    say("       ntn_tests.exe; the link step is covered by check 4 (functions shared across the boundary).")
    say("")

    # 1. Flags.
    say("== 1. Compile flags (our sources; -I/-D//Fo omitted) ==")
    for tgt, src, obj, flags, avx2, bdir in sorted(entries, key=lambda t: (t[5], t[0], t[1])):
        expect = True if tgt in AVX2_TARGETS else False if tgt in NO_AVX2_TARGETS else None
        status = "" if expect is None else ("ok" if avx2 == expect else "FAIL")
        if status == "FAIL":
            failures.append(f"flags: {src} ({tgt}) /arch:AVX2={avx2}, expected {expect}")
        say(f"[{bdir}] {tgt:14s} {src:34s} /arch:AVX2={'yes' if avx2 else 'no ':3s} {status}")
        say(f"    {flags}")
    for t in AVX2_TARGETS | NO_AVX2_TARGETS:
        if not any(e[0] == t for e in entries):
            failures.append(f"flags: target {t} not found in any compile database")
    say("")

    # 2/3. Disassembly per object.
    say("== 2-3. Disassembly per object (dumpbin /disasm) ==")
    say(f"{'build dir':26s} {'target':14s} {'object':34s} {'insns':>7s} {'VEX':>6s} {'ymm':>6s} {'zmm':>4s}")
    disasm: dict[Path, dict[str, list[tuple[str, str]]]] = {}
    for tgt, src, obj, flags, avx2, bdir in sorted(entries, key=lambda t: (t[5], t[0], t[1])):
        if not obj.exists():
            failures.append(f"missing object {obj}")
            continue
        disasm[obj] = disasm_by_function(obj)
        c = isa_counts([i for f in disasm[obj].values() for i in f])
        say(f"{bdir:26s} {tgt:14s} {Path(src).name + '.obj':34s} {c['insns']:7d} {c['vex']:6d} {c['ymm']:6d} {c['zmm']:4d}")
        if not avx2 and (c["vex"] or c["ymm"] or c["zmm"]):
            failures.append(f"VEX/ymm in non-AVX2 object {src} ({tgt}): {c}")
        if c["zmm"]:
            failures.append(f"zmm (AVX-512) in {src}")
        if not avx2 and c["insns"] == 0:
            failures.append(f"no machine code in {src} ({tgt}); /GL object? cannot inspect")
        if src in MUST_HAVE_YMM and c["ymm"] == 0:
            failures.append(f"AVX2 kernel {src} has no ymm instructions (not vectorised?)")
        if tgt == PROBE_TARGET and c["ymm"] == 0:
            failures.append("positive control: probe object has no ymm instructions")
    say("")

    # 4. COMDAT hazards, per link unit (each build directory links its own binary).
    say("== 4. Functions defined in both a non-AVX2 and an AVX2 object of one link unit ==")
    say("(the linker keeps one copy, chosen arbitrarily. FAIL if the AVX2 copy has VEX code, the")
    say(" non-AVX2 side is a scalar kernel or the bindings, and some object references the function;")
    say(" an unreferenced copy can never execute. Test/hello objects are reported as INFO.)")
    linked_avx2 = [e for e in entries if e[4] and e[0] != PROBE_TARGET and e[2] in disasm]
    for bdir in sorted({e[5] for e in entries}):
        plain_syms: dict[str, set[str]] = defaultdict(set)
        for e in entries:
            if e[5] == bdir and not e[4] and e[2] in disasm:
                for sym in external_functions(e[2]):
                    plain_syms[sym].add(f"{e[0]}:{Path(e[1]).name}")
        shared: dict[str, int] = {}  # symbol -> max VEX+ymm count over AVX2 copies
        where: dict[str, set[str]] = defaultdict(set)
        for e in linked_avx2:
            if e[5] != bdir:
                continue
            for sym in external_functions(e[2]):
                if sym in plain_syms:
                    c = isa_counts(disasm[e[2]].get(sym, []))
                    shared[sym] = max(shared.get(sym, 0), c["vex"] + c["ymm"])
                    if c["vex"] + c["ymm"]:
                        where[sym].add(f"{e[0]}:{Path(e[1]).name}")
        refs: set[str] = set()
        for e in entries:
            if e[5] == bdir and e[0] != PROBE_TARGET and e[2] in disasm:
                refs |= referenced_symbols(e[2])
        say(f"[{bdir}] relocation targets found: {len(refs)}")
        if not refs:
            failures.append(f"[{bdir}] no relocations parsed; the reference check would be vacuous")
        guarded = {sym for sym in shared
                   if any(o.split(":")[0] in NO_AVX2_TARGETS for o in plain_syms[sym])}
        n_guard_hz = 0
        for sym in sorted(shared):
            if not shared[sym]:
                continue
            referenced = sym in refs
            critical = sym in guarded and referenced
            n_guard_hz += critical
            tag = "FAIL" if critical else "INFO"
            say(f"[{bdir}] {tag} {sym}: VEX code in {', '.join(sorted(where[sym]))}; "
                f"also defined in {', '.join(sorted(plain_syms[sym]))}; "
                f"{'referenced' if referenced else 'referenced by no object in this link unit (dead copy)'}")
            if critical:
                failures.append(f"COMDAT hazard [{bdir}]: {sym}")
            else:
                info.append(f"[{bdir}] {sym} (VEX copy in {', '.join(sorted(where[sym]))}; "
                            f"{'referenced, test/hello only' if referenced else 'unreferenced'})")
        say(f"[{bdir}] shared functions: {len(shared)} (with ntn_scalar/_core: {len(guarded)}); "
            f"with VEX code in the AVX2 copy: {sum(1 for v in shared.values() if v)} "
            f"(referenced and shared with ntn_scalar/_core: {n_guard_hz})")
    say("")

    # 6 (printed after 5). Extension instruction families in AVX2 objects.
    ext_lines = []
    for e in sorted(linked_avx2 + [x for x in entries if x[0] == PROBE_TARGET and x[2] in disasm],
                    key=lambda t: (t[5], t[0], t[1])):
        counts = defaultdict(int)
        for insns in disasm[e[2]].values():
            for mnem, _ in insns:
                for fam, rx in EXT_FAMILIES.items():
                    if rx.match(mnem):
                        counts[fam] += 1
        untested = sorted(f for f in counts if f not in CPU_CHECK_TESTED)
        ext_lines.append(f"{e[5]:26s} {e[0]:14s} {Path(e[1]).name + '.obj':34s} "
                         + (" ".join(f"{k}={v}" for k, v in sorted(counts.items())) or "-")
                         + (f"  FAIL: not tested by cpu_check: {', '.join(untested)}" if untested else ""))
        if untested and e[0] != PROBE_TARGET:
            failures.append(f"{e[1]} ({e[5]}): instruction families not tested by cpu_check: {', '.join(untested)}")

    # 5. Dynamic initialisers.
    say("== 5. Dynamic initialisers (.CRT$XC*) in AVX2 objects ==")
    for e in linked_avx2:
        k = dynamic_initialisers(e[2])
        say(f"{e[5]:26s} {e[0]:14s} {Path(e[1]).name + '.obj':34s} {k}")
        if k:
            failures.append(f"dynamic initialiser in AVX2 object {e[1]}")
    say("")

    say("== 6. Instruction families beyond AVX/AVX2 in /arch:AVX2 objects ==")
    say(f"(cpu_check tests: AVX2, OS AVX state, {', '.join(sorted(CPU_CHECK_TESTED))}; any other family FAILs)")
    for ln in ext_lines:
        say(ln)
    say("")

    say("INFO list: shared functions whose AVX2 copy has VEX code, not counted as hazards:")
    for i in info:
        say("  - " + i)
    if not info:
        say("  (none)")
    say("")
    say("RESULT: " + ("PASS" if not failures else "FAIL"))
    for f in failures:
        say("  - " + f)

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
