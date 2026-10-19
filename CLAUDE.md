# NTN Channel Emulator
LEO NTN (S-band, 3GPP Rel-17) channel emulator: C++ engine + Python/SGP4
control plane. Portfolio project. Windows native (MSVC x64, Ninja).

## Build & test
Run from the repo root in PowerShell. The MSVC environment only exists inside
the `cmd /c` call, so configure, build, and ctest must all run inside it.
Always use the Build Tools install, amd64 target on an amd64 host.

```powershell
# Put vswhere on PATH (VsDevCmd.bat also calls it) and find the Build Tools
$env:PATH = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer;$env:PATH"
$vs = vswhere -latest -products Microsoft.VisualStudio.Product.BuildTools -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath

# Configure + build (Release, Ninja, x64)
cmd /c "`"$vs\Common7\Tools\VsDevCmd.bat`" -arch=amd64 -host_arch=amd64 -no_logo && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build"

# Sanity check: must print "Pointer size: 64 bit"
.\build\hello.exe
```

- C++ tests: append `&& ctest --test-dir build --output-on-failure` inside
  the same `cmd /c` string.
- Python tests: `.\.venv\Scripts\python.exe -m pytest python/tests`
- Always call Python as `.\.venv\Scripts\python.exe`.
- Benchmarks: Release builds only.
- `build/compile_commands.json` is generated for clangd.

## Physics conventions (wrong = silent bugs)
- One source of truth: range track rho(t). tau = rho/c;
  f_d = -f_c * rho_dot / c; carrier phase = -2*pi*f_c*tau(t).
  Never derive Doppler and delay independently.
- Positive Doppler = satellite approaching.
- Evaluate phase in float64 from segment polynomials; never accumulate
  phase in float32.
- Pipeline order: TDL fading first, then geometric delay/Doppler.
- Control plane passes C1-continuous polynomial segments via pybind11,
  never sampled vectors.
- SigMF: datatype cf32_le; core:frequency goes in captures, not global.
- SIMD via xsimd or Highway only, no raw intrinsics.

## Rules
- IMPORTANT: never edit tests/reference/ or the float64 golden model to
  make a test pass. Stop and tell me.
- Tests come from closed-form expectations (FFT tone peak, exact integer
  shifts, ideal-sinc reference), never from the implementation's output.
- No performance or accuracy claim in README without a measured number in
  docs/results.md.
- Each milestone ends with: tests green, numbers recorded, commit.
