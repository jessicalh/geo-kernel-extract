# h5-reader

Standalone Qt6/VTK reader for the per-TR-emits-its-own-group
`trajectory.h5` + 5-NPY topology sidecar produced by
`nmr_extract --trajectory` post-2026-05-13.

## What it does

Opens a per-protein analysis tree (`trajectory.h5` plus
`atoms_category_info.npy`, `residues.npy`, `bonds.npy`, `rings.npy`,
`ring_membership.npy`, and `extraction_manifest.json`), builds a typed
`QtProtein` / `QtConformation` / `QtFrame` object model, and presents
an interactive 3D viewer with per-frame animation. Calculator
contributions are visualisable — ring-current isosurfaces (BS and HM),
per-atom tensor glyphs, electric field arrows, water dipoles,
SASA surfaces, DSSP-coloured backbone, bond-order tubes, per-atom
inspection on click.

Read-only: never writes H5, never re-runs the extraction pipeline.

## Build matrix

| Platform     | Compiler         | Configs available                     | Default dev preset |
|--------------|------------------|---------------------------------------|--------------------|
| Linux x86_64 | GCC 13 (apt)     | `linux-debug` / `linux-rwdi` / `linux-release` | `linux-rwdi` |
| Linux ARM64  | GCC 13 (apt)     | same as Linux x86_64 (Ubuntu 24 LTS on DGX Spark) | `linux-rwdi` |
| macOS arm64  | AppleClang       | `mac-debug` / `mac-rwdi` / `mac-release` | `mac-rwdi` |
| Windows x64  | MSVC v143        | `win-debug` / `win-rwdi` / `win-release` | `win-rwdi` |

Build directory is always `build/<preset>/`, e.g. `build/linux-rwdi/`.

Per-platform setup, dependency paths, and target flags live in
`cmake/Platform-{Linux,Darwin,Windows}.cmake` — no `if(UNIX)/if(WIN32)`
branches in `CMakeLists.txt`. See
`notes/BUILD_LAYOUT_PLAN_2026-05-23.md` for the indirection rationale.

## Dependencies + pinned versions

| Dep      | Linux                       | macOS                        | Windows                       |
|----------|-----------------------------|------------------------------|-------------------------------|
| Qt6      | 6.4.x (apt)                 | Latest Qt Pro (installer)    | Latest Qt Pro (installer)     |
| VTK      | 9.5.x (from source, `~/VTK/`) | 9.x latest (from source)   | 9.x latest (from source)      |
| HDF5     | 1.10.x (apt)                | 1.14 (brew: `hdf5`)          | 1.14 (vcpkg: `hdf5`)          |
| Eigen3   | 3.4.x (apt)                 | 3.4.x (brew: `eigen`)        | 3.4.x (vcpkg: `eigen3`)       |
| HighFive | vendored at `../extern/HighFive/include/` | same             | same                          |
| CMake    | ≥ 3.21                      | ≥ 3.21                       | ≥ 3.21                        |
| Ninja    | ≥ 1.10                      | ≥ 1.10                       | ≥ 1.10                        |

The HDF5 1.10/1.14 split is deliberate: Linux stays on apt's 1.10
(matches the rest of the pipeline); Mac/Win get the current 1.14 that
their package managers ship. HighFive picks `H5Dvlen_reclaim` (≤1.12)
vs `H5Treclaim` (≥1.12) from the visible header; on Linux a
`BEFORE PRIVATE` include forces the system 1.10 header to win over any
transitive 1.14 path (handled in `cmake/Platform-Linux.cmake`).

If HighFive is missing from `../extern/HighFive/`, clone it:

```sh
git clone https://github.com/BlueBrain/HighFive ../extern/HighFive
```

or pass `-DHIGHFIVE_INCLUDE_DIR=/path/to/highfive/include` at configure time.

---

## Build — Linux (x86_64 or aarch64)

```sh
sudo apt install qt6-base-dev qt6-charts-dev libhdf5-dev libeigen3-dev \
                 cmake ninja-build g++
# VTK 9.5 — build from source following the VTK Build Guide.
# Default search location: $ENV{HOME}/VTK/lib/cmake/vtk-9.5/

cd h5-reader
cmake --preset linux-rwdi
cmake --build --preset linux-rwdi -j$(nproc)
./build/linux-rwdi/h5reader path/to/trajectory.h5
```

Alternative configs: `linux-debug` (sanitizer-friendly), `linux-release`
(stripped, LTO).

---

## Build — macOS (Apple Silicon, 14+)

```sh
# Qt Pro 6.10.x via the official Qt installer → ~/Qt/6.10.x/macos/
# VTK from source (see VTK Build Guide); default: ~/VTK/lib/cmake/vtk-9.x
brew install hdf5 eigen cmake ninja

cd h5-reader
cmake --preset mac-rwdi \
      -DH5READER_QT_DIR=$HOME/Qt/6.10.x/macos \
      -DH5READER_VTK_DIR=$HOME/VTK
cmake --build --preset mac-rwdi -j$(sysctl -n hw.ncpu)
./build/mac-rwdi/h5reader path/to/trajectory.h5
```

The `H5READER_*_DIR` hints are needed first time — once configured,
they're cached in `build/mac-rwdi/CMakeCache.txt` and you can drop them
from subsequent `--preset` calls.

`cmake/Platform-Darwin.cmake` is a stub at first sync — populate it
with the actual `H5READER_QT_DIR` / `H5READER_VTK_DIR` defaults your
Mac uses so future configures don't need the `-D` flags. See
[Troubleshooting → macOS](#macos-troubleshooting) for the common
failures and how to diagnose them.

---

## Build — Windows (11, MSVC 2022)

Open **x64 Native Tools Command Prompt for VS 2022**, then PowerShell:

```powershell
# Qt Pro 6.10.x via official installer → C:\Qt\6.10.x\msvc2022_64
# VTK from source (built per VTK Build Guide) → e.g. C:\VTK
# vcpkg per https://vcpkg.io/, then:
vcpkg install hdf5:x64-windows eigen3:x64-windows

cd h5-reader
cmake --preset win-rwdi `
      -DH5READER_QT_DIR="C:\Qt\6.10.x\msvc2022_64" `
      -DH5READER_VTK_DIR="C:\VTK" `
      -DH5READER_VCPKG_ROOT="C:\vcpkg"
cmake --build --preset win-rwdi
.\build\win-rwdi\h5reader.exe path\to\trajectory.h5
```

After build, deploy Qt + VTK DLLs alongside the exe:

```powershell
windeployqt.exe .\build\win-rwdi\h5reader.exe
xcopy /Y C:\VTK\bin\*.dll .\build\win-rwdi\
```

`cmake/Platform-Windows.cmake` is a stub at first sync — populate it
with your install paths so subsequent configures don't need the `-D`
flags. See [Troubleshooting → Windows](#windows-troubleshooting).

---

## Preset reference

| Name             | Build type        | Default? |
|------------------|-------------------|----------|
| `linux-debug`    | `Debug`           |          |
| `linux-rwdi`     | `RelWithDebInfo`  | ✓ dev    |
| `linux-release`  | `Release`         |          |
| `mac-debug`      | `Debug`           |          |
| `mac-rwdi`       | `RelWithDebInfo`  | ✓ dev    |
| `mac-release`    | `Release`         |          |
| `win-debug`      | `Debug`           |          |
| `win-rwdi`       | `RelWithDebInfo`  | ✓ dev    |
| `win-release`    | `Release`         |          |

Each preset writes to `build/<preset>/`. Presets with a `condition`
matching the host system are auto-hidden on other platforms by `cmake
--list-presets`.

## Hint cache variables

Override per environment with `-D<NAME>=<path>` at configure time, or
set them as defaults in `cmake/Platform-<OS>.cmake`.

| Variable               | Purpose                                                  |
|------------------------|----------------------------------------------------------|
| `H5READER_QT_DIR`      | Prepended to `CMAKE_PREFIX_PATH` before `find_package(Qt6)` |
| `H5READER_VTK_DIR`     | Prepended to `CMAKE_PREFIX_PATH` before `find_package(VTK)` |
| `H5READER_HDF5_DIR`    | Becomes `HDF5_ROOT` before `find_package(HDF5)`          |
| `H5READER_EIGEN_DIR`   | Becomes `Eigen3_DIR` (full path to the `cmake/` dir)     |
| `HIGHFIVE_INCLUDE_DIR` | Path to HighFive `include/` (default: `../extern/HighFive/include`) |

## Running the reader

Linux / macOS:

```sh
./build/<preset>/h5reader path/to/trajectory.h5
# or:
H5READER_PRESET=mac-rwdi ./launch_reader.sh path/to/trajectory.h5
```

Windows:

```powershell
.\build\win-rwdi\h5reader.exe path\to\trajectory.h5
# or:
.\launch_reader.ps1 path\to\trajectory.h5
```

Log stream: messages flow to stderr and to UDP port 9997 (unicast,
single subscriber on Linux). In another terminal:

```sh
python3 udp_listen.py
```

Do not run `udp_listen.py` and the reader at the same time on Linux —
they compete for the same UDP socket. The reader's own "Operations
Log" dock displays the stream live.

## Editor / language-server setup (clangd)

The build emits `build/<preset>/compile_commands.json`. Point your
editor's clangd at that directory — never symlink the file to the
project root (breaks on Windows, drifts as you switch presets).

**VS Code** (with the clangd extension), `.vscode/settings.json`:

```json
{
    "clangd.arguments": ["--compile-commands-dir=build/linux-rwdi"]
}
```

Replace `linux-rwdi` with whatever preset you build with.

**Neovim/Vim**: configure clangd with
`--compile-commands-dir=build/<preset>` in your LSP init.

**CLion / Qt Creator**: auto-detects the build directory.

**Command-line clangd**:

```sh
clangd --compile-commands-dir=build/<preset>
```

No `.clangd` file is committed — it would hardcode one preset and
silently break for anyone using a different build type or platform.

---

## Troubleshooting

### Linux troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `find_package(VTK)` fails | VTK not on `CMAKE_PREFIX_PATH` | `-DH5READER_VTK_DIR=$HOME/VTK` |
| `find_package(Qt6)` fails | apt Qt6 missing | `sudo apt install qt6-base-dev qt6-charts-dev` |
| HighFive `H5Treclaim` link error | System HDF5 1.10 vs micromamba/OpenBabel 1.14 transitive include | Already handled by `cmake/Platform-Linux.cmake` (`BEFORE PRIVATE /usr/include/hdf5/serial`); verify the include is present |
| Black window on launch | Missing VTK runtime libs | Check that `~/VTK/lib/*.so` is on `LD_LIBRARY_PATH` if VTK is not in a system location |
| `udp_listen.py` shows nothing | Reader holds the UDP socket | Close the reader first |

### macOS troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `find_package(Qt6)` fails | Qt installer not on `CMAKE_PREFIX_PATH` | `-DH5READER_QT_DIR=$HOME/Qt/6.10.x/macos` |
| `find_package(HDF5)` fails | Homebrew HDF5 not found by CMake | `-DH5READER_HDF5_DIR=$(brew --prefix hdf5)` |
| `find_package(Eigen3)` fails | Homebrew Eigen layout mismatch | `-DH5READER_EIGEN_DIR=$(brew --prefix eigen)/share/eigen3/cmake` |
| Linker error: `_main referenced from` | Qt deploy missing | Run `~/Qt/<ver>/macos/bin/macdeployqt build/mac-rwdi/h5reader.app` once |
| App bundle not generated | CMakeLists doesn't set `MACOSX_BUNDLE` on the target | (Open issue — currently the binary is a flat ELF/Mach-O; bundle support lands when the Mac sync starts) |
| `Library not loaded: @rpath/QtCore.framework` | rpath not set | `install_name_tool -add_rpath ~/Qt/6.10.x/macos/lib build/mac-rwdi/h5reader` |

### Windows troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `cl.exe` not found | Not running in VS dev prompt | Open "x64 Native Tools Command Prompt for VS 2022" first |
| `find_package(Qt6)` fails | Qt installer not on `CMAKE_PREFIX_PATH` | `-DH5READER_QT_DIR=C:\Qt\6.10.x\msvc2022_64` |
| `find_package(HDF5)` fails | vcpkg toolchain not used | Add `-DCMAKE_TOOLCHAIN_FILE=C:\vcpkg\scripts\buildsystems\vcpkg.cmake` (or set `H5READER_VCPKG_ROOT` in `Platform-Windows.cmake`) |
| App launches then dies with "Qt platform plugin could not be initialized" | Missing Qt DLLs next to exe | Run `windeployqt.exe .\build\win-rwdi\h5reader.exe` |
| App launches then dies with "vtkXXX.dll missing" | VTK DLLs not next to exe | `xcopy /Y C:\VTK\bin\*.dll .\build\win-rwdi\` |
| Minidumps in `%LOCALAPPDATA%\Temp` not written | Dbghelp.dll not linked | Already linked by `cmake/Platform-Windows.cmake`; check the exe imports list |

When a brand-new platform first builds, expect the first 2-3 configures
to fail until `H5READER_*_DIR` hints are tuned. Once a configure
succeeds, encode the working paths into
`cmake/Platform-<OS>.cmake` so future fresh checkouts find them
automatically.

---

## Platform status

| Platform | Build verified | Crash capture | Clean-quit signals |
|----------|----------------|---------------|--------------------|
| Linux x86_64 GCC | yes (Qt 6.4.2, VTK 9.5.2, HDF5 1.10.10) | sigaction + backtrace | self-pipe + SIGINT/SIGTERM → `aboutToQuit` → VTK finalise |
| Linux aarch64 GCC (DGX Spark) | preset present, not yet built | same as x86_64 | same as x86_64 |
| macOS arm64 Clang | preset present, not yet built | sigaction + backtrace | same POSIX path as Linux |
| Windows x64 MSVC | preset present, not yet built | stub (`MiniDumpWriteDump` pending wiring) | `SetConsoleCtrlHandler` + `QMetaObject::invokeMethod` quit |

Crash capture + clean-quit signals share a `CrashHandler::Install()` /
`InstallShutdownSignalHandlers()` surface — no caller writes `#ifdef
_WIN32`; the conditional code lives inside the handler classes.

## Scope and limits

See `notes/SCOPE.md` for the full statement. Short version:

- Single protein per H5. The reader does not load multiple H5s into
  one scene.
- Trajectory-animated. `QtConformation` is the trajectory; `QtFrame`
  is one sampled XTC frame.
- Rendering is honest per-frame: every frame reads its own data, runs
  any needed closed-form kernel re-evaluation (only BS and HM
  volumetric grids), and renders the result. No interpolation, no
  precomputed keyframes, no procedural fakery.
- Performance: on a fast workstation BS/HM re-evaluation is
  insignificant at normal ring counts. On slower hardware, enable the
  volumetric overlays selectively.
