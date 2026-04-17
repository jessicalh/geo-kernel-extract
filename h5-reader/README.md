# h5-reader

Standalone Qt6/VTK reader for analysis H5 files produced by
`nmr_extract --trajectory --analysis`.

## What it does

Opens an `{protein_id}_analysis.h5` file, loads all ~186 datasets into
memory, builds a typed protein/conformation/frame object model, and
presents an interactive 3D viewer with per-frame animation. Every
calculator's contributions are visualisable — ring-current isosurfaces
(BS and HM), per-atom tensor glyphs, electric field arrows, water
dipoles, H-bond vectors, SASA surfaces, DSSP-coloured backbone, bond-
order tubes, and per-atom inspection on click.

The reader is read-only: it never writes H5 files, never re-runs the
extraction pipeline, and never modifies its input.

## Dependencies

- **CMake** ≥ 3.21 (for presets v3)
- **Qt6** ≥ 6.4, commercial or open-source — Widgets, OpenGLWidgets, Network
- **VTK** ≥ 9.0, built with Qt support — modules listed in CMakeLists.txt
- **HDF5** C library (system install — libhdf5-dev or equivalent)
- **HighFive** headers — vendored in `../extern/HighFive/` relative to the
  reader. If you have a zip of `h5-reader/` alone, clone HighFive:
  `git clone https://github.com/BlueBrain/HighFive extern/HighFive`
  and point `-DHIGHFIVE_INCLUDE_DIR=` at `extern/HighFive/include`.
- **`../fileformat/analysis_file.{h,cpp}`** — the frozen H5 (de)serialiser,
  compiled in directly.

## Build

### Linux

```sh
cd h5-reader
cmake --preset linux-gcc
cmake --build --preset linux-gcc
./build/linux-gcc/h5reader --help
```

### macOS

```sh
cd h5-reader
cmake --preset mac-clang
cmake --build --preset mac-clang
./build/mac-clang/h5reader --help
```

Qt via Homebrew or the official installer. If VTK is not on the default
search path, pass `-DVTK_DIR=/path/to/VTK/lib/cmake/vtk-9.x`.

### Windows (MSVC)

From "x64 Native Tools Command Prompt for VS 2022":

```
cd h5-reader
cmake --preset win-msvc
cmake --build --preset win-msvc
.\build\win-msvc\h5reader.exe --help
```

After build, deploy the Qt and VTK DLLs next to the exe:

```
windeployqt.exe .\build\win-msvc\h5reader.exe
xcopy /Y C:\path\to\VTK\bin\*.dll .\build\win-msvc\
```

## Run

```
h5reader path/to/1B1V_4292_analysis.h5
```

With a separate UDP log tail in another terminal:

```
python3 udp_listen.py
```

Log messages (startup, H5 load, pick events, frame changes, render timings)
flow to both stderr and the UDP stream on port 9997. The reader's own
"Operations Log" dock also shows the stream in real time — but only one
socket can bind port 9997 on Linux at a time, so run `udp_listen.py`
ONLY when the reader is not running.

## Editor / language-server setup (clangd)

The build emits `build/<preset>/compile_commands.json`. Point your editor's
language server at that directory — do NOT symlink or copy the file to the
project root (both are fragile across platforms and presets).

**VS Code** (with the clangd extension): add to `.vscode/settings.json`,
```
"clangd.arguments": ["--compile-commands-dir=build/linux-gcc"]
```
replacing `linux-gcc` with your platform preset.

**Neovim / Vim** (coc-clangd, nvim-lspconfig, etc.): configure
`clangd` with `--compile-commands-dir=build/<preset>` in your LSP
client init.

**CLion / Qt Creator**: auto-detect the build directory; no
configuration needed.

**Command-line clangd**: invoke with
`clangd --compile-commands-dir=build/<preset>`.

Each developer configures their own editor. The project does not commit
a `.clangd` file because it would hardcode one preset name and silently
break for anyone using a different one (debug build, Windows preset,
macOS preset).

## Platform status

| Platform | Build verified | Crash capture | Clean-quit signals | Notes |
|----------|---------------|---------------|--------------------|-------|
| Linux GCC | yes (Qt 6.4.2, VTK 9.5.2, HDF5 1.10.10) | yes (sigaction + backtrace) | yes (self-pipe + SIGINT/SIGTERM → aboutToQuit → VTK finalise) | Primary development target |
| macOS Clang | preset present, not yet built on hardware | yes (sigaction + backtrace) | yes (same POSIX path as Linux) | Needs first build to verify |
| Windows MSVC | preset present, not yet built on hardware | stub (MiniDumpWriteDump pending) | yes (SetConsoleCtrlHandler + QMetaObject::invokeMethod quit) | Signal path compiled conditionally, untested on hardware |

Every platform exposes the same API — `CrashHandler::Install()` and
`InstallShutdownSignalHandlers()` — behind `#ifdef` internal branching.
No caller writes `#ifdef _WIN32`.

## Scope and limits

See `notes/SCOPE.md` for the full scope statement. Short version:

- Single protein per H5. The reader does not load multiple H5s into one
  scene.
- Trajectory-animated. `QtConformation` is the trajectory; `QtFrame` is
  one sampled XTC frame.
- Rendering is honest per-frame: every frame reads its own data, runs
  any needed closed-form kernel re-evaluation (only BS and HM volumetric
  grids), and renders the result. No interpolation, no precomputed
  keyframes, no procedural fakery.
- Performance: on a fast workstation BS/HM re-evaluation is insignificant
  at normal ring counts. On slower hardware, enable the volumetric
  overlays selectively.
