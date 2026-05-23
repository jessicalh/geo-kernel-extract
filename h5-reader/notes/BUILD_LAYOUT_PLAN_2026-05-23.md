# h5-reader Build System Layout Plan — 2026-05-23

**Status:** DRAFT, pending user review. No files change against this
plan until approved.

**Goal:** organise `h5-reader/`'s build so a single source tree builds
cleanly under all four runtime environments × two configurations,
with platform and config differences absorbed behind named indirection
rather than `if(UNIX)` chains in `CMakeLists.txt`.

---

## 1. Build targets

| Tag                | Host                       | Toolchain          | Use                              |
|--------------------|----------------------------|--------------------|----------------------------------|
| `linux-x86_64`     | This workstation, lab boxes| GCC 13 / apt       | Primary dev + adviser binary     |
| `linux-aarch64`    | DGX Spark                  | GCC 13 / apt       | ARM64 server build               |
| `mac-arm64`        | macOS 14+ Apple Silicon    | AppleClang         | Adviser binary; Mac primary      |
| `win-x86_64`       | Windows 11 MSVC 2022       | MSVC v143 / clang-cl | Adviser binary; Windows primary |

DGX Spark uses the same Linux presets — only `CMAKE_SYSTEM_PROCESSOR`
differs, which CMake auto-detects. No separate preset family needed.

## 2. Version pins per target

The pin matrix is **deliberately divergent** between Linux and
Mac/Win because the Linux environments are tied to system packaging
(apt + locally-built VTK) and the Mac/Win adviser builds should
follow the latest Qt commercial release plus VTK current.

| Dep      | Linux (x86_64 + aarch64)               | macOS (arm64)           | Windows (x86_64)        |
|----------|----------------------------------------|-------------------------|-------------------------|
| Qt6      | 6.4.x (apt: `qt6-base-dev`)            | Latest Qt Pro (6.10.x)  | Latest Qt Pro (6.10.x)  |
| VTK      | 9.5.x (built from source, `~/VTK/`)    | Latest 9.x              | Latest 9.x              |
| HDF5     | 1.10.x (apt: `libhdf5-dev`)            | Latest 1.14 (brew)      | Latest 1.14 (vcpkg)     |
| Eigen3   | 3.4.x (apt: `libeigen3-dev`)           | 3.4.x (brew)            | 3.4.x (vcpkg)           |
| HighFive | Vendored at `extern/HighFive/include/` | Same                    | Same                    |

The HDF5 split (1.10 on Linux, 1.14 on Mac/Win) is the load-bearing
divergence. HighFive picks `H5Dvlen_reclaim` (≤1.12) vs `H5Treclaim`
(≥1.12) from its visible header; Linux must keep the system 1.10
header path **first** on the include line, Mac/Win must not have any
such workaround. This lives entirely inside `Platform-*.cmake`.

## 3. Indirection structure

### 3.1 Preset tree (`CMakePresets.json`)

```
base                         no platform, no config — generator + compile flags only
├── linux-base               apt + ~/VTK/ defaults; HDF5 1.10
│   ├── linux-debug          CMAKE_BUILD_TYPE=Debug
│   ├── linux-rwdi           CMAKE_BUILD_TYPE=RelWithDebInfo   ← default dev
│   └── linux-release        CMAKE_BUILD_TYPE=Release
├── mac-base                 Homebrew prefix defaults; HDF5 1.14
│   ├── mac-debug
│   ├── mac-rwdi
│   └── mac-release
└── win-base                 vcpkg + Qt installer defaults; HDF5 1.14
    ├── win-debug
    ├── win-rwdi
    └── win-release
```

13 configure presets. Each leaf is two lines (`inherits` +
`CMAKE_BUILD_TYPE`). Matching `buildPresets` and `testPresets` 1:1.
A `packagePresets` block lands later when we wire CPack.

`condition: equals hostSystemName <PLATFORM>` filters per-host so a
listing on Linux only shows the four Linux presets.

### 3.2 `cmake/` helper modules

```
h5-reader/cmake/
├── Platform-Linux.cmake     HDF5 1.10 BEFORE-include workaround;
│                            VTK_DIR default = $ENV{HOME}/VTK/lib/cmake/vtk-9.5
│                            Qt6_DIR via apt default location.
├── Platform-Darwin.cmake    Homebrew prefix discovery (/opt/homebrew vs
│                            /usr/local); HDF5_ROOT = brew --prefix hdf5;
│                            CMAKE_PREFIX_PATH += $ENV{HOME}/Qt/6.10.x/macos.
├── Platform-Windows.cmake   vcpkg toolchain file inclusion (if available);
│                            Qt6_DIR = C:\Qt\6.10.x\msvc2022_64;
│                            HDF5 + Eigen via vcpkg.
├── BuildType-Debug.cmake    -O0 -g3 -fsanitize=address,undefined on UNIX;
│                            /fsanitize=address on MSVC; -DH5READER_DEBUG=1.
├── BuildType-Release.cmake  -O3 + LTO on UNIX; /O2 /GL on MSVC; NDEBUG.
│                            (RelWithDebInfo uses CMake's defaults — no
│                            module needed.)
└── DiscoverDeps.cmake       Single block of find_package(Qt6) +
                             find_package(VTK) + find_package(HDF5) +
                             find_package(Eigen3) + HighFive include
                             validation. Consumes H5READER_*_DIR hints.
```

`CMakeLists.txt` becomes thin orchestration:

```cmake
cmake_minimum_required(VERSION 3.21)
project(h5reader LANGUAGES CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(Platform-${CMAKE_SYSTEM_NAME})   # absorbs Linux/Darwin/Windows quirks
include(BuildType-${CMAKE_BUILD_TYPE} OPTIONAL)
include(DiscoverDeps)

add_executable(h5reader ${SOURCES})
target_link_libraries(h5reader PRIVATE ...)
```

No `if(UNIX)` / `if(WIN32)` blocks in the top file.

### 3.3 Hint cache variables

| Variable               | Purpose                                                  |
|------------------------|----------------------------------------------------------|
| `H5READER_QT_DIR`      | Prepended to `CMAKE_PREFIX_PATH` before `find_package(Qt6)` |
| `H5READER_VTK_DIR`     | Passed as `PATHS` to `find_package(VTK)`                 |
| `H5READER_HDF5_DIR`    | Becomes `HDF5_ROOT` before `find_package(HDF5)`          |
| `H5READER_EIGEN_DIR`   | Becomes `Eigen3_DIR` if set                              |
| `H5READER_VCPKG_ROOT`  | Windows only — sets `CMAKE_TOOLCHAIN_FILE`               |

Each `Platform-*.cmake` sets sane defaults; users override via
`cmake --preset <p> -DH5READER_QT_DIR=...`. Documented in README.

## 4. Vendoring policy

| Dep      | Vendor? | Reason                                                     |
|----------|---------|------------------------------------------------------------|
| HighFive | YES     | Header-only, ~80 KB. Already at `extern/HighFive/`.        |
| Eigen3   | NO      | Header-only but ~5 MB and ubiquitous; system is fine.       |
| Qt6      | NO      | Too large; commercial license install (Mac/Win) preferred. |
| VTK      | NO      | Built from source on every platform anyway.                |
| HDF5     | NO      | System binary; we pin the version, not the package source. |

If Eigen ever turns out to be unreliable across platforms we can vendor
it (`extern/Eigen/`) without disrupting anything else — just a fallback.

## 5. Migration sequence

Each step lands as one commit. Build is verified Linux-green before
moving to the next.

1. **Add `cmake/` modules without changing behaviour.** Create the
   four files; populate `Platform-Linux.cmake` and `DiscoverDeps.cmake`
   with the current logic verbatim. `CMakeLists.txt` adds the
   `include()` calls; the inline blocks they replace are deleted in
   the same commit. Linux build must stay green. Mac/Win modules are
   empty stubs at this stage.
2. **Refactor `CMakeLists.txt` to thin orchestration shape** (the
   block in §3.2). Source-file list stays. Verify Linux build.
3. **Add the missing 9 presets** (`*-debug`, `*-rwdi`, `*-release`
   for mac + win, and `linux-release`). Verify each Linux leaf
   configures and builds.
4. **Populate `Platform-Darwin.cmake`** based on the Homebrew install
   commands documented in README. Cannot verify here; user verifies
   first time on a Mac.
5. **Populate `Platform-Windows.cmake`** based on the Qt installer +
   vcpkg install commands documented in README. Cannot verify here.
6. **Add `BuildType-Debug.cmake` + `BuildType-Release.cmake`** with
   the sanitizer / LTO flags. Verify Linux-debug builds + runs.
7. **Add `launch_reader.ps1`** for Windows. Cross-check the existing
   `launch_reader.sh` works on Mac too (it should — Bash + DISPLAY
   stuff is Linux-only; the Mac path is `./build/mac-rwdi/h5reader.app`).
8. **README rewrite** — per-platform install + build sections,
   examples for each preset, cache-var override examples.

Steps 1–3 + 6 + 8 happen here on Linux. Steps 4, 5, 7 happen when the
user syncs to Mac and Windows respectively.

## 6. Per-platform install commands (target state for README)

### Linux (Ubuntu 24.04+)

```bash
sudo apt install qt6-base-dev qt6-charts-dev libhdf5-dev libeigen3-dev \
                 cmake ninja-build g++
# VTK: build from source per existing VTK_BUILD instructions
cmake --preset linux-rwdi
cmake --build --preset linux-rwdi -j$(nproc)
```

### macOS (Apple Silicon, 14+)

```bash
# Qt Pro 6.10.x via official installer to ~/Qt/
brew install hdf5 eigen cmake ninja vtk
cmake --preset mac-rwdi \
      -DH5READER_QT_DIR=$HOME/Qt/6.10.x/macos
cmake --build --preset mac-rwdi -j$(sysctl -n hw.ncpu)
```

### Windows (11 + MSVC 2022)

```powershell
# Qt Pro 6.10.x via official installer to C:\Qt\
# vcpkg via https://vcpkg.io/, then:
vcpkg install hdf5:x64-windows eigen3:x64-windows vtk:x64-windows
cmake --preset win-rwdi `
      -DH5READER_QT_DIR=C:\Qt\6.10.x\msvc2022_64 `
      -DH5READER_VCPKG_ROOT=C:\vcpkg
cmake --build --preset win-rwdi
```

## 7. Out of scope (explicitly)

- Cross-compilation between platforms.
- CPack / installer generation (lands after the cross-platform builds
  are green; separate session).
- Universal binaries on macOS (single-arch only — match the host).
- Static linking of Qt or VTK (we want shared, smaller binaries +
  easier license compliance on Mac/Win Pro licenses).
- Supporting both HDF5 1.10 and 1.14 on the *same* platform.
- Auto-fetching VTK or Qt via FetchContent — every platform's binary
  install dance is documented manually.

## 8. Open questions for the user

1. **Linux Qt pin.** Apt on Ubuntu 24.04 ships Qt 6.4.2. Acceptable, or
   should Linux also move to Qt 6.10 via the Qt installer / Qt OSS
   build? (Sticking with apt is simpler; moving aligns with Mac/Win.)
2. **HDF5 1.14 migration on Linux later.** Worth it now while we're
   restructuring, or defer to a separate decision?
3. **VTK on macOS.** Use Homebrew's `vtk` formula (currently 9.x) or
   build from source like we do on Linux? Homebrew is simpler but
   pulls a lot of transitive deps.
4. **Sanitizers in Debug.** `-fsanitize=address,undefined` is the
   default Debug config in §3.2 — OK, or just `-O0 -g3` (no
   sanitizer)? Sanitizers cost ~2× build time but catch real bugs.
5. **DGX Spark.** Confirm it can run `apt install qt6-base-dev` (i.e.
   that the Spark's Ubuntu has the Qt6 packages in its arm64 repos),
   or do we need a different Qt source for that environment?

Pick once and the rest of the structure falls out. Recommend
defaults: keep Linux Qt at 6.4 apt, defer HDF5 1.14, use Homebrew VTK
on macOS, enable sanitizers in Debug, verify Spark apt repos before
counting on apt-Qt6 there.
