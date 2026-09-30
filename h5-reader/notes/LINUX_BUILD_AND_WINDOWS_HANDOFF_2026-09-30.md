# Linux build and Windows handoff — 30 September 2026

## Scope and baseline

This work extends the current Qt Pro Reader with a Linux build using open source Qt, native CPU inference, multimedia, the published trajectory library, and an installed runtime bundle. The source baseline is `3c0303a32ade34eafcd6341591be337a811f1007` in `jessicalh/geo-kernel-extract`. The working checkout is `/home/jessica/src/h5-reader-linux`, branch `linux/oss-scope-20260930`; the standalone application is its `h5-reader` subdirectory. Read `LINUX_OPEN_SOURCE_SCOPE_2026-09-30.md` for the original investigation and acceptance matrix.

The actual binary/toolchain target is **Ubuntu 24.04, Linux x86_64**, GCC 13.3, CMake 3.28.3, Ninja 1.11.1, public Qt 6.10.2, VTK 9.5.2 built against that Qt, and PyTorch/LibTorch 2.7.1+cpu. HDF5 1.10.10 and Eigen 3.4 come from this host. A build on older Linux distributions, aarch64, or another operating system needs its own validation. The portable archive still depends on the host's glibc and graphics drivers.

Linux runtime deployment explicitly requires Qt **6.10 or newer**, because the plugin selection deployment API was introduced there. The shared application compilation requirement remains Qt >=6.8, including the existing Windows build.

The user has an active GPU workload. All work here uses CPU inference, Mesa software rendering, and software video encoding. **CUDA has not been integrated or validated by this work.** Wait for the user to release the GPU before any accelerator validation. Windows CPU/ROCm behavior is retained in the source; a fresh native Windows build and installed-runtime check remain necessary.

## Local toolchain and environment

The isolated toolchain root is `/home/jessica/opt/h5reader`:

| Component | Actual location |
| --- | --- |
| Open source Qt | `qt/6.10.2/gcc_64` |
| Matching VTK | `vtk-9.5.2-qt6.10.2` |
| Extracted Ubuntu LibArchive/patchelf packages | `deps/usr` |
| Test/export Python and CPU LibTorch | `cpu-venv` |
| Native helper build | `build/shielding-infer-cpu` |
| Model/helper/library staging | `ml-runtime` |
| Read-only data wrappers | `fixtures` |
| Build, acceptance, and packaging logs | `logs` |

The host shell originally supplied `/opt/orca/lib` and accelerator paths through `LD_LIBRARY_PATH`; `ldd` consequently selected ORCA's C++ runtime. Clear `LD_LIBRARY_PATH` and `LD_PRELOAD` for configuration, building, testing, and installation. Do not change the global shell or another application's environment.

Changing `H5READER_QT_DIR` or `H5READER_VTK_DIR` alone does not clear cached `Qt6_DIR`, individual Qt module directories, or `VTK_DIR`. **Use `cmake --fresh` when switching toolchain prefixes**, or start a new build directory. The previous VTK installation under `/home/jessica/VTK` used Qt 6.4.2. The dedicated VTK build now includes the selected Qt SDK's `lib` directory in its development runtime path, resolving the observed transitive QtOpenGL selection of the system Qt. Installed libraries receive paths relative to their own bundle locations.

Qt installation and the local dependency commands are recorded in `logs/toolchain-commands.sh`. The actual Qt installation command was:

```bash
/home/jessica/opt/h5reader/aqt-venv/bin/python -m aqt install-qt \
  linux desktop 6.10.2 linux_gcc_64 \
  -O /home/jessica/opt/h5reader/qt \
  -m qtcharts qthttpserver qtmultimedia qtwebsockets
```

The aqt environment contains `aqtinstall==3.3.0`. The development dependency payload contains Ubuntu `libarchive-dev` and `libarchive13t64` version `3.7.2-2ubuntu0.8`, plus `patchelf` version `0.18.0-1.1build1`, downloaded with `apt-get download` and extracted with `dpkg-deb -x` into `deps`. That avoided modifying global packages. The matching VTK configuration, selected modules, CUDA-off setting, build, and installation are reproducible through the existing local script:

```bash
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  bash /home/jessica/opt/h5reader/build-vtk-qt6102.sh
```

Its source input is `/home/jessica/builds/VTK-9.5.2`. This is a local prerequisite, not a file downloaded by the application build.

## CPU helper and unchanged model assets

Python is used to obtain LibTorch and run developer tests. The shipped inference process is the C++ `infer` executable; it does not start Python. The pinned developer dependency inventory is `python-cpu-requirements.txt`. For a fresh local environment:

```bash
reader_root=/home/jessica/opt/h5reader
reader_source=/home/jessica/src/h5-reader-linux/h5-reader
python3 -m venv "$reader_root/cpu-venv"
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  "$reader_root/cpu-venv/bin/python" -m pip install \
  --index-url https://download.pytorch.org/whl/cpu 'torch==2.7.1+cpu'
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  "$reader_root/cpu-venv/bin/python" -m pip install \
  --extra-index-url https://download.pytorch.org/whl/cpu \
  -r "$reader_root/python-cpu-requirements.txt"

reader_torch="$reader_root/cpu-venv/lib/python3.12/site-packages/torch"
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cmake --fresh -S "$reader_source/tools/experimental_shielding_ml" \
  -B "$reader_root/build/shielding-infer-cpu" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$reader_torch"
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cmake --build "$reader_root/build/shielding-infer-cpu" --parallel 4

mkdir -p "$reader_root/ml-runtime/lib" "$reader_root/ml-runtime/licenses/torch"
cp "$reader_root/build/shielding-infer-cpu/infer" "$reader_root/ml-runtime/infer"
cp -a "$reader_torch/lib/." "$reader_root/ml-runtime/lib/"
cp "$reader_torch/../torch-2.7.1+cpu.dist-info/LICENSE" \
   "$reader_torch/../torch-2.7.1+cpu.dist-info/NOTICE" \
   "$reader_root/ml-runtime/licenses/torch/"
cp /shared/2026Thesis/nmr-shielding-release/h5-reader/models/experimental_shielding_ml/model.ts \
   /shared/2026Thesis/nmr-shielding-release/h5-reader/models/experimental_shielding_ml/manifest.json \
   "$reader_root/ml-runtime/"
chmod +x "$reader_root/ml-runtime/infer"
sha256sum "$reader_root/ml-runtime/model.ts"
```

Expected model SHA256: `101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6`. The model is `F006-R007-v1-reader-runtime`, schema 3, helper protocol 2. The scientific model, feature contract, weights, and expected numerical values are unchanged. The Linux helper honors LibTorch ABI flags and resolves its libraries from `$ORIGIN/lib`.

Linux requires `infer`, `model.ts`, `manifest.json`, and `lib/{libc10.so,libtorch.so,libtorch_cpu.so,libtorch_global_deps.so}`, with the remaining runtime dependencies also staged. Both Linux `auto` and `cpu` use CPU directly, with no configured accelerator fallback. Explicit ROCm is unavailable on Linux. The Windows DLL layout and ROCm selection/fallback remain under their platform branches.

## Real scientific fixtures

`fixtures/prepare-fixtures.py` creates local `.LGS` wrappers over existing shared scientific data; it does not generate scientific arrays or alter the sources. Run it with the system Python used in this session, which has h5py and NumPy, or the configured CPU virtual environment:

```bash
/home/jessica/opt/h5reader/cpu-venv/bin/python \
  /home/jessica/opt/h5reader/fixtures/prepare-fixtures.py
```

| Purpose | Absolute `.LGS` path | Validated data |
| --- | --- | --- |
| Main REST, ring-current, export/shutdown | `/home/jessica/opt/h5reader/fixtures/1p9j/1p9j-july20260717.LGS` | 846 atoms, 751 frames, 751 DFT outputs |
| Native ML and alternate reload | `/home/jessica/opt/h5reader/fixtures/f006/full720-f006-A0A822IKI2.LGS` | A0A822IKI2_WT, 528 atoms, one pose; all 52 manifest-declared input NPY files present |
| Chignolin DFT | `/home/jessica/opt/h5reader/fixtures/chignolin/chignolin-july20260717.LGS` | 166 atoms, 5,001 frames, 5,001 DFT outputs |
| Trp-cage DFT | `/home/jessica/opt/h5reader/fixtures/trpcage/trpcage-july20260717.LGS` | 304 atoms, 1,501 frames, 1,501 DFT outputs |

All 7,253 DFT metadata records were checked for successful ORCA exit, a nonempty primary output, and agreement with the HDF5 original frame index and time. The wrappers use the July staging extraction and its copied ORCA jobs. The older Linux 1P9J default points at a different, 1,501-frame extraction; it is not the 751-frame Windows parity fixture. Pass actual `.LGS` files, because some REST tests read the supplied path as JSON.

**Use the multi-frame 1P9J trajectory for `H5READER_MODEL_INPUT_REST_FIXTURE`.** A single-pose export can finish before the shutdown request, making it unsuitable for testing cancellation. The dedicated export/shutdown rerun with 1P9J passed both tests. The static F006 member remains the successful ML inference fixture.

## Configure, build, and test

The following reproduces the configured Linux application build. It assumes the toolchain, helper/model staging, and fixture preparation above have completed.

```bash
cd /home/jessica/src/h5-reader-linux/h5-reader
reader_root=/home/jessica/opt/h5reader
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cmake --fresh -S . -B build/linux-oss -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON \
  -DCMAKE_INSTALL_PREFIX="$reader_root/package-stage" \
  -DCMAKE_PREFIX_PATH="$reader_root/deps/usr" \
  -DH5READER_QT_DIR="$reader_root/qt/6.10.2/gcc_64" \
  -DQt6_DIR="$reader_root/qt/6.10.2/gcc_64/lib/cmake/Qt6" \
  -DH5READER_VTK_DIR="$reader_root/vtk-9.5.2-qt6.10.2" \
  -DVTK_DIR="$reader_root/vtk-9.5.2-qt6.10.2/lib/cmake/vtk-9.5" \
  -DPython3_EXECUTABLE="$reader_root/cpu-venv/bin/python" \
  -DH5READER_PATCHELF_EXECUTABLE="$reader_root/deps/usr/bin/patchelf" \
  -DH5READER_THIRD_PARTY_NOTICES_DIR="$reader_root/license-material" \
  -DH5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML=ON \
  -DH5READER_EXPERIMENTAL_SHIELDING_ML_ROOT="$reader_root/ml-runtime" \
  -DH5READER_REST_FIXTURE_DEFAULT="$reader_root/fixtures/1p9j/1p9j-july20260717.LGS" \
  -DH5READER_REST_RELOAD_FIXTURE_DEFAULT="$reader_root/fixtures/f006/full720-f006-A0A822IKI2.LGS" \
  -DH5READER_F006_STATIC_FIXTURE_DEFAULT="$reader_root/fixtures/f006/full720-f006-A0A822IKI2.LGS" \
  -DH5READER_CHIGNOLIN_FIXTURE_DEFAULT="$reader_root/fixtures/chignolin/chignolin-july20260717.LGS" \
  -DH5READER_TRPCAGE_FIXTURE_DEFAULT="$reader_root/fixtures/trpcage/trpcage-july20260717.LGS" \
  -DH5READER_ENABLE_MODEL_INPUT_CONTRACT_TEST=ON \
  -DH5READER_MODEL_INPUT_REST_FIXTURE="$reader_root/fixtures/1p9j/1p9j-july20260717.LGS"
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cmake --build build/linux-oss --parallel 8
scripts/linux/test_cpu.sh build/linux-oss -j 2
```

`scripts/linux/test_cpu.sh` clears foreign runtime/plugin paths, masks CUDA/HIP/ROCr devices, selects CPU inference, forces Mesa/llvmpipe rendering and software Qt multimedia, and provides Xvfb for native GUI tests as well as REST tests. Keep `ffmpeg` available for the video decoding checks. Successful ML tests require runtime assets and the complete fixture: `H5READER_EXPECT_ML_SUCCESS=1` now fails missing prerequisites instead of allowing setup to skip.

For a focused rerun, append `-R 'h5reader_rest_ml_cpu|h5reader_rest_installed_ml_stale_dev_env|h5reader_model_input'` to the same test script. The generic publication test needs its own original/portable pair, described below; a default skip is not publication validation. Conditional Chignolin, Trp-cage, and missing-input tests must be interpreted in their dedicated fixture jobs.

## Install and build the archive

```bash
cd /home/jessica/src/h5-reader-linux/h5-reader
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cmake --install build/linux-oss \
  --prefix /home/jessica/opt/h5reader/package-stage
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
  cpack --config build/linux-oss/CPackConfig.cmake \
  -G TGZ -B /home/jessica/opt/h5reader/dist \
  -D CPACK_PACKAGE_FILE_NAME=h5reader-0.5.0-Linux-OSS-x86_64-20260930
```

The command above names the archive `h5reader-0.5.0-Linux-OSS-x86_64-20260930.tar.gz`; the unmodified CPack default is `h5reader-0.5.0-Linux-x86_64.tar.gz`. The install includes the executable, catalog, Qt/VTK/native dependencies, Qt platform/multimedia/TLS plugins, matching FFmpeg libraries, native ML runtime, `launch-h5reader`, and `launch-h5reader-cpu`. Deployment normalizes copied ELF runtime paths and fails unresolved or conflicting dependencies. It excludes the target OS's glibc and GPU driver implementations. The launcher clears foreign Qt/library overrides and starts the copy inside its own bundle. `H5READER_BUNDLED_DATA_DIR` can supply prepared local examples; this build does not copy the full remote catalog into the package.

The optional notices input above was prepared at `license-material` from the exact Qt 6.10.2 module sources, FFmpeg 7.1.2, ICU 73.2, VTK and available system-package copyrights. Its provenance and unresolved release decisions are recorded in the installed `licenses/upstream/` directory. The archive's adjacent `.manifest.json` records source/dependency hashes and bundled entries. These notices do not assign a license to the application or model.

Test an extracted copy outside the source/build tree, using `launch-h5reader` as its entry point. Repeat successful F006 CPU inference, stale-development fallback, real-data screenshots, and software video export from that relocated copy. Record the exact archive and hash below once final packaging and relocated acceptance finish.

The relocated local application is `/home/jessica/Applications/h5-reader-linux-20260930`. While the GPU remains occupied, start `launch-h5reader-cpu`, or use the application-menu entry **H5 Reader (Linux OSS, CPU)**, provided by `/home/jessica/.local/share/applications/h5reader-linux-oss.desktop`. The CPU launcher was independently checked through the process environment and confirmed to use Mesa llvmpipe. Its evidence is `review/installed-ml/cpu-launcher-check.json` under the toolchain root.

## Evidence already produced

- Native CPU inference and stale-development fallback CTest jobs passed using F006. Their assertions preserve the expected atom-16 isotropic value `29.76103401184082` and rank-2 tensor values within `0.001`, including propagation to the dashboard, scene tensor, and inspector.
- The final configured CTest run passed **31/31 tests, 0 failures, in 55.05 seconds**, including the real original/portable publication pair (`logs/ctest-final.log` and its XML report). The initial full run used the single-pose fixture for shutdown cancellation; the corrected 1P9J export/shutdown rerun passed both tests (`logs/model-input-trajectory-acceptance.log`).
- The relocated installed application's main REST suite passed **59 tests**, with **11 fixture-specific skips** covered by the separate applicable gates. Installed stale-development ML acceptance passed **3 tests**, with the negative-input-only case intentionally skipped. The separate negative-input job passed **2 tests**. Installed atom-16 isotropic output was exactly `29.76103401184082`; the maximum absolute rank-2 component difference from the reference was `1.66893e-6`. Installed screenshots are `review/installed-ml/window.png` and `review/installed-ml/atom16-closeup.png`.
- Missing required ML input was diagnosed through the REST manifest and error response: 2 passed; the 2 successful-inference cases were intentionally skipped for that incomplete-input job (`logs/ml-missing-input-acceptance.log`).
- Native trajectory-library acceptance downloaded exactly 382,396,716 bytes for the smallest catalog bundle, `bmr7057-20260919-files`, through Qt HTTPS; libarchive unpacked its declared 1,040,325,393 bytes and the application loaded 842 atoms across 100 frames. Download, extraction, and load took 120.83 seconds. Switching away and reopening used the cache in 0.23 seconds with exactly one network transfer overall. Reader exited with code 0.
- That downloaded package was compared against the independent original at `/mnt/expansion/traj_extract_s15_mopac/bmr7057/run_20260614T221950Z_1361196_20260727T201159Z.lgs`. Its manifest bytes match `provenance/source.LGS`. The native publication test compared all 100 frames, all 842 atom coordinates, times/indices, and every catalog field, with exact snapshot bytes; 263 arrays were loaded per frame. Result: **3 passed, 0 failed, 0 skipped**, 30.25 seconds.
- Built-in `/api/screenshot` scene/window captures of the downloaded trajectory were inspected: molecular geometry and the 100-frame UI timeline render under software OpenGL. Screenshots are observations at the captured state, not a claim that every possible visual state was examined.

Native library artifacts are under `/home/jessica/opt/h5reader/native-library-acceptance`: `report.json`, `reader.log`, `downloaded-ui.json`, `downloaded-samples.json`, `downloaded-scene.png`, `downloaded-window.png`, and `publication-comparison.{json,log}`. The reusable local REST harness is `fixtures/check-native-library.py`. It requires a new output directory and fetches only the selected package:

```bash
/home/jessica/opt/h5reader/cpu-venv/bin/python \
  /home/jessica/opt/h5reader/fixtures/check-native-library.py \
  --binary /home/jessica/src/h5-reader-linux/h5-reader/build/linux-oss/h5reader \
  --fixture /home/jessica/opt/h5reader/fixtures/f006/full720-f006-A0A822IKI2.LGS \
  --output /home/jessica/opt/h5reader/native-library-acceptance-rerun
```

Reproduce the exact publication check against the already downloaded cache:

```bash
cd /home/jessica/src/h5-reader-linux/h5-reader
H5READER_PUBLICATION_ORIGINAL=/mnt/expansion/traj_extract_s15_mopac/bmr7057/run_20260614T221950Z_1361196_20260727T201159Z.lgs \
H5READER_PUBLICATION_PORTABLE='/home/jessica/opt/h5reader/native-library-acceptance/cache/Beardsley Lab/h5 reader/trajectories/bmr7057-20260919-files/run.LGS' \
  scripts/linux/test_cpu.sh build/linux-oss -R '^h5reader_publication_tests$' -V
```

The final TGZ was extracted separately under `archive-check/`. Its CPU launcher passed the real F006 dashboard prediction test and a multi-frame BMR7057 check: frames 0/50/99 × atoms 16/421/841, with finite isotropic and five tensor components, correct scene frame/atom identity, CPU execution, and exact repeated frame-0 values. Four helper launches `[0, 50, 99, 0]` match the application's deliberate single-frame inference cache. Mesa llvmpipe and clean shutdown were verified. Reports: `logs/archive-ml-smoke.log` and `review/installed-trajectory-ml/acceptance.json`.

## Shared changes for the Windows side

| File | Change and Windows review point |
| --- | --- |
| `CMakeLists.txt` | Linux required-ML option, Linux fixture/test registrations, Linux ring-current fixture environment, Linux deployment include. Existing Windows packaging and test branches remain. |
| `CMakePresets.json` | Minimum CMake metadata corrected to 3.24, matching the project; Linux description updated. Confirm Windows uses CMake >=3.24. |
| `src/app/ReaderMainWindow.cpp` | Platform-specific runtime files/helper names; native Linux CPU policy and executable checks. Windows CPU DLLs, ROCm DLLs, auto selection, and fallback must be rerun on Windows. |
| `tests/rest/test_experimental_shielding_ml_runtime.py` | Portable Linux file/device expectations and strict required-acceptance prerequisites. Windows DLL/ROCm checks retained; missing required assets now fail acceptance. |
| `tools/experimental_shielding_ml/CMakeLists.txt` | Linux ABI flags and relative helper runtime path; existing MSVC/ROCm compatibility block retained. |
| `cmake/BuildType-Debug.cmake` | Non-MSVC sanitizer link options propagate to consumers of the instrumented static core. MSVC branch unchanged. |
| `cmake/Platform-Linux.cmake` | Selected HDF5 include directories and updated Linux toolchain guidance. |
| `cmake/DeployLinux.cmake`, `cmake/InstallLinuxRuntime.cmake.in`, `cmake/LaunchLinux.sh.in`, `cmake/LaunchLinuxCpu.sh.in` | New Linux deployment/runtime handling and CPU launcher; deployment requires Qt >=6.10. |
| `scripts/linux/test_cpu.sh` | New Linux CPU acceptance runner. |

Source review found no intended Windows product-behavior change; native Windows verification remains pending. The new files must travel with the shared patch. No model retraining, numerical algorithm change, or Windows installer redesign is part of this port.

Run the following from the same revision's `h5-reader` directory in the existing MSVC 2022 developer PowerShell environment (`Enter-VsDevShell`), with the established Qt Pro, VTK, vcpkg, model bundle, and scientific fixtures present:

```powershell
$cmake = 'C:/Qt/Tools/CMake_64/bin/cmake.exe'
$ctest = 'C:/Qt/Tools/CMake_64/bin/ctest.exe'
$cpack = 'C:/Qt/Tools/CMake_64/bin/cpack.exe'
$qt = 'C:/Qt/6.10.2/msvc2022_64'
$vtk = 'C:/Projects/VTK'
$vcpkg = 'C:/Projects/vcpkg'
$model = 'C:/Projects/reader-data/experimental-shielding-ml-bundle-20260722-F006-R007-v1'
$oneP9J = 'C:/projects/reader-data/1p9j-calibration-with-dft/1p9j-july20260717.LGS'
$env:PATH = "$qt/bin;$vtk/bin;$vcpkg/installed/x64-windows/bin;$env:PATH"
Remove-Item Env:H5READER_EXPERIMENTAL_SHIELDING_ML_DEVICE -ErrorAction SilentlyContinue
Remove-Item Env:HIP_VISIBLE_DEVICES -ErrorAction SilentlyContinue

& $cmake --fresh --preset win-rwdi `
  "-DCMAKE_TOOLCHAIN_FILE=$vcpkg/scripts/buildsystems/vcpkg.cmake" `
  "-DH5READER_QT_DIR=$qt" "-DH5READER_VTK_DIR=$vtk" `
  "-DH5READER_EXPERIMENTAL_SHIELDING_ML_ROOT=$model" `
  -DH5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML=ON `
  -DH5READER_REQUIRE_EXPERIMENTAL_SHIELDING_ML_ROCM=ON `
  -DBUILD_TESTING=ON -DH5READER_ENABLE_MODEL_INPUT_CONTRACT_TEST=ON `
  "-DH5READER_MODEL_INPUT_REST_FIXTURE=$oneP9J"
if ($LASTEXITCODE -ne 0) { throw 'Configure failed' }
& $cmake --build --preset win-rwdi
if ($LASTEXITCODE -ne 0) { throw 'Build failed' }
& $ctest --test-dir build/win-rwdi --output-on-failure
if ($LASTEXITCODE -ne 0) { throw 'CTest failed' }
& $ctest --test-dir build/win-rwdi --output-on-failure `
  -R 'installed_ml_stale_dev_env|ml_cpu_fallback|ml_rocm'
if ($LASTEXITCODE -ne 0) { throw 'Windows native ML verification failed' }
& $cpack --config build/win-rwdi/CPackConfig.cmake -G INNOSETUP
if ($LASTEXITCODE -ne 0) { throw 'Installer build failed' }
```

After installing, repeat the scientific fixture, CPU and ROCm inference, screenshot, download/cache, and video checks against the installed executable. Use the Windows machine's actual installed executable path for `H5READER_BINARY`; do not count a build-tree run as installed-package validation. The existing `launch_reader.ps1 -Package` explicitly requests NSIS while current CMake selects Inno Setup. That pre-existing discrepancy is documented in the scope note; the Linux work does not change it.

## Completion record and remaining technical checks

- `FINAL_CTEST_TOTAL`: **31/31 passed, 0 failures, 55.05 seconds; `logs/ctest-final.log` and XML report.**
- `FINAL_PACKAGE_PATH`: `/home/jessica/opt/h5reader/dist/h5reader-0.5.0-Linux-OSS-x86_64-20260930.tar.gz` (293,415,825 bytes).
- `FINAL_PACKAGE_SHA256`: `b6ac294741e54e0263a12770e15f29fd8f298be701532a4cc934d95fccf3fa82`.
- `FINAL_RELOCATED_PACKAGE_ACCEPTANCE`: **Relocated application at `/home/jessica/Applications/h5-reader-linux-20260930`; installed REST 59 passed, successful stale-development ML 3 passed, separate missing-input 2 passed. Both standard and CPU launchers verified.**
- Native Windows build/install/CPU/ROCm acceptance: **Pending on Windows.**
- CUDA: **Not validated; user GPU remains reserved.**

The user clarified that this package is purely for their advisor. The active scope is an advisor handoff, with Windows parity and eventual CUDA support as the remaining technical follow-up. Public-release licensing is background documentation for a possible future distribution project, outside the active checklist. Existing upstream notices are retained.

This note and the accompanying scope note are tracked for the advisor handoff. Local toolchains, scientific fixtures, test artifacts and the binary archive remain outside Git; their paths and reproduction commands are recorded above.

The Linux implementation and these notes are checked in on `linux/oss-scope-20260930`. The original reviewable source patch remains at `/home/jessica/opt/h5reader/dist/h5-reader-linux-oss.patch`, with its file list and SHA-256 in the adjacent `h5-reader-linux-oss.json`. That patch includes the 12 implementation files and passed `git apply --cached --check` against the exact baseline; it predates the documentation check-in. The Windows side can review the branch or apply the patch after `git apply --check` from the repository root. Merge any newer overlapping Windows work explicitly.
