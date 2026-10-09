# macOS build and validation — 9 October 2026

## Scope and status

This extends the existing Linux and Windows Reader with a native Apple Silicon build using the user's Qt Pro installation. Source baseline: `da6f25bb5f766cd338a6d7bbcdeb7e42eb99a139` in `jessicalh/geo-kernel-extract`. This note accompanies the local Git checkpoint authorized by the user on 9 October; nothing is pushed or published. Read `LINUX_BUILD_AND_WINDOWS_HANDOFF_2026-09-30.md` alongside this note for the established scientific acceptance gates.

The follow-up sections are chronological. Their later binary hashes and validation records supersede earlier ones; statements about an empty index or uncommitted work describe the state at the time of those checks. The final checkpoint includes source, tests, and these notes, without generated builds, SDKs, fixtures, or model weights. See `MACOS_DISTRIBUTION_READINESS_2026-10-09.md` for the separate assessment of external installation.

**The native application is installed and tested in `outputs/h5reader.app`.** The bounded 32-test CTest run, hardware picking, original F006 scientific reference, relocated application checks, and main 1P9J REST run have passed. The main REST run reports 66 passed and 22 fixture/opt-in skips. Two separate BMR prediction-video diagnostics passed with atom tracking, including a source-derived residue filter; the original fixed-camera test remains distinct. Final installation passed dependency and signature checks, with matching debug symbols, and reproduced the F006 reference using its bundled CPU helper. A configured deployment target is not evidence of execution on that older OS. Current execution host: macOS 27, arm64, AppleClang 21.0.0.21000334. The selected minimum deployment target is macOS 14.4. This build does not produce an Intel or universal application.

## Separate toolchains and local configuration

The current workspace is `/Users/jessicahansberry/Documents/Codex/2026-10-08/t`. Application source is `h5-reader-project/h5-reader`; generated files, dependencies, fixtures, and evidence live outside it under `work/`.

| Component | Selected version/location |
| --- | --- |
| Qt Pro | `~/QtPro/6.12.0/macos`, with Widgets, OpenGLWidgets, Network, Multimedia/FFmpeg, Charts, HttpServer, WebSockets, Test |
| Build tools | `~/QtPro/Tools/CMake/CMake.app/Contents/bin`, `~/QtPro/Tools/Ninja/ninja` |
| VTK | 9.5.2, `work/deps/vtk-9.5.2-qt6.12.0` |
| HDF5 and libaec | 1.14.6 and 1.1.7, `work/deps/hdf5-1.14.6` |
| LibArchive and liblzma | 3.8.9 and XZ 5.8.4, `work/deps/libarchive-3.8.9` |
| Eigen | 5.0.1 headers, `/opt/homebrew/opt/eigen/share/eigen3/cmake` |
| Boost.Math | 1.89.0 standalone headers, `work/deps/boost-math-1.89.0/include` |
| LibTorch | Native macOS arm64 2.7.1, `work/deps/torch-2.7.1/libtorch` |
| Test Python | Python 3.13, `work/test-venv`; exact packages in `work/logs/test-venv-requirements.txt` |
| ML runtime input | `work/ml-runtime`, containing only the model, manifest, helper, four required dylibs, and notices |
| Build/current install | `work/build/reader-mac-qt612`, `outputs` |

The ignored `CMakeUserPresets.json` defines `mac-qt612`, inheriting the tracked `mac-rwdi` preset. It supplies this machine's Qt, VTK, HDF5, Eigen, Boost.Math, LibArchive, Python, ML, fixture, and `work/notices` paths. The generic REST default is now the original 1P9J fixture. The test preset runs one job at a time to avoid overlapping desktop GUI sessions, prepends `work/test-venv/bin` to `PATH` for FFmpeg, and excludes only the unavailable dedicated Chignolin and Trp-cage fixture tests. Keep those paths local. Qt versions and commercial/open-source installations are selections of an SDK; they do not require different application source trees. A changed Qt prefix requires its own matching VTK prefix and either a new application build directory or `cmake --fresh` to clear cached package locations. Linux's Qt/VTK pair and deployment module remain independent.

## Native prerequisites

Build each dependency independently, with its own source, build, and installation prefix. All native dependency configurations here use Ninja, Release, `CMAKE_OSX_ARCHITECTURES=arm64`, `CMAKE_OSX_DEPLOYMENT_TARGET=14.4`, and shared libraries. No build writes into QtPro or changes the global toolchain.

Exact executed configure/build/install commands and archive checksums are retained in `work/logs/native-dependency-build-commands.txt`; the exact VTK configuration is in `work/logs/vtk-configure-command.txt`.

- Build libaec first into the HDF5 prefix, then HDF5 with `libaec_DIR` pointing there. Enable Deflate and SZIP; disable HDF5 tools, examples, tests, and unused language bindings.
- Build XZ's library into the LibArchive prefix, then LibArchive against that prefix. Ignore `/opt/homebrew` and `/usr/local` for this native library build to avoid inheriting a newer deployment floor. Enable XZ, zlib, bzip2, and expat; omit unused optional libraries and archive CLI programs. Installed LibArchive depends on the isolated liblzma plus macOS system libraries.
- Build VTK against the selected `H5READER_QT_DIR`. Disable broad module groups, wrapping, examples, tests, and remote modules; enable the components listed by `cmake/DiscoverDeps.cmake`. VTK enables their transitive dependencies. The recorded command uses the same fifteen explicitly requested Reader components.
- Eigen and Boost.Math are header-only inputs. Boost.Math supplies complete elliptic integrals where Apple libc++ lacks the standard functions; keep the discovered standard-library implementation on platforms that already provide it.

HDF5's Deflate/SZIP round trips and relocated CMake consumer passed. Native LibArchive XZ/gzip/bzip2/ZIP round trips passed from a relocated prefix. Evidence: `hdf5-{smoke,relocation,consumer}.log` and `libarchive-relocation-smoke.log`. Installed dependency architectures and minimum-OS records were inspected; actual execution has only been checked on macOS 27.

Locally recorded SHA-256 values for the downloaded VTK/LibTorch archives and the pinned model manifest:

| Input | SHA-256 |
| --- | --- |
| `VTK-9.5.2.tar.gz` | `cee64b98d270ff7302daf1ef13458dff5d5ac1ecb45d47723835f7f7d562c989` |
| `libtorch-macos-arm64-2.7.1.zip` | `aa89ac85b91c83d0f976f8d135330d51e38ab777b26ec24f312fd58d079314cb` |
| `manifest.json` | `5be2f99f2fb89cc0b56860b7783412eb72dba9bd243bb213f2b361c2c05e6886` |

The libaec/XZ/LibArchive archive hashes in the command log were compared with the installed Homebrew formula records. A locally recorded hash alone is not an upstream digest comparison.

## Native CPU helper and unchanged model

The [public source documentation](https://semantic.construction/thesis/software-design/sources.html) identifies the [source release at revision 8714b8ca6e859236f46d17e16c79dfe896780ff8](https://github.com/jessicalh/protein-extract-and-view-tools/tree/8714b8ca6e859236f46d17e16c79dfe896780ff8). Its `h5-reader/models/experimental_shielding_ml` directory supplies the original `model.ts` and `manifest.json`. The 7,376,595-byte model matches the existing Linux acceptance hash:

`101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6`

This is `F006-R007-v1-reader-runtime`, inference schema 3 and helper protocol 2. Do not retrain or regenerate it for the Mac port. Obtain the native [LibTorch 2.7.1 arm64 archive](https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-2.7.1.zip) separately. Python is used for developer tests; the shipped helper is C++ and does not start Python.

The following stages the current four-library CPU closure from an unpacked LibTorch SDK. Set the locations for the machine; the example matches this workspace. Run this block in Bash with errors stopping execution. It changes copied runtime files only.

```bash
set -e
reader_workspace="$HOME/Documents/Codex/2026-10-08/t"
reader_source="$reader_workspace/h5-reader-project/h5-reader"
reader_work="$reader_workspace/work"
reader_cmake="$HOME/QtPro/Tools/CMake/CMake.app/Contents/bin/cmake"
reader_ninja="$HOME/QtPro/Tools/Ninja/ninja"
reader_torch="$reader_work/deps/torch-2.7.1/libtorch"
reader_ml="$reader_work/ml-runtime"

"$reader_cmake" -S "$reader_source/tools/experimental_shielding_ml" \
  -B "$reader_work/deps/infer-build" -G Ninja \
  -DCMAKE_MAKE_PROGRAM="$reader_ninja" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=14.4 \
  -DCMAKE_PREFIX_PATH="$reader_torch"
"$reader_cmake" --build "$reader_work/deps/infer-build" --parallel 2
mkdir -p "$reader_ml/lib"
install -m 755 "$reader_work/deps/infer-build/infer" "$reader_ml/infer"
for reader_library in libc10 libtorch libtorch_cpu libomp; do
  install -m 755 "$reader_torch/lib/$reader_library.dylib" "$reader_ml/lib/"
done

# This exact upstream package refers to Homebrew libomp. Make the copied
# runtime resolve its own library; never patch the LibTorch or Qt SDK.
install_name_tool -change /opt/homebrew/opt/libomp/lib/libomp.dylib \
  '@rpath/libomp.dylib' "$reader_ml/lib/libtorch_cpu.dylib"
install_name_tool -id '@rpath/libomp.dylib' "$reader_ml/lib/libomp.dylib"
for reader_binary in "$reader_ml/infer" "$reader_ml"/lib/*.dylib; do
  codesign --force --sign - "$reader_binary"
  codesign --verify "$reader_binary"
done

reader_release=https://raw.githubusercontent.com/jessicalh/protein-extract-and-view-tools/8714b8ca6e859236f46d17e16c79dfe896780ff8/h5-reader/models/experimental_shielding_ml
curl --fail --location "$reader_release/model.ts" -o "$reader_ml/model.ts"
curl --fail --location "$reader_release/manifest.json" -o "$reader_ml/manifest.json"
reader_model_hash=$(shasum -a 256 "$reader_ml/model.ts")
test "${reader_model_hash%% *}" = \
  101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6
```

Keep upstream notices with the staged runtime under `licenses/` and supply additional prepared notices through `H5READER_THIRD_PARTY_NOTICES_DIR`. The helper's development runtime path is `@loader_path/lib`. The macOS staging module adds the bundle-relative Frameworks path to its copied helper and re-signs that copy. `auto` and `cpu` both use CPU; ROCm is unavailable on macOS. MPS acceleration is not implemented or claimed.

## Application build, tests, and installation

The existing ignored local preset reproduces this machine's configuration. Use a new build directory when changing architectures or dependency versions. Commands below reuse the variables above and Qt's installed tools:

```bash
cd "$reader_source"
export PATH="$reader_work/test-venv/bin:$HOME/QtPro/Tools/CMake/CMake.app/Contents/bin:$PATH"
cmake --preset mac-qt612
cmake --build --preset mac-qt612
ctest --preset mac-qt612
cmake --install "$reader_work/build/reader-mac-qt612" --component Runtime
cpack --config "$reader_work/build/reader-mac-qt612/CPackConfig.cmake" \
  -G TGZ -B "$reader_work/dist"
```

Run GUI/VTK tests in the logged-in macOS desktop session. Do not set `QT_QPA_PLATFORM=offscreen`: the embedded OpenGL widget needs a real context. Xvfb wrapping is Linux-only. Avoid inherited Qt/plugin/library overrides; select the intended SDK through CMake. The tests' FFmpeg command comes from the local Python environment, not an application runtime dependency.

The build installs `h5reader.app`, with models under `Contents/Resources/ml/experimental_shielding_ml`, helper under `Contents/Helpers/infer`, and shared libraries under `Contents/Frameworks`. Optional datasets use `Contents/Resources/datasets`. `DeployDarwin.cmake` collects native dependencies from the build executables while their development RPATHs can still resolve them, then invokes the selected SDK's Qt deployment API and checks the resulting dependency closure. It fails if a non-system dependency is unresolved or still resolves outside the bundle. Linux's ELF deployment is unchanged. The current local preset installs to `outputs`; the first successful install used `work/install/reader-mac-qt612`.

Preserve a matching `.dSYM` before removing build objects. For this RelWithDebInfo build, run `dsymutil` against `work/build/reader-mac-qt612/h5reader.app/Contents/MacOS/h5reader` and retain its output beside the archived package. `outputs/h5reader.app.dSYM` matches the installed executable's arm64 UUID `5632612A-B617-3821-BA76-57CC30FABDCA` (`work/logs/reader-final-uuids.log`); regenerate it if the executable is rebuilt. The deployment tool retains symbols, but an unstripped executable with a debug map is not a replacement for the `.dSYM`.

Extract the archive into a separate directory outside the source/build/install trees before acceptance. Launch the extracted executable directly for REST tests, or open the extracted `.app` normally. Successful launch from the build tree is not relocated-package validation. No Developer ID signing, notarization, public release, or older-macOS execution is claimed.

## Fixtures and validation record

The downloaded public BMR7057 trajectory contains 842 atoms and 100 frames, including reorientational dynamics and all manifest-declared numerical/mask source files. It has no attached DFT. Archive size: 382,396,716 bytes; expanded payload: 1,040,325,393 bytes. Local SHA-256: `5f43cf991769f437ced5b141034237e16a474154cfdc5cfa2ccff64911a0a37d`. Sizes match the HTTPS catalog and XZ integrity passed; the publisher supplied no digest. Provenance and extraction checks are in `work/fixtures/bmr7057-verification.json`.

The current preset's generic REST fixture is `work/fixtures/1p9j/1p9j-july20260717.LGS`; the separate F006 setting selects the original single-pose fixture below. BMR7057 remains available for no-ORCA prediction/video checks at `work/fixtures/bmr7057-20260919-files/run.LGS`. For actual BMR prediction acceptance, use `H5READER_EXPECT_AUTO_ML=1` with `test_selected_atom_predicts_without_orca`. Do not set `H5READER_EXPECT_ML_SUCCESS=1` for BMR7057: that flag's tests pin the separate 528-atom F006 fixture and its atom-16 reference value. Do not run DFT assertions against BMR7057.

The original Batcave acceptance fixtures are now available locally. They were fetched read-only through the existing `remoteadmin` account using `/Users/remoteadmin/.ssh/config_main_mesh` and the remote account `jessica@batcave`. Every transferred scientific payload was compared with its remote SHA-256. Original `source.LGS` wrappers are retained; only paths in the local working wrappers were rewritten, with optional raw-MD provenance omitted where unavailable.

- F006: `work/fixtures/f006/full720-f006-A0A822IKI2.LGS`, 528 atoms, the complete WT pose root: 289 files and 26,568,853 bytes. The ALA mutation directory was not needed or transferred. Evidence: `work/fixtures/f006/transfer-verification.json` and `complete-wt-inventory.json`.
- 1P9J: `work/fixtures/1p9j/1p9j-july20260717.LGS`, the unchanged 846-atom, 751-frame HDF5 trajectory, static topology, and all 751 original DFT metadata/primary-output pairs: 1,510 files and 2,906,247,615 bytes. Evidence: `work/fixtures/1p9j/transfer-verification.json`.
- 1P9J per-frame NPY snapshots are deliberately sparse: rows 0–30, 269, 270, 375, and 750, totaling 8,610 files and 373,574,676 bytes. The corresponding original indices are 0–60 in steps of two, 538, 540, 750, and 1500. Original indices match exactly and loaded NPY positions match HDF5 positions with maximum error 0. Optional missing files were also absent on Batcave. Evidence: `work/fixtures/1p9j/selected-npy-verification.json`. No scientific arrays were synthesized or edited. These snapshots support bounded checks; they do not establish all-frame metric or model-input-export coverage.

| Gate | Verified result or remaining limit |
| --- | --- |
| Native dependencies | HDF5 compression/relocation and LibArchive archive-format/relocation checks passed; logs above |
| Helper protocol | Seven synthetic native/Python cases passed, including malformed input/output and unavailable ROCm rejection; `work/logs/native-infer-protocol.json` |
| Original trained model, synthetic inputs | Python/native Torch 2.7.1 parity passed with maximum absolute error 0; `work/logs/f006-synthetic-runtime.json`; this is not scientific fixture acceptance |
| Reader bounded CTest run | 32/32 passed, including native tests, the 1P9J ring-current scientific case, F006 CPU inference, and installed-runtime fallback with a stale development path. Excludes `h5reader_rest_smoke`, `h5reader_rest_chignolin_fixture`, and `h5reader_rest_trpcage_fixture`; `work/logs/ctest-final.log` and `.xml` |
| Native hardware picking | 9/9 QtTest cases passed; `work/logs/hardware-picking.log` |
| Original F006 scientific reference, relocated app | 3 passed, 2 not-applicable skips; `work/logs/relocated-f006.log` and `.xml`. Tests confirm the original atom-16 reference of 29.76103401184082 ppm within 0.001 ppm, CPU execution, dashboard output, and tensor delivery to scene/inspector. Skips are the no-ORCA trajectory case and deliberately missing-input case; this complete fixture satisfies neither condition |
| Original 1P9J+DFT ring-current science | `h5reader_ring_current_case_tests` passed against the complete original 751-frame trajectory and DFT; included in `ctest-final`. Sparse NPY snapshots remain a limit for other checks |
| Earlier broad BMR7057 UI/REST run | 51 passed, 15 failed, 18 skipped, 3 deselected; `work/logs/bmr-rest.log` and `.xml`. Many failures require DFT or different fixture fields. This historical run is not a passing full REST result; the main suite subsequently used its original 1P9J fixture |
| Main 1P9J REST run | 66 passed, 22 skipped, 4 Pillow deprecation warnings, no failures, in 51.26 s; `work/logs/1p9j-rest.log` and `.xml`. Includes generic video export/decoding, requested frame range, partial-file finalization, and frame-interference recovery. Two stale test contracts were corrected before this run: canonical overlay names and optional raw-MD provenance paths; no production video code was changed |
| BMR predicted-video tracking diagnostic | 1 passed in 27.83 s; `work/logs/video-tracking.log` and `.xml`, with captures in `work/acceptance/video-tracking`. The diagnostic adds atom tracking and reuses the original test's pixel assertions unchanged: four predicted video frames match independently prepared references and each reference contains more than 100 tensor-colored pixels. This is distinct from passing the original fixed-camera test |
| BMR predicted video with source-derived residue filter | 1 passed in 29.33 s; `work/logs/video-tracking-source-residue.log` and `.xml`. Uses atom tracking and the selected atom's actual residue index from the source HDF5, while retaining the original pixel assertions. The captured view includes the real CG atom and its bonds; tensor arrows may be occluded by geometry |
| Relocated installed app | The first installed app was physically moved from `work/install/reader-mac-qt612/h5reader.app` to `outputs/h5reader.app`; the relocated F006 run used `outputs/h5reader.app/Contents/MacOS/h5reader`. Native launch and CPU inference passed after relocation; initial dependency closure is in `work/logs/reader-install-retry.log` |
| Final installation and signature | The latest build was installed to `outputs/h5reader.app`, including prepared notices. Full dependency closure passed (`work/logs/reader-final-install.log`), deep strict ad-hoc signature verification passed (`reader-final-codesign.log`), and the `.dSYM` UUID matches (`reader-final-uuids.log`) |
| Final installed smoke and scientific recheck | With inherited `H5READER*` environment overrides removed, the installed app used its bundled CPU helper and the unchanged original F006 scalar assertion passed. Captured state records exactly 29.76103401184082 ppm and a visible selected tensor. It then loaded the complete 751-frame 1P9J trajectory at frame 0 with no selection or filter. Evidence: `work/logs/final-installed-smoke.json`, `final-installed-ml-state.json`, retained `reader-live.log`, and `work/acceptance/visual/final-reader.png` |
| Diagnostic VTK/rendering run | VTK OpenGL diagnostic rebuild and install completed; `work/logs/vtk-gl-diagnostics-{configure,build,install}.log`. The successful main 1P9J REST and BMR tracking-video runs used this build. The REST harness deletes its application log on success, so the 1P9J result is not evidence of zero OpenGL error messages |
| Original-versus-published trajectory equality | Remote-to-local SHA-256 comparisons establish integrity of the fetched originals; a separate original-versus-published trajectory comparison is not recorded |
| Linux/Windows regression | Shared branches preserved in source; fresh native builds on those platforms remain separate validation |

The main 1P9J run's 22 skips are explicit: six Chignolin/Trp-cage fixture cases; two F006-specific reference cases; the no-ORCA prediction case; the deliberately missing-input case; two model-input-export opt-ins; two responsiveness/shutdown opt-ins; one combined DFT/prediction-pin export case; and seven no-ORCA, failing-helper, or separate-process tensor-video cases. Their exact reasons are retained in `1p9j-rest.log` and `.xml`. A skip does not count as validation of that behavior.

The original fixed-camera BMR video test failed because its reference had no visible tensor pixels. The tracking diagnostic demonstrates successful prediction/capture with usable framing. Source-data inspection also found that the selected BMR atom 16 belongs to residue index 1, while the original test filters residue index 0, hiding that atom but leaving its selection marker. The follow-up diagnostic using the actual HDF5 residue index also passed. Neither diagnostic changes production video code or establishes that the original fixed-camera test passed.

The initial validated installation's binary SHA-256 values are recorded in `work/logs/reader-final-binary-sha256.txt`. The subsequent trackpad change below supersedes the main executable:

| Installed file | SHA-256 |
| --- | --- |
| `Contents/MacOS/h5reader` | `02d1f40989a74247a067dab6531062ba2d80426cdb0038d8ff85d7b667b84e58` |
| `Contents/Helpers/infer` | `7dc919b8d9275ae3012be4871452509eca6f57f63f849c1655a2c71b551b1250` |
| `Contents/Resources/ml/experimental_shielding_ml/model.ts` | `101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6` |

The tensor-text report was checked against `work/acceptance/visual/tensor-window.png`, `tensor-tree.json`, and source. No Unicode conversion defect was found in those captures: the em dash and angstrom symbol render correctly, and data retains Greek letters and superscripts. The user subsequently identified a mark beside GLU. The reproduced panel (`glu-label-window.png` and `glu-label-tree.json`) contains a literal ASCII `GLU #2` in Atom Info and `A:GLU2:CG` in the tensor label; whether the reported mark was that hash is still unconfirmed. `sigma_iso` and related labels are also literal source text. The clear-selection button uses Qt's `SP_DialogCloseButton` artwork. No encoding workaround or platform font substitution was made.

The retained final `reader-live.log` contained no matches for the checked OpenGL/shader/framebuffer error patterns at the completion of the installed smoke check. This bounded log check does not establish the absence of every possible graphics issue. The passing original F006 reference test is distinct from synthetic parity. The main REST result covers its 66 executed tests, with the fixture and opt-in limits recorded above; fresh Windows and Linux builds and full cross-platform parity are not claimed. The original Git HEAD is unchanged, the index contains no staged changes, and the source diff whitespace check passed; nothing was committed.

## Trackpad zoom follow-up

The reader's camera event filter previously handled mouse dragging and wheel events only. Qt/VTK already registered pinch recognition, but the reader's deliberately quiet VTK interactor did not move the camera for those gestures. The first native-event regression run reproduced the missing zoom in free, atom-following, and plane-locked modes (`work/logs/pinch-before.log`).

`CameraInputFilter` now routes Qt pinch gestures and direct native zoom events through the existing `CameraComposer` dolly path. Incremental zoom is derived from the gesture's total scale, so interleaved rotation events and the final event do not repeat the last zoom. Direct native zoom also covers Qt platforms that deliver that event without a `QPinchGesture` recognizer. No OS-specific code or SDK changes were needed. `MoleculeScene::refreshCameraForInput` recomposes the current camera while paused without reloading positions or rebuilding overlays; mouse and wheel input use that same refresh path.

Focused verification: **12 QtTest cases passed, none failed or skipped**, including the three native pinch modes, repeated and reversed gestures, unchanged focal point and orientation, no extra zoom at gesture end or the next camera recomposition, wheel compatibility, direct native fallback, and the existing Retina/filtered hardware-picking cases. Evidence: `work/logs/pinch-and-picking.log`. These are injected Qt native events through the actual Qt recognizer and camera. The user subsequently confirmed that physical trackpad zoom works. Earlier broad scientific/REST results above predate this focused camera change and were not rerun unnecessarily.

The updated app is installed at the same `outputs/h5reader.app` path. Dependency closure and deep strict ad-hoc signature checks passed (`pinch-install.log`, `pinch-codesign.log`). Matching executable/dSYM UUID: `0020FA95-802F-3581-B8AE-E44481BD6375` (`pinch-uuids.log`). Main executable SHA-256: `27f62ceb4e3f89652b7d21f4f7043a4657bcf124008aac618107e0f8530cdc0e` (`pinch-binary-sha256.txt`). The installed app launched and restored the user's BMR7057 trajectory, frame, and selected atom; `pinch-installed-smoke.json`. `reader-live.log` is the current session log and changes on relaunch. All changes remain uncommitted.

## macOS toolbar text follow-up

The user reported white text on the checked toolbar actions. This is reproducible with an otherwise unmodified Qt 6.12 `QToolBar`/`QAction`, not a Reader-wide palette override: the forced light/Fusion palette in `main_reader.cpp` is Windows-only. The installed Qt source, `qtbase/src/plugins/styles/mac/qmacstyle_mac.mm`, paints checked toolbar backgrounds with a translucent neutral fill but uses the selected-menu `HighlightedText` role for text-only pressed/checked labels. In light appearance that produces white text on pale gray. Its disabled checked labels also fail to look consistently disabled in the comparison.

`MacWidgetStyle.h` (originally `MacToolbarStyle.h`, renamed during the checkbox follow-up below) applies a separate, toolbar-owned `QProxyStyle` to the Tools toolbar's existing text buttons, only when the active style is `macos`. It keeps QMacStyle's button/background/metrics and delegates text labels to Qt's common renderer using the button palette. There are no literal replacement colors, global palette changes, new control semantics, or installer-managed SDK edits. The existing `QAction` checked/enabled/triggered behavior remains intact. Linux/Windows and non-Mac styles do not enter this workaround. Recheck it on a Qt upgrade.

This follows [Qt's supported style-wrapper approach](https://doc.qt.io/qt-6/qproxystyle.html). [Apple's button guidance](https://developer.apple.com/design/human-interface-guidelines/buttons) permits short text labels when they communicate the action better than an icon. A broader toolbar redesign is separate; the current domain-specific labels remain useful. Do not switch this VTK window to Qt's unified title/toolbar option: [Qt documents that option as unsupported with QOpenGLWidget](https://doc.qt.io/qt-6/qmainwindow.html#unifiedTitleAndToolBarOnMac-prop).

The small standalone comparison (`work/toolbar-style-probe.cpp`) exercises the actual shared header with the installed Qt. Light/dark and checked/unchecked/pressed/disabled states were visually inspected, and two button clicks verified QAction synchronization and exactly two triggers. Evidence: `work/logs/toolbar-style-probe.log`, `work/acceptance/visual/toolbar-probe-{light,dark}.png`. It changes the probe application's appearance only, not system appearance. The Reader rebuilt, dependency closure and deep strict bundle signature passed, and executable/dSYM UUIDs match: `B893255C-DF6D-35EF-94CE-BFF349F85934`. Main executable SHA-256: `5bf0f5fb434fa1786be7d95b90862f3612e75f84b2867c48dd6434031f835ea9`. Evidence: `toolbar-build.log`, `toolbar-install.log`, `toolbar-codesign.log`, `toolbar-uuids.log`, and `toolbar-binary-sha256.txt` under `work/logs`. This supersedes the pinch-only executable above; the broader scientific tests were not rerun for this rendering-only change.

The installed Reader was reopened on BMR7057 at frame 0 with atom 649, its toolbar states, and its exact camera position/focal point/view-up restored (within `1e-10`). The ribbon action's checked state tracked an off/on round trip. The final application screenshot shows readable checked Ribbon/Rings labels, and the local source diff whitespace check passed with an empty index and unchanged HEAD. Evidence: `work/logs/toolbar-before-restart-state.json`, `toolbar-installed-smoke.json`, and `work/acceptance/visual/toolbar-after.png`. All changes remain uncommitted.

## Native tensor checkbox follow-up

The tensor panel's standard `QTreeWidget` check states and click handling were present, but the native checkbox indicators did not paint. An isolated Qt 6.12 comparison reproduced the issue with ordinary item-view checkboxes while standalone `QCheckBox` controls painted correctly. The installed QMacStyle source passes a nonzero item-view rectangle to AppKit as both the native view frame and drawing bounds. Translating the painter to the indicator and drawing that same native primitive at `(0, 0)` restores it on this macOS 27 host.

The existing small Mac style wrapper is now named `MacWidgetStyle.h` and also supplies this coordinate correction to the tensor tree only. It retains native Mac checkbox drawing, sizing, palette, clipping, checked/mixed/disabled states, and the standard delegate's hit rectangle and mouse/keyboard behavior. No replacement widget, manual checkmark artwork, alternate toggle API, application-wide style override, or SDK modification is involved. The existing model/view check-state and tensor show/hide connections are unchanged. Qt's [item-view style documentation](https://doc.qt.io/qt-6/style-reference-menus-views.html) describes the standard checkbox primitive used here.

Verification used the production header in `work/checkbox-style-probe.cpp`. Native checked, unchecked, mixed, disabled, and selected-row indicators were inspected in light and dark appearance. Mouse clicks in the style-reported checkbox rectangles and keyboard Space round trips produced exactly four model changes; clicks on disabled indicators produced none. Evidence: `work/logs/checkbox-style-probe.log` and `work/acceptance/visual/checkbox-probe-{light,light-selected,dark}.png`. The app rebuilt successfully, the deployment dependency audit and deep strict signature check passed, and executable/dSYM UUIDs match: `431F9BCB-7BE6-346A-A61E-EF2014578281`. Main executable SHA-256: `a8d38630ec0627a384d656d9d0896e0fd9babccde97728e177cee05f6d3d0cfe`. Build/install/signature/symbol evidence is in `work/logs/checkbox-{build,install,codesign,uuids}.log` and `checkbox-binary-sha256.txt`. The earlier broad scientific/REST tests were not repeated for this drawing-only correction.

The installed application's final capture (`work/acceptance/visual/tensor-checkbox-after.png`) shows both native blue checkboxes beside the Shielding and Bond orientation headings. Reopening restored BMR7057, atom 239, frame 8 (UI frame 9), camera pose, and tensor checked states. The residue filter was cleared at the user's explicit request. Evidence: `work/logs/checkbox-before-restart-state.json` and `checkbox-installed-smoke.json`. The source whitespace check passed and the index remains empty; nothing was committed.

## Deployment references

- [Qt deployment scripts](https://doc.qt.io/qt-6/qt-generate-deploy-script.html)
- [Qt runtime dependency deployment](https://doc.qt.io/qt-6/qt-deploy-runtime-dependencies.html)
- [Qt for macOS deployment and bundle layout](https://doc.qt.io/qt-6/macos-deployment.html)
- [Qt pinch gestures](https://doc.qt.io/qt-6/qpinchgesture.html)
- [Qt native gesture events and zoom increments](https://doc.qt.io/qt-6/qnativegestureevent.html)
