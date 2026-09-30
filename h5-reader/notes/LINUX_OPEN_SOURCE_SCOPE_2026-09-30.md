# Linux open source target and Windows preservation

## Objective and working baseline

Add a complete Linux target for the current Reader while preserving the deliberately developed Windows Qt Pro target. This is an extension of the current application. Keep molecular identities, calculations, model weights, manifests, interfaces and visualization behavior consistent across platforms.

The current development baseline is `jessicalh/geo-kernel-extract`, commit `3c0303a32ade34eafcd6341591be337a811f1007`, 30 September 2026 at 12:42 BST. Its message describes measurement geometry and visibility, safe reload, scrub-gap recovery, inspector navigation and CSA rotation fixes. It reports 176 sampled trajectory checks and 32 CTest targets passing on the Windows side, including a separate Trp-cage rerun. These are historical reports; this Linux investigation has not reproduced them yet.

The isolated working checkout is `/home/jessica/src/h5-reader-linux`, branch `linux/oss-scope-20260930`. The original `/shared/2026Thesis/nmr-shielding` checkout is a different producer-oriented state. The source release `/shared/2026Thesis/nmr-shielding-release` is commit `8714b8c` from 10 September and predates the trajectory library and today's fixes. The remote `windows-build` branch is older still, at `9cf5b53` from 30 August. Existing desktop launchers use a May binary.

User constraint: this Linux machine uses CPU or CUDA. Its GPU currently hosts an active large-model workload. Use CPU inference, software rendering and software encoding until the user explicitly says that workload has finished. Preserve Windows ROCm functionality. CUDA integration and validation remain a later hardware-specific work item.

## Evidence from the initial investigation

| Area | Finding | Consequence |
| --- | --- | --- |
| Qt | Current source requires Qt >=6.8 and Widgets, OpenGLWidgets, Network, Multimedia, Charts and HttpServer; requires QFFmpegMediaPlugin. | Retain features and install a modern OSS Qt. |
| Baseline configure | GCC 13.3 configures, then CMake rejects this host's Qt 6.4.2 at DiscoverDeps.cmake:31. | No application compilation or runtime success established by this probe. |
| Video | SceneVideoExporter uses QVideoFrameInput (introduced in Qt 6.8, FFmpeg only), H.264 or MPEG-4 in MP4. | Require a real encoder and successful decoded video, not merely the plugin target. |
| VTK | Local VTK 9.5.2 was built with Qt 6.4.2. | Build a dedicated VTK against the chosen Qt for reproducibility; old VTK is not proven ABI-incompatible. |
| Archive library | ScienceFilesQt requires LibArchive and uses C++20 privately. Host development package absent. | Provide headers and library with native tar/ZIP/xz support. |
| ML discovery | ReaderMainWindow requires infer.exe and Windows CPU/ROCm DLL lists on every platform, including development overrides. | Extend native runtime policy while keeping Windows rules. |
| ML package | Root CMake packages the helper/model/libraries only under WIN32. | Add Linux native helper, weights, manifest and shared-library deployment. |
| Runtime environment | Inherited LD_LIBRARY_PATH includes /opt/orca/lib and CUDA paths; ldd selected ORCA's libstdc++ and libgcc. | Use command-local clean build/test environments; no global shell changes. |
| Installed package | Linux currently installs executable/catalog/examples and makes a TGZ. | Add Qt, VTK, FFmpeg, TLS/platform/image plugins, native inference and relocatable library paths. |
| Tests | Windows-only ML registrations and hard-coded Windows runtime checks can skip Linux inference. Several other tests deliberately skip without fixtures. | Explicitly require successful CPU inference and record feature coverage and skips. |

The current catalog has 22 published bundles: 14,127,741,024 compressed bytes and 32,598,126,754 expanded bytes. Begin with existing local fixtures and one real download; do not download the full catalog as a prerequisite to compiling.

## Proposed Linux stack

Use Ubuntu 24.04 x86_64, GCC 13.3, CMake 3.28.3 and Ninja 1.11.1 already on this machine. Install public Qt 6.10.2 linux_gcc_64 with Charts, HttpServer, Multimedia and WebSockets into `/home/jessica/opt/h5reader/Qt`; build VTK 9.5.2 into a separate prefix under the same toolchain directory. Use the existing HDF5 1.10.10 and Eigen 3.4 initially, with consistent selected HDF5 headers and libraries. Provide LibArchive locally. Build native inference against a separate CPU-only LibTorch installation.

Public Qt 6.10.2 packages were confirmed in the official repository. Ubuntu 24.04/GCC13 is supported; the exact downloaded x86_64 archive names identify a RHEL build environment. These are distinct facts. Keep Qt's matching FFmpeg shared libraries: system FFmpeg 6.1 cannot automatically substitute for the Qt 6.10 package's 7.x library ABI.

Initial native inference layout:

```
ml/experimental_shielding_ml/
    infer
    model.ts
    manifest.json
    lib/libc10.so
    lib/libtorch.so
    lib/libtorch_cpu.so
    lib/libtorch_global_deps.so
    lib/<required transitive shared libraries>
```

The helper uses a relative runtime path to its lib directory. The Reader retains its existing process boundary and binary input/output protocol. CPU inference must work from the installed package with development environment variables absent or stale.

## Available model and fixture evidence

The source release contains `h5-reader/models/experimental_shielding_ml/model.ts` (7,376,595 bytes) and `manifest.json`. The model SHA256 was independently checked against the manifest:

`101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6`

The manifest is F006-R007-v1-reader-runtime, inference schema 3, helper protocol 2. It specifies native runtime without Python, six float32 output columns for isotropic shielding and rank-2 tensor components; rank 1 is neither predicted nor synthesized. The recorded export error is historical evidence, not a Linux numerical test. No native helper accompanies the located assets.

Existing potential fixtures include `/shared/2026Thesis/reader_staging/oneoffs/1p9j/extraction/run_20260717T222625Z.lgs` and `/shared/2026Thesis/reader_staging/trajectories/bmr4976/extraction/calcset_20260719T105114Z.lgs`. Their existence is verified; full compatibility and scientific completeness need testing. Chignolin and Trp-cage extraction data also exist, but current usable manifests with DFT references still need locating or constructing as local fixture wrappers. Keep source extraction files unchanged.

## Work sequence and Windows overlap

1. Establish the isolated OSS toolchain and compile the unchanged Reader as far as possible. Record concrete compiler issues rather than guessing at porting needs.
2. Extend native runtime discovery and the helper build. Windows `.exe`, DLL lists, CPU fallback, explicit ROCm behavior, status fields and bundled model protocol must remain intact.
3. Add Linux build/runtime packaging and strict CPU acceptance jobs. Preserve Windows packaging defaults, compiler flags and deployment script.
4. Run scientific, UI, REST, ML, publication and multimedia validation on real fixtures. Fix platform issues in the narrowest applicable layer.
5. Install and relocate the package; test it from a clean environment outside its source/build directories.
6. Supply the concrete shared patch and Windows verification commands to the Windows side. A fresh Windows build/test/package result is required to claim Windows preservation verified.

Confirmed shared overlap was reported to the user before implementation: `h5-reader/CMakeLists.txt`, `src/app/ReaderMainWindow.cpp`, `tests/rest/test_experimental_shielding_ml_runtime.py`, and the standalone ML helper CMake project. `ExperimentalShieldingMlStore.cpp` only needs touching if process/device handling requires it. Keep application calculations and visualization algorithms outside the porting changes unless an observed platform failure requires a focused fix.

Windows invariants include current Qt Pro discovery and user overrides, MSVC settings, VTK initialization, Debug console behavior, AVX2 and RWDI IPO, Dbghelp, CPU and ROCm required assets, stale development environment fallback, CTest DLL paths, windeployqt, vcpkg/MSVC dependency collection, and Inno Setup. The existing PowerShell launcher's explicit NSIS packaging argument differs from current CMake's Inno Setup setting; report this existing Windows issue separately rather than folding it into Linux work without coordination.

## Acceptance matrix

| Feature | Required evidence |
| --- | --- |
| Build | Clean Linux build with all required modules; no feature-disabling fallback. Relevant Debug/sanitizer linkage checked separately. |
| Data model | Existing loader, NPY, manifests, atom identities, frame map, field availability, physics, alignment and CSA tests pass. |
| Current UI fixes | Measurement geometry/visibility, reload, scrub gaps, frame detail, inspector navigation and rotated CSA tests execute successfully. |
| Scientific displays | Known ORCA tensors, ring-current displays, charts, selections, cameras, overlays and screenshots verified on actual fixtures. |
| Native ML | Installed CPU helper produces finite expected isotropic and tensor output; values reach chart, inspector and scene. Missing required inputs are diagnosed. Stale development settings recover through the installed runtime. |
| Learned activity | Load, radius, display rotations, sparse/missing frames, rejected input, reload and video interaction tests pass. Optional capture-generation Python tool tested separately. |
| Trajectory library | Catalog, HTTPS transfer, native xz extraction, open, cancel, cleanup, cache locations and reload protection work. |
| Publication | Exact all-frame source-versus-portable comparison configured on a known pair; corpus checker uses explicit expected count and requires predictions. |
| Model export | Enable fixed model-input export and shutdown contracts with a complete fixture. |
| Video | Complete and decode MP4 with exact requested, changing frames; stop produces finalized partial video; playback frame and UI state recover. Use software encoding while GPU is occupied. |
| Package | Relocated installed app opens data, predicts, renders in software, serves REST, downloads and exports screenshots/video with no build-tree or developer-shell dependency. |
| Windows | Windows side builds the same shared revision, runs existing CTest plus CPU/ROCm inference checks, builds installer and checks installed app. |

Full CTest runs need an outer Xvfb/software-rendering environment for C++ GUI tests too. Set actual `.LGS` fixture files rather than directory defaults; learned-activity tests open the fixture as JSON. Use a different reload fixture. Supply ffmpeg for video decode verification. Configure optional publication and model-input tests explicitly. A green metric sweep can include unavailable ML descriptors, so it cannot replace successful inference tests.

The corpus checker samples first/middle/last frames and a small tensor atom subset. Its reported coverage must remain explicit; it is not exhaustive scientific verification of every frame. Every required conditional test should run in an appropriate matrix job; record and explain remaining skips.

## Advisor handoff scope

The user clarified that the package is purely for their advisor. Public-release preparation is outside the current task; the licensing investigation below is retained as background for any later public distribution. It does not hold up the advisor handoff.

## Open source distribution background

The required Qt modules are available as open source. Charts and HttpServer offer GPLv3 rather than an LGPL-only route. Retaining those modules avoids a UI/server rewrite. This repository lacks a tracked project license and has incomplete local licensing evidence for the embedded science-files code and model assets. Record the owner's distribution/license decision before publishing an OSS release; do not silently assign a license or publish the repository. Local implementation and validation can proceed while this decision is resolved.

## Primary external references

- [Qt 6.10 supported platforms](https://doc.qt.io/qt-6.10/supported-platforms.html)
- [Official public Qt 6.10.2 Linux package repository](https://download.qt.io/online/qtsdkrepository/linux_x64/desktop/qt6_6102/qt6_6102/)
- [QVideoFrameInput API and FFmpeg requirement](https://doc.qt.io/qt-6.10/qvideoframeinput.html)
- [Qt Multimedia backends and deployment](https://doc.qt.io/qt-6.10/qtmultimedia-index.html)
- [Qt Charts licensing](https://doc.qt.io/qt-6.10/qtcharts-index.html#licenses)
- [Qt HTTP Server licensing](https://doc.qt.io/qt-6.10/qthttpserver-index.html#licenses)
- [Qt runtime deployment API](https://doc.qt.io/qt-6.10/qt-deploy-runtime-dependencies.html)

## Implementation and acceptance update

The initial findings above record the starting state. The Linux CPU implementation now builds and runs from a relocated installation at `/home/jessica/Applications/h5-reader-linux-20260930`. The application menu entry **H5 Reader (Linux OSS, CPU)** invokes `launch-h5reader-cpu`, which disables accelerator inference, hardware video codecs, and hardware OpenGL rendering for this session's occupied-GPU constraint.

- Toolchain: public Qt 6.10.2 archives at `/home/jessica/opt/h5reader/qt/6.10.2/gcc_64`, VTK 9.5.2 compiled against that Qt, GCC 13.3, system HDF5 1.10.10 and Eigen 3.4, isolated LibArchive 3.7.2. Linux deployment explicitly requires Qt 6.10 APIs; Windows dependency discovery is unchanged.
- Native ML: CPU Torch 2.7.1, original model SHA-256 unchanged, no CUDA/HIP build. Python/native protocol parity produced identical outputs in three graph cases. Real F006 inference matched the existing isotropic and five tensor-component reference values in the dashboard, inspector, and scene. An omitted required input was correctly rejected with REST 409 and its exact missing filename.
- Final configured CTest run: **31 of 31 passed**, including all-frame publication comparison, model-input export and shutdown, dense Chignolin and Trp-cage DFT, CPU ML and stale-development-runtime fallback. Log: `/home/jessica/opt/h5reader/logs/ctest-final.log`; machine-readable result: `ctest-final.xml` in the same directory.
- The relocated application's main REST suite passed **59 tests**. Its 11 skips are fixture-specific ML, model-input, Chignolin and Trp-cage tests exercised in separate jobs; missing-input rejection was also exercised separately. Video tests created finalized MP4 files and decoded the expected changing frame sequence. Optional learned-activity capture tests passed separately (2 tests).
- Real native catalog download: BMR7057, 382,396,716 compressed bytes, 842 atoms and 100 frames. HTTPS transfer, native extraction, opening and cache reuse passed. Every frame, atom coordinate and catalog field matched the independent original through the native publication test. This is exhaustive for that dataset, not for all 22 catalog entries.
- The installed app passed its REST and CPU ML checks with ORCA and system Qt library/plugin paths deliberately present in the parent environment. Copied ELF files use relative runtime paths, and the launcher isolates the application runtime. Screenshots came from the app's own `/api/screenshot` interface; logs report Mesa llvmpipe.
- Debug static-library sanitizer link propagation was corrected and verified with a minimal consumer using the real CMake module. A complete Debug sanitizer suite has not been claimed.
- Windows dispatch behavior matched baseline across nine branch-policy scenarios; a native Windows Qt Pro/MSVC build remains the Windows side's acceptance step. No CUDA workload was run. No public release or new project/model license has been declared.

Reproducible commands and Windows overlap are recorded in `LINUX_BUILD_AND_WINDOWS_HANDOFF_2026-09-30.md` in this directory. Both notes are tracked with the Linux implementation for the advisor handoff; local toolchains, fixtures and binary artifacts remain outside Git.
