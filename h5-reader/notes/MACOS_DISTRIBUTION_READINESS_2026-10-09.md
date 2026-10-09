# macOS distribution readiness — 9 October 2026

## Recommendation and scope

Ship the first public Mac release as an **Apple Silicon application in a Developer ID signed, notarized DMG**, with a drag-to-Applications installation. The existing bundle already carries Qt, VTK, HDF5, archive support, the native CPU inference helper, and the original model. End users should not need Qt, Homebrew, Python, or a compiler. A clean-machine test must confirm that conclusion for the final signed package.

This is an assessment accompanying the local source checkpoint, not a public-release declaration. No release signing, notarization submission, upload, or remote Git push was performed. The current development application remains available at `outputs/h5reader.app`, relative to the workspace containing `h5-reader-project` and `work`.

Keep release packaging in the Darwin deployment layer. Retain the independent Windows installer and Linux deployment paths, and keep machine-specific SDK paths in ignored local presets. Commercial versus open-source Qt is a toolchain/licensing choice, not a separate application source tree. An Intel build would require its own compatible dependency and ML runtime set and acceptance run; some universal Qt frameworks do not make this application universal.

## Evidence from the current application

| Item | Observed state |
| --- | --- |
| Application | Version `0.5.0`, native arm64 executable; approximately 443 MiB installed |
| Native dependencies | 125 Mach-O files inspected; every file supports arm64. The installed deployment audit resolves all non-system runtime dependencies inside the bundle |
| Load commands | No absolute non-system dependency loads or absolute install IDs remain. Two Torch libraries retain an unused upstream build RPATH; remove it from release copies before final signing |
| OS floor | Bundle and highest dependency minimum are macOS 14.4. Execution has only been checked on this macOS 27 host |
| Identity | `CFBundleIdentifier=jlh-test.h5reader`; empty icon and copyright fields. These need release values |
| Signature | Deep, strict ad-hoc signature verification passes. Gatekeeper execution assessment returns `rejected` (exit 3); this is not a distribution signature |
| Signing capability | One valid `Developer ID Application: Jessica Hansberry (5P76G9955J)` identity is present. Xcode's `notarytool` and `stapler` are installed. Notarization authentication has not been verified |
| Package generator | The non-Windows CPack default is currently `TGZ`; no Mac release DMG workflow is configured |
| Debug recovery | Matching executable/dSYM UUID `431F9BCB-7BE6-346A-A61E-EF2014578281`; dSYM retained separately, approximately 110 MiB |
| User data | Reader instance data and download caches use Qt user locations; inference intermediates use a temporary directory; crash reports use the user's data directory. Final read-only-install validation remains to be done |
| Included data | Model and manifest are bundled. Large scientific fixtures are outside the application; no sample dataset directory is configured |

The inspected executable SHA-256 is `a8d38630ec0627a384d656d9d0896e0fd9babccde97728e177cee05f6d3d0cfe`. The model SHA-256 is `101eba14f6a6e891fdf558f5126411b71ca042a9cfe09c67b98f86aa3ee79ff6`. The machine-readable audit is retained locally in `work/logs/macos-distribution-audit.json`. The two unused Torch RPATHs are in `libtorch.dylib` and `libtorch_cpu.dylib`, both referring to `/Users/runner/work/pytorch/pytorch/pytorch/build/lib`; this is upstream build metadata, not a currently resolved external dependency.

Existing validation includes 32/32 bounded CTest cases, 66 passing main REST cases with 22 explicit skips, the original F006 CPU scientific reference, relocated launch/inference, 12 passing pinch/picking cases, and native toolbar/checkbox checks. The broad suites preceded the focused camera and drawing follow-ups; those follow-ups have their own validation. These checks do not substitute for testing the final signed download on a second Mac. Exact coverage and unresolved observations are in `MACOS_BUILD_AND_VALIDATION_2026-10-09.md`.

## Testing preview packaged at the user's request

A limited tester handoff can precede the public-release work below. After an initial blocked launch, a tester can normally approve this ad-hoc-signed app through **System Settings > Privacy & Security > Open Anyway**. This is an exception for the individual app, without disabling Gatekeeper generally; managed Macs may restrict it. [Apple's instructions](https://support.apple.com/en-gb/102445).

The preview is `outputs/H5-Reader-0.5.0-preview-ef98b7c-macOS-arm64.dmg`, **130,678,128 bytes** (about 125 MiB), with a sibling `.sha256` file. SHA-256: `a5f5d3fb835d292849f16ba42fe88f87017255a872c215276e0b221edc746a9b`. It contains the unchanged app from source checkpoint `ef98b7c`, an Applications shortcut, and a short `Read Me.txt` covering installation, developer approval, opening data, and reporting problems. Datasets and debug symbols are not included. No upload was performed.

CPack's `DragNDrop` generator packaged the already installed application without rebuilding or re-signing it. All 534 file/symlink inventory entries match the original bundle, including regular-file contents and executable permissions. Disk-image integrity, SHA-256, and the packaged app's deep strict signature check passed. The three existing F006 runtime, scientific scalar, and scene/inspector tests passed in 3.21 seconds directly from the read-only mounted image, in a path containing spaces. The native process used a system-only PATH with development runtime overrides cleared; additional assertions confirmed the installed CPU runtime. This remains a same-Mac check, not clean-machine or older-OS acceptance.

Packaging configuration, manifest, and test evidence are in `work/logs/macos-preview-CPackConfig.cmake`, `macos-preview-manifest.json`, `macos-preview-f006.log`, `.xml`, and `macos-preview-installed-runtime.json`. The first test attempt encountered the expected single-instance guard while the existing Reader was open. After saving and closing that instance, the package tests passed and the original app was reopened with its dataset, frame, selection, overlays, and camera restored. Matching debug symbols remain in `outputs/h5reader.app.dSYM`.

## Work required before a public release

### 1. Establish release identity and metadata

Choose a stable publisher-owned reverse-DNS bundle identifier, retain it across updates, add a proper `.icns` application icon, and set version/build and copyright metadata in the Mac bundle configuration. Keep Qt application/organization names stable unless user-data migration is intentional. Give the download an explicit version and architecture, for example `H5-Reader-<version>-macOS-arm64.dmg`.

Qt 6.12 lists macOS 14.4 as its supported minimum. Treat that as a candidate support floor until the application and all dependencies pass on that OS; otherwise state a narrower tested range. [Qt supported platforms](https://doc.qt.io/qt-6/supported-platforms.html).

### 2. Add a repeatable release-signing stage

Keep ordinary local builds ad-hoc signed. For an explicit release build, accept the signing identity as local/CI configuration and finish all resource copies, dependency fixups, RPATH cleanup, and any stripping before signing. Preserve the matching dSYM and a manifest of versions and hashes separately.

Sign the nested native code and the complete application using Developer ID, hardened runtime, and a secure timestamp. Qt's `macdeployqt` supports `-sign-for-notarization=<identity>`, which enables those signing options; use the selected Qt SDK through the existing deployment API and verify the helper, third-party libraries, plugins, and outer bundle afterward. Start with minimal entitlements and exercise VTK, video export, and TorchScript under the signed runtime before adding any demonstrated exception. [Qt macOS deployment](https://doc.qt.io/qt-6/macos-deployment.html).

The current identity means obtaining an application-signing certificate is probably not the blocker. Notarization credentials and account access still need establishing in a local Keychain or protected CI configuration, outside Git. Apple requires Developer ID signing and hardened runtime for this distribution route. Submit the final distribution artifact with `notarytool`, inspect its result, staple the accepted ticket, and validate it. An ad-hoc signature or a successful local launch is insufficient. [Apple notarization workflow](https://developer.apple.com/documentation/security/notarizing-macos-software-before-distribution), [Apple notarization troubleshooting](https://developer.apple.com/documentation/security/resolving-common-notarization-issues).

### 3. Produce a Mac installation package

Add an Apple-only CPack `DragNDrop` configuration around the existing installed `.app`, using its Applications link. Keep Linux's `TGZ` and Windows' `INNOSETUP` choices intact. A package that merely copies this self-contained app does not need privileged installer scripts. [CMake DragNDrop generator](https://cmake.org/cmake/help/latest/cpack_gen/dmg.html).

Make the pipeline order explicit: install to a fresh staging directory, finalize and sign the application, create and sign the DMG, notarize, staple, then test and record the final download hash. If CPack performs installation again, signing must run on that final staged bundle; do not accidentally package a newly ad-hoc-signed copy. Perform no bundle edits after its final signature. Retain one matching symbol archive per release, outside the user download.

### 4. Finish redistribution materials

Prepared upstream notices are already installed under `Contents/Resources/licenses`, including Qt, VTK, HDF5, HighFive, Eigen, Boost.Math, LibArchive/liblzma, libaec, nanoflann, and PyTorch. The current CMake input makes the additional notices directory optional; make complete notices a release gate and preserve the dependency inventory and source/build provenance alongside it.

The tracked tree does not contain a Reader-level license, and the model manifest does not state distribution terms. Record the intended application and model terms before publishing. Confirm coverage under the actual commercial Qt agreement for the modules shipped, including Charts, whose alternatives are commercial licensing or GPLv3. [Qt Charts licensing](https://doc.qt.io/qt-6/qtcharts-index.html).

Qt's commercial license does not settle the separate terms of all bundled third-party code. The installed Qt inventory identifies FFmpeg 9.0.1; verify its exact build configuration and fulfill the applicable notices and source/relinking obligations for the binaries actually shipped. Keep a precise record rather than assuming an aggregate notice file completes that work. [Qt third-party license inventory](https://doc.qt.io/qt-6/licenses-used-in-qt.html).

### 5. Accept the actual download on another Mac

Use a clean standard user account on a Mac without the development dependencies, including the oldest OS claimed. Download the final artifact through a browser so normal quarantine/Gatekeeper checks apply. Confirm normal opening and copying to Applications without security bypass instructions.

Exercise first launch, opening local `.LGS` data through the application, collection downloads and caches, CPU prediction against the unchanged F006 reference, video export, Retina picking, pinch zoom, and tensor/toolbar controls in light and dark appearance. Check application operation when its bundle is not writable, user-selected paths containing spaces/non-ASCII text, and reopening after an app replacement. Collect crash/log evidence and match any crash to the retained symbols. Check ticket validation and launch without network after installation.

Run fresh Linux and Windows builds/tests before presenting the checkpoint as a validated cross-platform release; preserving their source branches is not fresh platform execution evidence. Mac signing and packaging changes should remain scoped to Darwin.

## Useful follow-ups that need not block the first download

- **Finder document opening:** no `.LGS` document declaration or `QFileOpenEvent` handling was found. File > Open already provides the initial route; add and test Finder double-click/open-with support if it is promised in installation instructions.
- **Approachable first use:** provide brief installation/opening instructions and either a small authorized sample or a clear collection link. Do not include the multi-gigabyte validation experiments in the app.
- **Size:** the installed application is about 443 MiB, with Torch CPU the largest component. Audit deployed plugin/framework requirements before trimming; do not manually prune the Qt SDK. Compression will determine the eventual download size.
- **Updates:** start with replacement of the application while preserving user data. Automatic updating and an Intel build can be separate work if needed.

The remaining work is chiefly a controlled release process and external acceptance, rather than rebuilding the application's portability architecture. The local checkpoint is a recoverable development milestone; the gates above define what would justify calling its Mac package ready to distribute.
