# Mac Qt accessibility crash — 11 October 2026

## Cause and conservative correction

The recorded Reader crash (`h5reader-2026-10-10-225634.ips`, PID 4845) occurred
in `-[QMacAccessibilityElement accessibilitySelectedChildren]` while macOS
accessibility inspected the run-selection table. Qt's synthesized table rows,
columns, and placeholder cells share their parent's accessibility identifier.
Two cleanup sites incorrectly deleted that parent interface, leaving selected
cell references invalid.

The proposed upstream [QTBUG-149612](https://qt-project.atlassian.net/browse/QTBUG-149612)
[change 765434](https://codereview.qt-project.org/c/qt/qtbase/+/765434)
adds the existing `isManagedByParent` ownership guard at those two sites and
extends Qt's native table regression. The change was still **NEW**, not merged,
when inspected on 11 October. This is a local, version-pinned backport; reassess
it when Qt supplies a corrected SDK.

`tools/macos/qt-cocoa` builds only the Cocoa plugin and native regression from
copies of the selected Qt 6.12.0 source. It checks the exact Qt version and source
hash and refuses build/install paths inside the SDK. It does not modify the
installed Qt plugin or sources. Reader's optional `H5READER_MACOS_COCOA_PLUGIN`
deployment input selects this independently built dependency in a clean staging
prefix. The Linux and Windows deployment paths remain unchanged. Standard Qt
widgets and accessibility remain enabled.

## Evidence on this Mac

Apple M3, macOS 27.0.1, arm64, Qt Pro 6.12.0, Xcode 27. The diagnostic build
instruments the Cocoa plugin and native Qt regression with ASan and UBSan;
it does not instrument all of the vendor Qt frameworks or all of Reader.

- Original SDK plugin: the extended native table regression fails because
  `QAccessible::accessibleInterface(tableId)` becomes null.
- Corrected plugin: all **11 native Mac accessibility tests pass**, including
  table identity, selected cells, column changes, trees, notifications,
  checkboxes, tabs, and windows.
- Separate ASan/UBSan plugin/test build: the same **11 tests pass**, with no
  reported sanitizer errors. Logs verify that the corrected plugin was loaded.
- Reader trajectory-library QtTest suite: **71 passed, 0 failed, 1 skipped**.
  The skip requires a separately supplied published-package fixture. A new
  `QAbstractItemModelTester` check covers filtering, empty selection, catalogue
  replacement, and persistent-index invalidation.
- Developer ID signed, deployed Reader: real accessibility scans of the
  176-entry public catalogue, filters for 4976 and 5292, a no-match filter,
  mouse selection of Ubiquitin, keyboard selection of Profilin, catalogue
  reopening/reloading, included-starter opening, and frame advancement pass.
  The supervised process (PID 14904) exits normally with code **0**, with no
  new crash reports during this run. The loaded plugin path is inside the app.
- The replacement installed in `/Applications` launches from the existing
  desktop Finder alias (PID 15021), loads Chignolin with all 5,001 frames, and
  survives the live catalogue accessibility scan. It is left open with Chignolin.
- The deployment closure check and deep strict signature verification pass.
  The deployed plugin UUID matches the corrected build, not the SDK plugin.

## Qt quality tools

Qt Creator's installed Clang-Tidy (LLVM 21.1.2) and Clazy 1.16 were used; no
extra tool installation was needed.

- Clang-Tidy's analyzer/lifetime checks on `ReaderCollectionDialog.cpp` report
  no enabled findings.
- Clazy level 1 on the picker completes. It flags an existing helper signal
  argument that should be fully qualified (`sciencefiles::Download::Result`),
  plus compiler warnings about Qt 6.12 variadic logging macros in C++17 mode.
- An optional additional Clazy pass over the entire large trajectory-library
  test translation unit was stopped after ten minutes and is **incomplete**.
  Its partial output also flags the existing broad `<QtTest>` module include.
  This does not affect the completed picker analysis or the executed QtTests.
- Clang-Tidy on the affected Cocoa implementation reports three nullability
  warnings and one possible `NSAttributedString` retain-count leak. A matching
  run on the original SDK source produces the **same four diagnostics**.
  These are recorded as upstream follow-ups; they are not silently called clean
  or folded into this ownership backport.

The saved evidence is under
`work/acceptance/macos-qt-accessibility-fix/`: original/fixed native regression
logs and XML, sanitizer results, analyzer logs and baseline comparison,
Reader tests, native picker AX snapshots/screenshot, deployment UUIDs,
signing evidence, normal exit, and crash delta. Matching Cocoa symbols are
retained at `outputs/libqcocoa-Qt6.12.0-QTBUG-149612.dSYM`. Reader's unchanged
main executable still matches the retained `outputs/h5reader.app.dSYM`.

## Distribution boundary

The earlier notarized `db081f8` DMG contains the original Qt plugin and is
superseded for new testing by this correction. Its notarization does **not**
cover the changed bundle. A replacement DMG needs fresh Developer ID signing,
Apple notarization, ticket stapling, and final package checks before being
presented as the distributable candidate. No Git push or public release is
part of this repair. This Mac's results are not evidence of testing on another
Mac, older supported macOS, Linux, or Windows.
