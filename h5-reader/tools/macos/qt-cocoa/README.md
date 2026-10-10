# Qt 6.12.0 Cocoa accessibility ownership backport

The Mac platform plugin can delete its parent table accessibility interface
when a synthesized row or cell is released. Reader's run picker exposed this
as a crash in `accessibilitySelectedChildren`. The proposed upstream fix is
[QTBUG-149612](https://qt-project.atlassian.net/browse/QTBUG-149612),
[change 765434](https://codereview.qt-project.org/c/qt/qtbase/+/765434).
It guards the two cleanup sites with the existing `isManagedByParent` predicate.
Accessibility and the standard Qt widgets remain enabled.

This is an optional standalone dependency build, separate from Reader, VTK,
HDF5, and the Windows/Linux builds. It reads the selected Qt 6.12.0 SDK and its
installer-provided sources, copies only Cocoa and the native accessibility test
into the build directory, then applies the proposed upstream patch there.
The source hash and exact Qt version are checked. No installed SDK file is
edited. Reassess/remove the backport when adopting a vendor build containing
the fix; do not reuse the binary with a different Qt version.

For example, with machine-specific paths supplied by the developer:

```sh
cmake -S tools/macos/qt-cocoa -B /path/to/cocoa-build -G Ninja \
  -DCMAKE_PREFIX_PATH=/path/to/Qt/6.12.0/macos \
  -DH5READER_QTBASE_SOURCE_DIR=/path/to/Qt/6.12.0/Src/qtbase \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=14.4
cmake --build /path/to/cocoa-build
ctest --test-dir /path/to/cocoa-build --output-on-failure
```

Run the native test from a logged-in Mac GUI session. It covers normal widgets,
notifications, selected table cells, table identity, and column-count changes.
`-DH5READER_COCOA_SANITIZERS=ON` instruments the plugin and the native regression
with AddressSanitizer and UndefinedBehaviorSanitizer for a separate diagnostic
build; package the ordinary build.

Supply the result to Reader's local configure preset:

```sh
cmake --preset your-mac-preset \
  -DH5READER_MACOS_COCOA_PLUGIN=/path/to/cocoa-build/plugins/platforms/libqcocoa.dylib
cmake --build --preset your-mac-preset
cmake --install /path/to/reader-build --prefix /path/to/new-empty-staging-prefix
```

Use a **clean staging prefix** for this deployment. Qt's deployment tool keeps
the supplied plugin and handles the remaining Qt frameworks and signing.
The plugin uses the bundle-relative Frameworks path. Check the deployed
plugin UUID against the build product, run the real Reader picker accessibility
test, then perform final Developer ID signing and notarization of the changed
package. The old notarization does not cover the changed bundle.
