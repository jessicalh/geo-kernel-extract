# tomlplusplus (vendored)

Single-header [toml++](https://github.com/marzer/tomlplusplus) v3.4.0, used by
`src/io/ReaderInputManifest` to parse `h5reader_manifest.toml`.

* Source: <https://raw.githubusercontent.com/marzer/tomlplusplus/v3.4.0/toml.hpp>
* SHA256: `6b5172ad4dd6519aec67b919181fa7a38a2234131e5b2afa232dfe444819783e`
* License: MIT (see `LICENSE`)

Vendored as a single-header drop-in — same pattern as `extern/HighFive/` at the
repo root. CMake adds `extern/tomlplusplus/include` to the h5reader target's
include path; no separate compile step.

## Updating

Replace `include/toml.hpp` with a newer single-header release, update the
URL / SHA256 above, rebuild. The TOML spec is stable; minor version bumps
have not broken `ReaderInputManifest`'s call surface in this repo's history.
