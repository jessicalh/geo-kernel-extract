# h5-reader — session handoff (2026-05-27): buffered strip chart + DFT shielding channel

This session rearchitected the strip chart into a **real, value-buffered** instrument
and added a **per-frame DFT shielding channel** — the advisor demo: pick atoms, play,
watch DFT σ move with structure. Built green on `linux-rwdi`; the data path is verified
in-app on `:1` (store parsed + validated a real frame). Full detail in memory
`project_h5reader_killer_app_multiatom_compare_20260526` (the 2026-05-27-later section).

---

## ⇩ PASTE PROMPT (next session)

```
Resume the h5-reader killer app — a CONTINUE. Deep-dive first; living state is the memory
store + the code. Work ONLY in /shared/2026Thesis/nmr-shielding/h5-reader/.

LANDED + committed: the value-buffer strip chart + DFT shielding channel.
- src/model/StripChartChannel.h — ChannelBuffer (std::vector we own; append-on-play, past
  kept, future never sampled, NO decimate, NaN-gaps) + ChannelSource (functor). Plain data.
- src/app/StripChartDock.{h,cpp} — two named Tracks (structural geometry + DFT shielding),
  setFrame APPENDS (extendTo: forward backfills, backward no-op), Qt Charts fed only the
  visible-window slice, gate on frameCount()>1. FFT flat-period (F6) + leading-zero (F7) fixed.
- src/model/DftShieldingStore.{h,cpp} — QObject citizen: lazy ORCA parse keyed by ORIGINAL
  frame index, meta.json->out_primary (not glob), STRICT validate (846 atoms + no holes +
  total≈dia+para; verified 0.001ppm). sample() cheap/cached; requestFrame() parses.
- Conformation::originalFrameIndex (base identity; TrajectoryConformation -> h5 frame_indices).
- QtLoadResult.runPath + ReaderMainWindow::locateDftJobsDir (dft/ sibling of extract/).
- notes/RUN_DESCRIPTOR_TOML.md — TOML format SPEC (parser NOT built; token-light).

DEMO dataset: /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft
— loads via extract/trajectory.h5 (846 atoms, 751 frames); dft/jobs has ~500 of 751 ORCA
frames (gaps real). Desktop launcher "h5-reader — 1P9J + DFT shielding".

NEXT (not started): descriptor PARSER (toml++); more DFT scalars (dia/para, |T2|, T1) as a
channel selector; PROMOTE DFT parse to a worker thread (frameReady is already the seam) so
continuous play does not hitch; the richer/VTK strip widget (reads the SAME ChannelBuffer).

RULES: never link/read the nmr_shielding library; Qt/VTK discipline (CENSUS_REGISTER, ACONNECT,
ASSERT_THREAD, QPointer); charts are Qt Charts (not a 2nd VTK surface); new types born honest
(no Qt prefix) but FULL CITIZENS stay QObjects; run on :1 (never offscreen); clangd
"not found / private member" diagnostics are editor noise — the cmake build is truth; commit
scoped `git add h5-reader/` (isolates the main-tree Codex); UDP 9997 has one consumer.
```

---

## Landed this session

- **`StripChartChannel.h`** — the value-buffer seam (`ChannelBuffer` + `ChannelSource`), the
  durable half that a richer widget will reuse. Plain data, born honest, non-QObject.
- **`StripChartDock` rewrite** — append-on-play over owned buffers; Qt Charts fed the bounded
  visible-window slice (display-strides only at µs scale; **buffer never decimated**); two
  panels (structural geometry + DFT shielding) sharing the playhead; FFT over the collected
  prefix. Fixed F6 (flat-spectrum false period) and F7 (leading-zero).
- **`DftShieldingStore`** — QObject citizen mirroring the `Conformation` facade: lazy parse of
  `dft/jobs/…_fNNNNNN_…/*_nmr.out` via `meta.json→out_primary`, strict validation, `frameReady`.
- **Plumbing** — `Conformation::originalFrameIndex`, `QtLoadResult.runPath`,
  `ReaderMainWindow::locateDftJobsDir`, the DFT panel follows the focus atom (σ total T₀).
- **Docs** — `notes/RUN_DESCRIPTOR_TOML.md` (format spec); this handoff; desktop launcher.

## Deferred / known (not this work)
- **Descriptor parser** — format written; reading it is future (loads work by convention now).
- **DFT parse is synchronous** on the GUI thread (fine at 5 fps; `frameReady` is the worker seam).
- **RESIDUAL_RENDER_DROP** — pre-existing intermittent ball-and-stick drop (`notes/RESIDUAL_RENDER_DROP.md`).
