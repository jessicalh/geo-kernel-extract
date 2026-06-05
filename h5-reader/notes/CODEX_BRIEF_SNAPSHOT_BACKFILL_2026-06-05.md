# Codex — fix: strip backfills NPY snapshot values on `snapshotReady` (complete the intended pattern)

## FIRST — read the qt6-cpp skill (open the files; full fs access)
`/home/jessica/.claude/skills/qt6-cpp/SKILL.md` + `references/architecture.md` (signal/slot, the async-readiness
reaction pattern). `ACONNECT`/`ASSERT_THREAD`, no new `QTimer`, portable. Root
`/shared/2026Thesis/nmr-shielding/h5-reader`, scope `h5-reader/src` only, no library link, no H5 write.
**Verify from the code first.**

## The bug (proven by a continuity test — not a guess)
Strip series for NPY snapshot sources sample the per-frame value *before* the async snapshot is resident, so
they write a NaN gap that is never revisited. A traced `npy:pos` atom-50 z-series came back almost all `nan`,
except frame 3 (`27.7406`, dead-on the NPY ground truth `27.741`) which only landed because its snapshot
happened to be resident. The data and pipeline are correct; only the *timing* is wrong.

The architecture already declares the right pattern: `Conformation` (`src/model/Conformation.h:84-99`) — a
**single resident-snapshot slot**; `requestSnapshot(frame)` is async; "consumers cope with `snapshot()==null`
and **react to `snapshotReady(frame)`**". `snapshotReady` carries the frame index. The strip half-implemented
this — it `requestSnapshot`s (`DashboardDisplayController.cpp:2538`) then samples immediately
(`:2546`) without ever reacting to `snapshotReady`. (DFT is already synchronous — `DftShieldingStore`
`requestFrame` v1 — so it has NO race; do not touch the DFT path.)

Because there is ONE resident slot, sync-loading a past frame would evict the display's frame — wrong. The fix
is to fill each frame's value the instant THAT frame is the resident one, i.e. in the `snapshotReady(frame)`
handler. Keep the buffer's contiguous-from-0 invariant: the initial gap stays as a placeholder; we overwrite
it in place when the snapshot lands.

## Read first
- `model::SignalBuffer` (`src/model/SignalTimeSeries.{h,cpp}`) — wraps `channel` (a `ChannelBuffer`) +
  `statuses` + `gapReasons`. `append(FrameSignalSample)` pushes status/gapReason and `channel.append(value or
  nullopt)`. `lastFrame()` delegates to `channel`.
- `model::ChannelBuffer` (`src/model/StripChartChannel.h`) — `values`/`valid` vectors + `yMin/yMax/hasRange`.
  `append` updates the range only on valid samples (leading gaps don't pin the range — preserve that).
- `DashboardDisplayController`: `ActiveSeries` (`.h:153`) has `buffer` (SignalBuffer), `sample`
  (functor `FrameSignalSample(frame)`), `needsFrameSnapshot`. `conformation_` is set at `:1195`. The repaint
  signal is `stripTracksChanged()` (already emitted at `:1286`, `:1302`, `:1728`).

## Deliverable (additive; ~30 lines; do NOT change extendToFrame or the DFT path)
1. **`ChannelBuffer::set(std::size_t i, std::optional<double> v)`** — overwrite `values[i]`/`valid[i]` in place
   (precondition `i < size()`). On a finite `v`: set value, `valid[i]=1`, and extend the range exactly like
   `append` (init `yMin/yMax/hasRange` if first valid, else min/max). On `nullopt`/non-finite: NaN + `valid[i]=0`
   (don't shrink the range). (We only ever call it to fill a gap, so it never removes a valid sample.)
2. **`SignalBuffer::backfill(std::size_t i, FrameSignalSample sample)`** — overwrite `statuses[i]`,
   `gapReasons[i]`, and `channel.set(i, sample.value-or-nullopt-by-status)` in place. Precondition `i < size()`.
3. **`DashboardDisplayController::onSnapshotReady(std::size_t frame)`** (`ASSERT_THREAD(this)`):
   ```
   bool changed = false;
   for (ActiveSeries& s : series_) {
       if (!s.needsFrameSnapshot || !s.sample) continue;
       if (static_cast<long long>(frame) > s.buffer.lastFrame()) continue;   // not sampled yet
       if (s.buffer.isValidAt(frame)) continue;                              // already filled — add a tiny accessor if needed
       s.buffer.backfill(frame, s.sample(frame));                            // snapshot is resident now → real value
       changed = true;
   }
   if (changed) emit stripTracksChanged();
   ```
   Add a minimal `SignalBuffer::isValidAt(i)` (reads `channel.valid[i]`) or reuse an existing accessor — your
   call, keep it small.
4. **Wire it**: `ACONNECT(conformation_, &model::Conformation::snapshotReady, this,
   &DashboardDisplayController::onSnapshotReady);` where `conformation_` is set (`:1195`). Guard null.

Keep the temporary `event=strip_sample` diagnostic log in `extendToFrame` — I'll use it to verify the fix
(expect `finite=1` values after a snapshot lands), then a follow-up removes it.

## Output
Build green (`cmake --build .../build/linux-rwdi --target h5reader -j$(nproc)`). Do NOT launch, do NOT git.
Append a `notes/INSTRUMENT_DASHBOARD_2026-06-05.md` line. Note any code that contradicts this brief (e.g. if a
valid-at accessor or the range-extend already exists, use it and say so).
