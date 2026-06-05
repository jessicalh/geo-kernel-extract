# Codex — strip history persists across rebuilds (the real fix) + revert the moot backfill

## FIRST — read the qt6-cpp skill (open the files; full fs access)
`/home/jessica/.claude/skills/qt6-cpp/SKILL.md` + `references/architecture.md` + `references/model-view.md`.
`ACONNECT`/`ASSERT_THREAD`, no new `QTimer`, portable. Root `/shared/2026Thesis/nmr-shielding/h5-reader`,
scope `h5-reader/src`, no library link, no H5 write. **Verify from the code first.**

## What we proved (so you fix the right thing)
A continuity test traced an NPY value (`npy:pos` atom-50 z) through to the strip. The data is correct when it
lands; the strips are empty/sparse because **`DashboardDisplayController::rebuild()` discards all strip history
on every call**: `series_ = std::move(next)` (`DashboardDisplayController.cpp:1734`) with `next` built fresh
(`buildGenericTracks`, `:1716`, makes new empty `SignalBuffer`s). `rebuild()` fires on signal add/remove,
mode-toggle, **atom-focus change**, and active-panel change — so accumulated history is thrown away on nearly
every interaction, and the re-fill from frame 0 can't recover the past (the single resident-snapshot slot only
holds the current frame). `requestSnapshot` is **synchronous** (not an async race) and the frame mapping is
clean (751↔751) — neither is the bug.

The earlier "backfill on snapshotReady" idea was based on a WRONG (async) diagnosis. **Revert it.**

## Part 1 — REVERT the moot snapshot-backfill you added
Remove what the snapshot-backfill brief added (it's a no-op given synchronous `requestSnapshot`, and its
`snapshotReady→onSnapshotReady` connection is dead weight):
- `DashboardDisplayController::onSnapshotReady` (`.cpp:~1313` + the `.h` declaration) and the
  `snapshotReady → onSnapshotReady` ACONNECT(s) (`~:1196-1201`).
- `SignalBuffer::backfill` (`src/model/SignalTimeSeries.{h,cpp}`).
- `ChannelBuffer::set` (`src/model/StripChartChannel.h`).
**KEEP** the temporary `event=strip_sample` diagnostic log in `extendToFrame` (I use it to verify; a later
follow-up removes it). If any of those methods are now referenced elsewhere, tell me rather than leaving
half-wired code.

## Part 2 — persist strip history across rebuilds
In `rebuild()`, immediately BEFORE `series_ = std::move(next)` (`:1734`), carry each existing buffer into its
matching new entry so history survives the rebuild:
- Key each `ActiveSeries` by the stable triple **(`signal.id`, `channel.id`, `displayModeId`)**.
- Build a lookup of the OLD `series_` by that key (a `QHash<QString, int>` over indices, or similar — Qt
  idiom, your call).
- For each entry in `next`, if its key matches an old entry, **`std::move` the old `buffer` into the new
  entry** (`next[i].buffer = std::move(series_[j].buffer);`). Genuinely new metrics keep their fresh empty
  buffer.
- Then `series_ = std::move(next)` as before. `extendToFrame(frame_)` on the next line will extend a carried
  buffer only if needed (no-op if already at `frame_`); it must NOT shrink or reset it.

Effect: a selected strip accumulates as you play and **survives** add/remove/focus/panel rebuilds; because you
no longer reconstruct the past, the single-slot eviction that produced the gaps never happens on the
forward-play path (the current frame's snapshot is resident when sampled). Do NOT change `extendToFrame`'s
sampling, the active-panel render filter, or the DFT path.

Note: a metric added mid-playback still can't fill frames BEFORE it was added (those snapshots are gone) — that
is acceptable and expected ("fills from when added"); do not try to backfill them.

## Output
Build green (`cmake --build .../build/linux-rwdi --target h5reader -j$(nproc)`). Do NOT launch, do NOT git.
Append a `notes/INSTRUMENT_DASHBOARD_2026-06-05.md` line (revert + persistence). Note any code that
contradicts this brief.
