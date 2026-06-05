# Codex — strip-series dump route (continuity-test instrument)

## FIRST — read the qt6-cpp skill (open the files; full fs access)
`/home/jessica/.claude/skills/qt6-cpp/SKILL.md` + `references/architecture.md`. ACONNECT/ASSERT_THREAD, no
new QTimer, portable. Working root `/shared/2026Thesis/nmr-shielding/h5-reader`, scope `h5-reader/src` only,
no library link, no H5 write. **Verify from the code first.**

## Why
We're running an end-to-end continuity test: an independently-read NPY value (atom A, field F, per frame) must
be shown to reach the strip's `ChannelBuffer` and then the painted strip. Checkpoint 2 is "what the buffer
actually holds" — and that's not exposed anywhere. Add a read-only REST route that dumps the sampled series so
we can compare it numerically against the NPY ground truth. Pure observability, no behavior change.

## Read first
- `model::ChannelBuffer` (`src/model/StripChartChannel.h`) — `values` (std::vector<double>, NaN at gaps),
  `valid` (uint8_t), `lastFrame()`, `size()`. Series is indexed by frame from 0 (values[f] = value at frame f).
- How `DashboardDisplayController` stores active strip tracks + their `ChannelBuffer`s (the `stripTracks`
  it pushes to `DashboardStripDock`/`StripStackWidget`). Find the live track list + each track's
  signal id / descriptor id / display mode / channel id / buffer.

## Deliverable — GET `/dashboard/strip/series`
Return one object per active strip track:
```json
{ "tracks": [ { "signal_id","descriptor_id","mode","channel","label",
                "last_frame": int, "count": int,
                "values": [double...], "valid": [0|1...] } ] }
```
`values`/`valid` are the raw `ChannelBuffer` contents in frame order (so `values[k]` = the sample at frame k).
Emit NaN as `null` in JSON (so a gap is visible). Add a `const` accessor on `DashboardDisplayController` (and
wherever the track holds its buffer) if needed; RestServer just reads it. `ASSERT_THREAD(this)`, guard nulls.
(Optional: accept `?signal_id=` to filter, but the full dump is fine.)

## Output
Build green (`cmake --build .../build/linux-rwdi --target h5reader -j$(nproc)`). Do NOT launch, do NOT git.
Append a short `notes/INSTRUMENT_DASHBOARD_2026-06-05.md` line. Note any code that contradicts this brief.
