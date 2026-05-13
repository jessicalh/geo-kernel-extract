# Fix-tests register — open test-shape items

Living register of test-shape problems noticed in passing. Items land
in real fixes when scope and timing allow; this file is the queue in
the meantime. New items append at the bottom; closed items move to the
"Closed" section with the resolving commit.

---

## 1. Narrow `RunConfiguration` in trajectory tests slides the suite away from production

**First noted:** 2026-04-25.

**Symptom.** `tests/test_gromacs_streaming.cpp` (and likely siblings)
constructs a one-off `RunConfiguration` per test, enabling only the
narrow `ConformationResult` set the TR under test depends on, and
explicitly setting `skip_apbs = true`, `skip_coulomb = true`,
`skip_dssp = true`, with no AIMNet2 requirement. Comments in the test
say things like `// narrow config; no AIMNet2 needed.` See e.g.
`test_gromacs_streaming.cpp:144-159, 165` for the
`BsWelfordAttachAndFinalizeTest` shape.

**Why this is wrong.** Production runs use
`RunConfiguration::PerFrameExtractionSet()` which exercises the full
classical stack + APBS + AIMNet2. Tests that exercise narrow configs
catch single-TR plumbing but never surface integration issues that
production hits — including AIMNet2 startup
(`reference_nvrtc_rpath_fix` memory) and APBS paths. The
`feedback_no_convenience_exclusions` rule is explicit: *"APBS +
AIMNet2 belong in the baseline. Excluding core calculators for
test-speed slides the suite away from production and compounds into
AI-driven partiality."* This is exactly that.

**Why it landed.** Prior AI sessions chose the narrow-config shape for
test speed and to side-step the AIMNet2 setup hassle (the recurring
NVRTC rpath issue). User report 2026-04-25: *"this got AI'd under my
nose and I have been told doesn't happen by other yous."* The pattern
is in the tree at every TR test now and removing it is a real
refactor, not a one-test edit.

**What right looks like.**

- Tests that exercise the trajectory plumbing run against
  `PerFrameExtractionSet()` (or `FullFatFrameExtraction()` where MOPAC
  participation is the point). Narrow shapes are reserved for
  per-TR discipline tests (Frame0Semantics, FinalizeIdempotency,
  H5RoundTrip), which test invariants of one TR — those don't need
  the full pipeline.
- AIMNet2 model is supplied through a CTest fixture: a wrapper
  script `scripts/run_with_cuda_env.sh` that exports the bundled-cu13
  `LD_LIBRARY_PATH` (per `reference_nvrtc_rpath_fix`) and execs the
  test, registered via `set_tests_properties(... ENVIRONMENT ...)`,
  OR patchelf the bundled `libnvrtc.so.13` to add `$ORIGIN` so it
  resolves its sibling builtins without env-var help. Either lands
  in version control once and stops drifting.
- Test runtime budget is paid honestly. APBS adds ~30s/frame, AIMNet2
  adds startup. Pipeline-shaped tests run on stride-99999 fixtures
  (frame 0 only) so the cost stays bounded; full-trajectory cost is
  reserved for `fleet_smoke_tests` and equivalents.

**Sequencing.** Land alongside the AIMNet2-CTest infrastructure, not
as a standalone refactor. `cold_read_review_20260424.md §5` already
flagged the inline-vs-separate-file test sprawl as an unresolved
follow-up; consolidating that and removing the narrow-config shape
are both bigger than a single sitting and want batching.

**Scope of debt.** Every per-TR test in `test_gromacs_streaming.cpp`
post-2026-04-23 uses a narrow `RunConfiguration`. The
`BondLengthStatsEndToEnd / Frame0Semantics / FinalizeIdempotency /
H5RoundTrip` quartet is the canonical shape; six other tests follow
suit. Refactor target: ~700 lines of test rewiring + the launcher
infrastructure decision.

---

## 2. AIMNet2 / NVRTC test invocation — RESOLVED 2026-04-25

**The procedural failure.** Across many sessions the canonical
test-invocation motion was `./build/<binary>` directly. AIMNet2 /
libtorch's CUDA JIT requires the bundled cu13 nvrtc-builtins on
`LD_LIBRARY_PATH`; without it, the binary fails at first inference
with `nvrtc: error: failed to open libnvrtc-builtins.so.13.0`. The
recurring fix was an `LD_LIBRARY_PATH` export in the invoking shell
that drifted out across sessions, so the same NVRTC error appeared
50+ times. Worse, individual tests with narrow configs avoided AIMNet2
entirely and passed, normalising the one AIMNet2-using test
(`GromacsStreaming.TrajectoryRunDrivesLoop`) as "expected failure" —
sliding the calculator out of the customary test path even when
operators were nominally diligent.

**The fix that landed.**

1. `scripts/run_with_cuda_env.sh` — universal launcher. Resolves the
   bundled cu13 directory via `python3 -c "import nvidia.cu13"` and
   prepends `${path}/lib` to `LD_LIBRARY_PATH` before exec'ing its
   args. One source of truth for the path; works on any machine with
   the cu13 pip package installed.
2. `CMakeLists.txt` — every `gtest_add_tests` switched to
   `gtest_discover_tests` with `PROPERTIES ENVIRONMENT
   "LD_LIBRARY_PATH=..."`. Path resolved at configure time via the
   same Python lookup. CMake errors out at configure if `nvidia.cu13`
   is not importable, surfacing a deployment problem early instead of
   at first AIMNet2 use.
3. `spec/TEST_FRAMEWORK.md` documented invocation switched from
   `./build/<binary>` to `ctest`. Manual invocation when needed
   routes through the launcher script. No `./build/<binary>` in the
   documented path means the procedure can't auto-skip the env
   setup.
4. Memory entry `feedback_test_invocation_via_ctest` records the
   rule for future sessions.

**Verification.** `GromacsStreaming.TrajectoryRunDrivesLoop` —
broken across at least 2 days (prior memory entry from
2026-04-23) — passes via `ctest -R 'TrajectoryRunDrivesLoop'`
in 904 s on batcave (full 25 ns trajectory, AIMNet2 every frame).

---

## 3. Item 1 (narrow-config test sprawl) is still open

The narrow-`RunConfiguration` issue documented in §1 is unrelated to
the AIMNet2 invocation fix in §2 — narrow configs sidestep AIMNet2
*by attaching nothing that needs it*, not by failing to set up
the environment. Now that `ctest` is the documented path and AIMNet2
runs cleanly through it, the next refactor that revisits the test
shape can move trajectory-plumbing tests onto `PerFrameExtractionSet()`
without dreading the rpath issue.

---

(extend with new items as they surface)
