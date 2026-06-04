# Strip-chart analysis methods — done, planned, deprioritized

> **Historical idea bank — not current truth (trued 2026-06-04).** Useful
> provenance, but old "planned" language must not drive current work without a
> fresh check against `UI_STATE_OVERVIEW_2026-06-04.md`.

The strip-chart dock (`src/app/StripChartDock`) turns the selected geometric
observable into views that reveal MOTION. This note tracks which analysis
methods are in, which are worthy future work, and which were tried and pulled —
so the queue is never lost (the project keeps every proposed method on
principle).

## Centerpiece — DONE (2026-05-27)

**Dihedral / angle / distance time series + FFT power spectrum.**
- Time-domain strip: the observable vs frame (scrolling window, science grids,
  playhead, digital readout).
- Frequency-domain: `model::ComputePowerSpectrum` (mean-subtract + Hann window,
  Eigen FFT), dominant-period readout. The dihedral is phase-unwrapped before
  the transform so a ±180°-straddling oscillation reads as one clean peak.

Why this is the headline (user, 2026-05-27): seeing periodicity in a dihedral is
**directly observable and assumption-free** — the most compelling, defensible
demonstration. Angular observables carry the periodicity (rotamer hops, ring
flips, methyl spins, librations); distances are mostly diffusive.

## Worthy future methods (user-endorsed, ranked)

Each slots in as a new tab below the always-visible time strip (retab the dock
when the first of these lands).

1. **PCA / dihedral-PCA collective modes** *(user: "would also have been
   worthy").* Coordinate-covariance PCA → dominant concerted 3-D motion; or
   circular dPCA (cos/sin embedding, Altis/Stock) → coupled angular mode; FFT
   the projection. The data reveals the periodic *direction*; finds coupled
   multi-dihedral periodicity rather than transforming each in isolation.
2. **Recurrence plot / RQA** *(user: "and so would 4").* Image
   R(t,t′)=[‖x(t)−x(t′)‖<ε] over the 3-D state: periodic = diagonal stripes,
   quasi-periodic = broken diagonals, diffusive = none. Non-parametric, visual,
   vetting-viewer-friendly.
3. **Spectrogram / wavelet** (time-frequency). A global FFT assumes a constant
   period; MD periodicity is often intermittent (a ring flips, then sits). A
   short-time FFT / CWT shows *when* each frequency is active.

Lighter upgrades to the existing FFT panel: **multitaper/Welch** PSD for cleaner
low-variance peaks; **cross-spectral coherence** between two channels for
coupled periodicity; a **dominant-component overlay** (reconstruct the peak
frequency's cosine and draw it over the time series — the literal "fft fit").

## Deprioritized — prototyped and pulled (2026-05-27)

**Orientational ACF → S² / τc / J(ω)** (the NMR relaxation statistic). Track the
selection's 3-D orientation (bond vector / plane normal), P₂ autocorrelation
C(τ), order parameter + correlation time; FT is the spectral density J(ω).

Pulled because: computed in the **lab frame without removing global tumbling**,
S²/τc may not validly apply to our data — presenting a relaxation statistic that
may not hold is an overclaim risk at a viva. Dihedral periodicity is the honest,
more impressive demonstration (user). The math (`ComputeOrientationalAcf`) was
straightforward and is recoverable from git/this note **if** an aligned /
tumbling-removed path is built first; only then would S²/τc/J(ω) be defensible.
