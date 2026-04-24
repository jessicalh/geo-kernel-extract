# TICKLERS.md

Tentative observations surfaced while summarising reference material. Each entry is a specific claim from a specific paper that **might** bear on our code, statistics, or thesis argument — but where we do not yet have the framing to decide.

**What this file is for:** a parking lot. Review at the start of relevant coding or stats sessions to see whether an entry has become actionable. Entries are explicitly not load-bearing.

**What this file is NOT for:**
- Primary citations we plan to use in arguments (those become evidence cards once we have the PDF in hand and the argument shape is sharp)
- Claims we have validated against our own data (those become feedback or project memories)
- Claims we have rejected (those get resolved and removed, with a one-line note in the journal below)
- Papers to acquire (those go to `../references/PENDING_ACQUISITIONS.md`)

**Entry format:**
```
## YYYY-MM-DD — [TENTATIVE] — Short title

**Source:** CiteKey → path to summary

**Claim (paraphrased):** specific, narrow

**Why it might matter:** which code surface or statistical analysis would be affected

**How to validate:** what evidence would confirm or reject this

**Status:** tentative | validated | rejected | superseded — one line on what changed
```

Status transitions: tentative → validated (promoted to memory / evidence card / applied in code) OR tentative → rejected (with reason) OR tentative → superseded (a better source replaced this). Either way, the entry moves out of the "active" section to the journal at the bottom.

---

# Active tickers

## 2026-04-23 — [TENTATIVE] — Johnson-Bovey two-loop separation may be unnecessary

**Source:** agarwal-1977-canjchem → `agarwal-1977-ring-currents-local-anisotropy-paracyclophane-summary.txt`

**Claim (paraphrased):** Agarwal et al. fit [10]-paracyclophane out-of-plane proton shifts against several values of the JB two-loop separation S. S=0 (single loop at ring plane) gave r=0.9928; S=1.28 Å gave r=0.8883; S=1.63 Å gave r=0.7634. They conclude the two-loop structure is numerically unnecessary once local anisotropic contributions are handled separately. The result rests on 1977 parameterisations of ¹³C shielding tensors and Allinger mechanics geometries — we do not yet have framing that generalises it to protein BS kernel evaluation.

**Why it might matter:** `BiotSavartResult.cpp` likely assumes a non-zero loop separation (Case 1995 used 0.64 Å; Cross-Wright 1985 used 1.28 Å). If S is effectively irrelevant when local-anisotropy kernels are present in the same model, then its value is carrying no information in our Stage 1 ridge regression — which would change how we interpret BS-kernel weights relative to McConnell bond-anisotropy kernels.

**How to validate:** Run the Stage 1 ridge with BS kernels computed at S=0 vs. S=0.64 vs. S=1.28 and check whether per-element R² and kernel weights move. If they move substantially, the separation IS load-bearing for us; if they don't, Agarwal's point generalises.

**Status:** tentative — 1977 paper, no protein-scale framing yet.

## 2026-04-23 — [TENTATIVE] — Ring normal vectors from 3 atoms are unstable under MD fluctuations

**Source:** sahakyan-vendruscolo-2013-jpcb → `sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-summary.txt`

**Claim (paraphrased):** Sahakyan & Vendruscolo report that computing an aromatic ring's normal vector from three ring atoms produces unstable directions when those atoms undergo out-of-plane fluctuations (e.g. in MD). They propose averaging two independently computed normals from disjoint atom subsets for stability. Their context is MD-restrained chemical-shift simulations; unclear whether the instability magnitude matters for our Stage 2 trajectory averaging.

**Why it might matter:** `BiotSavartResult.cpp` and `HaighMallionResult.cpp` both rely on a ring normal for their geometric factors. If our implementation uses a fixed 3-atom definition and the normal wobbles materially frame-to-frame in a 25 ns trajectory, T2 predictions could carry a bias that disappears (or worsens) depending on the choice of atoms used. For flat rings (PHE, TYR, HIS, TRP aromatic ring) the out-of-plane wobble is small but real.

**How to validate:** For one protein's 600-frame trajectory, log the ring-normal vector per frame (per ring) using the current 3-atom method, compute its angular variance. If variance is < 2° the choice probably doesn't matter; if variance reaches 5–10° or more, adopting the Sahakyan-V averaging is worth doing.

**Status:** tentative — pending a diagnostic run.

## 2026-04-23 — [TENTATIVE] — Case 1995 ring current intensities differ from Giessner-Prettre especially for HIS, TRP-5mem

**Source:** case-1995-jbnmr → `case-1995-ring-current-calibration-summary.txt`

**Claim (paraphrased):** Case 1995's DFT-calibrated ring-current intensity factors are larger than the 1969/1987 Giessner-Prettre & Pullman values, with the largest gaps at histidine (0.53 → 1.35–1.40, ~2.6×) and the 5-membered ring of tryptophan (0.56 → 1.02–1.32). Phe and Tyr move less (1.00 → 1.27–1.46, 0.94 → 1.10–1.24). Whichever set of intensities our library uses will change the per-atom kernel magnitudes, especially in pockets near HIS / TRP — and both appear frequently in proteins as active-site participants.

**Why it might matter:** If `data/calculator_params.toml` carries Giessner-Prettre era intensities, our BS kernel output for HIS/TRP-5mem neighbours is systematically low. If it carries Case 1995 values, we are aligned with the DFT-calibrated literature. Either way, thesis methods chapter needs to state which set we use and why.

**How to validate:** Read `data/calculator_params.toml` and compare. One-line fix if mismatched — but do not change values without checking whether the Stage 1 regression was already trained on the current ones.

**Status:** tentative — requires a peek at params file before coding.

## 2026-04-23 — [TENTATIVE] — Fast bond-length vibrations average out within ~200 fs

**Source:** de-gortari-2010-jacs → `de-gortari-2010-time-averaging-nmr-chemical-shifts-mlf-peptide-summary.txt`

**Claim (paraphrased):** In CP-MD of MLF peptide, C=O, Cα–N and other bond lengths oscillate with ~20 fs period. ¹³C shifts computed along this trajectory fluctuate by up to 20 ppm instantaneously, but the fluctuations average out within ~200 fs. This justifies their use of ps-spaced snapshots without explicit zero-point vibrational correction.

**Why it might matter:** If we sample our 25 ns trajectory at 1 ps spacing (or wider), we are above the 200 fs averaging window for fast bond vibrations — so those don't contaminate our per-frame kernel values. Useful as a defense of our sampling rate; suggests that ZPVC is not the dominant remaining error.

**How to validate:** Not actionable right now — more an argument-ready fact to park for the thesis methods chapter or a committee question. If Stage 2 residuals turn out larger than expected and we rule out other causes, revisiting ZPVC per Vaara et al. 1998 (Facelli's ref) may be warranted.

**Status:** tentative / argument-ready.

## 2026-04-23 — [TENTATIVE] — NICS_zz is more accurate than NICS_iso for aromaticity, by analogy supports T2 > T0

**Source:** gershoni-poranne-stanger-2015-csr → `gershoni-poranne-stanger-2015-magnetic-criteria-aromaticity-summary.txt`

**Claim (paraphrased):** The aromaticity-criteria literature has progressively moved from NICS_iso toward NICS_zz and NICS_πzz because the z-component of the shielding tensor at a ring-plane probe isolates π-ring-current effects, while the isotropic trace averages over σ-framework contamination. This is a parallel argument to ours — that T2 carries information the trace T0 does not — coming from an entirely different subfield.

**Why it might matter:** If a reviewer or committee asks "why should we care about T2 separately from T0?", the NICS_zz story is independent prior art for the general "tensor component beats scalar trace" move. Useful in the thesis introduction or defense, not in code.

**How to validate:** Not actionable for code. Would be worth checking the specific NICS_zz validation papers (Gershoni-Poranne refs 91–97) before citing; they specify quantitatively how much better NICS_zz is, which would sharpen the analogy.

**Status:** tentative / argument-ready.

## 2026-04-24 — [TENTATIVE] — Stage 2 validation: shift-improvement regions should coincide with R1/R2 ms-exchange regions

**Source:** markwick-2010-jacs-commun → `markwick-2010-enhanced-conformational-sampling-chemical-shifts-summary.txt`

**Claim (paraphrased):** Markwick et al. observed that the per-residue improvement in chemical-shift prediction (X-ray → MD-averaged) correlates spatially with regions of ms-timescale exchange dynamics independently measured by R1/R2/hetNOE relaxation on IκBα (Cervantes et al. 2009). Ubiquitin, with more uniform dynamics, showed smaller and less-patterned improvement. Their argument: if MD averaging is recovering real dynamics physics rather than noise, the improvement should spatially track the independent dynamics signal.

**Why it might matter:** This is the cross-validation template for our Stage 2 shift-tensor predictions. If our 25 ns MD averaging is picking up real dynamics (rather than just sampling noise that smooths the static prediction), the per-residue improvement vs. static DFT should spatially coincide with regions showing relaxation-derived mobility — measurable via BMRB R1/R2 entries or via internal analyses of the trajectory itself. Failing to see this pattern would be evidence that our averaging is a statistical artifact, not physics.

**How to validate:** For our 10 calibration proteins (Stage 2 sample), pull per-residue R1/R2/hetNOE from BMRB where available and correlate with per-residue improvement in tensor-prediction residuals (Stage 2 predicted minus static reference vs. experimental). Ubiquitin is in the set already; see if it behaves similarly to Markwick's (less improvement in uniform regions, more in flexible loops). This is a diagnostic, not a training target — we should not optimise for it.

**Status:** tentative — pending Stage 2 completion and analysis. If the pattern holds, this becomes an argument-ready defense of the averaging approach.

---

## 2026-04-24 — [TENTATIVE] — Stage 2 calibration set should be audited for dynamical diversity

**Source:** markwick-2010-jacs-commun → `markwick-2010-enhanced-conformational-sampling-chemical-shifts-summary.txt`

**Claim (paraphrased):** In Markwick's demonstration, the magnitude of MD-averaging improvement scaled with the system's dynamical heterogeneity — IκBα (heterogeneous) showed large improvement, ubiquitin (more uniform) showed smaller improvement. A system with no dynamical interest would presumably show negligible improvement.

**Why it might matter:** Our 10-protein Stage 2 calibration set needs enough dynamical diversity to produce a signal we can measure. If the current 10 proteins are all similarly well-ordered, our Stage 2 residuals will be dominated by other sources of error (force-field, functional choice, sampling), and we may conclude the kernels don't benefit from MD averaging when the real answer is "not on this set." Need to check the dynamical profile of the 10.

**How to validate:** Look at BMRB R1/R2 entries for the 10 calibration proteins; compute S²  distributions or equivalent. If dynamical range across the 10 is narrow, either supplement with a more dynamic protein or reframe Stage 2 claims around the subset that does show diversity.

**Status:** tentative — coding-session diagnostic.

---

## 2026-04-23 — [TENTATIVE] — Classical ring-current models cap around 0.22 ppm RMS on protein ¹H shifts

**Source:** moyna-1998-jcics → `moyna-zauhar-williams-1998-comparison-ring-current-methods-structure-refinement-summary.txt`

**Claim (paraphrased):** Moyna et al. compared HM, JB, and reparametrized point-dipole models against 2992 protein ¹H assignments and found RMS errors of 0.220 / 0.228 / 0.226 ppm respectively. This is the scalar-shift accuracy floor of classical (non-DFT-calibrated) ring-current kernels against experimental data as of 1998.

**Why it might matter:** Provides a benchmark floor against which our DFT-calibrated tensor kernels can be reported. If our Stage 1 scalar-shift residual (from the T0 channel alone) is comparable or better than 0.22 ppm, we are at the state-of-the-art classical level despite the much larger protein set and more kernels.

**How to validate:** Pull the T0-channel residuals from Stage 1 ridge regression on protein ¹H subset. Compare to 0.22 ppm. Not a one-line lookup — needs a small script — but not load-bearing either, just framing.

**Status:** tentative / framing comparator.

---

# Journal (resolved entries)

Entries move here when validated, rejected, or superseded. Keep one line on what changed so future readers see the trail.

(Empty as of first write, 2026-04-23.)
