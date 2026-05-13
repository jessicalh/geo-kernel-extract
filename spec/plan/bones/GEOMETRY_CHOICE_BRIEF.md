# GeometryChoice: Agent Brief

A GeometryChoice is a runtime record of one geometric decision made by a
calculator. Created during evaluation, deposited on the conformation.
Every geometry choice in a calculator gets one. Must be drawable.

## What a GeometryChoice holds

- **What**: the decision pattern (radial threshold, topological exclusion, etc.)
- **Value**: the constant used (15.0 Å, factor 0.5, -12.0 nA, etc.)
- **Why**: physics text with citation at point of use
- **Which objects**: references to model objects touched (Ring, Bond, atom indices)
- **What happened**: how many evaluations, how many accepted/rejected
- **DomainKind**: drives visualization (sphere, highlighted atoms, shell, etc.)

## Agent workflow (per calculator)

1. Read this brief
2. Read learn/bones/CalculationArea.h (pattern vocabulary, DomainKind enum)
3. Read your calculator's section of GEOMETRIC_KERNEL_CATALOGUE.md
4. Read your calculator's .cpp file
5. Read KernelEvaluationFilter.h (the filters your calculator uses)
6. For every geometric decision in the code, create a GeometryChoice

Do NOT read OBJECT_MODEL.md or CONSTITUTION.md. Those are session docs.

## Rules

- Binary identical .npy output. No behavior changes.
- Every constant, threshold, exclusion, intensity, and cutoff gets a GeometryChoice.
- If it can't be drawn in a viewer given a conformation, it's incomplete.
- Object references use the typed indices from the model (ring_index, bond_index, atom_index).
- Physics text is the actual justification, not a label. "Ring current B-field
  decays as 1/r^3; at 15 Å contribution <0.03% of 3 Å value" not "ring cutoff."

## Existing pattern vocabulary (from learn/bones/CalculationArea.h)

| Type | DomainKind | Holds | Draws as |
|------|-----------|-------|----------|
| RadialThreshold | Spatial/Shell | radius + sense (include/exclude) | sphere at source |
| SourceRelativeExclusion | SourceRelative | factor of source extent | scaled sphere |
| RingBondedExclusion | Topological | per-ring excluded atom sets | highlighted atoms |
| SelfSourceExclusion | Topological | (identity only) | source atom highlight |
| SequenceGate | Sequence | min residue separation | backbone coloring |
| ShellBoundary | Shell | bin edge distance | concentric shell |
| SwitchingFunction | Switching | onset + cutoff | two shells + gradient |
| DecayFunction | Decay | decay length + unit | opacity gradient |
| RingCurrent | RingMagnitude | ring type + intensity (nA) | ring color by I |
| LobeOffset | RingGeometry | ring type + offset (Å) | lobe visualization |
| NumericalAccuracy | Numerical | threshold distance | wireframe sphere |
| ValueClamp | Numerical | max value + unit | (diagnostic only) |
| ValueGate | ValueThreshold | floor value + unit | (diagnostic only) |
| SentinelValue | Sentinel | marker value | (diagnostic only) |
| WholeDomain | WholeDomain | (identity only) | whole protein highlight |

## Per-calculator manifest

Line ranges are in GEOMETRIC_KERNEL_CATALOGUE.md for the physics.
File sizes are the .cpp line counts.

### 1. BiotSavartResult (442 lines) — kernel catalogue lines 119-295

Choices:
- RingHorizon 15.0 Å (IncludeWithin) — per ring iterated
- SingularityGuard 0.1 Å (ExcludeWithin)
- MultipoleInnerBoundary factor=0.5 — per ring, uses ring diameter
- RingBondedExclusion — per ring, topology walk
- PheRingCurrent -12.0 nA — per PHE ring
- TyrRingCurrent -11.28 nA — per TYR ring
- Trp6RingCurrent -12.48 nA — per TRP6 ring
- Trp5RingCurrent -6.72 nA — per TRP5 ring
- Trp9RingCurrent -19.20 nA — per TRP9 ring
- HisRingCurrent -5.16 nA — per HIS ring
- SixMemberedLobeHeight 0.64 Å — PHE/TYR/TRP6
- Trp5LobeHeight 0.52 Å — TRP pyrrole
- Trp9LobeHeight 0.60 Å — TRP indole perimeter
- HisLobeHeight 0.50 Å — HIS imidazole
- RingShellStrong 3.0 Å — counting
- RingShellModerate 5.0 Å — counting
- RingShellFarField 8.0 Å — counting
- RingShellBackground 12.0 Å — counting
- RingProximityDecay tau=4.0 Å — exp weighting

### 2. HaighMallionResult (390 lines) — catalogue lines 191-295

Shares ring-level choices with BiotSavart, plus:
- HmRefineNear 2.0 Å — adaptive quadrature L0→L1
- HmRefineClose 1.0 Å — adaptive quadrature L1→L2

### 3. McConnellResult (379 lines) — catalogue lines 298-377

Choices:
- BondAnisotropyHorizon 10.0 Å (IncludeWithin) — per bond iterated
- SingularityGuard 0.1 Å (ExcludeWithin)
- MultipoleInnerBoundary factor=0.5 — per bond, uses bond length
- SelfSourceExclusion — atom != bond endpoints
- NoDataMarker 99.0 Å — sentinel for nearest distances

### 4. CoulombResult (391 lines) — catalogue lines 380-437

Choices:
- CoulombWholeProtein — no cutoff, N² sum
- SelfSourceExclusion — i != j
- ChargeNoiseFloor 1e-15 e — skip negligible charges
- EFieldSanityClamp 100 V/Å — ceiling on E magnitude
- NearZeroNorm 1e-10 Å — skip degenerate directions

### 5. HBondResult (391 lines) — catalogue lines 562-592

Choices:
- HBondMaxReach 50.0 Å (IncludeWithin) — per H-bond partner
- HBondProximityShell 3.5 Å — counting within shell
- SelfSourceExclusion — atom != donor/acceptor
- MultipoleInnerBoundary factor=0.5 — uses N...O distance
- HBondSequenceGate min_sep=2 — through-bond exclusion
- NoDataMarker 99.0 Å — sentinel

### 6. PiQuadrupoleResult (277 lines) — catalogue lines 440-504

Choices:
- RingHorizon 15.0 Å — per ring iterated
- SingularityGuard 0.1 Å
- MultipoleInnerBoundary factor=0.5 — per ring
- RingBondedExclusion — per ring

### 7. RingSusceptibilityResult (255 lines) — catalogue lines 508-538

Choices:
- RingHorizon 15.0 Å — per ring iterated
- SingularityGuard 0.1 Å
- MultipoleInnerBoundary factor=0.5 — per ring
- RingBondedExclusion — per ring

### 8. DispersionResult (378 lines) — catalogue lines 542-560

Choices:
- DispersionTaper onset=4.3 Å, cutoff=5.0 Å — CHARMM switching
- RingBondedExclusion — per ring (reimplemented independently)
- MultipoleInnerBoundary factor=0.5

### 9. MopacCoulombResult (319 lines) — catalogue lines 595-673

Choices:
- CoulombWholeProtein — no cutoff
- SelfSourceExclusion — i != j
- ChargeNoiseFloor 1e-15 e
- EFieldSanityClamp 100 V/Å
- NearZeroNorm 1e-10 Å

### 10. MopacMcConnellResult (312 lines) — catalogue lines 676-757

Choices:
- MopacBondAnisHorizon 10.0 Å (IncludeWithin)
- SingularityGuard 0.1 Å
- MultipoleInnerBoundary factor=0.5
- SelfSourceExclusion
- MopacBondOrderFloor 1e-6 — skip negligible bonds

## Execution order

Interactive first (establish the pattern):
1. BiotSavartResult — most choices, most complex, sets the template
2. McConnellResult — different source type (bonds vs rings), validates generality

Then agents (one per calculator, check in between):
3. CoulombResult
4. HBondResult
5. HaighMallionResult (shares ring choices with BS, adds quadrature)
6. PiQuadrupoleResult
7. RingSusceptibilityResult
8. DispersionResult
9. MopacCoulombResult
10. MopacMcConnellResult

## Token budget for an agent

| What | Tokens | Required |
|------|--------|----------|
| This brief | ~2k | yes |
| learn/bones/CalculationArea.h | ~5k | yes |
| Calculator .cpp | ~3-6k | yes |
| Kernel catalogue section | ~2-3k | yes |
| KernelEvaluationFilter.h | ~3k | yes |
| **Total** | **~15-19k** | |

Leaves ~30k+ for the agent to think and write code.
