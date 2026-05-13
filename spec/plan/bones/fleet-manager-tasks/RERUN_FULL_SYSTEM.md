# Fleet Task: Re-run Production with Full-System Trajectory Output

**Priority: BLOCKING for solvent calculators (WaterFieldResult,
HydrationShellResult). Code is built and waiting.**

## The problem

The original production runs used this in `prod.mdp`:

```
compressed-x-grps = Protein
```

This discarded water and ion coordinates during simulation. The .xtc
files in results/ contain only protein atoms (e.g., 4876 of 154,944
for 1Q8K). The .gro has the full system (one frame). The water
positions needed for explicit-solvent NMR feature extraction were
never saved.

## The fix

Change one line in the .mdp and restart. The simulations are
equilibrated — PLUMED biasing is converged, the FES is built.
This is not a new run, it's a replay to capture what was thrown away.

```
compressed-x-grps = System    # was: Protein
```

## Disk budget

Full-system .xtc sizes (estimated):

| Scenario | Atoms | Per frame | 1ns (100 frames) | 5 walkers | 685 proteins |
|----------|-------|-----------|-------------------|-----------|-------------|
| Full system | ~155K | ~1 MB | ~100 MB | 500 MB | **340 GB** |
| Protein + water within 15A | ~25K | ~170 KB | ~17 MB | 85 MB | **58 GB** |
| Protein only (current) | ~5K | ~35 KB | ~3.5 MB | 17 MB | 12 GB |

**Decision: Save full system, post-process on batcave.**

10ns full-system output = 679 GB per machine against 1.2 TB free
(57%). Comfortable. Don't be clever on the fleet — save everything,
filter locally after harvest. The 15A water shell can be extracted
on batcave with `gmx trjconv` if storage becomes tight.

Measured disk (2026-04-12): scan2 and scan3 each have 1.2 TB free.
scan1 not reachable (VPN). Verify before starting.

## How much trajectory do we need?

We don't need 10ns again. The biasing is converged. We need enough
frames to:
1. Sample the thermally accessible conformational range (already mapped by FES)
2. Capture water packing variation around each atom
3. Feed the frame selection algorithm (TBD — clustering, diversity)

**10ns of production — same as original run.** The goal is variety
in the conformational ensemble, not speed. 10ns at 10ps stride
gives 1000 frames per walker with full water positions. Frame
selection happens downstream; save everything now.

Wall time (median 364 ns/day): 10ns takes ~40 minutes per protein.

Total wall time for 685 proteins on 3 machines:
- 685 × 40 min / 3 machines ≈ 6.3 days
- With overhead (restarts, I/O, scheduling): ~8 days

Disk per machine: ~679 GB (228 proteins × 3 GB each).
Available: 1.2 TB. Headroom: 521 GB.

## Step-by-step procedure

### 1. Verify scan machine disk space

```bash
ssh scan1 "df -h /opt/"
ssh scan2 "df -h /opt/"
# scan3 when recovered
```

Need ~40 GB free per machine for protein+near-water, or ~115 GB
for full system.

### 2. Prepare modified .mdp

On batcave, create `prod_fullsystem.mdp`:
```bash
cp /shared/fleet/clean_fleet/results/1Q8K_10023/run_params/prod.mdp \
   /shared/fleet/clean_fleet/prod_fullsystem.mdp

# Edit: change compressed-x-grps
sed -i 's/compressed-x-grps  = Protein/compressed-x-grps  = System/' \
    /shared/fleet/clean_fleet/prod_fullsystem.mdp

# Edit: reduce nsteps for 2ns (current is 5000000 = 10ns)
sed -i 's/nsteps.*/nsteps = 1000000/' \
    /shared/fleet/clean_fleet/prod_fullsystem.mdp
```

### 3. Deploy to scan machines

```bash
for machine in scan1 scan2; do
    scp /shared/fleet/clean_fleet/prod_fullsystem.mdp \
        ${machine}:/opt/fleet-reprep1/prod_fullsystem.mdp
done
```

### 4. Generate new .tpr files and run

For each protein on each machine:

```bash
PROTEIN=1Q8K_10023
STAGING=/opt/fleet-reprep1/staging/${PROTEIN}
RESULTS=/opt/fleet-reprep1/staging/${PROTEIN}  # or results path

for w in 0 1 2 3 4; do
    WDIR=${RESULTS}/walker_${w}

    # Generate new tpr from modified mdp + existing state
    gmx_mpi grompp \
        -f /opt/fleet-reprep1/prod_fullsystem.mdp \
        -c ${WDIR}/md.gro \
        -p ${STAGING}/${PROTEIN}_topol.top \
        -r ${WDIR}/md.gro \
        -o ${WDIR}/md_fullsys.tpr \
        -maxwarn 1

    # Run 2ns production with full-system output
    cd ${WDIR}
    gmx_mpi mdrun -s md_fullsys.tpr \
        -deffnm md_fullsys \
        -plumed ../plumed_mpi.dat \
        -ntomp $(nproc) \
        -gpu_id 0
done
```

Output: `md_fullsys.xtc` (full system), `md_fullsys.gro`, `md_fullsys.edr`.
Original `md.xtc` (protein-only) is untouched.

### 5. Harvest full-system trajectories

Add `md_fullsys.xtc` and `md_fullsys.edr` to the harvest script's
WALKER_FILES list, or harvest manually:

```bash
PROTEIN=1Q8K_10023
MACHINE=scan1
DEST=/shared/fleet/clean_fleet/results/${PROTEIN}

for w in 0 1 2 3 4; do
    rsync -avz --partial \
        ${MACHINE}:/opt/fleet-reprep1/staging/${PROTEIN}/walker_${w}/md_fullsys.xtc \
        ${DEST}/walker_${w}/
    rsync -avz --partial \
        ${MACHINE}:/opt/fleet-reprep1/staging/${PROTEIN}/walker_${w}/md_fullsys.edr \
        ${DEST}/walker_${w}/
done
```

### 6. Verify

```python
import struct
with open(f"{DEST}/walker_0/md_fullsys.xtc", "rb") as f:
    f.read(4)
    n = struct.unpack(">i", f.read(4))[0]
    assert n > 100000, f"Expected full system, got {n} atoms"
    print(f"OK: {n} atoms (full system)")
```

## What to do first (test on one protein)

Pick a small protein on scan2 (scan1 VPN status unknown).
Run steps 2-6 for just that one. Verify the .xtc has full-system
atoms. Copy to `tests/data/` on batcave. The extraction code can
then be tested end-to-end with water.

Check which small proteins are on scan2:
```bash
ssh scan2 "head -5 /opt/fleet-reprep1/scan2_reprep_results.txt"
```

## scan3 recovery

scan3 has 228 proteins. When recovered:
- Verify data integrity (md.gro, COLVAR, HILLS exist for all walkers)
- Same procedure as above for the re-run
- If disk is damaged, the .gro from results/ on batcave has the full
  system restart state — deploy back to scan3 and re-run from there

## Timeline

- Test one protein: ~1 hour (grompp + 10ns run + verify + harvest)
- Full fleet on 2 machines: ~8 days of wall time
- Full fleet on 3 machines (with scan3): ~6 days
- Harvest to batcave: ~2-3 days (network bound, ~680 GB per machine)

Total: ~10-12 days from start to all data on batcave.
Fits in the 40-day window with 28+ days to spare for extraction
and analysis.

## Post-harvest: filter to 15A shell on batcave (optional)

If batcave disk gets tight, strip distant water locally:

```bash
# For each protein/walker:
gmx_mpi select -s md_fullsys.tpr \
    -selrpos atom \
    -select "group Protein or (resname SOL and within 1.5 of group Protein)" \
    -os near_water.ndx

gmx_mpi trjconv -f md_fullsys.xtc -s md_fullsys.tpr \
    -n near_water.ndx -o md_nearwater.xtc
```

This reduces ~3 GB/protein to ~500 MB/protein. Do this AFTER
verifying the full-system extraction works end-to-end.
