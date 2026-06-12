# BUILD_STATUS

Status: GREEN.

Built:
- Consolidated FP-native emit in `h5reader_extract --case consolidated_emit` / `--case fp_emit`.
- One-pass 720 static + 1P9J trajectory emission with global append-order `row_id`.
- `rows.csv`, `column_support.csv`, `sidecar_support.csv`, target sidecars, native-axis sidecars/maps, topology CSVs, e3nn sidecars, and JSON manifests.

Canonical inputs used:
- 720 root: `/shared/2026Thesis/shielding-calcsets/data/720-mutants/full720_20260609_180327`
- 1P9J LGS: `/shared/2026Thesis/shielding-calcsets/data/1p9j-work/extract_stride1_mopac_20260609_180327/md_20260610T103128Z.lgs`
- Code verification: canonical constants at `h5-reader/src/rediscover/CanonicalSpineInputs.h:33` and `h5-reader/src/rediscover/CanonicalSpineInputs.h:35`; hard root/pose guards at `h5-reader/src/rediscover/CanonicalSpineGuard.cpp:36`, `h5-reader/src/rediscover/CanonicalSpineGuard.cpp:53`, and `h5-reader/src/rediscover/CanonicalSpineGuard.cpp:103`; consolidated runner invokes them at `h5-reader/src/rediscover/ConsolidatedEmit.cpp:559` and `h5-reader/src/rediscover/ConsolidatedEmit.cpp:632`.

Run:
```bash
./build/linux-rwdi/h5reader_extract --case consolidated_emit \
  --root720 /shared/2026Thesis/shielding-calcsets/data/720-mutants/full720_20260609_180327 \
  --run /shared/2026Thesis/shielding-calcsets/data/1p9j-work/extract_stride1_mopac_20260609_180327/md_20260610T103128Z.lgs \
  --out /shared/2026Thesis/nmr-shielding-emit-build/output/consolidated_fp_emit_20260612T000000Z
```

Survival gate:
- PASS. The command exited 0 with `rows=1744962` and `scoped_fields=122`.
- Gate implementation/call: `h5-reader/src/rediscover/ConsolidatedEmit.cpp:1050` and `h5-reader/src/rediscover/ConsolidatedEmit.cpp:689`.
- Output rows: `rows.csv` has 1,744,963 lines including header.
- Support files: `column_support.csv` has 110 lines; `sidecar_support.csv` has 367 lines.
- Final sidecars: 248 `.npy` files, 0 `.rawtmp` files, output size 30G.
- 720 `not-produced-in-dataset` fields with 1P9J support: exactly 20.

The 20 1P9J-only fields reported `not-produced-in-dataset` on 720:
`coulomb_efg`, `coulomb_E`, `coulomb_E_backbone`, `coulomb_E_sidechain`, `coulomb_E_aromatic`, `coulomb_efg_backbone`, `coulomb_efg_sidechain`, `coulomb_efg_aromatic`, `coulomb_scalars`, `coulomb_aromatic_E_proj`, `coulomb_aromatic_n_src`, `water_efield`, `water_efield_first`, `water_efg`, `water_efg_first`, `water_shell_counts`, `hydration_shell`, `water_polarization`, `gromacs_energy`, `bonded_energy`.

Commits:
- `523aa29` Add consolidated emit scaffold
- `22e7d5d` Wire consolidated FP emit runner
