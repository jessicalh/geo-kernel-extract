Piece 3 postmortem, 2026-06-03
Run dir: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-piece3-loop3
Commit: HEAD (this atomic commit; exact hash printed by Codex after commit)
Drift: landed additive C++ Catalog/per_atom_substrate emit for missing classical mechanism channels.
Loader added: none; reused existing QtTrajectoryH5 H5 buffers and added guarded buffer accessors only.
Absent H5 slabs: hbond_nearest_dist, hbond_nearest_dir, hbond_is_donor, hbond_is_acceptor, eeq_coordination_number.
Fail-loud absent behavior: present=0 and NaN payloads; absent slab list recorded in manifest.
Build gate: PASS, cmake --build build/linux-gcc --target h5reader_extract h5reader_rediscover_tests.
Unit gate: PASS, h5reader_rediscover_tests 6/6.
Extraction gate: PASS, rows=558360, atoms=846, DFT frames=660.
R/uniqueness gate: PASS, dense row_id and unique (atom_index, original_frame_index).
Sidecar shape gate: PASS, row sidecars first dim 558360; driver_by_atom first dim 846.
New-channel flags gate: PASS, binary flags; present rows finite by C++ audit.
Axis gate: PASS, EEQ charge constant per atom; frame channels vary by nonzero ranges.
Oracle parity gate: PASS, old CSV prefix identical and unchanged sidecars byte-identical.
Classical parity note: features_classical grew 45->89; old 45-column prefix byte-identical.
Query parity gate: PASS, ring_identity and default query_results unchanged.
Backbone regression gate: PASS, backbone_audit.npy byte-identical.
H-bond shielding T2: present=558360, range=[-0.704945, 0.729606]; count range=[0, 9].
H-bond geometry/flags: present=0 for nearest distance, nearest direction, donor flag, acceptor flag.
Pi-quadrupole T2: present=558360, range=[-1.30225, 0.999791].
Dispersion T2: present=558360, range=[-0.0988106, 0.1163].
Haigh-Mallion T2: present=558360, range=[-1.11584, 1.40612].
Ring-susceptibility T2: present=558360, range=[-0.490185, 0.348028].
Water E-field: present=558360, range=[-3.55169, 3.38441].
Water shell counts: n_first present=558360 range=[0, 9]; n_second present=558360 range=[0, 24].
Hydration scalars: half-shell range=[0, 1]; dipole-cos range=[-0.999891, 0.999994].
SASA: present=558360, range=[0, 53.573]; normal present=558360, range=[-0.999998, 0.999999].
EEQ: charge present=558360 range=[-1.6913, 1.10947]; coordination_number present=0.
