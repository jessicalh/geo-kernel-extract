Piece3B final run: /tmp/rediscover-runs/2026-06-03-per-atom-substrate-piece3b-final
Commits: chunk1 9965c2b target split + ring paths; chunk2 1d0f56a method paths; chunk3 5a39288 hbond + conditioners.
Source location pass: read ../python/nmr_extract/_catalog.py/QtFieldCatalog stems; DSSP raw partners from /trajectory/dssp8_time_series; ORCA targets from BuildTarget raw total/dia/para.
No producer, CMake, or ctest edits. C++ emitted all values; per-source bond/ring data folded in memory and discarded.
Rows/axis gate: 846 atoms x 660 DFT frames = 558360 rows; dense row_id and unique (atom_index,original_frame_index): PASS.
Size gate: 431 appended float64 columns = 1.793 GiB append; final run 3235370959 bytes. PASS, no hidden source/pair axis.
Oracle parity: prior rows, old NPYs, target_T0/T2, ring_identity, query_results, and backbone_audit byte-identical at each chunk. PASS.
T1 frame: frame_verified; max_angle_deg=0.00023711862716899507, max_rmsd_A=0.0005132701081822394.
Absent located legacy dense slabs: hbond_nearest_dist, hbond_nearest_dir, hbond_is_donor, hbond_is_acceptor. Replacement source hbond_scalars is present.
Excluded intentionally: coulomb_*, delta_*/wt_/mut_*, mopac_global, tripeptide_bb_shielding.
target_T1 558360[-85.7885,114.0645]; target_dia_T0 558360[-795.166,1752.0627]; target_para_T0 558360[-1495.5123,918.8117].
target_dia_T1 558360[-1176.28,1075.0165]; target_dia_T2 558360[-1770.4165,1726.0222].
target_para_T1 558360[-1082.467,1199.4835]; target_para_T2 558360[-2157.9944,2135.6672].
bs_shielding 558360[-1.3517,1.6906]; bs_per_type_T0 558360[-0.4026,0.06564]; bs_per_type_T2 558360[-0.6961,0.8392].
bs_total_B 558360[-2.0626e-06,2.3180e-06]; bs_ring_counts 558360[0,10].
hm_shielding 558360[-1.1158,1.4061]; hm_per_type_T0 558360[-0.3332,0.07793]; hm_per_type_T2 558360[-0.5699,0.6734].
ringchi_shielding_full 558360[-0.4902,0.4783]; pq_per_type_T0 558360[-0.01317,0.1772]; pq_per_type_T2 558360[-1.0187,0.68895].
disp_per_type_T0 558360[0,0.05723]; disp_per_type_T2 558360[-0.07904,0.08578].
mc_category_T2 558360[-1.2935,1.5022]; mc_scalars 558360[-0.4097,9.4350].
mopac_mc_category_T2 558360[-1.5509,1.4959]; mopac_mc_scalars 558360[-0.5756,9.4350].
mopac_bond_orders_aggregate 291720[0.188,10]; mopac_scalars 558360[-0.946882,5.07024].
mopac_coulomb_E 558360[-8.8195,8.7364]; mopac_coulomb_efg_backbone 558360[-21.9263,20.7756].
mopac_coulomb_efg_aromatic 558360[-14.1757,12.9593]; mopac_coulomb_scalars 558360[-9.3043,9.3117].
aimnet2_efg 558360[-19.8528,17.4840]; aimnet2_efg_aromatic 558360[-11.0695,9.6624]; aimnet2_efg_backbone 558360[-13.9893,12.2718].
water_efg 558360[-4.9581,6.2133]; water_efield_first 558360[-2.9659,3.1243]; eeq_cn 558360[0.82618,3.68238].
larsen_1pHB_T2 558360[-98.1275,107.1526]; larsen_2pHB_T2 558360[-88.1157,97.9731].
larsen_1pHaB_T2 558360[-28.1396,30.5940]; larsen_2pHaB_T2 558360[-91.0316,87.8694]; larsen_water 558360[0,2.07].
hbond_scalars 558360[-0.4177,11.4324]; dssp_chemical_flags 558360[0,1]; dssp_hbond_energy 558360[-3.968,0].
dssp_raw_hbond_backup 141240[-3.968,53]; dssp_ss8 558360[0,1]; dssp_chi 558360[-1,1].
omega_actual 545160[-3.14157,3.14157]; pyramidalization 558360[-0.27929,0.31701]; ring_geometry 45540[-0.99985,85.67966].
