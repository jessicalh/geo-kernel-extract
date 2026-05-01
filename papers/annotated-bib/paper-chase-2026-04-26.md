# Paper-chase log — 2026-04-26

A chase session is one acquisition pass for a specific need: thesis section, bib slot, methods justification, etc. Each pass records scope, what was searched, what got in, what stayed out.

Convention going forward:
- One file per chase date: `paper-chase-YYYY-MM-DD.md` in `papers/annotated-bib/` (or wherever the chase served).
- Citations in Harvard (Cite Them Right) format.
- "Got in" entries name the local basename so the file maps cleanly to `references/<basename>.pdf`.
- Blocked or deferred entries name the DOI so they can be retried via institutional access.

---

## Scope of this pass

Triggered by the annotated-bibliography submission. Three target topics:

1. **r²SCAN effectiveness for NMR shielding in an experimental context** — methods-chapter justification for the functional choice; the corpus had no paper doing this directly.
2. **Geometric kernel and equivariant ML applications with full shielding tensor (T2 / SO(3) irrep structure)** — the corpus had Kondor 2025 in theoretical form but no primary research demonstrating tensor-output ML.
3. **MD + NMR experimental validation papers, especially protein-scale** — to complement Markwick 2010 with other primary research in the same lane.

Plus secondary items: experimental CSA tensor anchors, modern protein shift predictors (SHIFTX2 / UCBShift / ProCS15 / LEGOLAS), and a Sahakyan-Vendruscolo group methyl-shift companion paper.

## Searches run

- `molecular dynamics chemical shift prediction protein NMR experimental validation 2024 2025`
- `r2SCAN rSCAN meta-GGA NMR shielding chemical shift benchmark protein`
- `equivariant neural network NMR chemical shielding tensor anisotropy SO(3)`
- `protein chemical shielding tensor anisotropy CSA experimental measurement solid state 2024`
- `protein mutation chemical shift prediction BMRB machine learning ridge regression`
- `LEGOLAS protein NMR chemical shift prediction machine learning 2025`
- `"Robustelli" OR "Vendruscolo" molecular dynamics NMR chemical shift validation force field`

## Got in

All ingested into `references/` with chunks in `references-text/` and page renders in `references-images/`. Summary + keyword files (`references-meta/<basename>-{summary,keywords}.txt`) are not auto-generated — manual step per the scholarship workflow.

### r²SCAN for NMR

- **Yates, J.R. and Bartók, A.P. (2024)** 'Accurate predictions of chemical shifts with the rSCAN and r²SCAN mGGA exchange–correlation functionals', *Faraday Discussions*, 255, pp. 192–202. doi:10.1039/D4FD00142G. → `yates-bartok-2024-rscan-r2scan-chemical-shifts`. Open access (CC-BY).

### Equivariant / tensor-output ML for NMR shielding

- **Venetos, M.C., Wen, M. and Persson, K.A. (2023)** — *blocked here, see Out of reach.*
- **Ben Mahmoud, C., Rosset, L.A.M., Yates, J.R. and Deringer, V.L. (2024)** 'Graph-neural-network predictions of solid-state NMR parameters from spherical tensor decomposition', *Journal of Chemical Physics*, 163, 024118. doi:10.1063/5.0274240. → `ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors` (arXiv preprint version ingested; JCP published version in `_backup_2026-04-26/`).
- **Kellner, M., Holmes, J.B., Rodriguez-Madrid, R., Viscosi, F., Zhang, Y., Emsley, L. and Ceriotti, M. (2025)** 'A deep learning model for chemical shieldings in molecular organic solids including anisotropy' (ShiftML3.0), *Journal of Physical Chemistry Letters*. arXiv:2506.13146. → `kellner-2025-shiftml3-anisotropy-organic-solids`.
- **Bånkestad, M., Dorst, K.M., Widmalm, G. and Rönnols, J. (2023)** 'Carbohydrate NMR chemical shift predictions using E(3) equivariant graph neural networks', arXiv:2311.12657. → `bankestad-2023-carbohydrate-equivariant-gnn-shifts`.

### MD + NMR experimental, protein-scale

- **Robustelli, P., Stafford, K.A. and Palmer, A.G. III (2012)** 'Interpreting protein structural dynamics from NMR chemical shifts', *Journal of the American Chemical Society*, 134(14), pp. 6045–6056. doi:10.1021/ja300265w. → `robustelli-stafford-palmer-2012-protein-dynamics-shifts`.
- **Robustelli, P., Kohlhoff, K., Cavalli, A. and Vendruscolo, M. (2010)** 'Using NMR chemical shifts as structural restraints in molecular dynamics simulations of proteins', *Structure*, 18(8), pp. 923–933. doi:10.1016/j.str.2010.04.016. → `robustelli-cavalli-vendruscolo-2010-shifts-md-restraints`.
- **Li, D.-W. and Brüschweiler, R. (2010)** 'Certification of molecular dynamics trajectories with NMR chemical shifts', *Journal of Physical Chemistry Letters*, 1(1), pp. 246–248. doi:10.1021/jz9001345. → `li-bruschweiler-2010-md-trajectory-certification-shifts`.
- **Li, D.-W. and Brüschweiler, R. (2012)** 'PPM: a side-chain and backbone chemical shift predictor for the assessment of protein conformational ensembles', *Journal of Biomolecular NMR*, 54(3), pp. 257–265. doi:10.1007/s10858-012-9668-8. → `li-bruschweiler-2012-ppm-shift-predictor-ensembles`.
- **Karp, J.M., Erylimaz, E. and Cowburn, D. (2014)** 'Correlation of chemical shifts predicted by molecular dynamics simulations for partially disordered proteins', *Journal of Biomolecular NMR*, 61(1), pp. 35–45. doi:10.1007/s10858-014-9879-2. → `karp-2014-md-shifts-disordered-proteins`.

### Experimental CSA tensor anchors

- **Hartman, J.D., Day, G.M. and Beran, G.J.O. (2018)** 'Accurate ¹³C and ¹⁵N molecular crystal chemical shielding tensors from fragment-based electronic structure theory', *Solid State Nuclear Magnetic Resonance*, 96, pp. 10–18. doi:10.1016/j.ssnmr.2018.09.003. → `hartman-day-beran-2018-fragment-shielding-tensors-crystals`.

### Modern protein shift predictors (queue from earlier supplemental search)

- **Han, B., Liu, Y., Ginzinger, S.W. and Wishart, D.S. (2011)** 'SHIFTX2: significantly improved protein chemical shift prediction', *Journal of Biomolecular NMR*, 50(1), pp. 43–57. doi:10.1007/s10858-011-9478-4. → `han-2011-shiftx2-protein-chemical-shift-prediction`.
- **Larsen, A.S., Bratholm, L.A., Christensen, A.S., Channir, M. and Jensen, J.H. (2015)** 'ProCS15: a DFT-based chemical shift predictor for backbone and Cβ atoms in proteins', *PeerJ*, 3, e1344. doi:10.7717/peerj.1344. → `larsen-2015-procs15-dft-chemical-shift-predictor`.
- **Bratholm, L.A., Christensen, A.S., Hamelryck, T. and Jensen, J.H. (2017)** 'Protein structure refinement using a quantum mechanics-based chemical shielding predictor', *Chemical Science*, 8(3), pp. 2061–2072. doi:10.1039/C6SC04344E. → `bratholm-2017-procs15-qm-chemical-shielding-refinement`.
- **Li, J., Bennett, K.C., Liu, Y., Martin, M.V. and Head-Gordon, T. (2020)** 'Accurate prediction of chemical shifts for aqueous protein structure on "real-world" data', *Chemical Science*, 11(12), pp. 3180–3191. doi:10.1039/C9SC06561J. → `li-2020-ucbshift-real-world-data`.
- **Ptaszek, A.L., Li, J., Konrat, R., Platzer, G. and Head-Gordon, T. (2024)** 'UCBShift 2.0: bridging the gap from backbone to side-chain protein chemical shift prediction for protein structures', *Journal of the American Chemical Society*, 146(46), pp. 31733–31745. doi:10.1021/jacs.4c10474. → `ptaszek-2024-ucbshift-2-side-chain-prediction`.
- **Yang, Z., Chakraborty, M. and White, A.D. (2021)** 'Predicting chemical shifts with graph neural networks', *Chemical Science*, 12(32), pp. 10802–10809. doi:10.1039/D1SC01895G. → `yang-white-2021-gnn-chemical-shifts`.
- **Zhu, H., Hu, L., Yang, Y. and Chen, Z. (2025)** 'A novel approach to protein chemical shift prediction from sequences using a protein language model', *Digital Discovery*, 4(1), pp. 331–337. doi:10.1039/D4DD00367E. → `zhu-2025-plm-cs-sequence-chemical-shifts`.
- **Klukowski, P., Riek, R. and Güntert, P. (2022)** 'Rapid protein assignments and structures from raw NMR spectra with the deep learning technique ARTINA', *Nature Communications*, 13, 6151. doi:10.1038/s41467-022-33879-5. → `klukowski-2022-artina-raw-nmr-deep-learning`.
- **Paruzzo, F.M., Hofstetter, A., Musil, F., De, S., Ceriotti, M. and Emsley, L. (2018)** 'Chemical shifts in molecular solids by machine learning', *Nature Communications*, 9, 4501. doi:10.1038/s41467-018-06972-x. → `paruzzo-2018-shiftml-molecular-solids`.

### Sahakyan-Vendruscolo methyl shifts companion

- **Sahakyan, A.B., Vranken, W.F., Cavalli, A. and Vendruscolo, M. (2011)** 'Structure-based prediction of methyl chemical shifts in proteins', *Journal of Biomolecular NMR*, 50(4), pp. 331–346. doi:10.1007/s10858-011-9524-2. → `sahakyan-2011-methyl-chemical-shifts-proteins`.

## Out of reach (this session)

PMC anti-bot and ChemRxiv direct-PDF rejected the bot user-agent from this environment. Try via Birkbeck institutional access.

- **Venetos, M.C., Wen, M. and Persson, K.A. (2023)** 'Machine learning full NMR chemical shift tensors of silicon oxides with equivariant graph neural networks', *Journal of Physical Chemistry A*, 127(12), pp. 2388–2398. doi:10.1021/acs.jpca.2c07530. PMC10026072.
- **Darrows, M.Y., Kodituwakku, D., Xue, J., Pickering, I., Terrel, N.S. and Roitberg, A.E. (2025)** 'LEGOLAS: a machine learning method for rapid and accurate predictions of protein NMR chemical shifts', *Journal of Chemical Theory and Computation*, 21(8), pp. 4266–4275. doi:10.1021/acs.jctc.5c00026. ChemRxiv mirror at https://chemrxiv.org/engage/chemrxiv/article-details/6776ea28fa469535b991c521.

## Backups (for review)

`references/incoming/_backup_2026-04-26/` holds the two duplicates that did not need re-ingesting:

- `024118_1_5.0274240__JCP_published_version_of_ben-mahmoud-2024.pdf` — JCP final version of the spherical-tensor-decomposition paper. Currently the arXiv preprint is what's ingested. If you'd rather cite the JCP version as canonical, swap by re-ingesting from this backup with the same basename.
- `d4fd00142g__duplicate_of_yates-bartok-2024.pdf` — exact duplicate of the already-ingested r²SCAN paper.

## Notes for the next chase

- PMC PDFs reject this environment's bot user-agent with a "Preparing to download" interstitial. Europe PMC mirror (`europepmc.org/backend/ptpmcrender.fcgi?accid=PMC<ID>&blobtype=pdf`) also rejects. Birkbeck institutional access is the practical path; or have Jessica fetch via browser and drop into `references/incoming/`.
- ArXiv direct PDFs (`arxiv.org/pdf/<id>`) work cleanly.
- RSC direct PDFs work with a browser-like User-Agent and a Referer header.
- Background curl jobs in a chained `cd && curl &` only respect the `cd` for the foreground process; subsequent backgrounded `curl &` calls run from the original cwd. Use `cd; curl & curl & wait` (separate statements) or pre-build full paths.
- The ingest script (`scripts/references/ingest_pdf.sh`) ignores the input filename and writes outputs under the basename argument; files don't need pre-renaming, just the right basename on call.
