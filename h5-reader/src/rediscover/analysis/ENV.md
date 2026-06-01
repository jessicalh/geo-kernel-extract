# Environment for the equivariant T2 fitter (e3nn) + the change-of-basis test

The equivariant fitter (`equiv_t2_e3nn.py`) and the pinned change-of-basis test
(`test_change_of_basis.py`) need **torch + e3nn**. Per the lead's decision
(MODEL_PLACEMENT_PROPOSAL.md §7 Q1; brief decision 1): run + validate in the
SYSTEM Python where e3nn 0.6.0 / torch 2.11.0+cu130 are already installed; PIN
the versions (`requirements-e3nn.txt`); do not burn time on a fresh torch build.

## The env used / validated

- Python: `/usr/bin/python3` (3.12), with user-site `~/.local/lib/python3.12/site-packages`.
- torch 2.11.0+cu130, e3nn 0.6.0 (declared equivariant dep). numpy 1.26.4,
  pandas 2.1.4, scipy 1.11.4, matplotlib 3.6.3 (system).
- The PySR scalar-distillation venv (`analysis/venv`, julia/PySR) is SEPARATE and
  unchanged; it has no torch/e3nn. Only `pysr_distill.py` uses it.

## The LD_LIBRARY_PATH gotcha (torch cu130 segfaults without it)

The cu130 torch wheel dlopens CUDA 13 shared libs from the `nvidia/*` wheels.
Without them on `LD_LIBRARY_PATH`, `import torch` SEGFAULTS (memory
`feedback_cuda_ld_path`). The libs live under the user site-packages:

    NV=$HOME/.local/lib/python3.12/site-packages/nvidia
    TORCHLIB=$HOME/.local/lib/python3.12/site-packages/torch/lib
    export LD_LIBRARY_PATH="$NV/cu13/lib:$NV/cudnn/lib:$NV/nccl/lib:$NV/cusparselt/lib:$NV/nvshmem/lib:$TORCHLIB:$LD_LIBRARY_PATH"

(`$NV/cu13/lib` carries libcudart.so.13, libcublas.so.13, libnvrtc.so.13, etc.;
confirmed present on this machine.)

## Run commands

Change-of-basis fixture test (derives + pins C, asserts orthogonality + the
Wigner-D equivariance round-trip vs the C++ library tensor):

    cd src/rediscover/analysis
    <export LD_LIBRARY_PATH as above>
    /usr/bin/python3 test_change_of_basis.py
    # or: /usr/bin/python3 -m pytest test_change_of_basis.py -q

Equivariant T2 fitter (gate: reproduce/beat frame-split T2 R²=0.467, |T2| r=0.756):

    REDISCOVER_OUT=/tmp/rediscover-out-v2 \
    /usr/bin/python3 equiv_t2_e3nn.py /tmp/rediscover-out-v2 --cross both

The fitter REQUIRES the emitted local-frame target NPY sidecar
`rediscover_ring_current_sources_target_local_T2.npy` in the out_dir (it fails
loud if absent — there is no Python projection fallback). The canonical
`/tmp/rediscover-out-v2/` currently has only the CSVs; point at a dir that also
has the NPY sidecars (e.g. a fresh extractor run, or `/tmp/rediscover-rebuild-npy/`
/ `/tmp/rdc-composed/` which carry the matching `*_target_local_T2.npy`), OR
re-run the extractor into v2 so the sidecar is regenerated alongside the CSVs.

After the first successful run, freeze the change-of-basis constant: run
`python3 change_of_basis.py`, paste the printed `_C_FROZEN` array into
`change_of_basis.py` so the model path never needs e3nn merely to load the matrix.

## Why a fresh agent could not run this here

The Claude Agent-tool sandbox in this environment denies `import torch` (and even
`import numpy`) — the same restriction that denies the compiler (STATE.md: the
author-with-Claude / build-with-codex split). So this pass AUTHORED the e3nn
model + the change-of-basis test + the script cleanups and proved the boundary by
grep; the RUN + the constant-freeze + the commit are the lead/codex half, with the
exact commands above.
