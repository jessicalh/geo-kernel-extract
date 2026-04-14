Trajectory Data
===============

Load the H5 master file produced by ``--trajectory`` mode.

::

   from nmr_extract import load_trajectory
   traj = load_trajectory("output/trajectory.h5")

   traj.positions          # (T, N, 3) per-frame xyz
   traj.rollup.bs_T0.mean  # (N,) ring current isotropic mean
   traj.bonds.length_mean  # (B,) mean bond length

.. autofunction:: nmr_extract.load_trajectory

.. autoclass:: nmr_extract.TrajectoryData
   :members:
   :undoc-members:

.. autoclass:: nmr_extract.TrajectoryRollup
   :members:
   :undoc-members:

.. autoclass:: nmr_extract.BondRollup
   :members:
   :undoc-members:
