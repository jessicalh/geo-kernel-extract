Per-Ring Data
=============

Sparse (P, 57) table of per-(atom, ring) pair contributions, plus a
(R, 10) ring geometry reference table.

RingContributions
-----------------

.. autoclass:: nmr_extract.RingContributions
   :members:
   :undoc-members:
   :exclude-members: COLS

Column layout::

   [0]     atom_index
   [1]     ring_index
   [2]     ring_type           0-7
   [3]     distance            Angstroms
   [4]     rho                 Angstroms
   [5]     z                   Angstroms (signed)
   [6]     theta               radians
   [7]     mcconnell_factor    (3cos^2 theta - 1) / r^3
   [8]     exp_decay           exp(-distance / 4.0)
   [9:18]  bs_G                SphericalTensor -- BS shielding kernel
   [18:27] hm_H                SphericalTensor -- HM raw integral (pure T2)
   [27:36] hm_G                SphericalTensor -- HM shielding kernel
   [36:45] pq_G                SphericalTensor
   [45:54] chi_G               SphericalTensor
   [54]    disp_scalar         1/r^6
   [55]    disp_contacts       vertex contact count
   [56]    gaussian_density    (placeholder)

RingGeometry
------------

.. autoclass:: nmr_extract.RingGeometry
   :members:
   :undoc-members:
   :exclude-members: COLS
