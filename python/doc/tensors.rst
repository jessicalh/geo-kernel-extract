Tensor Types
============

All tensor wrappers expose:

- ``.data`` -- ``numpy.ndarray``
- ``.torch()`` -- ``torch.Tensor`` (zero-copy from numpy)
- ``.irreps`` -- ``e3nn.o3.Irreps``

SphericalTensor
---------------

.. autoclass:: nmr_extract.SphericalTensor
   :members:
   :undoc-members:
   :exclude-members: IRREPS

.. autoclass:: nmr_extract.ShieldingTensor
   :show-inheritance:

.. autoclass:: nmr_extract.EFGTensor
   :show-inheritance:

VectorField
-----------

.. autoclass:: nmr_extract.VectorField
   :members:
   :undoc-members:
   :exclude-members: IRREPS

Per-type decompositions
-----------------------

.. autoclass:: nmr_extract.PerRingTypeT0
   :members:
   :undoc-members:
   :exclude-members: IRREPS

.. autoclass:: nmr_extract.PerRingTypeT2
   :members:
   :undoc-members:
   :exclude-members: IRREPS

.. autoclass:: nmr_extract.PerBondCategoryT2
   :members:
   :undoc-members:
   :exclude-members: IRREPS

MOPAC scalars
-------------

.. autoclass:: nmr_extract.MopacScalars
   :members:
   :undoc-members:

.. autoclass:: nmr_extract.MopacGlobal
   :members:
   :undoc-members:
