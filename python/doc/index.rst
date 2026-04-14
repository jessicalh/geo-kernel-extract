nmr_extract
===========

Typed reader for NMR shielding tensor extraction data.  Loads NPY files
produced by the C++ extractor into Python objects with ``e3nn.o3.Irreps``,
``torch.Tensor``, and ``numpy`` access.

Install::

   pip install -e /path/to/nmr-shielding/python

Quick start::

   from nmr_extract import load

   p = load("path/to/extraction/directory")
   p.biot_savart.shielding.T2      # ndarray (N, 5)
   p.biot_savart.shielding.irreps  # Irreps("1x0e+1x1o+1x2e")
   p.biot_savart.shielding.torch() # torch.Tensor (N, 9)
   p.ring_contributions.bs.T2      # ndarray (P, 5), per-ring

.. toctree::
   :maxdepth: 2

   loading
   trajectory
   tensors
   ring
   groups
   enums
   catalog
