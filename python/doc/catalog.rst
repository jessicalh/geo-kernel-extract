Catalog
=======

The format contract between C++ and Python.  One entry per NPY file
the extractor can produce.

.. autoclass:: nmr_extract.ArraySpec
   :members:

``CATALOG`` is a ``dict[str, ArraySpec]`` mapping every NPY filename
stem to its metadata.  74 entries covering identity, ring calculators,
bond calculators, MOPAC, APBS, Orca DFT, AIMNet2, water field,
hydration geometry, EEQ, and mutation delta arrays.

::

   from nmr_extract import CATALOG
   for stem, spec in CATALOG.items():
       print(f"{stem:30s} {spec.group:20s} required={spec.required}")
