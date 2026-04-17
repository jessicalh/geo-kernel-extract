# analysis-speculative/

Speculative analysis planning for the fleet ensemble investigation.

Everything here is ideas under active assumption-testing. Nothing here
is settled. The extraction engine, the calibration, and the 28
realities ARE settled (see learn/). This directory is where we think
about what to do with 685 proteins x 5005 frames of geometry-only
extraction data, experimental RefDB shifts, and three months.

## Constraints (real, not aspirational)

- 4x RTX 5090 + 128 GB CPU RAM (primary machine)
- Possible short-term access to large-memory R machine (~2 TB)
- Python + R + LaTeX (no notebooks)
- nmr-extract SDK: geometry-only mode, ~1.3s/frame, ~10 MB/frame
- RefDB shifts for all 685 fleet proteins
- Must produce a working GNN model for the thesis
- Must produce analytical physics findings that the GNN corroborates

## What this is NOT

- Not a model-building plan (that's nmr-training/)
- Not calibration (that's learn/, it's done)
- Not extraction engineering (that's the C++ codebase)

## What this IS

Using the extraction as a physics observatory on 685 proteins worth
of conformational ensemble data. Finding signals. Testing whether
the physics the calibration validated on mechanical mutants shows up
in experimental shifts from real proteins. Discovering things we
didn't predict.
