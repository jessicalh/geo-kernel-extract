"""
Analysis trajectory reader for NMR shielding H5 files.

Provides a typed Python object model mirroring the C++ Protein/Atom/Residue/Ring/Bond
hierarchy, populated from the analysis H5 written by nmr_extract --trajectory --analysis.

Usage:
    from trajectory_reader import load_analysis, Element

    protein = load_analysis("path/to/md_analysis.h5")
    for atom in protein.atoms:
        if atom.element == Element.C:
            print(atom.pdb_name, atom.ring_current.bs_T0_per_type[:, 0])
"""

from .enums import (
    Element,
    AtomRole,
    Hybridisation,
    BondOrder,
    BondCategory,
    RingTypeIndex,
    AminoAcid,
)
from .model import (
    AnalysisProtein,
    AnalysisAtom,
    AnalysisResidue,
    AnalysisRing,
    AnalysisBond,
    TimeSeries,
    RingCurrentData,
    EfgData,
    BondAnisoData,
    QuadrupoleData,
    DispersionData,
    HBondData,
    SasaData,
    WaterData,
    ChargesData,
    AimNet2EmbeddingData,
    PredictionsData,
    ProjectionsData,
    DihedralData,
    DsspData,
    FrameView,
)
from .loader import load_analysis

__all__ = [
    "load_analysis",
    # Enums
    "Element",
    "AtomRole",
    "Hybridisation",
    "BondOrder",
    "BondCategory",
    "RingTypeIndex",
    "AminoAcid",
    # Model
    "AnalysisProtein",
    "AnalysisAtom",
    "AnalysisResidue",
    "AnalysisRing",
    "AnalysisBond",
    "TimeSeries",
    "RingCurrentData",
    "EfgData",
    "BondAnisoData",
    "QuadrupoleData",
    "DispersionData",
    "HBondData",
    "SasaData",
    "WaterData",
    "ChargesData",
    "AimNet2EmbeddingData",
    "PredictionsData",
    "ProjectionsData",
    "DihedralData",
    "DsspData",
    "FrameView",
]
