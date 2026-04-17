"""
Trajectory analysis H5 reader SDK.

Usage:
    from trajectory_reader_v2 import load_analysis
    protein = load_analysis("path/to/analysis.h5")
"""

from .loader import load_analysis
from .model import AnalysisProtein
from .enums import (
    Element, AtomRole, Hybridisation, BondOrder, BondCategory,
    RingTypeIndex, AminoAcid,
)

__all__ = [
    "load_analysis",
    "AnalysisProtein",
    "Element", "AtomRole", "Hybridisation", "BondOrder", "BondCategory",
    "RingTypeIndex", "AminoAcid",
]
