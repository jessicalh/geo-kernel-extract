"""
Typed enums matching the C++ enum values in src/Types.h.

Element uses ATOMIC NUMBERS as stored in the H5 (atoms/element).
All other enums use the C++ ordinal values.
"""

from enum import IntEnum


class Element(IntEnum):
    """Protein NMR elements. Values are atomic numbers (H5 storage format)."""
    H = 1
    C = 6
    N = 7
    O = 8
    S = 16

    @property
    def symbol(self) -> str:
        return self.name

    @property
    def covalent_radius(self) -> float:
        return _COVALENT_RADII[self]

    @property
    def electronegativity(self) -> float:
        return _ELECTRONEGATIVITIES[self]


_COVALENT_RADII = {
    Element.H: 0.31,
    Element.C: 0.76,
    Element.N: 0.71,
    Element.O: 0.66,
    Element.S: 1.05,
}

_ELECTRONEGATIVITIES = {
    Element.H: 2.20,
    Element.C: 2.55,
    Element.N: 3.04,
    Element.O: 3.44,
    Element.S: 2.58,
}


class AtomRole(IntEnum):
    """NMR-relevant atom classification. Values match C++ enum class AtomRole."""
    BackboneN = 0
    BackboneCA = 1
    BackboneC = 2
    BackboneO = 3
    SidechainC = 4
    SidechainN = 5
    SidechainO = 6
    SidechainS = 7
    AromaticC = 8
    AromaticN = 9
    AmideH = 10
    AlphaH = 11
    MethylH = 12
    AromaticH = 13
    HydroxylH = 14
    OtherH = 15
    Unknown = 16


class Hybridisation(IntEnum):
    """Values match C++ enum class Hybridisation."""
    sp = 0
    sp2 = 1
    sp3 = 2
    Unassigned = 3


class BondOrder(IntEnum):
    """Values match C++ enum class BondOrder."""
    Single = 0
    Double = 1
    Triple = 2
    Aromatic = 3
    Peptide = 4
    Unknown = 5


class BondCategory(IntEnum):
    """Values match C++ enum class BondCategory."""
    PeptideCO = 0
    PeptideCN = 1
    BackboneOther = 2
    SidechainCO = 3
    Aromatic = 4
    Disulfide = 5
    SidechainOther = 6
    Unknown = 7


class RingTypeIndex(IntEnum):
    """
    8 aromatic ring types. Values match C++ enum class RingTypeIndex.
    """
    PheBenzene = 0
    TyrPhenol = 1
    TrpBenzene = 2
    TrpPyrrole = 3
    TrpPerimeter = 4
    HisImidazole = 5
    HidImidazole = 6
    HieImidazole = 7

    @property
    def short_name(self) -> str:
        return _RING_TYPE_NAMES[self]

    @property
    def ring_size(self) -> int:
        return _RING_SIZES[self]


_RING_TYPE_NAMES = {
    RingTypeIndex.PheBenzene: "PHE",
    RingTypeIndex.TyrPhenol: "TYR",
    RingTypeIndex.TrpBenzene: "TRP6",
    RingTypeIndex.TrpPyrrole: "TRP5",
    RingTypeIndex.TrpPerimeter: "TRP9",
    RingTypeIndex.HisImidazole: "HIS",
    RingTypeIndex.HidImidazole: "HID",
    RingTypeIndex.HieImidazole: "HIE",
}

_RING_SIZES = {
    RingTypeIndex.PheBenzene: 6,
    RingTypeIndex.TyrPhenol: 6,
    RingTypeIndex.TrpBenzene: 6,
    RingTypeIndex.TrpPyrrole: 5,
    RingTypeIndex.TrpPerimeter: 9,
    RingTypeIndex.HisImidazole: 5,
    RingTypeIndex.HidImidazole: 5,
    RingTypeIndex.HieImidazole: 5,
}


class AminoAcid(IntEnum):
    """
    20 standard amino acids. Values match C++ enum class AminoAcid.
    """
    ALA = 0
    ARG = 1
    ASN = 2
    ASP = 3
    CYS = 4
    GLN = 5
    GLU = 6
    GLY = 7
    HIS = 8
    ILE = 9
    LEU = 10
    LYS = 11
    MET = 12
    PHE = 13
    PRO = 14
    SER = 15
    THR = 16
    TRP = 17
    TYR = 18
    VAL = 19
    Unknown = 20

    @property
    def three_letter_code(self) -> str:
        return self.name


# Lookup from three-letter code string to AminoAcid enum.
_AMINO_ACID_FROM_CODE = {aa.name: aa for aa in AminoAcid}


def amino_acid_from_code(code: str) -> AminoAcid:
    """Convert a three-letter residue name (e.g. 'ALA') to AminoAcid enum."""
    return _AMINO_ACID_FROM_CODE.get(code, AminoAcid.Unknown)
