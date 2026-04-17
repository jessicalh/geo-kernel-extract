"""
Typed enums mirroring the C++ Types.h definitions exactly.

Values are the C++ enum underlying values. These are what the H5 stores.

IMPORTANT: Element is stored in the H5 as enum index (H=0, C=1, N=2, O=3, S=4),
NOT as atomic number, despite the spec text saying "1=H, 6=C, 7=N, 8=O".
The actual data uses enum ordinal values.
"""

import enum


class Element(enum.IntEnum):
    """C++ enum class Element { H, C, N, O, S, Unknown }"""
    H = 0
    C = 1
    N = 2
    O = 3
    S = 4
    Unknown = 5

    @property
    def symbol(self) -> str:
        _symbols = {0: "H", 1: "C", 2: "N", 3: "O", 4: "S", 5: "?"}
        return _symbols[self.value]

    @property
    def atomic_number(self) -> int:
        _z = {0: 1, 1: 6, 2: 7, 3: 8, 4: 16, 5: 0}
        return _z[self.value]

    @property
    def covalent_radius(self) -> float:
        _r = {0: 0.31, 1: 0.76, 2: 0.71, 3: 0.66, 4: 1.05, 5: 0.0}
        return _r[self.value]

    @property
    def electronegativity(self) -> float:
        _en = {0: 2.20, 1: 2.55, 2: 3.04, 3: 3.44, 4: 2.58, 5: 0.0}
        return _en[self.value]


class AtomRole(enum.IntEnum):
    """
    C++ enum class AtomRole — NMR-relevant classification.

    Heavy backbone: 0-3
    Heavy sidechain: 4-9
    Hydrogen: 10-15
    Unknown: 16
    """
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

    @property
    def is_backbone(self) -> bool:
        return self.value <= 3

    @property
    def is_hydrogen(self) -> bool:
        return 10 <= self.value <= 15


class Hybridisation(enum.IntEnum):
    """C++ enum class Hybridisation { sp, sp2, sp3, Unassigned }"""
    sp = 0
    sp2 = 1
    sp3 = 2
    Unassigned = 3


class BondOrder(enum.IntEnum):
    """C++ enum class BondOrder"""
    Single = 0
    Double = 1
    Triple = 2
    Aromatic = 3
    Peptide = 4
    Unknown = 5


class BondCategory(enum.IntEnum):
    """C++ enum class BondCategory — finer than order, needed for McConnell"""
    PeptideCO = 0
    PeptideCN = 1
    BackboneOther = 2
    SidechainCO = 3
    Aromatic = 4
    Disulfide = 5
    SidechainOther = 6
    Unknown = 7


class RingTypeIndex(enum.IntEnum):
    """
    C++ enum class RingTypeIndex — 8 ring types.

    The per-type arrays in ring_current/ use these as the last axis index.
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
    def type_name(self) -> str:
        _names = {
            0: "PHE", 1: "TYR", 2: "TRP6", 3: "TRP5",
            4: "TRP9", 5: "HIS", 6: "HID", 7: "HIE",
        }
        return _names[self.value]

    @property
    def ring_size(self) -> int:
        if self.value in (0, 1, 2):
            return 6
        if self.value in (3, 5, 6, 7):
            return 5
        if self.value == 4:
            return 9
        return 0


class AminoAcid(enum.IntEnum):
    """
    C++ enum class AminoAcid — the 20 standard amino acids.

    Not stored directly in the H5 (residues/ uses three-letter string names),
    but needed for typed identity in the Python model.
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

    @property
    def is_aromatic(self) -> bool:
        return self in (AminoAcid.PHE, AminoAcid.TYR, AminoAcid.TRP, AminoAcid.HIS)


# Lookup table: three-letter code string -> AminoAcid enum
_AA_FROM_CODE: dict[str, AminoAcid] = {}
for _aa in AminoAcid:
    _AA_FROM_CODE[_aa.name] = _aa
# Add common protonation variant names
_AA_FROM_CODE["HID"] = AminoAcid.HIS
_AA_FROM_CODE["HIE"] = AminoAcid.HIS
_AA_FROM_CODE["HIP"] = AminoAcid.HIS
_AA_FROM_CODE["ASH"] = AminoAcid.ASP
_AA_FROM_CODE["GLH"] = AminoAcid.GLU
_AA_FROM_CODE["CYX"] = AminoAcid.CYS
_AA_FROM_CODE["CYM"] = AminoAcid.CYS
_AA_FROM_CODE["LYN"] = AminoAcid.LYS


def amino_acid_from_code(code: str) -> AminoAcid:
    """Resolve a three-letter residue name to AminoAcid enum."""
    return _AA_FROM_CODE.get(code, AminoAcid.Unknown)
