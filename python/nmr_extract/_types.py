"""Physical enums matching the C++ RingTypeIndex and BondCategory."""

from enum import IntEnum


class RingType(IntEnum):
    PHE = 0            # Phe benzene
    TYR = 1            # Tyr phenol
    TRP_benzene = 2    # Trp 6-ring
    TRP_pyrrole = 3    # Trp 5-ring
    TRP_perimeter = 4  # Trp indole 9-atom perimeter
    HIS = 5            # His imidazole (neutral, N-delta)
    HID = 6            # Hid imidazole (doubly protonated)
    HIE = 7            # Hie imidazole (neutral, N-epsilon)


N_RING_TYPES = 8


class BondCategory(IntEnum):
    backbone_total = 0
    sidechain_total = 1
    aromatic_total = 2
    CO_nearest = 3
    CN_nearest = 4


N_BOND_CATEGORIES = 5
