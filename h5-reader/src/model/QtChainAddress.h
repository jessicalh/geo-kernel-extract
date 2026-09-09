// QtChainAddress — typed wrapper for the (chain_id, residue_number,
// insertion_code) addressing triple.
//
// Addressing strings are genuine projection data: PDB chain identifiers can
// be multi-character and insertion codes are free-form identifiers. They must
// not be compared as chemistry properties.
//
// The wrapper deliberately DISABLES operator==. Equality requires the
// explicit IsSameAddress() predicate so a caller writing
// `if (residue_a.address == residue_b.address)` gets a compile error
// telling them to think about what comparison they want — by chain?
// by number? by full triple? The type system enforces the distinction.

#pragma once

#include <QString>

namespace h5reader::model {

struct QtChainAddress {
    QString chainId;        // free-form chain identifier; can be multi-char
    int residueNumber = 0;  // PDB residue sequence number
    QString insertionCode;  // PDB column 27 ("" or "A", "B", ...)

    // The ONE typed equality predicate. Explicit by design — no
    // operator== silent comparison.
    bool IsSameAddress(const QtChainAddress& other) const {
        return chainId == other.chainId && residueNumber == other.residueNumber && insertionCode == other.insertionCode;
    }

    // operator== is DELETED on purpose — catches "snobol-style"
    // misuse at compile time. Any code wanting equality must call
    // IsSameAddress() and thus actively chose to compare addresses.
    bool operator==(const QtChainAddress&) const = delete;
    bool operator!=(const QtChainAddress&) const = delete;
};

}  // namespace h5reader::model
