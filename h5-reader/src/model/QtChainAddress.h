// QtChainAddress — typed wrapper for the (chain_id, residue_number,
// insertion_code) addressing triple.
//
// Per the no-strings discipline (notes/H5_READER_REWRITE_DESIGN_2026-05-23.md
// §2): addressing strings are genuine projection data (free-form PDB
// convention, multimeric chain identifiers like "AA"/"AB", insertion
// codes per PDB column 27), but they must NOT be compared directly
// like chemistry properties.
//
// The wrapper deliberately DISABLES operator==. Equality requires the
// explicit IsSameAddress() predicate so a future agent writing
// `if (residue_a.address == residue_b.address)` gets a compile error
// telling them to think about what comparison they want — by chain?
// by number? by full triple? The discipline is enforced by the
// language, not by code review.
//
// User decision §11.A (2026-05-23): "DELETED" — strongest enforcement.

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
