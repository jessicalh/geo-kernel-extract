#pragma once
//
// IupacAtomName: a typed wrapper around an IUPAC atom-name string.
//
// The string is implementation detail. The object is the model. Comparison
// happens through the typed `operator==`, never through bare `strcmp` or
// `std::string ==` in calculator / result / loader / stamper code. This
// extends the typed-boundary discipline (PATTERNS.md): we still cross from
// PDB-string to typed identity at the loader, but the typed identity is now
// a first-class value type with its own equality, hash, and ordering, not
// a raw string member.
//
// Source convention: Markley et al. 1998, J. Biomol. NMR 12:1-23.
// Names are the IUPAC-IUBMB-IUPAB recommendations (with PDB/BMRB's "H" still
// accepted as a synonym for "HN" at the PDB-loading boundary).
//
// Storage is std::string — same memory cost as the previous bare-string
// field, plus typed equality and ordering. const char* literals convert
// implicitly so the AminoAcidType table populates without IupacAtomName{...}
// noise.
//

#include <functional>
#include <ostream>
#include <string>
#include <utility>

namespace nmr {

class IupacAtomName {
public:
    IupacAtomName() = default;

    // Implicit from string literals — keeps AminoAcidType table entries terse:
    //   { "HB2", E::H, false, ... }
    /*implicit*/ IupacAtomName(const char* name)
        : name_(name ? name : "") {}

    // Implicit from std::string — for runtime parsing at the loader boundary:
    //   atom->iupac_name = parsed_pdb_name;
    /*implicit*/ IupacAtomName(std::string name)
        : name_(std::move(name)) {}

    // Underlying string accessor for diagnostics / serialisation. Calculator
    // code SHOULD NOT call this and then re-do string compares — use the
    // typed `==` instead.
    const std::string& AsString() const { return name_; }
    const char*        c_str()    const { return name_.c_str(); }
    bool               empty()    const { return name_.empty(); }
    size_t             size()     const { return name_.size(); }

    // Typed equality. The string compare lives here — once, not scattered.
    bool operator==(const IupacAtomName& other) const {
        return name_ == other.name_;
    }
    bool operator!=(const IupacAtomName& other) const {
        return !(*this == other);
    }

    // Ordering for std::map keys.
    bool operator<(const IupacAtomName& other) const {
        return name_ < other.name_;
    }

private:
    std::string name_;
};

// Stream output for diagnostics, gtest EXPECT_EQ failure messages, etc.
inline std::ostream& operator<<(std::ostream& os, const IupacAtomName& n) {
    return os << n.AsString();
}

}  // namespace nmr

namespace std {
template <>
struct hash<nmr::IupacAtomName> {
    size_t operator()(const nmr::IupacAtomName& n) const noexcept {
        return std::hash<std::string>{}(n.AsString());
    }
};
}  // namespace std
