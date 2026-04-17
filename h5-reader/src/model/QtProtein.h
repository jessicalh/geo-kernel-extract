// QtProtein — identity and topology of one protein, independent of
// geometry. Owns the rings as unique_ptr<QtRing> so the polymorphic
// hierarchy works. Atoms, residues, bonds are value-type vectors.
//
// Constructed by io::QtProteinLoader at H5 load time. Non-copyable,
// non-movable — QtFrames and QtConformation hold raw const pointers
// back at the protein for topology lookups; the pointers must not
// dangle. Same pattern as nmr::Protein in src/Protein.h.

#pragma once

#include "QtAtom.h"
#include "QtBond.h"
#include "QtResidue.h"
#include "QtRing.h"

#include <QString>
#include <cstddef>
#include <memory>
#include <vector>

// The loader needs access to the protein's private vectors to populate
// them during H5 decode. Forward-declare here so the friend line below
// can name the type without pulling the whole header in.
namespace h5reader::io { class QtProteinLoader; }

namespace h5reader::model {

class QtProtein {
public:
    QtProtein() = default;
    ~QtProtein() = default;

    QtProtein(const QtProtein&)            = delete;
    QtProtein& operator=(const QtProtein&) = delete;
    QtProtein(QtProtein&&)                 = delete;
    QtProtein& operator=(QtProtein&&)      = delete;

    // ----- Identity -----
    const QString& proteinId() const { return proteinId_; }
    void setProteinId(const QString& id) { proteinId_ = id; }

    // ----- Atoms -----
    size_t atomCount() const { return atoms_.size(); }
    const QtAtom& atom(size_t i) const { return atoms_[i]; }
    const std::vector<QtAtom>& atoms() const { return atoms_; }

    // ----- Residues -----
    size_t residueCount() const { return residues_.size(); }
    const QtResidue& residue(size_t i) const { return residues_[i]; }
    const std::vector<QtResidue>& residues() const { return residues_; }

    // ----- Bonds -----
    size_t bondCount() const { return bonds_.size(); }
    const QtBond& bond(size_t i) const { return bonds_[i]; }
    const std::vector<QtBond>& bonds() const { return bonds_; }

    // ----- Rings -----
    size_t ringCount() const { return rings_.size(); }
    const QtRing& ring(size_t i) const { return *rings_[i]; }
    const std::vector<std::unique_ptr<QtRing>>& rings() const { return rings_; }

    // Ring counts by type — useful for inventory printing.
    int ringCountByType(RingTypeIndex t) const;

private:
    // QtProteinLoader populates these private vectors at H5 load time.
    // Everyone else reads via the const accessors above.
    friend class ::h5reader::io::QtProteinLoader;

    QString                               proteinId_;
    std::vector<QtAtom>                   atoms_;
    std::vector<QtResidue>                residues_;
    std::vector<QtBond>                   bonds_;
    std::vector<std::unique_ptr<QtRing>>  rings_;
};

}  // namespace h5reader::model
