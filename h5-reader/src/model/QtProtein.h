// QtProtein — identity + topology of one protein, mirroring
// `nmr::Protein`. Non-copyable + non-movable (back-pointers from
// a Conformation must never dangle).
//
// Owns:
//   - QtAtom[] (substrate-typed identity per atom)
//   - QtAtomNames[] (parallel projection layer for display names)
//   - QtResidue[] (substrate-typed identity per residue)
//   - QtTopology (wraps bonds + rings + ring_membership; provides
//     reverse-index caches for per-atom / per-ring lookups)
//
// Does NOT own:
//   - Conformation (the trajectory or single pose) — held separately on
//     the QtLoadResult and accessed via the loader's surface
//
// The library's nmr::Protein owns its conformations; the reader keeps
// them separate so the existing ReaderMainWindow's ownership pattern
// (loader returns protein + conformation, MainWindow holds both) stays
// intact through the rewrite without invasive MainWindow changes.

#pragma once

#include "QtAtom.h"
#include "QtAtomNames.h"
#include "QtBond.h"
#include "QtResidue.h"
#include "QtResidueNames.h"
#include "QtRing.h"
#include "QtRingMembership.h"
#include "QtTopology.h"

#include <QString>
#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::io {
class QtProteinLoader;
}

namespace h5reader::model {

class QtProtein {
public:
    QtProtein() = default;
    ~QtProtein() = default;

    QtProtein(const QtProtein&) = delete;
    QtProtein& operator=(const QtProtein&) = delete;
    QtProtein(QtProtein&&) = delete;
    QtProtein& operator=(QtProtein&&) = delete;

    // ----- Identity -----
    const QString& proteinId() const { return proteinId_; }

    // ----- Atoms (substrate-typed identity) -----
    std::size_t atomCount() const { return atoms_.size(); }
    const QtAtom& atom(std::size_t i) const { return atoms_[i]; }
    const std::vector<QtAtom>& atoms() const { return atoms_; }
    void setAtomPartialCharge(std::size_t i, double charge) {
        atoms_[i].partialCharge = charge;
        atoms_[i].hasPartialCharge = true;
    }

    // ----- Atom names (projection layer; explicit accessor) -----
    const QtAtomNames& atomNames(std::size_t i) const { return atomNames_[i]; }

    // ----- Residues (substrate-typed) -----
    std::size_t residueCount() const { return residues_.size(); }
    const QtResidue& residue(std::size_t i) const { return residues_[i]; }
    const std::vector<QtResidue>& residues() const { return residues_; }

    // ----- Name projections (selectable; labels are projections, NOT
    // identity — the reader is not label-driven). atomNames() above is the
    // verbatim atom projection; residueNames() is the verbatim residue
    // projection. The *Label() helpers pick a convention (and, for
    // residues, a source: Verbatim from the sidecar vs Derived from the
    // typed AminoAcid). -----
    const QtResidueNames& residueNames(std::size_t i) const { return residueNames_[i]; }
    QString atomLabel(std::size_t i, NamingConvention conv) const;
    QString residueLabel(std::size_t i, NamingConvention conv, NamingSource src) const;

    // ----- Topology (bonds + rings + ring_membership) -----
    const QtTopology& topology() const { return *topology_; }
    QtTopology& mutableTopology() { return *topology_; }

    // ----- Bonds (delegated through topology) -----
    std::size_t bondCount() const { return topology_ ? topology_->bondCount() : 0; }
    const QtBond& bond(std::size_t i) const { return topology_->bondAt(i); }

    // ----- Rings (delegated through topology) -----
    std::size_t ringCount() const { return topology_ ? topology_->ringCount() : 0; }
    const QtRing& ring(std::size_t i) const { return topology_->ringAt(i); }

    // ----- Ring membership (delegated through topology) -----
    std::size_t ringMembershipCount() const { return topology_ ? topology_->ringMembershipCount() : 0; }

    // ----- Counts by typed enum (cached at load) -----
    int ringCountByType(RingTypeIndex t) const;

private:
    friend class ::h5reader::io::QtProteinLoader;

    QString proteinId_;
    std::vector<QtAtom> atoms_;
    std::vector<QtAtomNames> atomNames_;
    std::vector<QtResidueNames> residueNames_;  // verbatim projection, parallel to residues_
    std::vector<QtResidue> residues_;
    std::unique_ptr<QtTopology> topology_;
};

}  // namespace h5reader::model
