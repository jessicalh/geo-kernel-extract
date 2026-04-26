#include "Protein.h"
#include "AminoAcidType.h"
#include <algorithm>
#include <map>
#include <set>
#include <cstdio>

namespace nmr {

size_t Protein::AddAtom(std::unique_ptr<Atom> atom) {
    size_t idx = atoms_.size();
    atoms_.push_back(std::move(atom));
    return idx;
}

size_t Protein::AddResidue(Residue residue) {
    size_t idx = residues_.size();
    residues_.push_back(std::move(residue));
    return idx;
}


// ============================================================================
// Conformation factory methods
// ============================================================================

ProteinConformation& Protein::AddConformation(
    std::vector<Vec3> positions,
    std::string description)
{
    auto conf = std::make_unique<ProteinConformation>(
        this, std::move(positions), std::move(description));
    conformations_.push_back(std::move(conf));
    return *conformations_.back();
}

CrystalConformation& Protein::AddCrystalConformation(
    std::vector<Vec3> positions,
    double resolution, double r_factor,
    double temperature, std::string pdb_id)
{
    auto conf = std::make_unique<CrystalConformation>(
        this, std::move(positions), resolution, r_factor,
        temperature, std::move(pdb_id));
    crystal_index_ = conformations_.size();
    conformations_.push_back(std::move(conf));
    return static_cast<CrystalConformation&>(*conformations_.back());
}

PredictionConformation& Protein::AddPrediction(
    std::vector<Vec3> positions,
    std::string method,
    double confidence)
{
    auto conf = std::make_unique<PredictionConformation>(
        this, std::move(positions), std::move(method), confidence);
    prediction_indices_.push_back(conformations_.size());
    conformations_.push_back(std::move(conf));
    return static_cast<PredictionConformation&>(*conformations_.back());
}

MDFrameConformation& Protein::AddMDFrame(
    std::vector<Vec3> positions,
    int walker, double time_ps, double weight,
    double rmsd_nm, double rg_nm)
{
    auto conf = std::make_unique<MDFrameConformation>(
        this, std::move(positions), walker, time_ps, weight, rmsd_nm, rg_nm);
    md_frame_indices_.push_back(conformations_.size());
    conformations_.push_back(std::move(conf));
    return static_cast<MDFrameConformation&>(*conformations_.back());
}

PredictionConformation& Protein::PredictionAt(size_t i) {
    return static_cast<PredictionConformation&>(*conformations_[prediction_indices_[i]]);
}

const PredictionConformation& Protein::PredictionAt(size_t i) const {
    return static_cast<const PredictionConformation&>(*conformations_[prediction_indices_[i]]);
}

MDFrameConformation& Protein::MDFrameAt(size_t i) {
    return static_cast<MDFrameConformation&>(*conformations_[md_frame_indices_[i]]);
}

const MDFrameConformation& Protein::MDFrameAt(size_t i) const {
    return static_cast<const MDFrameConformation&>(*conformations_[md_frame_indices_[i]]);
}

DerivedConformation& Protein::AddDerived(
    const ProteinConformation& /*parent*/,
    std::string description,
    std::vector<Vec3> positions)
{
    // parent is accepted for future use (tracking derivation lineage)
    auto conf = std::make_unique<DerivedConformation>(
        this, std::move(positions), std::move(description));
    conformations_.push_back(std::move(conf));
    return static_cast<DerivedConformation&>(*conformations_.back());
}

ProteinConformation& Protein::Conformation() {
    if (conformations_.empty()) {
        fprintf(stderr, "FATAL: Protein::Conformation() -- no conformations.\n");
        std::abort();
    }
    return *conformations_[0];
}

const ProteinConformation& Protein::Conformation() const {
    if (conformations_.empty()) {
        fprintf(stderr, "FATAL: Protein::Conformation() -- no conformations.\n");
        std::abort();
    }
    return *conformations_[0];
}

CrystalConformation& Protein::CrystalConf() {
    if (crystal_index_ == SIZE_MAX) {
        fprintf(stderr, "FATAL: Protein::CrystalConf() -- no crystal conformation.\n");
        std::abort();
    }
    return static_cast<CrystalConformation&>(*conformations_[crystal_index_]);
}

const CrystalConformation& Protein::CrystalConf() const {
    if (crystal_index_ == SIZE_MAX) {
        fprintf(stderr, "FATAL: Protein::CrystalConf() -- no crystal conformation.\n");
        std::abort();
    }
    return static_cast<const CrystalConformation&>(*conformations_[crystal_index_]);
}


// ============================================================================
// FinalizeConstruction: every loader must call this after adding all atoms
// and residues. Ensures backbone indices, bonds, and rings are all detected.
// ============================================================================

void Protein::FinalizeConstruction(const std::vector<Vec3>& positions,
                                    double bond_tolerance) {
    // Layer 2: symbolic topology (no geometry needed)
    CacheResidueBackboneIndices();
    DetectAromaticRings();

    // Layer 3: geometric topology (the geometry→topology boundary)
    topology_ = CovalentTopology::Resolve(atoms_, rings_, residues_,
                                           positions, bond_tolerance);

    // Copy connectivity back to Atom objects for convenient calculator access.
    // Calculators read atom.bond_indices and atom.parent_atom_index directly.
    // CovalentTopology is the authority; these are copies for access convenience.
    for (size_t i = 0; i < atoms_.size(); ++i) {
        atoms_[i]->bond_indices = topology_->BondIndicesFor(i);
        atoms_[i]->parent_atom_index = topology_->HydrogenParentOf(i);
    }

    // Layer 2.5: IUPAC symbolic topology (Markley 1998).
    // Stamps typed per-atom identity from the AminoAcidType table.
    // Additive to the existing layers; calculators ignoring it see no change.
    StampAtomTopology();
    StampRingMembership();
    StampResidueContext();
    StampPartnerAtomIndex();
    StampChiPositions();
}


// ============================================================================
// IUPAC symbolic topology layer (Markley 1998, J Biomol NMR 12:1-23).
//
// Source PDF: references/markley-1998-iupac-nmr-nomenclature-recommendations.pdf
//
// PDB LOADING BOUNDARY (continued): the PDB atom-name string is used here for
// one final lookup against the AminoAcidType template, after which all
// downstream identity work uses the typed AtomTopology payload.
// ============================================================================

void Protein::StampAtomTopology() {
    std::map<std::string, int> unstamped_by_residue;
    int stamped_count = 0;
    int unstamped_count = 0;

    for (const auto& res : residues_) {
        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (size_t ai : res.atom_indices) {
            Atom& atom = *atoms_[ai];
            const AminoAcidAtom* templ = aatype.FindAtom(atom.iupac_name.c_str());

            // PDB/BMRB use "H" for the backbone amide H; IUPAC recommends "HN".
            // The table uses "H"; normalise the PDB-side name on lookup.
            if (!templ && atom.iupac_name == "HN") {
                templ = aatype.FindAtom("H");
            }

            if (templ && templ->topology.stamped) {
                atom.topology = templ->topology;
                ++stamped_count;
            } else {
                unstamped_by_residue[aatype.three_letter_code]++;
                ++unstamped_count;
            }
        }
    }

    if (unstamped_count > 0) {
        std::string types_str;
        for (const auto& kv : unstamped_by_residue) {
            if (!types_str.empty()) types_str += ", ";
            types_str += kv.first + "(" + std::to_string(kv.second) + ")";
        }
        fprintf(stderr,
                "StampAtomTopology: %d atoms stamped, %d unstamped "
                "[residue counts: %s] — table population in progress.\n",
                stamped_count, unstamped_count, types_str.c_str());
    }
}


void Protein::StampRingMembership() {
    for (auto& atom : atoms_) atom->ring_indices.clear();

    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        const Ring& ring = *rings_[ri];
        for (size_t ai : ring.atom_indices) {
            atoms_[ai]->ring_indices.push_back(ri);
        }
    }
}


void Protein::StampResidueContext() {
    const size_t n = residues_.size();
    for (size_t ri = 0; ri < n; ++ri) {
        Residue& res = residues_[ri];

        // Previous in same chain
        if (ri > 0 && residues_[ri - 1].chain_id == res.chain_id) {
            res.prev_residue_type = residues_[ri - 1].type;
            res.is_n_terminal = false;
        } else {
            res.prev_residue_type = AminoAcid::Unknown;
            res.is_n_terminal = true;
        }

        // Next in same chain
        if (ri + 1 < n && residues_[ri + 1].chain_id == res.chain_id) {
            res.next_residue_type = residues_[ri + 1].type;
            res.is_c_terminal = false;
        } else {
            res.next_residue_type = AminoAcid::Unknown;
            res.is_c_terminal = true;
        }
    }
}


void Protein::StampPartnerAtomIndex() {
    for (const auto& res : residues_) {
        for (size_t ai : res.atom_indices) {
            Atom& atom_a = *atoms_[ai];
            atom_a.partner_atom_index = SIZE_MAX;

            const AtomTopology& ta = atom_a.topology;
            if (ta.diastereotopic_index == DiastereotopicIndex::None) continue;

            DiastereotopicIndex partner_idx =
                (ta.diastereotopic_index == DiastereotopicIndex::Two)
                    ? DiastereotopicIndex::Three
                    : DiastereotopicIndex::Two;

            for (size_t bi : res.atom_indices) {
                if (bi == ai) continue;
                const AtomTopology& tb = atoms_[bi]->topology;
                if (tb.locant == ta.locant
                    && tb.branch_index == ta.branch_index
                    && tb.diastereotopic_index == partner_idx) {
                    atom_a.partner_atom_index = bi;
                    break;
                }
            }
        }
    }
}


void Protein::StampChiPositions() {
    for (const auto& res : residues_) {
        const AminoAcidType& aatype = res.AminoAcidInfo();
        const size_t n_chi = aatype.chi_angles.size();
        for (size_t chi_idx = 0; chi_idx < n_chi && chi_idx < 4; ++chi_idx) {
            const ChiAngleDef& chi = aatype.chi_angles[chi_idx];
            for (int pos = 0; pos < 4; ++pos) {
                const char* atom_name = chi.atoms[pos];
                for (size_t ai : res.atom_indices) {
                    if (atoms_[ai]->iupac_name == atom_name) {
                        atoms_[ai]->topology.chi_position[chi_idx] =
                            static_cast<int8_t>(pos);
                        break;
                    }
                }
            }
        }
    }
}


// ============================================================================
// DetectAromaticRings -- from residue types and atom presence
//
// PDB LOADING BOUNDARY: string comparisons are used here to match PDB atom
// names against ring definitions from the AminoAcidType table. This is the
// translation from PDB naming to typed ring objects. After this function,
// rings are typed objects (PheBenzeneRing, HisImidazoleRing, etc.) with
// atom indices -- no further string work.
//
// HIS tautomer detection: checks for "HD1" and "HE2" hydrogen atoms to
// determine protonation state (HID/HIE/HIP). This is at the loading
// boundary where the PDB hydrogen names are the authoritative source.
// ============================================================================

void Protein::DetectAromaticRings() {
    rings_.clear();

    // For TRP: track the 5-ring and 6-ring indices to set fused partners
    struct TrpRingIndices {
        size_t benzene_idx = SIZE_MAX;
        size_t pyrrole_idx = SIZE_MAX;
        size_t perimeter_idx = SIZE_MAX;
    };
    std::map<size_t, TrpRingIndices> trp_rings;

    for (size_t ri = 0; ri < residues_.size(); ++ri) {
        const Residue& res = residues_[ri];
        const AminoAcidType& aatype = res.AminoAcidInfo();

        if (!aatype.is_aromatic) continue;

        // Build a name->atom_index map for this residue
        std::map<IupacAtomName, size_t> name_to_idx;
        for (size_t ai : res.atom_indices) {
            name_to_idx[atoms_[ai]->iupac_name] = ai;
        }

        for (const auto& ring_def : aatype.rings) {
            // Check all ring atoms are present
            std::vector<size_t> atom_indices;
            bool all_present = true;
            for (const char* aname : ring_def.atom_names) {
                auto it = name_to_idx.find(IupacAtomName(aname));
                if (it == name_to_idx.end()) {
                    all_present = false;
                    break;
                }
                atom_indices.push_back(it->second);
            }
            if (!all_present) continue;

            // Determine ring type -- for HIS, use protonation_variant_index
            // if ProtonationDetectionResult has run (preferred), otherwise
            // fall back to string check at the PDB loading boundary.
            RingTypeIndex effective_type = ring_def.type_index;
            if (res.type == AminoAcid::HIS) {
                if (res.protonation_variant_index >= 0) {
                    // Typed path: read from protonation_variant_index
                    // (set by ProtonationDetectionResult)
                    // 0 = HID, 1 = HIE, 2 = HIP
                    switch (res.protonation_variant_index) {
                        case 0: effective_type = RingTypeIndex::HidImidazole; break;
                        case 1: effective_type = RingTypeIndex::HieImidazole; break;
                        case 2: effective_type = RingTypeIndex::HisImidazole; break;
                        default: break;  // keep HisImidazole (ambiguous)
                    }
                } else {
                    // PDB LOADING BOUNDARY fallback: protonation not yet
                    // detected (e.g., crystal structure with no H). Check
                    // hydrogen atom names to determine tautomer.
                    bool has_HD1 = name_to_idx.find(IupacAtomName("HD1")) != name_to_idx.end();
                    bool has_HE2 = name_to_idx.find(IupacAtomName("HE2")) != name_to_idx.end();
                    if (has_HD1 && has_HE2) {
                        effective_type = RingTypeIndex::HisImidazole;
                    } else if (has_HD1) {
                        effective_type = RingTypeIndex::HidImidazole;
                    } else if (has_HE2) {
                        effective_type = RingTypeIndex::HieImidazole;
                    }
                    // else: no H on either N -- keep HisImidazole (ambiguous)
                }
            }

            auto ring = CreateRing(effective_type);
            ring->atom_indices = std::move(atom_indices);
            ring->parent_residue_index = ri;
            ring->parent_residue_number = res.sequence_number;

            size_t ring_idx = rings_.size();
            rings_.push_back(std::move(ring));

            // Track TRP rings for fused partner assignment
            if (res.type == AminoAcid::TRP) {
                auto& trp = trp_rings[ri];
                if (effective_type == RingTypeIndex::TrpBenzene)
                    trp.benzene_idx = ring_idx;
                else if (effective_type == RingTypeIndex::TrpPyrrole)
                    trp.pyrrole_idx = ring_idx;
                else if (effective_type == RingTypeIndex::TrpPerimeter)
                    trp.perimeter_idx = ring_idx;
            }
        }
    }

    // Set fused partner indices for TRP rings
    for (auto& kv : trp_rings) {
        auto& trp = kv.second;
        if (trp.benzene_idx != SIZE_MAX && trp.pyrrole_idx != SIZE_MAX) {
            rings_[trp.benzene_idx]->fused_partner_index = trp.pyrrole_idx;
            rings_[trp.pyrrole_idx]->fused_partner_index = trp.benzene_idx;
        }
    }
}


// ============================================================================
// CacheResidueBackboneIndices
//
// PDB LOADING BOUNDARY: string comparisons are used here to match PDB atom
// names ("N", "CA", "C", "O", "H", "HA", "CB") against each residue's atoms.
// This populates typed backbone index fields (res.N, res.CA, etc.) which are
// atom INDICES, not strings. After this function, backbone identity is
// determined entirely by these typed indices -- no further string work.
// ============================================================================

void Protein::CacheResidueBackboneIndices() {
    for (auto& res : residues_) {
        for (size_t ai : res.atom_indices) {
            const IupacAtomName& name = atoms_[ai]->iupac_name;
            if      (name == "N")   res.N  = ai;
            else if (name == "CA")  res.CA = ai;
            else if (name == "C")   res.C  = ai;
            else if (name == "O")   res.O  = ai;
            else if (name == "H" || name == "HN")  res.H  = ai;
            else if (name == "HA" || name == "HA2") res.HA = ai;
            else if (name == "CB")  res.CB = ai;
        }

        // PDB LOADING BOUNDARY (continued): match chi angle atom names
        // against the residue's atoms. After this, chi angles are typed
        // index tuples — no further string work.
        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (int ci = 0; ci < aatype.chi_angle_count && ci < 4; ++ci) {
            const ChiAngleDef& def = aatype.chi_angles[ci];
            for (int j = 0; j < 4; ++j) {
                for (size_t ai : res.atom_indices) {
                    if (atoms_[ai]->iupac_name == def.atoms[j]) {
                        res.chi[ci].a[j] = ai;
                        break;
                    }
                }
            }
        }
    }
}


}  // namespace nmr
