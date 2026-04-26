#include "ProtonationDetectionResult.h"
#include "Protein.h"
#include "OperationLog.h"
#include <map>
#include <set>
#include <string>

namespace nmr {

std::unique_ptr<ProtonationDetectionResult> ProtonationDetectionResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("ProtonationDetectionResult::Compute",
        "residues=" + std::to_string(conf.ProteinRef().ResidueCount()));

    auto result = std::make_unique<ProtonationDetectionResult>();
    result->conf_ = &conf;

    // Non-const access to protein for setting protonation_variant_index.
    // This is the ONLY place that writes this field. The Protein is the
    // owner of residues; we access through const_cast at the loading boundary
    // because protonation detection is a one-time topology annotation that
    // logically belongs on the Protein, not the conformation.
    Protein& protein = const_cast<Protein&>(conf.ProteinRef());
    size_t n_residues = protein.ResidueCount();
    result->variant_names_.resize(n_residues);

    // Build disulfide bond set: pairs of SG atom indices that are bonded
    std::set<size_t> disulfide_sg;
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if (bond.category == BondCategory::Disulfide) {
            disulfide_sg.insert(bond.atom_index_a);
            disulfide_sg.insert(bond.atom_index_b);
        }
    }

    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue& res = protein.MutableResidueAt(ri);
        const AminoAcidType& aatype = res.AminoAcidInfo();

        if (!aatype.is_titratable) continue;
        if (aatype.variants.empty()) continue;

        // PDB LOADING BOUNDARY: build name set from atom name strings.
        // This is the ONLY place protonation detection reads PDB atom names.
        std::map<IupacAtomName, size_t> name_to_idx;
        for (size_t ai : res.atom_indices) {
            // PDB LOADING BOUNDARY: string comparison for hydrogen detection
            name_to_idx[protein.AtomAt(ai).iupac_name] = ai;
        }

        // Check whether any hydrogen atoms are present in this residue
        bool has_any_H = false;
        for (size_t ai : res.atom_indices) {
            if (protein.AtomAt(ai).element == Element::H) {
                has_any_H = true;
                break;
            }
        }

        std::string variant_name;
        int variant_idx = -1;

        if (res.type == AminoAcid::HIS) {
            // PDB LOADING BOUNDARY: check "HD1" and "HE2" hydrogen names
            bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();
            bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();

            if (has_HD1 && has_HE2) {
                // Doubly protonated: HIP (variant index 2 in HIS variants)
                variant_name = "HIP";
                variant_idx = 2;
            } else if (has_HD1) {
                // Delta-protonated: HID (variant index 0)
                variant_name = "HID";
                variant_idx = 0;
            } else if (has_HE2) {
                // Epsilon-protonated: HIE (variant index 1)
                variant_name = "HIE";
                variant_idx = 1;
            } else if (!has_any_H) {
                // No hydrogens at all: unresolved (crystal structure without H)
                result->unresolved_count_++;
                continue;
            }
            // else: has other H but not HD1/HE2 — ambiguous, leave unresolved
        }
        else if (res.type == AminoAcid::ASP) {
            // PDB LOADING BOUNDARY: check "HD2" hydrogen name
            bool has_HD2 = name_to_idx.find("HD2") != name_to_idx.end();
            if (has_HD2) {
                variant_name = "ASH";
                variant_idx = 0;  // ASH is variants[0] for ASP
            } else if (!has_any_H) {
                result->unresolved_count_++;
                continue;
            }
            // else: has H but not HD2 — charged ASP (default, no variant)
        }
        else if (res.type == AminoAcid::GLU) {
            // PDB LOADING BOUNDARY: check "HE2" hydrogen name
            bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();
            if (has_HE2) {
                variant_name = "GLH";
                variant_idx = 0;  // GLH is variants[0] for GLU
            } else if (!has_any_H) {
                result->unresolved_count_++;
                continue;
            }
            // else: has H but not HE2 — charged GLU (default, no variant)
        }
        else if (res.type == AminoAcid::CYS) {
            // Check for disulfide: SG bonded to another SG
            // PDB LOADING BOUNDARY: look up "SG" by name
            auto sg_it = name_to_idx.find("SG");
            if (sg_it != name_to_idx.end() &&
                disulfide_sg.count(sg_it->second) > 0) {
                variant_name = "CYX";
                variant_idx = 0;  // CYX is variants[0] for CYS
            } else {
                // PDB LOADING BOUNDARY: check "HG" hydrogen name
                bool has_HG = name_to_idx.find("HG") != name_to_idx.end();
                if (!has_HG && !has_any_H) {
                    result->unresolved_count_++;
                    continue;
                }
                // has_HG means free CYS (default), no variant needed
            }
        }
        else if (res.type == AminoAcid::LYS) {
            // PDB LOADING BOUNDARY: check "HZ1", "HZ2", "HZ3" hydrogen names
            bool has_HZ1 = name_to_idx.find("HZ1") != name_to_idx.end();
            bool has_HZ2 = name_to_idx.find("HZ2") != name_to_idx.end();
            bool has_HZ3 = name_to_idx.find("HZ3") != name_to_idx.end();

            if (has_HZ1 && has_HZ2 && has_HZ3) {
                // Fully protonated LYS (charged) — default, no variant
            } else if (has_HZ1 || has_HZ2) {
                // Missing HZ3 but has some: LYN (deprotonated)
                variant_name = "LYN";
                variant_idx = 0;  // LYN is variants[0] for LYS
            } else if (!has_any_H) {
                result->unresolved_count_++;
                continue;
            }
        }
        else if (res.type == AminoAcid::TYR) {
            // PDB LOADING BOUNDARY: check "HH" hydrogen name
            bool has_HH = name_to_idx.find("HH") != name_to_idx.end();
            if (!has_HH && has_any_H) {
                // Has H but not HH: deprotonated tyrosinate TYM
                variant_name = "TYM";
                variant_idx = 0;  // TYM is variants[0] for TYR
            } else if (!has_any_H) {
                result->unresolved_count_++;
                continue;
            }
        }

        if (variant_idx >= 0) {
            res.protonation_variant_index = variant_idx;
            result->variant_names_[ri] = variant_name;
            result->assigned_count_++;

            OperationLog::Log(OperationLog::Level::Info, LogProtonation,
                "ProtonationDetectionResult",
                "residue " + std::to_string(res.sequence_number) +
                " " + ThreeLetterCodeForAminoAcid(res.type) +
                " -> " + variant_name);
        }
    }

    OperationLog::Log(OperationLog::Level::Info, LogProtonation,
        "ProtonationDetectionResult::Compute",
        "assigned=" + std::to_string(result->assigned_count_) +
        " unresolved=" + std::to_string(result->unresolved_count_));

    return result;
}


std::string ProtonationDetectionResult::VariantNameAt(size_t residue_index) const {
    if (residue_index >= variant_names_.size()) return "";
    return variant_names_[residue_index];
}

}  // namespace nmr
