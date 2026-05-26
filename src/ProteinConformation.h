#pragma once
//
// ProteinConformation: one geometric instance of a protein — its
// positions and everything computed from them. Positions are const
// after construction; new geometry means a new conformation.
//
// Computed results attach as typed ConformationResult singletons,
// reached by conf.Result<T>() / conf.HasResult<T>().
//

#include "Types.h"
#include "ConformationAtom.h"
#include "ConformationResult.h"
#include "GeometryChoice.h"
#include "Ring.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <typeindex>
#include <memory>
#include <string>
#include <sstream>

namespace nmr {

class Protein;

class ProteinConformation {
public:
    ProteinConformation(const Protein* protein,
                        std::vector<Vec3> positions,
                        std::string description = "");
    virtual ~ProteinConformation() = default;

    // The owning Protein (borrowed — the caller must ensure it outlives
    // this conformation).
    const Protein& ProteinRef() const { return *protein_; }

    const std::string& Description() const { return description_; }

    size_t AtomCount() const { return atoms_.size(); }

    Vec3 PositionAt(size_t i) const { return atoms_[i].Position(); }
    const std::vector<Vec3>& Positions() const { return positions_; }

    const ConformationAtom& AtomAt(size_t i) const { return atoms_[i]; }
    ConformationAtom& MutableAtomAt(size_t i) { return atoms_[i]; }

    // Attach a result. Checks dependencies. Returns true on success.
    bool AttachResult(std::unique_ptr<ConformationResult> result);

    // Test-only: attach or replace a result with no dependency or
    // duplicate checks; production code must use AttachResult. For
    // accumulation-path tests of TrajectoryResults whose source result
    // has a deep dependency tree (e.g. the AIMNet2 family).
    void ForceAttachResultForTesting(std::unique_ptr<ConformationResult> result);

    // Type-safe result access. Result<T>() aborts if T is not attached;
    // HasResult<T>() reports presence.
    template<typename T>
    T& Result() {
        auto it = results_.find(std::type_index(typeid(T)));
        if (it == results_.end()) {
            std::ostringstream msg;
            msg << "Result<" << typeid(T).name() << "> not attached. Attached: [";
            bool first = true;
            for (const auto& kv : results_) {
                if (!first) msg << ", ";
                msg << kv.second->Name();
                first = false;
            }
            msg << "]";
            // Calling Result<T> for an unattached T is a programming
            // error: abort with a diagnostic listing what is attached.
            fprintf(stderr, "FATAL: %s\n", msg.str().c_str());
            std::abort();
        }
        return static_cast<T&>(*it->second);
    }

    template<typename T>
    const T& Result() const {
        auto it = results_.find(std::type_index(typeid(T)));
        if (it == results_.end()) {
            fprintf(stderr, "FATAL: Result<%s> not attached.\n", typeid(T).name());
            std::abort();
        }
        return static_cast<const T&>(*it->second);
    }

    template<typename T>
    bool HasResult() const {
        return results_.find(std::type_index(typeid(T))) != results_.end();
    }

    const std::unordered_map<std::type_index, std::unique_ptr<ConformationResult>>&
    AllResults() const { return results_; }

    // Ring geometry storage (set by GeometryResult)
    std::vector<RingGeometry> ring_geometries;

    // Ring pair properties (set by GeometryResult)
    struct RingPair {
        size_t ring_a = 0;
        size_t ring_b = 0;
        double center_distance = 0.0;
        double normal_dot = 0.0;
        double normal_cross_mag = 0.0;
        bool is_fused = false;
    };
    std::vector<RingPair> ring_pairs;

    // Bond geometry (set by GeometryResult)
    std::vector<double> bond_lengths;
    std::vector<Vec3> bond_directions;
    std::vector<Vec3> bond_midpoints;

    // Global geometry (set by GeometryResult)
    Vec3 bounding_min = Vec3::Zero();
    Vec3 bounding_max = Vec3::Zero();
    Vec3 center_of_geometry = Vec3::Zero();
    double radius_of_gyration = 0.0;

    // Pre-built collections (set by GeometryResult)
    std::map<RingTypeIndex, std::vector<size_t>> rings_by_type;
    std::map<BondCategory, std::vector<size_t>> bonds_by_category;
    std::map<AminoAcid, std::vector<size_t>> residues_by_type;

    // Geometric decisions recorded by calculators during Compute(),
    // via the GeometryChoiceBuilder factory.
    std::vector<GeometryChoice> geometry_choices;

protected:
    const Protein* protein_;
    std::string description_;
    std::vector<Vec3> positions_;
    std::vector<ConformationAtom> atoms_;
    std::unordered_map<std::type_index, std::unique_ptr<ConformationResult>> results_;
};


class CrystalConformation : public ProteinConformation {
public:
    CrystalConformation(const Protein* protein,
                        std::vector<Vec3> positions,
                        double resolution,
                        double r_factor,
                        double temperature,
                        std::string pdb_id);

    double ResolutionAngstroms() const { return resolution_angstroms_; }
    double RFactor() const { return r_factor_; }
    double TemperatureKelvin() const { return temperature_kelvin_; }
    const std::string& PdbId() const { return pdb_id_; }

private:
    double resolution_angstroms_;
    double r_factor_;
    double temperature_kelvin_;
    std::string pdb_id_;
};


class PredictionConformation : public ProteinConformation {
public:
    PredictionConformation(const Protein* protein,
                           std::vector<Vec3> positions,
                           std::string method,
                           double confidence = std::nan(""));

    const std::string& PredictionMethod() const { return method_; }
    double Confidence() const { return confidence_; }

private:
    std::string method_;
    double confidence_;
};


class MDFrameConformation : public ProteinConformation {
public:
    MDFrameConformation(const Protein* protein,
                        std::vector<Vec3> positions,
                        int walker,
                        double time_ps,
                        double weight,
                        double rmsd_nm,
                        double rg_nm);

    int Walker() const { return walker_; }
    double TimePicoseconds() const { return time_ps_; }
    double BoltzmannWeight() const { return weight_; }
    double RmsdNanometres() const { return rmsd_nm_; }
    double RadiusOfGyrationNm() const { return rg_nm_; }

private:
    int walker_;
    double time_ps_;
    double weight_;
    double rmsd_nm_;
    double rg_nm_;
};


class DerivedConformation : public ProteinConformation {
public:
    DerivedConformation(const Protein* protein,
                        std::vector<Vec3> positions,
                        std::string derivation_description);

    const std::string& DerivationDescription() const { return derivation_description_; }

private:
    std::string derivation_description_;
};

}  // namespace nmr
