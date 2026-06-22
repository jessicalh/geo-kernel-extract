// AnalysisElement -- stepwise analysis-object model.
//
// One long-lived object exists per atom, ring, and residue.  The object owns its
// retained series and is called once for every trajectory step by the
// analysis-object driver.  The shared Body/model/catalog/indexes are borrowed
// read-only; the only shared runtime helper here is a cadence/context surface
// that answers frame identity and exposes per-step source snapshots.

#pragma once

#include "AnalysisBody.h"
#include "PerAtomSubstrate.h"

#include <QJsonObject>
#include <QString>
#include <QStringList>

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace h5reader::model {
class QtConformationSnapshot;
}

namespace h5reader::rediscover {

struct MopacFrameBond {
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double wiberg_order = 0.0;
};

class AnalysisCadence {
public:
    explicit AnalysisCadence(const FrameMap& frameMap);

    std::size_t stepCount() const { return original_by_step_.size(); }
    std::size_t sigmaStepCount() const { return sigma_rows_.size(); }
    std::size_t originalIndex(std::size_t step) const { return original_by_step_[step]; }
    bool sigmaPresent(std::size_t step) const { return sigma_present_[step]; }
    const std::vector<bool>& sigmaMask() const { return sigma_present_; }
    const std::vector<std::size_t>& sigmaRows() const { return sigma_rows_; }
    std::optional<std::size_t> nominalStride() const { return nominal_stride_; }

private:
    std::vector<std::size_t> original_by_step_;
    std::vector<bool> sigma_present_;
    std::vector<std::size_t> sigma_rows_;
    std::optional<std::size_t> nominal_stride_;
};

class AnalysisObjectContext {
public:
    AnalysisObjectContext(const Body& body, const AnalysisCadence& cadence);

    const Body& body() const { return body_; }
    const AnalysisCadence& cadence() const { return cadence_; }

    std::shared_ptr<const model::QtConformationSnapshot> snapshot(std::size_t step) const;
    const std::vector<MopacFrameBond>& mopacBonds(std::size_t step) const;
    std::vector<MopacFrameBond> mopacBondsForAtom(std::size_t step, std::size_t atom) const;
    std::optional<double> mopacWiberg(std::size_t step, std::size_t atomA, std::size_t atomB) const;
    std::vector<double> mopacWibergSeries(std::size_t atomA, std::size_t atomB) const;

private:
    void rebuildMopacCache(std::size_t step) const;

    const Body& body_;
    const AnalysisCadence& cadence_;
    mutable std::size_t cached_step_ = static_cast<std::size_t>(-1);
    mutable std::shared_ptr<const model::QtConformationSnapshot> cached_snapshot_;
    mutable std::vector<MopacFrameBond> cached_mopac_bonds_;
    mutable std::map<std::uint64_t, std::vector<double>> mopac_wiberg_by_pair_;
};

class AnalysisElement {
public:
    AnalysisElement(const AnalysisObjectContext& context,
                    QString objectType,
                    std::size_t modelIndex);
    virtual ~AnalysisElement() = default;

    AnalysisElement(const AnalysisElement&) = delete;
    AnalysisElement& operator=(const AnalysisElement&) = delete;

    virtual void Calculate(std::size_t step) = 0;
    virtual QJsonObject Truth() const = 0;

    const QString& objectType() const { return object_type_; }
    std::size_t modelIndex() const { return model_index_; }

protected:
    const AnalysisObjectContext& context_;
    const Body& body_;
    const AnalysisCadence& cadence_;
    QString object_type_;
    std::size_t model_index_ = 0;
};

class AnalysisStructure : public AnalysisElement {
public:
    AnalysisStructure(const AnalysisObjectContext& context,
                      QString objectType,
                      std::size_t modelIndex);
};

class AnalysisRing final : public AnalysisStructure {
public:
    AnalysisRing(const AnalysisObjectContext& context,
                 std::size_t ringId,
                 double nearCenterCutoffA);

    void Calculate(std::size_t step) override;
    QJsonObject Truth() const override;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

class AnalysisResidue final : public AnalysisStructure {
public:
    AnalysisResidue(const AnalysisObjectContext& context, std::size_t residueIndex);

    void Calculate(std::size_t step) override;
    QJsonObject Truth() const override;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

struct AnalysisObjectPassConfig {
    PerAtomSubstrateConfig per_atom;
    bool ndjson = false;
    // The small-run emit selector (ACCUMULATOR_RESPEC work-item 9). When
    // non-empty, ONLY the listed atoms emit an accumulator/object JSON; the
    // FULL protein still runs every step as the field/EFG/ring source
    // environment for every emitted atom, and every emitted atom runs the FULL
    // trajectory. Each entry is "residue_number:atom_name" (e.g. "ASP7:CG"),
    // resolved against the loaded protein. Empty => all atoms emit (production:
    // byte-identical path).
    QStringList only_atoms;
};

struct AnalysisObjectPassDiagnostics {
    std::size_t atom_count = 0;
    std::size_t ring_count = 0;
    std::size_t residue_count = 0;
    std::size_t step_count = 0;
    std::size_t sigma_step_count = 0;
    std::size_t calculate_calls = 0;
    std::size_t atom_sigma_folds = 0;
    std::size_t atom_relationships = 0;
    std::size_t ring_near_center_hits = 0;
    std::size_t residue_frame_folds = 0;
    std::size_t mapped_bonds = 0;
    std::size_t mismatch_events = 0;
    std::size_t accumulator_responses = 0;
    std::size_t accumulator_contexts = 0;
    std::size_t bounded_sigma_rows = 0;
    std::size_t bounded_sigma_atoms = 0;
    qint64 bounded_sigma_bytes = 0;
    std::size_t classical_source_rows = 0;
    qint64 classical_source_bytes = 0;
    std::size_t classical_source_leaf_rows = 0;
    qint64 classical_source_leaf_bytes = 0;
    std::size_t source_family_matrix_rows = 0;
    qint64 source_family_matrix_bytes = 0;
    std::size_t subspace_overlap_rows = 0;
    qint64 subspace_overlap_bytes = 0;
    std::size_t eta2_by_well_rows = 0;
    qint64 eta2_by_well_bytes = 0;
    bool sigma_mask_recorded = false;
    bool field_vectors_retained = false;
    bool full_sigma_tensors_retained = false;
    bool oxygen_gate_passed = false;
    QString bounded_sigma_path;
    QString classical_source_path;
    QString classical_source_leaf_path;
    QString source_family_matrix_path;
    QString subspace_overlap_path;
    QString eta2_by_well_path;
    QString manifest_path;
};

bool RunAnalysisObjectPass(const Body& body,
                           const QString& outDir,
                           const AnalysisObjectPassConfig& config,
                           AnalysisObjectPassDiagnostics* diagnostics,
                           QString* errOut);

}  // namespace h5reader::rediscover
