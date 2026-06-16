// AnalysisAtom -- long-lived per-atom first-pass collector.
//
// The object owns its accumulated state and only borrows the immutable
// rediscover Body. This pass deliberately exposes diagnostics only: enough to
// prove every atom folded every DFT-bearing frame.

#pragma once

#include "AnalysisBody.h"
#include "PerAtomSubstrate.h"

#include <QString>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace h5reader::rediscover {

struct AnalysisAtomDiagnostics {
    std::size_t atom_count = 0;
    std::size_t dft_rows = 0;
    std::size_t frame_events = 0;
    std::size_t relationship_folds = 0;
    std::size_t self_property_folds = 0;
    std::size_t relationship_organs = 0;
    std::size_t self_property_organs = 0;
    std::size_t dft_present = 0;
    std::size_t sigma_folds = 0;
    std::size_t dihedral_folds = 0;
    std::size_t mopac_scalar_folds = 0;
    std::size_t efg_folds = 0;
    double sample_value = 0.0;
    QString sample_label;
    bool has_sample = false;
};

struct AnalysisAtomTruth {
    std::size_t frame_events = 0;
    std::size_t relationship_folds = 0;
    std::size_t self_property_folds = 0;
    std::size_t relationship_organs = 0;
    std::size_t self_property_organs = 0;
    std::size_t dft_present = 0;
    std::size_t sigma_folds = 0;
    std::size_t dihedral_folds = 0;
    std::size_t mopac_scalar_folds = 0;
    std::size_t efg_folds = 0;
    double sample_value = 0.0;
    QString sample_label;
    bool has_sample = false;
};

class AnalysisAtom {
public:
    AnalysisAtom(const Body& body, std::size_t atom_index);
    virtual ~AnalysisAtom() = default;

    AnalysisAtom(const AnalysisAtom&) = delete;
    AnalysisAtom& operator=(const AnalysisAtom&) = delete;

    std::size_t atomIndex() const { return atom_index_; }

    virtual void observeFrame(std::size_t h5_row) = 0;
    virtual AnalysisAtomTruth diagnostics() const = 0;

protected:
    const Body& body_;
    const std::size_t atom_index_;
};

class NeighbourhoodAccumulatingAnalysisAtom final : public AnalysisAtom {
public:
    NeighbourhoodAccumulatingAnalysisAtom(const Body& body,
                                          std::size_t atom_index,
                                          PerAtomSubstrateConfig config);

    void observeFrame(std::size_t h5_row) override;
    AnalysisAtomTruth diagnostics() const override;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

AnalysisAtomDiagnostics RunAnalysisAtomFirstPass(const Body& body,
                                                 const PerAtomSubstrateConfig& config);

}  // namespace h5reader::rediscover
