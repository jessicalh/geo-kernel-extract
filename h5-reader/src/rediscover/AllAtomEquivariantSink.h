// AllAtomEquivariantSink -- lab-frame all-atom equivariant geometry carrier.
//
// This is deliberately separate from BroadBackboneSink. BroadBackbone is a
// backbone-local-frame scalar/ridge substrate; this sink writes the corrected
// e3nn-style substrate in the single H5/ORCA-aligned molecular frame. No
// per-atom local frame is emitted or implied.

#pragma once

#include "RediscoverTypes.h"

#include <QSaveFile>
#include <QString>
#include <QStringList>
#include <QTextStream>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

struct AllAtomEquivariantTargetRecord {
    int32_t atom_index = -1;
    int32_t residue_index = -1;
    int32_t residue_number = 0;
    int amino_acid = -1;
    int element = -1;
    QString atom_name;
    int32_t h5_row = -1;
    int32_t original_index = -1;
    double time_ps = 0.0;

    DftTarget target;

    bool apbs_efield_present = false;
    Vec3 apbs_efield_lab = Vec3::Zero();
    bool apbs_efg_present = false;
    std::array<double, 5> apbs_efg_T2 = {};

    bool aimnet2_charge_present = false;
    double aimnet2_charge = 0.0;
    bool aimnet2_crg_present = false;
    double aimnet2_crg_scalar = 0.0;
    Vec3 aimnet2_crg_lab = Vec3::Zero();
    bool aimnet2_embedding_present = false;
    const float* aimnet2_embedding = nullptr;
    std::size_t aimnet2_embedding_dims = 0;
};

struct AllAtomEquivariantSourceRecord {
    int32_t target_atom_index = -1;
    int32_t target_residue_index = -1;
    int32_t h5_row = -1;
    int32_t original_index = -1;
    double time_ps = 0.0;

    QString mechanism;
    QString source_kind;
    QString category;
    int category_ord = -1;
    int32_t source_id = -1;
    int32_t source_atom_index = -1;
    int32_t source_residue_index = -1;
    int32_t source_residue_number = 0;
    int source_element = -1;
    QString source_atom_name;

    Vec3 disp = Vec3::Zero();
    double r = 0.0;
    double inv_r3 = 0.0;
    double cos_theta = 0.0;
    double dipolar = 0.0;
    Vec3 orientation_a = Vec3::Zero();
    Vec3 orientation_b = Vec3::Zero();

    double source_value = 0.0;
    double source_value_2 = 0.0;
    QString source_units;
    bool source_present = true;

    Vec3 value_vec = Vec3::Zero();
    double value_vec_mag = 0.0;
    std::array<double, 5> tensor_T2 = {};
    bool tensor_present = false;

    int32_t ring_index = -1;
    int ring_type_index = -1;
    int ring_aromaticity = -1;
    int ring_size = 0;
    int ring_nitrogen = 0;
    bool ring_fused = false;
    double ring_intensity = 0.0;
    double ring_jb_offset = 0.0;

    int32_t bond_index = -1;
    int bond_category = -1;
    int bond_order = -1;
    int32_t bond_atom_a = -1;
    int32_t bond_atom_b = -1;
    int bond_elem_a = -1;
    int bond_elem_b = -1;
    double bond_length = 0.0;
    bool source_is_self_or_bonded = false;

    QString charge_source;
    double source_q_e = 0.0;
    double q_over_r3 = 0.0;

    bool aimnet2_embedding_present = false;
    std::size_t aimnet2_embedding_dims = 0;
};

class AllAtomEquivariantSink {
public:
    AllAtomEquivariantSink(const QString& outDir, const QString& caseName,
                           std::size_t embeddingDims = 256);
    ~AllAtomEquivariantSink();

    bool Ok() const { return ok_; }
    int64_t WriteTarget(const AllAtomEquivariantTargetRecord& row);
    void WriteSource(int64_t rowId, const AllAtomEquivariantSourceRecord& row);
    bool Commit();

    std::size_t targetRowsWritten() const { return targetRows_; }
    std::size_t sourceRowsWritten() const { return sourceRows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }
    std::size_t embeddingDims() const { return embeddingDims_; }

private:
    QString caseName_;
    QString sourcesPath_;
    QString targetsPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> sourcesFile_;
    std::unique_ptr<QSaveFile> targetsFile_;
    std::unique_ptr<QTextStream> sourcesOut_;
    std::unique_ptr<QTextStream> targetsOut_;
    bool ok_ = false;
    int64_t nextRowId_ = 0;
    std::size_t targetRows_ = 0;
    std::size_t sourceRows_ = 0;
    std::size_t embeddingDims_ = 256;

    std::vector<double> targetT2_;             // (targetRows, 5), lab frame
    std::vector<double> targetSigmaIso_;       // (targetRows,)
    std::vector<double> targetRaw_;            // (targetRows, 9), lab frame
    std::vector<double> apbsEfield_;           // (targetRows, 3), lab frame
    std::vector<double> apbsEfgT2_;            // (targetRows, 5), lab frame
    std::vector<double> aimnet2Charge_;        // (targetRows,)
    std::vector<double> aimnet2Crg_;           // (targetRows, 3), lab frame
    std::vector<double> aimnet2CrgScalar_;     // (targetRows,)
    std::vector<float> aimnet2Embedding_;      // (targetRows, embeddingDims), f32
};

}  // namespace h5reader::rediscover
