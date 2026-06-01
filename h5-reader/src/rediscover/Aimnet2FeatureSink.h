// Aimnet2FeatureSink -- per-atom AIMNet2 feature carrier.
//
// Emits one row per backbone (atom, frame), with target identity in CSV and
// row-aligned NPY sidecars for the wide payloads. The charge-response-gradient
// is deliberately labelled as CRG, not polarizability: it is the AIMNet2
// d(sum(q^2))/dr scalar/vector feature.

#pragma once

#include "../model/Types.h"

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

using model::Mat3;
using model::SphericalTensor;
using model::Vec3;

struct Aimnet2FeatureRow {
    int32_t atom_index = -1;
    int32_t residue_index = -1;
    int32_t residue_number = 0;
    int amino_acid = -1;
    int element = -1;
    QString atom_name;
    int frame_variant = 0;
    bool frame_valid = false;
    int32_t frame_anchor_atom_index = -1;
    Vec3 frame_z = Vec3::Zero();
    Vec3 frame_x = Vec3::Zero();
    Vec3 frame_y = Vec3::Zero();
    int32_t h5_row = -1;
    int32_t original_index = -1;
    double time_ps = 0.0;

    bool dft_present = false;
    Mat3 dft_total_raw = Mat3::Zero();
    SphericalTensor dft_total_decomp;
    std::array<double, 5> dft_target_local_T2 = {};
    bool dft_local_frame_valid = false;

    bool aimnet2_charge_present = false;
    double aimnet2_charge = 0.0;

    bool crg_present = false;
    double crg_scalar = 0.0;
    Vec3 crg_vector_lab = Vec3::Zero();
    Vec3 crg_vector_local = Vec3::Zero();

    bool embedding_present = false;
    const float* embedding = nullptr;
    std::size_t embedding_dims = 0;
};

struct Aimnet2FeatureStats {
    std::size_t rows = 0;
    std::size_t dft_present = 0;
    std::size_t frame_valid = 0;
    std::size_t charge_present = 0;
    std::size_t crg_present = 0;
    std::size_t embedding_present = 0;
    std::size_t embedding_dims = 256;
};

bool FiniteAimnetVec3(const Vec3& v);

class Aimnet2FeatureSink {
public:
    Aimnet2FeatureSink(const QString& outDir, const QString& caseName,
                       std::size_t embeddingDims = 256);
    ~Aimnet2FeatureSink();

    bool Ok() const { return ok_; }
    void Write(const Aimnet2FeatureRow& row);
    bool Commit();

    std::size_t rowsWritten() const { return rows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }
    const QString& aggregatedPath() const { return aggregatedPath_; }
    std::size_t embeddingDims() const { return embeddingDims_; }

private:
    QString caseName_;
    QString aggregatedPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> aggregatedFile_;
    std::unique_ptr<QTextStream> aggregatedOut_;
    bool ok_ = false;
    std::size_t rows_ = 0;
    int64_t nextRowId_ = 0;
    std::size_t embeddingDims_ = 256;

    std::vector<double> crgVectorLocal_;     // (rows, 3)
    std::vector<double> crgVectorLab_;       // (rows, 3), audit payload
    std::vector<double> crgScalar_;          // (rows,)
    std::vector<float> embedding_;           // (rows, embeddingDims), float32
    std::vector<double> targetLocalT2_;      // (rows, 5)
    std::vector<double> targetLabT2_;        // (rows, 5)
};

}  // namespace h5reader::rediscover
