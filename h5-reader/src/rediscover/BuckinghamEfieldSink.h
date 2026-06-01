// BuckinghamEfieldSink — per_atom_feature carrier for APBS E-field -> DFT T0.
//
// One row per backbone DFT-present (atom, frame). The scalar fit reads only
// the producer-emitted E_proj, |E|, and sigma_iso columns. The row also carries
// the full DFT total tensor plus T1/T2 sidecars for audit completeness; T1 is
// emitted but convention-unverified and must not be fitted.

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

struct BuckinghamEfieldRow {
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
    bool apbs_efield_present = false;
    Vec3 efield_local = Vec3::Zero();
    double e_proj = 0.0;  // local z component
    double e_mag = 0.0;
    QString efield_units;

    Mat3 dft_total_raw = Mat3::Zero();
    SphericalTensor dft_total_decomp;
};

struct BuckinghamEfieldStats {
    std::size_t rows = 0;
    std::size_t dft_present = 0;
    std::size_t frame_valid = 0;
    std::size_t apbs_efield_present = 0;
    std::size_t finite_efield = 0;
    double min_efield_magnitude = 0.0;
    double max_efield_magnitude = 0.0;
};

bool FiniteVec3(const Vec3& v);

class BuckinghamEfieldSink {
public:
    BuckinghamEfieldSink(const QString& outDir, const QString& caseName);
    ~BuckinghamEfieldSink();

    bool Ok() const { return ok_; }
    void Write(const BuckinghamEfieldRow& row);
    bool Commit();

    std::size_t rowsWritten() const { return rows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }
    const QString& aggregatedPath() const { return aggregatedPath_; }

private:
    QString caseName_;
    QString aggregatedPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> aggregatedFile_;
    std::unique_ptr<QTextStream> aggregatedOut_;
    bool ok_ = false;
    std::size_t rows_ = 0;
    int64_t nextRowId_ = 0;

    std::vector<double> featureFieldLocal_;  // (rows, 3), local APBS E field
    std::vector<double> targetT1_;           // (rows, 3), emitted; unverified convention
    std::vector<double> targetT2_;           // (rows, 5), emitted for completeness
};

}  // namespace h5reader::rediscover
