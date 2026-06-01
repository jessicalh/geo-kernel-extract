// EfgFeatureSink — per_atom_feature carrier for the APBS EFG -> DFT T2 cell.
//
// This is intentionally small and direct: one aggregated row per DFT-present
// (atom, frame), plus corrected local-frame T2 sidecars and lab-frame audit
// sidecars. It is the second non-source_sum output shape after broad_backbone;
// #29 should unify these sibling carriers, but this spike keeps the ring/mc
// oracle path untouched.

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
using model::Vec3;

struct EfgFeatureRow {
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
    bool apbs_efg_present = false;
    std::array<double, 5> efg_feature_T2 = {};      // local frame
    std::array<double, 5> dft_target_T2 = {};       // local frame
    std::array<double, 5> efg_feature_lab_T2 = {};  // lab-frame audit payload
    std::array<double, 5> dft_target_lab_T2 = {};   // lab-frame audit payload
    QString efg_units;
};

struct EfgFeatureStats {
    std::size_t rows = 0;
    std::size_t dft_present = 0;
    std::size_t frame_valid = 0;
    std::size_t apbs_efg_present = 0;
    std::size_t finite_efg = 0;
    double min_efg_magnitude = 0.0;
    double max_efg_magnitude = 0.0;
};

std::array<double, 5> DecomposeEfgLibraryT2(const Mat3& efg);
Mat3 ReconstructLibraryT2(const std::array<double, 5>& t2);
double T2Magnitude(const std::array<double, 5>& t2);
bool FiniteT2(const std::array<double, 5>& t2);

class EfgFeatureSink {
public:
    EfgFeatureSink(const QString& outDir, const QString& caseName);
    ~EfgFeatureSink();

    bool Ok() const { return ok_; }
    void Write(const EfgFeatureRow& row);
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

    std::vector<double> featureT2_;     // (rows, 5), APBS EFG local-frame T2
    std::vector<double> targetT2_;      // (rows, 5), DFT target local-frame T2
    std::vector<double> featureLabT2_;  // (rows, 5), APBS EFG lab-frame T2
    std::vector<double> targetLabT2_;   // (rows, 5), DFT target lab-frame T2
};

}  // namespace h5reader::rediscover
