#pragma once

#include "../model/Types.h"

#include <QString>
#include <QStringList>

#include <array>
#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

enum class RowColType { String, Int, Double, Bool };
enum class RowNativeAxis { Row, Protein, Frame, Atom, Residue, Ring, Bond };
enum class RowDatasetScope { Both, OneP9J, Static720 };

enum class ResidueClass : int {
    Hydrophobic = 0,
    Polar = 1,
    Aromatic = 2,
    PosCharged = 3,
    NegCharged = 4,
    Glycine = 5,
    Proline = 6,
    Cys = 7,
    Other = 8,
};

struct RowColumnSpec {
    QString name;
    RowColType type = RowColType::Double;
    QString unit;
    QString irrep;
    RowNativeAxis nativeAxis = RowNativeAxis::Row;
    bool timeVarying = false;
    RowDatasetScope scope = RowDatasetScope::Both;
    double nullSentinel = 0.0;
    bool signFlipLegal = false;
    bool isFeature = true;
};

struct RowDesignRow {
    std::vector<QString> values;
    std::array<double, 5> targetT2 = {};
    bool targetT2Present = false;
    std::array<double, 162> ringTensors = {};
    bool ringTensorsPresent = false;
    const float* embedding = nullptr;
    std::size_t embeddingDims = 0;
};

ResidueClass ClassifyResidue(model::AminoAcid aa);
const char* NameForResidueClass(ResidueClass c);

const std::vector<RowColumnSpec>& RowDesignSchema();
QStringList RowDesignHeader();
int RowDesignColumnIndex(const QString& name);

}  // namespace h5reader::rediscover
