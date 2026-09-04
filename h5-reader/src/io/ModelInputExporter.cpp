#include "ModelInputExporter.h"

#include "QtFieldCatalog.gen.h"
#include "QtNpyWriter.h"

#include "../model/Conformation.h"
#include "../model/QtConformationSnapshot.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLoggingCategory>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <initializer_list>
#include <limits>
#include <utility>

namespace h5reader::io {

namespace {

Q_LOGGING_CATEGORY(cModelInput, "h5reader.modelinput")

constexpr std::size_t kTorsionComponentCount = 6;
constexpr std::size_t kDsspStateCount = 8;
constexpr std::size_t kDsspHbondCount = 4;

enum class ScalarChannel : std::size_t {
    BackboneAngleFirst = 0,
    TorsionCosFirst = 5,
    TorsionSinFirst = 11,
    OmegaCos = 17,
    OmegaSin = 18,
    OmegaIsXpro = 19,
    Pyramidalization = 20,
    DsspStateFirst = 21,
    DsspHbondFirst = 29,
    DsspPpii = 33,
    TotalFormalCharge = 34,
    Count = 35,
};

enum class ApplicabilityChannel : std::size_t {
    CbetaResidual = 0,
    TorsionFirst = 1,
    Pyramidalization = 7,
    Count = 8,
};

constexpr std::size_t channelIndex(ScalarChannel channel) {
    return static_cast<std::size_t>(channel);
}

constexpr std::size_t channelIndex(ApplicabilityChannel channel) {
    return static_cast<std::size_t>(channel);
}

constexpr std::size_t kScalarWidth = channelIndex(ScalarChannel::Count);
constexpr std::size_t kApplicabilityWidth = channelIndex(ApplicabilityChannel::Count);
constexpr std::size_t kL1Width = 1;

enum class OutputFile : std::size_t {
    Positions,
    Scalars,
    ScalarValidity,
    Applicability,
    L1,
    Rotations,
    Translations,
    FitStatus,
    FitRmsd,
    Atoms,
    Bonds,
    Frames,
    Count,
};

constexpr std::size_t kOutputFileCount = static_cast<std::size_t>(OutputFile::Count);

const std::array<FieldKind, 5> kAngleFields{
    FieldKind::TauNCAC,
    FieldKind::AngleNCACB,
    FieldKind::AngleCBCAC,
    FieldKind::AngleCprevNCA,
    FieldKind::AngleCACNnext,
};

const std::array<FieldKind, 5> kAngleValidFields{
    FieldKind::TauNCACValid,
    FieldKind::AngleNCACBValid,
    FieldKind::AngleCBCACValid,
    FieldKind::AngleCprevNCAValid,
    FieldKind::AngleCACNnextValid,
};

constexpr std::array<const char*, kOutputFileCount> kOutputFiles{
    "positions.npy",
    "scalars.npy",
    "scalar_valid.npy",
    "applicability.npy",
    "l1.npy",
    "rotations.npy",
    "translations.npy",
    "fit_status.npy",
    "fit_rmsd_A.npy",
    "atoms.csv",
    "bonds.csv",
    "frames.csv",
};

bool checkedProduct(std::initializer_list<std::size_t> factors, std::size_t* result) {
    std::size_t product = 1;
    for (const std::size_t factor : factors) {
        if (factor != 0 && product > std::numeric_limits<std::size_t>::max() / factor) {
            return false;
        }
        product *= factor;
    }
    *result = product;
    return true;
}

QString fieldName(FieldKind kind) {
    return QString::fromLatin1(FieldSpecFor(kind).stem);
}

const model::NpyColumn* requireColumn(const model::QtConformationSnapshot& snapshot,
                                      FieldKind kind,
                                      std::size_t expectedRows,
                                      int expectedColumns,
                                      QString* error) {
    if (!snapshot.has(kind)) {
        *error = QStringLiteral("required field is absent: %1").arg(fieldName(kind));
        return nullptr;
    }

    const model::NpyColumn& column = snapshot.column(kind);
    if (column.rows < 0 || column.cols < 0 || static_cast<std::size_t>(column.rows) != expectedRows
        || column.cols != expectedColumns || column.data.size() != expectedRows * static_cast<std::size_t>(expectedColumns)) {
        *error = QStringLiteral("%1 has shape (%2,%3), expected (%4,%5)")
                     .arg(fieldName(kind))
                     .arg(column.rows)
                     .arg(column.cols)
                     .arg(static_cast<qulonglong>(expectedRows))
                     .arg(expectedColumns);
        return nullptr;
    }
    return &column;
}

bool readMask(const model::NpyColumn& column,
              std::size_t row,
              int component,
              const QString& name,
              std::uint8_t* value,
              QString* error) {
    const double raw = column.row(row)[component];
    if (raw != 0.0 && raw != 1.0) {
        *error = QStringLiteral("%1 contains a non-binary value at row %2").arg(name).arg(static_cast<qulonglong>(row));
        return false;
    }
    *value = raw == 1.0 ? 1 : 0;
    return true;
}

bool readFiniteFloat(const model::NpyColumn& column,
                     std::size_t row,
                     int component,
                     const QString& name,
                     float* value,
                     QString* error) {
    const double raw = column.row(row)[component];
    const float narrowed = static_cast<float>(raw);
    if (!std::isfinite(raw) || !std::isfinite(narrowed)) {
        *error = QStringLiteral("%1 contains a nonfinite value at row %2").arg(name).arg(static_cast<qulonglong>(row));
        return false;
    }
    *value = narrowed;
    return true;
}

bool readMaskedFloat(const model::NpyColumn& values,
                     std::size_t valueRow,
                     int valueComponent,
                     const model::NpyColumn& validity,
                     std::size_t validityRow,
                     int validityComponent,
                     const QString& valueName,
                     const QString& validityName,
                     float* value,
                     std::uint8_t* valid,
                     QString* error) {
    if (!readMask(validity, validityRow, validityComponent, validityName, valid, error)) {
        return false;
    }
    if (*valid == 0) {
        *value = 0.0f;
        return true;
    }
    return readFiniteFloat(values, valueRow, valueComponent, valueName, value, error);
}

bool readInteger(const model::NpyColumn& column, std::size_t row, const QString& name, int* value, QString* error) {
    const double raw = column.row(row)[0];
    if (!std::isfinite(raw) || std::trunc(raw) != raw || raw < static_cast<double>(std::numeric_limits<int>::min())
        || raw > static_cast<double>(std::numeric_limits<int>::max())) {
        *error = QStringLiteral("%1 contains a non-integer value at row %2").arg(name).arg(static_cast<qulonglong>(row));
        return false;
    }
    *value = static_cast<int>(raw);
    return true;
}

QString csvField(QString value) {
    value = value.trimmed();
    if (!value.contains(QLatin1Char(',')) && !value.contains(QLatin1Char('"')) && !value.contains(QLatin1Char('\n'))
        && !value.contains(QLatin1Char('\r'))) {
        return value;
    }
    value.replace(QStringLiteral("\""), QStringLiteral("\"\""));
    return QLatin1Char('"') + value + QLatin1Char('"');
}

bool writeTextFile(const QString& path, const QString& text, QString* error) {
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        *error = QStringLiteral("could not open %1: %2").arg(path, file.errorString());
        return false;
    }
    const QByteArray bytes = text.toUtf8();
    if (file.write(bytes) != bytes.size()) {
        *error = QStringLiteral("could not write %1: %2").arg(path, file.errorString());
        file.cancelWriting();
        return false;
    }
    if (!file.commit()) {
        *error = QStringLiteral("could not commit %1: %2").arg(path, file.errorString());
        return false;
    }
    return true;
}

template <typename Enum>
QString enumNumber(Enum value) {
    return QString::number(static_cast<int>(value));
}

QString atomsCsv(const model::QtProtein& protein, const std::vector<int>& atomRoles, const std::vector<int>& hybridisations) {
    QString csv;
    QTextStream stream(&csv);
    stream << "atom_index,chain_id,residue_number,insertion_code,residue_index,"
              "previous_residue_index,next_residue_index,amber_atom_name,bmrb_atom_name,"
              "parent_atom_index,element,iupac_atom,iupac_residue,amber_residue,"
              "terminal_state,residue_variant,formal_charge,polar_h_kind,atom_role,"
              "hybridisation,backbone_role,planar_group,locant,ring_position,aromatic,"
              "exchangeable\n";

    for (std::size_t index = 0; index < protein.atomCount(); ++index) {
        const model::QtAtom& atom = protein.atom(index);
        const model::QtResidue& residue = protein.residue(static_cast<std::size_t>(atom.residueIndex));
        const model::QtAtomNames& atomNames = protein.atomNames(index);
        const model::QtResidueNames& residueNames = protein.residueNames(static_cast<std::size_t>(atom.residueIndex));

        const QStringList fields{
            QString::number(atom.atomIndex),
            csvField(residue.address.chainId),
            QString::number(residue.address.residueNumber),
            csvField(residue.address.insertionCode),
            QString::number(atom.residueIndex),
            QString::number(residue.prevResidueIndex),
            QString::number(residue.nextResidueIndex),
            csvField(atomNames.amber),
            csvField(atomNames.bmrb),
            QString::number(atom.parentAtomIndex),
            QString::number(atom.AtomicNumber()),
            csvField(atomNames.iupac),
            csvField(residueNames.iupac),
            csvField(residueNames.amber),
            enumNumber(residue.terminalState),
            QString::number(residue.protonationVariantIndex),
            QString::number(atom.formalCharge),
            enumNumber(atom.polarH),
            QString::number(atomRoles[index]),
            QString::number(hybridisations[index]),
            enumNumber(atom.backboneRole),
            enumNumber(atom.planarGroup),
            enumNumber(atom.locant),
            enumNumber(atom.ringPositionPrimary),
            atom.aromatic ? QStringLiteral("1") : QStringLiteral("0"),
            atom.isExchangeable ? QStringLiteral("1") : QStringLiteral("0"),
        };
        stream << fields.join(QLatin1Char(',')) << QLatin1Char('\n');
    }
    return csv;
}

QString bondsCsv(const model::QtProtein& protein) {
    QString csv;
    QTextStream stream(&csv);
    stream << "bond_index,atom_a,atom_b,order,category\n";
    for (std::size_t index = 0; index < protein.bondCount(); ++index) {
        const model::QtBond& bond = protein.bond(index);
        stream << bond.bondIndex << QLatin1Char(',') << bond.atomIndexA << QLatin1Char(',') << bond.atomIndexB
               << QLatin1Char(',') << static_cast<int>(bond.order) << QLatin1Char(',') << static_cast<int>(bond.category)
               << QLatin1Char('\n');
    }
    return csv;
}

QString
framesCsv(bool trajectory, const std::vector<std::uint64_t>& originalFrames, const std::vector<double>& timesPicoseconds) {
    QString csv;
    QTextStream stream(&csv);
    stream << "frame_row,frame_id,original_frame_index,time_ps\n";
    for (std::size_t frame = 0; frame < originalFrames.size(); ++frame) {
        stream << static_cast<qulonglong>(frame) << QLatin1Char(',');
        if (!trajectory) {
            stream << "static,,\n";
            continue;
        }
        stream << QStringLiteral("frame_%1").arg(static_cast<qulonglong>(originalFrames[frame]), 6, 10, QLatin1Char('0'))
               << QLatin1Char(',') << static_cast<qulonglong>(originalFrames[frame]) << QLatin1Char(',')
               << QString::number(timesPicoseconds[frame], 'g', 17) << QLatin1Char('\n');
    }
    return csv;
}

}  // namespace

ModelInputExporter::ModelInputExporter(const model::QtProtein* protein,
                                       model::Conformation* conformation,
                                       QString outputDirectory)
    : protein_(protein)
    , conformation_(conformation)
    , outputDirectory_(std::move(outputDirectory)) {}

ModelInputExporter::~ModelInputExporter() {
    if (phase_ != Phase::Finished)
        removeWrittenFiles();
}

bool ModelInputExporter::prepare(QString* error) {
    if (prepared_) {
        *error = QStringLiteral("model-input export was already prepared");
        return false;
    }
    if (!protein_ || !conformation_) {
        *error = QStringLiteral("no complete Reader load is available");
        return false;
    }

    const QFileInfo outputInfo(outputDirectory_);
    if (!outputInfo.isAbsolute() || !outputInfo.exists() || !outputInfo.isDir()) {
        *error = QStringLiteral("output_directory must name an existing absolute directory");
        return false;
    }
    QDir output(outputInfo.canonicalFilePath());
    if (output.path().isEmpty()) {
        *error = QStringLiteral("output_directory could not be resolved");
        return false;
    }
    const QDir::Filters entries = QDir::AllEntries | QDir::NoDotAndDotDot | QDir::Hidden | QDir::System;
    if (!output.entryList(entries).isEmpty()) {
        *error = QStringLiteral("output_directory must be empty");
        return false;
    }
    outputDirectory_ = output.path();

    frameCount_ = conformation_->frameCount();
    atomCount_ = protein_->atomCount();
    bondCount_ = protein_->bondCount();
    trajectory_ = conformation_->asTrajectory() != nullptr;
    if (frameCount_ == 0 || atomCount_ == 0 || conformation_->protein() != protein_) {
        *error = QStringLiteral("loaded protein and conformation axes are incomplete");
        return false;
    }
    if (!trajectory_) {
        conformation_->requestSnapshot(0);
        const std::shared_ptr<const model::QtConformationSnapshot> snapshot = conformation_->snapshot(0);
        if (!snapshot) {
            *error = QStringLiteral("static position snapshot is unavailable");
            return false;
        }
        if (!requireColumn(*snapshot, FieldKind::Pos, atomCount_, 3, error)) {
            return false;
        }
    }

    residueForAtom_.resize(atomCount_);
    std::int64_t formalCharge = 0;
    for (std::size_t index = 0; index < atomCount_; ++index) {
        const model::QtAtom& atom = protein_->atom(index);
        if (atom.atomIndex != static_cast<std::int32_t>(index) || atom.residueIndex < 0
            || static_cast<std::size_t>(atom.residueIndex) >= protein_->residueCount()) {
            *error = QStringLiteral("atom identity axis is inconsistent at row %1").arg(static_cast<qulonglong>(index));
            return false;
        }
        residueForAtom_[index] = static_cast<std::size_t>(atom.residueIndex);
        formalCharge += atom.formalCharge;
        if (atom.IsBackboneNitrogen() || atom.IsBackboneAlphaCarbon() || atom.IsBackboneCarbonylCarbon()
            || atom.IsBackboneCarbonylOxygen()) {
            fitAtomIndices_.push_back(index);
        }
    }
    totalFormalCharge_ = static_cast<float>(formalCharge);

    for (std::size_t index = 0; index < protein_->residueCount(); ++index) {
        if (protein_->residue(index).residueIndex != static_cast<std::int32_t>(index)) {
            *error = QStringLiteral("residue identity axis is inconsistent at row %1").arg(static_cast<qulonglong>(index));
            return false;
        }
    }
    for (std::size_t index = 0; index < bondCount_; ++index) {
        const model::QtBond& bond = protein_->bond(index);
        if (bond.bondIndex != static_cast<std::int32_t>(index) || bond.atomIndexA < 0 || bond.atomIndexB < 0
            || static_cast<std::size_t>(bond.atomIndexA) >= atomCount_
            || static_cast<std::size_t>(bond.atomIndexB) >= atomCount_) {
            *error = QStringLiteral("bond axis is inconsistent at row %1").arg(static_cast<qulonglong>(index));
            return false;
        }
    }
    if (trajectory_ && fitAtomIndices_.size() < 3) {
        *error = QStringLiteral("trajectory alignment needs at least three typed NCACO atoms");
        return false;
    }

    std::size_t positionValues = 0;
    std::size_t scalarValues = 0;
    std::size_t applicabilityValues = 0;
    if (!checkedProduct({frameCount_, atomCount_, 3}, &positionValues)
        || !checkedProduct({frameCount_, atomCount_, kScalarWidth}, &scalarValues)
        || !checkedProduct({frameCount_, atomCount_, kApplicabilityWidth}, &applicabilityValues)) {
        *error = QStringLiteral("model-input array dimensions overflow");
        return false;
    }

    sourcePositions_.frameCount = frameCount_;
    sourcePositions_.atomCount = atomCount_;
    sourcePositions_.values.resize(frameCount_ * atomCount_);
    originalFrameIndices_.resize(frameCount_);
    frameTimesPicoseconds_.resize(frameCount_);
    positions_.resize(positionValues);
    scalars_.assign(scalarValues, 0.0f);
    scalarValid_.assign(scalarValues, 0);
    applicability_.assign(applicabilityValues, 0.0f);
    l1_.assign(positionValues * kL1Width, 0.0f);
    rotations_.resize(frameCount_ * 9);
    translations_.resize(frameCount_ * 3);
    fitStatus_.resize(frameCount_);
    fitRmsdAngstrom_.resize(frameCount_);
    atomRoles_.resize(atomCount_);
    hybridisations_.resize(atomCount_);

    prepared_ = true;
    qCInfo(cModelInput).noquote() << "model-input export prepared"
                                  << "| frames=" << static_cast<qulonglong>(frameCount_)
                                  << "| atoms=" << static_cast<qulonglong>(atomCount_)
                                  << "| bonds=" << static_cast<qulonglong>(bondCount_) << "| output=" << outputDirectory_;
    return true;
}

std::optional<ModelInputExporter::Result> ModelInputExporter::advance() {
    if (!prepared_)
        return finishFailure(QStringLiteral("model-input export was not prepared"));
    if (phase_ == Phase::Finished)
        return finishFailure(QStringLiteral("model-input export is already finished"));
    if (cancelRequested_)
        return finishFailure(QStringLiteral("model-input export cancelled"), true);
    if (!conformation_)
        return finishFailure(QStringLiteral("loaded conformation was released"));

    try {
        QString error;
        switch (phase_) {
        case Phase::Positions:
            if (!collectPositionFrame(nextFrame_, &error))
                return finishFailure(error);
            ++nextFrame_;
            if (nextFrame_ == frameCount_) {
                nextFrame_ = 0;
                phase_ = Phase::Alignment;
            }
            return std::nullopt;

        case Phase::Alignment:
            if (!buildAlignment(&error))
                return finishFailure(error);
            phase_ = Phase::Features;
            return std::nullopt;

        case Phase::Features:
            if (!projectFeatureFrame(nextFrame_, &error))
                return finishFailure(error);
            ++nextFrame_;
            if (nextFrame_ == frameCount_)
                phase_ = Phase::Files;
            return std::nullopt;

        case Phase::Files:
            if (!writeNextFile(&error))
                return finishFailure(error);
            ++nextFile_;
            if (nextFile_ == kOutputFileCount)
                return finishSuccess();
            return std::nullopt;

        case Phase::Finished:
            break;
        }
    } catch (const std::exception& error) {
        return finishFailure(QStringLiteral("model-input export failed: %1").arg(QString::fromLocal8Bit(error.what())));
    } catch (...) {
        return finishFailure(QStringLiteral("model-input export failed"));
    }
    return finishFailure(QStringLiteral("model-input export entered an invalid phase"));
}

bool ModelInputExporter::collectPositionFrame(std::size_t frame, QString* error) {
    if (frame >= frameCount_) {
        *error = QStringLiteral("position frame is outside the loaded trajectory");
        return false;
    }
    originalFrameIndices_[frame] = static_cast<std::uint64_t>(conformation_->originalFrameIndex(frame));
    frameTimesPicoseconds_[frame] = conformation_->timePicoseconds(frame);
    if (trajectory_ && !std::isfinite(frameTimesPicoseconds_[frame])) {
        *error = QStringLiteral("frame %1 has a nonfinite time").arg(static_cast<qulonglong>(frame));
        return false;
    }

    for (std::size_t atom = 0; atom < atomCount_; ++atom) {
        const model::Vec3 position = conformation_->atomPosition(frame, atom);
        if (!position.allFinite()) {
            *error = QStringLiteral("position is nonfinite at frame %1 atom %2")
                         .arg(static_cast<qulonglong>(frame))
                         .arg(static_cast<qulonglong>(atom));
            return false;
        }
        sourcePositions_.values[frame * atomCount_ + atom] = position;
    }
    return true;
}

bool ModelInputExporter::buildAlignment(QString* error) {
    if (trajectory_) {
        alignment_ = model::BuildScientificAlignment(sourcePositions_, fitAtomIndices_);
        if (!alignment_.ok) {
            *error = QStringLiteral("scientific NCACO alignment failed: %1").arg(alignment_.error);
            return false;
        }
    } else {
        alignment_.ok = true;
        alignment_.fitAtomIndices = fitAtomIndices_;
        alignment_.mean.converged = true;
        alignment_.mean.iterations = 0;
        alignment_.referencePositions.reserve(fitAtomIndices_.size());
        for (const std::size_t atom : fitAtomIndices_)
            alignment_.referencePositions.push_back(sourcePositions_.at(0, atom));
        model::ScientificFrameFit fit;
        fit.status = model::ScientificFitStatus::Valid;
        fit.rotation = model::Mat3::Identity();
        fit.translation = model::Vec3::Zero();
        fit.rmsdAngstrom = 0.0;
        alignment_.frameFits.push_back(fit);
    }

    if (alignment_.frameFits.size() != frameCount_) {
        *error = QStringLiteral("alignment transform count differs from frame count");
        return false;
    }
    for (std::size_t frame = 0; frame < frameCount_; ++frame) {
        const model::ScientificFrameFit& fit = alignment_.frameFits[frame];
        if (!fit.valid() || !fit.rotation.allFinite() || !fit.translation.allFinite() || !std::isfinite(fit.rmsdAngstrom)) {
            *error = QStringLiteral("alignment is invalid at frame %1").arg(static_cast<qulonglong>(frame));
            return false;
        }
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                rotations_[frame * 9 + static_cast<std::size_t>(row * 3 + column)] = fit.rotation(row, column);
            }
        }
        translations_[frame * 3] = fit.translation.x();
        translations_[frame * 3 + 1] = fit.translation.y();
        translations_[frame * 3 + 2] = fit.translation.z();
        fitStatus_[frame] = static_cast<std::uint8_t>(fit.status);
        fitRmsdAngstrom_[frame] = fit.rmsdAngstrom;
    }

    qCInfo(cModelInput).noquote() << (trajectory_ ? "scientific NCACO alignment complete"
                                                  : "static identity transform selected")
                                  << "| frames=" << static_cast<qulonglong>(frameCount_)
                                  << "| fit_atoms=" << static_cast<qulonglong>(fitAtomIndices_.size());
    return true;
}

bool ModelInputExporter::projectFeatureFrame(std::size_t frame, QString* error) {
    conformation_->requestSnapshot(frame);
    const std::shared_ptr<const model::QtConformationSnapshot> snapshot = conformation_->snapshot(frame);
    if (!snapshot) {
        *error = QStringLiteral("calculator snapshot is unavailable at frame %1").arg(static_cast<qulonglong>(frame));
        return false;
    }
    if (trajectory_ && snapshot->frameIndex() != originalFrameIndices_[frame]) {
        *error = QStringLiteral("snapshot identity differs at frame row %1").arg(static_cast<qulonglong>(frame));
        return false;
    }

    const std::size_t residueCount = protein_->residueCount();
    std::array<const model::NpyColumn*, 5> angles{};
    std::array<const model::NpyColumn*, 5> angleValid{};
    for (std::size_t index = 0; index < angles.size(); ++index) {
        angles[index] = requireColumn(*snapshot, kAngleFields[index], residueCount, 1, error);
        if (!angles[index])
            return false;
        angleValid[index] = requireColumn(*snapshot, kAngleValidFields[index], residueCount, 1, error);
        if (!angleValid[index])
            return false;
    }

    const model::NpyColumn* torsionCos =
        requireColumn(*snapshot, FieldKind::DSSPTorsionCos, residueCount, kTorsionComponentCount, error);
    const model::NpyColumn* torsionSin =
        requireColumn(*snapshot, FieldKind::DSSPTorsionSin, residueCount, kTorsionComponentCount, error);
    const model::NpyColumn* torsionValid =
        requireColumn(*snapshot, FieldKind::DSSPTorsionValid, residueCount, kTorsionComponentCount, error);
    const model::NpyColumn* omegaCos = requireColumn(*snapshot, FieldKind::OmegaCos, residueCount, 1, error);
    const model::NpyColumn* omegaSin = requireColumn(*snapshot, FieldKind::OmegaSin, residueCount, 1, error);
    const model::NpyColumn* omegaIsXpro = requireColumn(*snapshot, FieldKind::OmegaIsXpro, residueCount, 1, error);
    const model::NpyColumn* omegaValid = requireColumn(*snapshot, FieldKind::OmegaValid, residueCount, 1, error);
    const model::NpyColumn* pyramidalization = requireColumn(*snapshot, FieldKind::Pyramidalization, atomCount_, 1, error);
    const model::NpyColumn* pyramidalizationValid =
        requireColumn(*snapshot, FieldKind::PyramidalizationValid, atomCount_, 1, error);
    const model::NpyColumn* dsspObserved = requireColumn(*snapshot, FieldKind::DSSPObserved, atomCount_, 1, error);
    const model::NpyColumn* dsspSs8 = requireColumn(*snapshot, FieldKind::DSSPSs8, atomCount_, kDsspStateCount, error);
    const model::NpyColumn* dsspHbond =
        requireColumn(*snapshot, FieldKind::DSSPHBondEnergy, atomCount_, kDsspHbondCount, error);
    const model::NpyColumn* dsspPpii = requireColumn(*snapshot, FieldKind::DSSPPpii, atomCount_, 1, error);
    const model::NpyColumn* enrichmentRole = requireColumn(*snapshot, FieldKind::EnrichmentRole, atomCount_, 1, error);
    const model::NpyColumn* enrichmentHybridisation =
        requireColumn(*snapshot, FieldKind::EnrichmentHybridisationClass, atomCount_, 1, error);
    if (!torsionCos || !torsionSin || !torsionValid || !omegaCos || !omegaSin || !omegaIsXpro || !omegaValid
        || !pyramidalization || !pyramidalizationValid || !dsspObserved || !dsspSs8 || !dsspHbond || !dsspPpii
        || !enrichmentRole || !enrichmentHybridisation) {
        return false;
    }

    const model::ScientificFrameFit& fit = alignment_.frameFits[frame];
    model::Vec3 centroid = model::Vec3::Zero();
    for (std::size_t atom = 0; atom < atomCount_; ++atom) {
        centroid += fit.rotation * sourcePositions_.at(frame, atom) + fit.translation;
    }
    centroid /= static_cast<double>(atomCount_);
    for (std::size_t atom = 0; atom < atomCount_; ++atom) {
        const model::Vec3 centered = fit.rotation * sourcePositions_.at(frame, atom) + fit.translation - centroid;
        const std::size_t offset = (frame * atomCount_ + atom) * 3;
        for (int component = 0; component < 3; ++component) {
            const float value = static_cast<float>(centered(component));
            if (!std::isfinite(value)) {
                *error = QStringLiteral("centred position is outside float32 range at frame %1 atom %2")
                             .arg(static_cast<qulonglong>(frame))
                             .arg(static_cast<qulonglong>(atom));
                return false;
            }
            positions_[offset + static_cast<std::size_t>(component)] = value;
        }
    }

    std::vector<model::Vec3> cbResidual(residueCount, model::Vec3::Zero());
    std::vector<std::uint8_t> cbValid(residueCount, 0);
    for (std::size_t residueIndex = 0; residueIndex < residueCount; ++residueIndex) {
        const model::QtResidue& residue = protein_->residue(residueIndex);
        if (!residue.HasN() || !residue.HasCA() || !residue.HasC() || !residue.HasCB()) {
            continue;
        }
        const std::array<std::int32_t, 4> atomIndices{residue.N, residue.CA, residue.C, residue.CB};
        bool indicesValid = true;
        for (const std::int32_t atom : atomIndices) {
            indicesValid = indicesValid && atom >= 0 && static_cast<std::size_t>(atom) < atomCount_;
        }
        if (!indicesValid) {
            *error =
                QStringLiteral("residue %1 has an out-of-range backbone atom cache").arg(static_cast<qulonglong>(residueIndex));
            return false;
        }

        const model::Vec3& nPosition = sourcePositions_.at(frame, static_cast<std::size_t>(residue.N));
        const model::Vec3& caPosition = sourcePositions_.at(frame, static_cast<std::size_t>(residue.CA));
        const model::Vec3& cPosition = sourcePositions_.at(frame, static_cast<std::size_t>(residue.C));
        const model::Vec3& observedCb = sourcePositions_.at(frame, static_cast<std::size_t>(residue.CB));
        const model::Vec3 n = nPosition - caPosition;
        const model::Vec3 c = cPosition - caPosition;
        const model::Vec3 idealCb = caPosition - 0.58273431 * c.cross(n) - 0.56802827 * n - 0.54067466 * c;
        cbResidual[residueIndex] = fit.rotation * (observedCb - idealCb);
        if (!cbResidual[residueIndex].allFinite()) {
            *error = QStringLiteral("C-beta residual is nonfinite at residue %1").arg(static_cast<qulonglong>(residueIndex));
            return false;
        }
        cbValid[residueIndex] = 1;
    }

    for (std::size_t atom = 0; atom < atomCount_; ++atom) {
        int atomRole = 0;
        int hybridisation = 0;
        if (!readInteger(*enrichmentRole, atom, fieldName(FieldKind::EnrichmentRole), &atomRole, error)
            || !readInteger(*enrichmentHybridisation,
                            atom,
                            fieldName(FieldKind::EnrichmentHybridisationClass),
                            &hybridisation,
                            error)) {
            return false;
        }
        if (frame == 0) {
            atomRoles_[atom] = atomRole;
            hybridisations_[atom] = hybridisation;
        } else if (atomRoles_[atom] != atomRole || hybridisations_[atom] != hybridisation) {
            *error = QStringLiteral("fixed atom categories change at frame %1 atom %2")
                         .arg(static_cast<qulonglong>(frame))
                         .arg(static_cast<qulonglong>(atom));
            return false;
        }

        const std::size_t residue = residueForAtom_[atom];
        const std::size_t scalarBase = (frame * atomCount_ + atom) * kScalarWidth;
        const auto masked = [&](std::size_t outputChannel,
                                const model::NpyColumn& values,
                                std::size_t valueRow,
                                int valueComponent,
                                const model::NpyColumn& validity,
                                std::size_t validityRow,
                                int validityComponent,
                                FieldKind valueKind,
                                FieldKind validityKind) {
            return readMaskedFloat(values,
                                   valueRow,
                                   valueComponent,
                                   validity,
                                   validityRow,
                                   validityComponent,
                                   fieldName(valueKind),
                                   fieldName(validityKind),
                                   &scalars_[scalarBase + outputChannel],
                                   &scalarValid_[scalarBase + outputChannel],
                                   error);
        };
        const auto required =
            [&](std::size_t outputChannel, const model::NpyColumn& values, std::size_t row, int component, FieldKind kind) {
                const bool ok =
                    readFiniteFloat(values, row, component, fieldName(kind), &scalars_[scalarBase + outputChannel], error);
                if (ok) {
                    scalarValid_[scalarBase + outputChannel] = 1;
                }
                return ok;
            };

        for (std::size_t index = 0; index < angles.size(); ++index) {
            if (!masked(channelIndex(ScalarChannel::BackboneAngleFirst) + index,
                        *angles[index],
                        residue,
                        0,
                        *angleValid[index],
                        residue,
                        0,
                        kAngleFields[index],
                        kAngleValidFields[index])) {
                return false;
            }
        }
        for (std::size_t component = 0; component < kTorsionComponentCount; ++component) {
            if (!masked(channelIndex(ScalarChannel::TorsionCosFirst) + component,
                        *torsionCos,
                        residue,
                        static_cast<int>(component),
                        *torsionValid,
                        residue,
                        static_cast<int>(component),
                        FieldKind::DSSPTorsionCos,
                        FieldKind::DSSPTorsionValid)) {
                return false;
            }
        }
        for (std::size_t component = 0; component < kTorsionComponentCount; ++component) {
            if (!masked(channelIndex(ScalarChannel::TorsionSinFirst) + component,
                        *torsionSin,
                        residue,
                        static_cast<int>(component),
                        *torsionValid,
                        residue,
                        static_cast<int>(component),
                        FieldKind::DSSPTorsionSin,
                        FieldKind::DSSPTorsionValid)) {
                return false;
            }
        }
        if (!masked(channelIndex(ScalarChannel::OmegaCos),
                    *omegaCos,
                    residue,
                    0,
                    *omegaValid,
                    residue,
                    0,
                    FieldKind::OmegaCos,
                    FieldKind::OmegaValid)
            || !masked(channelIndex(ScalarChannel::OmegaSin),
                       *omegaSin,
                       residue,
                       0,
                       *omegaValid,
                       residue,
                       0,
                       FieldKind::OmegaSin,
                       FieldKind::OmegaValid)
            || !required(channelIndex(ScalarChannel::OmegaIsXpro), *omegaIsXpro, residue, 0, FieldKind::OmegaIsXpro)
            || !masked(channelIndex(ScalarChannel::Pyramidalization),
                       *pyramidalization,
                       atom,
                       0,
                       *pyramidalizationValid,
                       atom,
                       0,
                       FieldKind::Pyramidalization,
                       FieldKind::PyramidalizationValid)) {
            return false;
        }
        for (std::size_t component = 0; component < kDsspStateCount; ++component) {
            if (!required(channelIndex(ScalarChannel::DsspStateFirst) + component,
                          *dsspSs8,
                          atom,
                          static_cast<int>(component),
                          FieldKind::DSSPSs8)) {
                return false;
            }
        }
        for (std::size_t component = 0; component < kDsspHbondCount; ++component) {
            if (!masked(channelIndex(ScalarChannel::DsspHbondFirst) + component,
                        *dsspHbond,
                        atom,
                        static_cast<int>(component),
                        *dsspObserved,
                        atom,
                        0,
                        FieldKind::DSSPHBondEnergy,
                        FieldKind::DSSPObserved)) {
                return false;
            }
        }
        if (!masked(channelIndex(ScalarChannel::DsspPpii),
                    *dsspPpii,
                    atom,
                    0,
                    *dsspObserved,
                    atom,
                    0,
                    FieldKind::DSSPPpii,
                    FieldKind::DSSPObserved)) {
            return false;
        }
        const std::size_t chargeChannel = channelIndex(ScalarChannel::TotalFormalCharge);
        scalars_[scalarBase + chargeChannel] = totalFormalCharge_;
        scalarValid_[scalarBase + chargeChannel] = 1;

        const std::size_t applicabilityBase = (frame * atomCount_ + atom) * kApplicabilityWidth;
        applicability_[applicabilityBase + channelIndex(ApplicabilityChannel::CbetaResidual)] =
            static_cast<float>(cbValid[residue]);
        for (std::size_t component = 0; component < kTorsionComponentCount; ++component) {
            applicability_[applicabilityBase + channelIndex(ApplicabilityChannel::TorsionFirst) + component] =
                static_cast<float>(scalarValid_[scalarBase + channelIndex(ScalarChannel::TorsionCosFirst) + component]);
        }
        applicability_[applicabilityBase + channelIndex(ApplicabilityChannel::Pyramidalization)] =
            static_cast<float>(scalarValid_[scalarBase + channelIndex(ScalarChannel::Pyramidalization)]);

        const std::size_t l1Base = (frame * atomCount_ + atom) * kL1Width * 3;
        if (cbValid[residue] != 0) {
            for (int component = 0; component < 3; ++component) {
                const float value = static_cast<float>(cbResidual[residue](component));
                if (!std::isfinite(value)) {
                    *error = QStringLiteral("C-beta residual is outside float32 range at frame %1 residue %2")
                                 .arg(static_cast<qulonglong>(frame))
                                 .arg(static_cast<qulonglong>(residue));
                    return false;
                }
                l1_[l1Base + static_cast<std::size_t>(component)] = value;
            }
        }
    }

    const std::size_t interval = std::max<std::size_t>(1, frameCount_ / 4);
    if (frame + 1 == frameCount_ || (frame + 1) % interval == 0) {
        qCInfo(cModelInput).noquote() << "model-input projection"
                                      << "| frames=" << static_cast<qulonglong>(frame + 1) << "/"
                                      << static_cast<qulonglong>(frameCount_);
    }
    return true;
}

bool ModelInputExporter::writeNextFile(QString* error) {
    if (nextFile_ >= kOutputFiles.size()) {
        *error = QStringLiteral("model-input file index is outside the contract");
        return false;
    }
    const QString path = QDir(outputDirectory_).filePath(QString::fromLatin1(kOutputFiles[nextFile_]));
    if (QFileInfo::exists(path)) {
        *error = QStringLiteral("refusing to replace existing output file: %1").arg(path);
        return false;
    }

    QtNpyWriter::Result written;
    bool ok = false;
    switch (static_cast<OutputFile>(nextFile_)) {
    case OutputFile::Positions:
        written = QtNpyWriter::WriteFloat32(path, {frameCount_, atomCount_, 3}, positions_);
        ok = written.ok;
        break;
    case OutputFile::Scalars:
        written = QtNpyWriter::WriteFloat32(path, {frameCount_, atomCount_, kScalarWidth}, scalars_);
        ok = written.ok;
        break;
    case OutputFile::ScalarValidity:
        written = QtNpyWriter::WriteUInt8(path, {frameCount_, atomCount_, kScalarWidth}, scalarValid_);
        ok = written.ok;
        break;
    case OutputFile::Applicability:
        written = QtNpyWriter::WriteFloat32(path, {frameCount_, atomCount_, kApplicabilityWidth}, applicability_);
        ok = written.ok;
        break;
    case OutputFile::L1:
        written = QtNpyWriter::WriteFloat32(path, {frameCount_, atomCount_, kL1Width, 3}, l1_);
        ok = written.ok;
        break;
    case OutputFile::Rotations:
        written = QtNpyWriter::WriteFloat64(path, {frameCount_, 3, 3}, rotations_);
        ok = written.ok;
        break;
    case OutputFile::Translations:
        written = QtNpyWriter::WriteFloat64(path, {frameCount_, 3}, translations_);
        ok = written.ok;
        break;
    case OutputFile::FitStatus:
        written = QtNpyWriter::WriteUInt8(path, {frameCount_}, fitStatus_);
        ok = written.ok;
        break;
    case OutputFile::FitRmsd:
        written = QtNpyWriter::WriteFloat64(path, {frameCount_}, fitRmsdAngstrom_);
        ok = written.ok;
        break;
    case OutputFile::Atoms:
        ok = writeTextFile(path, atomsCsv(*protein_, atomRoles_, hybridisations_), error);
        break;
    case OutputFile::Bonds:
        ok = writeTextFile(path, bondsCsv(*protein_), error);
        break;
    case OutputFile::Frames:
        ok = writeTextFile(path, framesCsv(trajectory_, originalFrameIndices_, frameTimesPicoseconds_), error);
        break;
    case OutputFile::Count:
        break;
    }

    if (!ok) {
        if (error->isEmpty())
            *error = written.error;
        return false;
    }
    writtenFiles_.push_back(path);
    return true;
}

ModelInputExporter::Result ModelInputExporter::finishSuccess() {
    phase_ = Phase::Finished;
    writtenFiles_.clear();
    Result result;
    result.ok = true;
    result.frameCount = frameCount_;
    result.atomCount = atomCount_;
    result.bondCount = bondCount_;
    result.trajectory = trajectory_;
    result.alignment = alignment_;
    qCInfo(cModelInput).noquote() << "model-input export complete"
                                  << "| frames=" << static_cast<qulonglong>(frameCount_)
                                  << "| atoms=" << static_cast<qulonglong>(atomCount_) << "| output=" << outputDirectory_;
    return result;
}

ModelInputExporter::Result ModelInputExporter::finishFailure(const QString& error, bool cancelled) {
    removeWrittenFiles();
    phase_ = Phase::Finished;
    Result result;
    result.cancelled = cancelled;
    result.error = error;
    result.frameCount = frameCount_;
    result.atomCount = atomCount_;
    result.bondCount = bondCount_;
    result.trajectory = trajectory_;
    return result;
}

void ModelInputExporter::removeWrittenFiles() {
    for (const QString& path : std::as_const(writtenFiles_)) {
        if (!QFile::remove(path)) {
            qCWarning(cModelInput).noquote() << "could not remove incomplete model-input output"
                                             << "| path=" << path;
        }
    }
    writtenFiles_.clear();
}

}  // namespace h5reader::io
