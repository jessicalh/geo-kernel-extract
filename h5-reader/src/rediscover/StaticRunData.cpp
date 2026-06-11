#include "StaticRunData.h"

#include "Catalog.h"
#include "CanonicalSpineGuard.h"
#include "SphericalBasis.h"

#include "../io/DftShieldingLoader.h"
#include "../io/QtFieldCatalog.gen.h"
#include "../io/QtNpyReader.h"
#include "../io/QtProteinLoader.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"

#include <QDir>
#include <QFileInfo>
#include <QLoggingCategory>

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cStaticRun, "h5reader.rediscover.static_run")

using ProducerPathMap = std::unordered_map<int, QString>;

QString fieldStem(const io::FieldSpec& spec) {
    return QString::fromUtf8(spec.stem.data(), static_cast<qsizetype>(spec.stem.size()));
}

ProducerPathMap enumerateProducerNpys(const QString& dir) {
    ProducerPathMap out;
    QDir d(dir);
    const QStringList npys = d.entryList(QStringList{QStringLiteral("*.npy")}, QDir::Files, QDir::Name);
    for (const QString& fname : npys) {
        const std::string stem = QFileInfo(fname).completeBaseName().toStdString();
        const std::optional<io::FieldKind> kind = io::FindFieldByStem(stem);
        if (!kind) continue;
        const io::FieldSpec& spec = io::FieldSpecFor(*kind);
        if (spec.group == io::FieldGroup::Topology || *kind == io::FieldKind::AtomsCategoryInfo)
            continue;
        out[static_cast<int>(*kind)] = d.filePath(fname);
    }
    return out;
}

bool loadStaticFieldArray(const ProducerPathMap& paths,
                          io::FieldKind kind,
                          bool required,
                          bool keepFloat,
                          const model::QtProtein& protein,
                          StaticNpyArray* out,
                          bool* present_out,
                          QString* err_out) {
    if (present_out) *present_out = false;
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    const auto pathIt = paths.find(static_cast<int>(kind));
    if (pathIt == paths.end()) {
        if (required && err_out) {
            *err_out = QStringLiteral("required static producer NPY is missing: %1.npy")
                           .arg(fieldStem(spec));
        }
        return !required;
    }
    const QString& path = pathIt->second;
    auto arr = io::QtNpyReader::ReadArrayWidened(path);
    if (!arr.ok) {
        if (err_out) *err_out = arr.error;
        return false;
    }

    StaticNpyArray s;
    s.stem = fieldStem(spec);
    s.path = path;
    s.rows = arr.rows;
    s.cols = arr.cols;
    if (spec.axis == io::NativeAxis::Protein && arr.cols == 1
        && spec.cols > 0 && arr.rows == static_cast<std::size_t>(spec.cols)) {
        s.rows = 1;
        s.cols = arr.rows;
    }
    s.dtype_descr = QString::fromStdString(arr.descr);
    s.values = std::move(arr.data);
    if (keepFloat) {
        s.floatValues.reserve(s.values.size());
        for (double v : s.values) s.floatValues.push_back(static_cast<float>(v));
    }

    if (spec.cols != -1 && s.cols != static_cast<std::size_t>(spec.cols)) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 columns; catalog %3 expects %4")
                           .arg(path)
                           .arg(s.cols)
                           .arg(fieldStem(spec))
                           .arg(spec.cols);
        }
        return false;
    }

    if (spec.axis == io::NativeAxis::Atom && s.rows != protein.atomCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 atoms")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.atomCount());
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::Residue && s.rows != protein.residueCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 residues")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.residueCount());
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::AromaticRing
        && s.rows != protein.topology().aromaticRingCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 aromatic rings")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.topology().aromaticRingCount());
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::SaturatedRing
        && s.rows != protein.topology().saturatedRingCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 saturated rings")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.topology().saturatedRingCount());
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::Protein && s.rows != 1) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; protein-axis fields must have 1 row")
                           .arg(path)
                           .arg(s.rows);
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::Ring && s.rows != protein.topology().ringCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 rings")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.topology().ringCount());
        }
        return false;
    }
    if (spec.axis == io::NativeAxis::RingMembership
        && s.rows != protein.topology().ringMembershipCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 ring memberships")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.topology().ringMembershipCount());
        }
        return false;
    }
    if (kind == io::FieldKind::MOPACTopologyBondOrdersFull
        && s.rows != protein.topology().bondCount()) {
        if (err_out) {
            *err_out = QStringLiteral("%1 has %2 rows; topology has %3 bonds")
                           .arg(path)
                           .arg(s.rows)
                           .arg(protein.topology().bondCount());
        }
        return false;
    }
    *out = std::move(s);
    if (present_out) *present_out = true;
    return true;
}

struct StaticProducerRequest {
    ArrayId id;
    bool required = false;
    bool keepFloat = false;
};

bool loadStaticArray(const ProducerPathMap& paths,
                     const StaticProducerRequest& request,
                     const model::QtProtein& protein,
                     StaticNpyArray* out,
                     bool* present_out,
                     QString* err_out) {
    const std::optional<io::FieldKind> kind = ProducerFieldFor(request.id);
    if (!kind) {
        if (err_out) *err_out = QStringLiteral("row_design ArrayId has no producer FieldKind");
        return false;
    }
    return loadStaticFieldArray(paths, *kind, request.required, request.keepFloat,
                                protein, out, present_out, err_out);
}

Mat3 unpackMat3(const StaticNpyArray& a, std::size_t atom) {
    Mat3 m = Mat3::Zero();
    if (atom >= a.rows || a.cols < 9) return m;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            m(r, c) = a.value(atom, static_cast<std::size_t>(r * 3 + c));
    return m;
}

bool loadStaticDft(const QString& poseDir,
                   const ProducerPathMap& paths,
                   const model::QtProtein& protein,
                   std::shared_ptr<const model::DftShieldingFrame>* out,
                   QString* err_out) {
    StaticNpyArray total, dia, para;
    bool present = false;
    const StaticProducerRequest totalReq{ArrayId::DftTotalRaw, true, false};
    const StaticProducerRequest diaReq{ArrayId::DftDiaRaw, true, false};
    const StaticProducerRequest paraReq{ArrayId::DftParaRaw, true, false};
    if (!loadStaticArray(paths, totalReq, protein, &total, &present, err_out)) return false;
    if (!loadStaticArray(paths, diaReq, protein, &dia, &present, err_out)) return false;
    if (!loadStaticArray(paths, paraReq, protein, &para, &present, err_out)) return false;
    if (total.cols < 9 || dia.cols < 9 || para.cols < 9) {
        if (err_out) *err_out = QStringLiteral("static ORCA NPYs must be Nx9 tensors");
        return false;
    }

    model::DftShieldingFrame frame;
    frame.valid = true;
    frame.atoms.resize(protein.atomCount());
    for (std::size_t atom = 0; atom < protein.atomCount(); ++atom) {
        model::DftAtomShielding& dst = frame.atoms[atom];
        dst.element = protein.atom(atom).element;
        dst.total_raw = unpackMat3(total, atom);
        dst.dia_raw = unpackMat3(dia, atom);
        dst.para_raw = unpackMat3(para, atom);
        dst.total = DecomposeLibrary(dst.total_raw);
        dst.dia = DecomposeLibrary(dst.dia_raw);
        dst.para = DecomposeLibrary(dst.para_raw);
    }

    QString validateErr;
    if (!io::ValidateDftFrame(frame, &protein, poseDir, &validateErr)) {
        if (err_out) *err_out = validateErr;
        return false;
    }
    *out = std::make_shared<const model::DftShieldingFrame>(std::move(frame));
    return true;
}

}  // namespace

std::optional<RunData> StaticRunData::Load(const QString& poseDir, QString* err_out) {
    if (!ValidateCanonical720Pose(poseDir, err_out)) return std::nullopt;

    const QFileInfo poseInfo(poseDir);
    if (!poseInfo.isDir()) {
        if (err_out) *err_out = QStringLiteral("static pose dir does not exist: %1").arg(poseDir);
        return std::nullopt;
    }
    const QString extractionManifest = QDir(poseDir).filePath(QStringLiteral("extraction_manifest.json"));
    if (!QFileInfo::exists(extractionManifest)) {
        if (err_out)
            *err_out = QStringLiteral("static pose dir is missing extraction_manifest.json: %1").arg(poseDir);
        return std::nullopt;
    }

    io::QtLoadResult loaded = io::QtProteinLoader::LoadPose(poseDir, extractionManifest);
    if (!loaded.ok) {
        if (err_out) *err_out = QStringLiteral("static pose load failed: %1").arg(loaded.error);
        return std::nullopt;
    }

    RunData run;
    run.protein = std::move(loaded.protein);
    run.conformation = std::move(loaded.conformation);
    run.manifest.kind = io::CalcsetManifest::Kind::SinglePose;
    run.manifest.dataset_id = QStringLiteral("720_static");
    run.manifest.protein_id = loaded.proteinId;
    run.manifest.calcset_root_abspath = QFileInfo(poseDir).absoluteFilePath();
    run.manifest.single_pose = io::CalcsetManifest::SinglePose{
        io::CalcsetManifest::PoseKind::Orca,
        QFileInfo(poseDir).absoluteFilePath(),
        QFileInfo(extractionManifest).absoluteFilePath()};

    const ProducerPathMap producerPaths = enumerateProducerNpys(poseDir);
    for (const auto& item : producerPaths) {
        const io::FieldKind kind = static_cast<io::FieldKind>(item.first);
        StaticNpyArray arr;
        QString arrErr;
        bool present = false;
        const bool keepFloat = kind == io::FieldKind::AIMNet2Aim;
        if (loadStaticFieldArray(producerPaths, kind, false, keepFloat, *run.protein,
                                 &arr, &present, &arrErr)) {
            if (present) run.producerArrays[static_cast<int>(kind)] = std::move(arr);
        } else if (!arrErr.isEmpty()) {
            if (err_out) *err_out = arrErr;
            return std::nullopt;
        }
    }
    const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
    if (!run.producerArray(io::FieldKind::Pos)) {
        if (err_out) *err_out = QStringLiteral("static pose is missing required pos.npy");
        return std::nullopt;
    }
    if (const std::optional<io::FieldKind> ff = ProducerFieldFor(ArrayId::Ff14sbCharge)) {
        if (const StaticNpyArray* q = run.producerArray(*ff)) {
            for (std::size_t atom = 0; atom < atomCount; ++atom)
                run.protein->setAtomPartialCharge(atom, q->value(atom));
        }
    }

    std::shared_ptr<const model::DftShieldingFrame> dftFrame;
    if (!loadStaticDft(poseDir, producerPaths, *run.protein, &dftFrame, err_out)) return std::nullopt;
    run.dft.Insert(0, std::move(dftFrame));
    run.frameMap = FrameMap::Static(0, true);

    qCInfo(cStaticRun).noquote() << "StaticRunData ready | protein=" << run.manifest.protein_id
                                 << "| atoms=" << atomCount
                                 << "| producer_arrays=" << run.producerArrays.size();
    return run;
}

}  // namespace h5reader::rediscover
