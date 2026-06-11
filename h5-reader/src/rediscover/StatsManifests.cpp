#include "StatsManifests.h"

#include "../io/QtFieldCatalog.gen.h"

#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QStringList>
#include <QTextStream>

namespace h5reader::rediscover {

namespace {

bool writeText(const QString& path, const QByteArray& bytes, QString* err_out) {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(path);
        return false;
    }
    if (f.write(bytes) != bytes.size()) {
        if (err_out) *err_out = QStringLiteral("short write to %1").arg(path);
        return false;
    }
    return true;
}

QString typeName(RowColType t) {
    switch (t) {
    case RowColType::String: return QStringLiteral("string");
    case RowColType::Int: return QStringLiteral("int");
    case RowColType::Double: return QStringLiteral("double");
    case RowColType::Bool: return QStringLiteral("bool");
    }
    return QStringLiteral("unknown");
}

QString axisName(RowNativeAxis a) {
    switch (a) {
    case RowNativeAxis::Row: return QStringLiteral("row");
    case RowNativeAxis::Protein: return QStringLiteral("protein");
    case RowNativeAxis::Frame: return QStringLiteral("frame");
    case RowNativeAxis::Atom: return QStringLiteral("atom");
    case RowNativeAxis::Residue: return QStringLiteral("residue");
    case RowNativeAxis::Ring: return QStringLiteral("ring");
    case RowNativeAxis::Bond: return QStringLiteral("bond");
    }
    return QStringLiteral("unknown");
}

std::size_t populatedCountAt(const RowDesignStats& stats, std::size_t idx) {
    return idx < stats.populatedCounts.size() ? stats.populatedCounts[idx] : 0;
}

QString fieldStem(io::FieldKind kind) {
    const io::FieldSpec& spec = io::FieldSpecFor(kind);
    return QString::fromUtf8(spec.stem.data(), static_cast<qsizetype>(spec.stem.size()));
}

QString npyName(io::FieldKind kind) {
    return fieldStem(kind) + QStringLiteral(".npy");
}

QString joinedNpyNames(std::initializer_list<io::FieldKind> kinds) {
    QStringList names;
    for (io::FieldKind kind : kinds) names << npyName(kind);
    return names.join(QStringLiteral(" + "));
}

QJsonArray schemaJson(const std::vector<RowColumnSpec>& schema, const RowDesignStats& stats) {
    QJsonArray cols;
    for (std::size_t i = 0; i < schema.size(); ++i) {
        const RowColumnSpec& c = schema[i];
        QJsonObject o;
        o.insert(QStringLiteral("name"), c.name);
        o.insert(QStringLiteral("type"), typeName(c.type));
        o.insert(QStringLiteral("unit"), c.unit);
        o.insert(QStringLiteral("irrep"), c.irrep);
        o.insert(QStringLiteral("native_axis"), axisName(c.nativeAxis));
        o.insert(QStringLiteral("time_varying"), c.timeVarying);
        o.insert(QStringLiteral("is_feature"), c.isFeature);
        o.insert(QStringLiteral("populated_count"), static_cast<double>(populatedCountAt(stats, i)));
        cols.push_back(o);
    }
    return cols;
}

QJsonArray columnProvenanceJson(const std::vector<RowColumnSpec>& schema, const RowDesignStats& stats) {
    QJsonArray cols;
    for (std::size_t i = 0; i < schema.size(); ++i) {
        const RowColumnSpec& c = schema[i];
        QJsonObject o;
        o.insert(QStringLiteral("name"), c.name);
        o.insert(QStringLiteral("native_axis"), axisName(c.nativeAxis));
        o.insert(QStringLiteral("time_varying"), c.timeVarying);
        o.insert(QStringLiteral("populated_count"), static_cast<double>(populatedCountAt(stats, i)));

        if (c.name == QStringLiteral("ring_bs_T0")
            || c.name == QStringLiteral("ring_bs_absT2")
            || c.name == QStringLiteral("ring_hm_T0")
            || c.name == QStringLiteral("ring_hm_absT2")) {
            o.insert(QStringLiteral("source"), c.name.startsWith(QStringLiteral("ring_bs"))
                                               ? npyName(io::FieldKind::BSShielding)
                                               : npyName(io::FieldKind::HMShielding));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("ring_geometry"));
            o.insert(QStringLiteral("physical_quantity"),
                     QStringLiteral("unit-current geometric l<=2 multipole; not shielding sigma"));
        } else if (c.name.startsWith(QStringLiteral("ring_jb_T0_"))) {
            o.insert(QStringLiteral("source"), npyName(io::FieldKind::BSPerTypeT0)
                                              + QStringLiteral(" + literature ring-current intensities"));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("ring_jb"));
            o.insert(QStringLiteral("physical_quantity"),
                     QStringLiteral("per-aromatic-type intensity-scaled T0 shielding attribution"));
        } else if (c.name.startsWith(QStringLiteral("ring_jb_"))) {
            o.insert(QStringLiteral("source"), joinedNpyNames({io::FieldKind::BSPerTypeT0,
                                                               io::FieldKind::BSPerTypeT1,
                                                               io::FieldKind::BSPerTypeT2})
                                              + QStringLiteral(" + literature ring-current intensities"));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("ring_jb"));
            o.insert(QStringLiteral("physical_quantity"),
                     QStringLiteral("intensity-scaled ring-current shielding sigma"));
        } else if (c.name.startsWith(QStringLiteral("ring_hm_jb_"))) {
            o.insert(QStringLiteral("source"), joinedNpyNames({io::FieldKind::HMPerTypeT0,
                                                               io::FieldKind::HMPerTypeT1,
                                                               io::FieldKind::HMPerTypeT2})
                                              + QStringLiteral(" + literature ring-current intensities"));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("ring_hm_jb"));
            o.insert(QStringLiteral("physical_quantity"),
                     QStringLiteral("intensity-scaled Haigh-Mallion ring-current shielding sigma"));
        } else if (c.name == QStringLiteral("ring_cos_phi")
                   || c.name == QStringLiteral("ring_sin_phi")
                   || c.name == QStringLiteral("ring_cyl_z")
                   || c.name == QStringLiteral("ring_cyl_rho")
                   || c.name == QStringLiteral("ring_angle_to_normal")
                   || c.name == QStringLiteral("nearest_ring_dist")
                   || c.name == QStringLiteral("nearest_ring_identity_ord")) {
            o.insert(QStringLiteral("source"), QStringLiteral("ring geometry cache / nearest ring center"));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("ring_geometry"));
            o.insert(QStringLiteral("physical_quantity"), QStringLiteral("geometric descriptor"));
        } else {
            o.insert(QStringLiteral("source"), QStringLiteral("row_design emitter"));
            o.insert(QStringLiteral("mechanism"), QStringLiteral("mixed"));
            o.insert(QStringLiteral("physical_quantity"), c.unit.isEmpty() ? QStringLiteral("dimensionless")
                                                                           : c.unit);
        }
        cols.push_back(o);
    }
    return cols;
}

}  // namespace

bool WriteRowDesignManifests(const QString& outDir,
                             const std::vector<RowColumnSpec>& schema,
                             const RowDesignStats& stats,
                             const RunData& run,
                             const ConditioningSpec& spec,
                             const RowDesignManifestContext& context,
                             QString* err_out) {
    QDir().mkpath(outDir);
    const QString datasetId = !context.datasetId.isEmpty()
                                  ? context.datasetId
                                  : (run.manifest.dataset_id.isEmpty() ? QStringLiteral("unknown")
                                                                       : run.manifest.dataset_id);
    QJsonObject base;
    base.insert(QStringLiteral("emitter"), QStringLiteral("row_design"));
    base.insert(QStringLiteral("schema_version"), 1);
    base.insert(QStringLiteral("dataset_id"), datasetId);
    base.insert(QStringLiteral("protein_id"), run.protein ? run.protein->proteinId() : QString());
    base.insert(QStringLiteral("pose_kind"), run.poseKind() == PoseKind::Static ? QStringLiteral("static")
                                                                                : QStringLiteral("trajectory"));
    base.insert(QStringLiteral("row_count"), static_cast<double>(stats.rows));
    base.insert(QStringLiteral("atom_count"), static_cast<double>(stats.atoms));
    base.insert(QStringLiteral("dft_rows"), static_cast<double>(stats.dftRows));
    base.insert(QStringLiteral("dft_present"), static_cast<double>(stats.dftPresent));
    base.insert(QStringLiteral("phi_psi_present"), static_cast<double>(stats.phiPsiPresent));
    base.insert(QStringLiteral("phi_psi_eligible"), static_cast<double>(stats.phiPsiEligible));
    base.insert(QStringLiteral("phi_psi_finite_eligible"), static_cast<double>(stats.phiPsiFiniteEligible));
    base.insert(QStringLiteral("embedding_present"), static_cast<double>(stats.embeddingPresent));
    base.insert(QStringLiteral("fixture"), context.fixture);
    base.insert(QStringLiteral("subsampling"), QStringLiteral("none"));
    base.insert(QStringLiteral("grain"), QStringLiteral("absolute shielding; sigma_iso T0 primary; T2 sidecar"));
    base.insert(QStringLiteral("conditioning_config_hash"), spec.configHash());
    base.insert(QStringLiteral("root_path"), context.rootPath);
    base.insert(QStringLiteral("lgs_path"), run.manifest.lgs_path_abspath);

    QJsonObject audit = base;
    audit.insert(QStringLiteral("columns"), schemaJson(schema, stats));
    if (!writeText(QDir(outDir).filePath(QStringLiteral("schema_audit.json")),
                   QJsonDocument(audit).toJson(QJsonDocument::Indented), err_out)) return false;

    QJsonObject prov = base;
    prov.insert(QStringLiteral("source_stem_audit_hash"), QStringLiteral("not_computed"));
    prov.insert(QStringLiteral("required_target_stems"), QJsonArray{
        fieldStem(io::FieldKind::OrcaTotal),
        fieldStem(io::FieldKind::OrcaDiamagnetic),
        fieldStem(io::FieldKind::OrcaParamagnetic),
    });
    prov.insert(QStringLiteral("columns"), columnProvenanceJson(schema, stats));
    prov.insert(QStringLiteral("ring_restore_note"),
                QStringLiteral("ring_bs_* and ring_hm_* are raw unit-current geometric multipoles; ring_jb_* columns are literature-intensity-scaled shielding sigma. Per-type order verified from RingTypeIndex: phe,tyr,trp6,trp5,trp9,his,hid,hie; pro is saturated and zero."));
    prov.insert(QStringLiteral("sidecars"), QJsonArray{});
    if (!writeText(QDir(outDir).filePath(QStringLiteral("column_provenance.json")),
                   QJsonDocument(prov).toJson(QJsonDocument::Indented), err_out)) return false;

    QJsonObject irrep = base;
    QJsonArray irreps;
    for (const RowColumnSpec& c : schema) {
        QJsonObject o;
        o.insert(QStringLiteral("name"), c.name);
        o.insert(QStringLiteral("irrep"), c.irrep);
        o.insert(QStringLiteral("sign_flip_legal"), c.signFlipLegal);
        irreps.push_back(o);
    }
    irrep.insert(QStringLiteral("columns"), irreps);
    if (!writeText(QDir(outDir).filePath(QStringLiteral("column_irrep_schema.json")),
                   QJsonDocument(irrep).toJson(QJsonDocument::Indented), err_out)) return false;

    QJsonObject nullSpec = base;
    nullSpec.insert(QStringLiteral("double_null"), QStringLiteral("NaN"));
    nullSpec.insert(QStringLiteral("int_null"), -1);
    nullSpec.insert(QStringLiteral("bool_null"), 0);
    nullSpec.insert(QStringLiteral("embedding_null"), QStringLiteral("NaN"));
    if (!writeText(QDir(outDir).filePath(QStringLiteral("null_spec.json")),
                   QJsonDocument(nullSpec).toJson(QJsonDocument::Indented), err_out)) return false;

    QJsonObject region = base;
    region.insert(QStringLiteral("rama_region_classifier"), QStringLiteral("row_design_v1_basin_grid"));
    region.insert(QStringLiteral("angle_convention"), QStringLiteral("IUPAC signed radians"));
    if (!writeText(QDir(outDir).filePath(QStringLiteral("region_spec.json")),
                   QJsonDocument(region).toJson(QJsonDocument::Indented), err_out)) return false;

    QFile columnSupport(QDir(outDir).filePath(QStringLiteral("column_support.csv")));
    if (!columnSupport.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write column_support.csv");
        return false;
    }
    QTextStream cs(&columnSupport);
    cs << "column,populated_count,row_count\n";
    for (std::size_t i = 0; i < schema.size(); ++i)
        cs << schema[i].name << "," << populatedCountAt(stats, i) << "," << stats.rows << "\n";

    QFile support(QDir(outDir).filePath(QStringLiteral("support_table.csv")));
    if (!support.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write support_table.csv");
        return false;
    }
    QTextStream s(&support);
    s << "metric,count\n";
    s << "rows," << stats.rows << "\n";
    s << "atoms," << stats.atoms << "\n";
    s << "dft_rows," << stats.dftRows << "\n";
    s << "dft_present," << stats.dftPresent << "\n";
    s << "phi_psi_present," << stats.phiPsiPresent << "\n";
    s << "phi_psi_eligible," << stats.phiPsiEligible << "\n";
    s << "phi_psi_finite_eligible," << stats.phiPsiFiniteEligible << "\n";
    s << "dssp_present," << stats.dsspPresent << "\n";
    s << "embedding_present," << stats.embeddingPresent << "\n";
    return true;
}

}  // namespace h5reader::rediscover
