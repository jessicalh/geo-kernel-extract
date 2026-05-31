#include "OutputManifest.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>

namespace h5reader::rediscover {

bool WriteOutputManifest(const QString& outDir, const std::vector<OutputEntry>& outputs,
                         const DftFrameAlignment& alignment, int rc, QString* err_out) {
    QJsonObject root;
    root.insert(QStringLiteral("substrate_version"), QStringLiteral("0.2.0-rediscover-surface"));
    root.insert(QStringLiteral("rc"), rc);

    QJsonObject align;
    align.insert(QStringLiteral("n_frames"), alignment.n_frames);
    align.insert(QStringLiteral("n_atoms_used"), alignment.n_atoms_used);
    align.insert(QStringLiteral("mean_angle_deg"), alignment.mean_angle_deg);
    align.insert(QStringLiteral("max_angle_deg"), alignment.max_angle_deg);
    align.insert(QStringLiteral("mean_rmsd_A"), alignment.mean_rmsd_A);
    align.insert(QStringLiteral("max_rmsd_A"), alignment.max_rmsd_A);
    align.insert(QStringLiteral("t2_components"),
                 alignment.max_angle_deg < 1.0
                     ? QStringLiteral("FRAME-ALIGNED (comparable as emitted)")
                     : QStringLiteral("ROTATED (T0/|T2| only until tensor rotation)"));
    root.insert(QStringLiteral("dft_frame_alignment"), align);

    QJsonObject conventions;
    QJsonObject sh;
    sh.insert(QStringLiteral("convention"), QStringLiteral("library_isometric_real_sh"));
    sh.insert(QStringLiteral("basis_ordering_ref"), QStringLiteral("src/Types.cpp::SphericalTensor::Decompose"));
    sh.insert(QStringLiteral("normalization"), QStringLiteral("isometric_real_sh"));
    conventions.insert(QStringLiteral("spherical_harmonics"), sh);
    conventions.insert(QStringLiteral("pbc_mode"), QStringLiteral("none_protein_whole_upstream"));
    conventions.insert(QStringLiteral("output_carrier"), QStringLiteral("per_relationship_schema"));
    conventions.insert(QStringLiteral("wide_array_payloads"),
                       QStringLiteral("sidecar_npy_entries_documented_in_python/nmr_extract/_catalog.py"));
    QJsonObject target;
    target.insert(QStringLiteral("T0"), QStringLiteral("dft_sigma_iso"));
    target.insert(QStringLiteral("T1_status"), QStringLiteral("unverified_emitted_not_discarded"));
    target.insert(QStringLiteral("T2_status"), QStringLiteral("emitted_with_frame_alignment_diagnostic"));
    conventions.insert(QStringLiteral("dft_target"), target);
    QJsonObject multipoles;
    multipoles.insert(QStringLiteral("origin"), QStringLiteral("target_atom"));
    multipoles.insert(QStringLiteral("charge_source_required"), true);
    multipoles.insert(QStringLiteral("exclude_residue_supported"), true);
    conventions.insert(QStringLiteral("charge_multipoles"), multipoles);
    root.insert(QStringLiteral("conventions"), conventions);

    QJsonArray arr;
    for (const OutputEntry& e : outputs) {
        QJsonObject o;
        o.insert(QStringLiteral("relationship"), e.relationship);
        o.insert(QStringLiteral("relationship_kind"), e.relationshipKind);
        QJsonObject files;
        files.insert(QStringLiteral("sources_csv"), e.sourcesCsv);
        files.insert(QStringLiteral("aggregated_csv"), e.aggregatedCsv);
        QJsonArray npys;
        for (const QString& name : e.sidecarNpys) npys.append(name);
        files.insert(QStringLiteral("array_payload_npys"), npys);
        o.insert(QStringLiteral("outputs"), files);
        QJsonObject counts;
        counts.insert(QStringLiteral("cases"), static_cast<qint64>(e.cases));
        counts.insert(QStringLiteral("source_rows"), static_cast<qint64>(e.sourceRows));
        counts.insert(QStringLiteral("aggregated_rows"), static_cast<qint64>(e.aggregatedRows));
        o.insert(QStringLiteral("counts"), counts);
        arr.append(o);
    }
    root.insert(QStringLiteral("relationships"), arr);

    QSaveFile f(QStringLiteral("%1/manifest.json").arg(outDir));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot open manifest.json in %1").arg(outDir);
        return false;
    }
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    if (!f.commit()) {
        if (err_out) *err_out = QStringLiteral("cannot commit manifest.json in %1").arg(outDir);
        return false;
    }
    return true;
}

}  // namespace h5reader::rediscover
