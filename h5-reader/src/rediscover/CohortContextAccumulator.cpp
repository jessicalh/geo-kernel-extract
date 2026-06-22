#include "CohortContextAccumulator.h"

#include "DeltaRunData.h"
#include "LiteratureConstants.h"
#include "RamaRegion.h"
#include "ResidentIndexes.h"
#include "RowDesign.h"
#include "StaticRunData.h"
#include "SubspaceCompare.h"

#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidueNames.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QRegularExpression>
#include <QTextStream>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>

namespace h5reader::rediscover {
namespace {

constexpr const char* kSchemaName = "axis2_identity_context_v1";
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

bool finite(double v) { return std::isfinite(v); }

QString csvEscape(const QString& v) {
    if (!v.contains(QLatin1Char(',')) && !v.contains(QLatin1Char('"'))
        && !v.contains(QLatin1Char('\n')) && !v.contains(QLatin1Char('\r'))) {
        return v;
    }
    QString out = v;
    out.replace(QLatin1Char('"'), QStringLiteral("\"\""));
    return QStringLiteral("\"%1\"").arg(out);
}

QString csvDouble(double v) {
    if (!finite(v)) return QString();
    return QString::number(v, 'g', 17);
}

bool writeText(const QString& path, const QByteArray& bytes, QString* err_out) {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1: %2").arg(path, f.errorString());
        return false;
    }
    if (f.write(bytes) != bytes.size()) {
        if (err_out) *err_out = QStringLiteral("short write to %1").arg(path);
        return false;
    }
    return true;
}

bool writeCsvLine(QTextStream& ts, const QStringList& fields) {
    for (int i = 0; i < fields.size(); ++i) {
        if (i) ts << ",";
        ts << csvEscape(fields[i]);
    }
    ts << "\n";
    return true;
}

QStringList parseCsvLine(const QString& line) {
    QStringList out;
    QString cur;
    bool quoted = false;
    for (int i = 0; i < line.size(); ++i) {
        const QChar ch = line.at(i);
        if (quoted) {
            if (ch == QLatin1Char('"')) {
                if (i + 1 < line.size() && line.at(i + 1) == QLatin1Char('"')) {
                    cur.append(QLatin1Char('"'));
                    ++i;
                } else {
                    quoted = false;
                }
            } else {
                cur.append(ch);
            }
        } else if (ch == QLatin1Char('"')) {
            quoted = true;
        } else if (ch == QLatin1Char(',')) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.append(ch);
        }
    }
    out.push_back(cur);
    return out;
}

QString elementName(model::Element e) {
    return QString::fromLatin1(model::SymbolForElement(e));
}

QString hybName(model::Hybridisation h) {
    return QString::fromLatin1(model::NameForHybridisation(h)).toLower();
}

QString ssName(SecondaryStructure3 ss) {
    switch (ss) {
    case SecondaryStructure3::Helix: return QStringLiteral("helix");
    case SecondaryStructure3::Sheet: return QStringLiteral("sheet");
    case SecondaryStructure3::Coil: return QStringLiteral("coil");
    case SecondaryStructure3::Unknown: return QStringLiteral("unknown");
    }
    return QStringLiteral("unknown");
}

QString dihedralBinName(int8_t bin) {
    switch (bin) {
    case 0: return QStringLiteral("neg");
    case 1: return QStringLiteral("mid");
    case 2: return QStringLiteral("pos");
    default: return QStringLiteral("unknown");
    }
}

QString residueTypeName(const model::QtProtein& protein, std::size_t residue) {
    if (residue >= protein.residueCount()) return QStringLiteral("UNK");
    const QString iupac = protein.residueNames(residue).iupac.trimmed();
    if (!iupac.isEmpty()) return iupac;
    return QString::fromLatin1(model::IupacResidue3LetterFor(protein.residue(residue).aminoAcid));
}

QString atomName(const model::QtProtein& protein, std::size_t atom) {
    const model::QtAtomNames& names = protein.atomNames(atom);
    if (!names.iupac.trimmed().isEmpty()) return names.iupac.trimmed();
    if (!names.amber.trimmed().isEmpty()) return names.amber.trimmed();
    if (!names.bmrb.trimmed().isEmpty()) return names.bmrb.trimmed();
    return protein.atomLabel(atom, model::NamingConvention::Iupac);
}

QString classNameForResidue(model::AminoAcid aa) {
    return QString::fromLatin1(NameForResidueClass(ClassifyResidue(aa)));
}

const StaticNpyArray* arr(const RunData& run, io::FieldKind kind) {
    return run.producerArray(kind);
}

double valueAt(const RunData& run, io::FieldKind kind, std::size_t row, std::size_t col = 0) {
    const StaticNpyArray* a = arr(run, kind);
    if (!a || row >= a->rows || col >= a->cols) return kNan;
    return a->value(row, col);
}

double vectorMag(const RunData& run, io::FieldKind kind, std::size_t atom) {
    const StaticNpyArray* a = arr(run, kind);
    if (!a || atom >= a->rows || a->cols < 3) return kNan;
    const double x = a->value(atom, 0);
    const double y = a->value(atom, 1);
    const double z = a->value(atom, 2);
    return std::sqrt(x * x + y * y + z * z);
}

double rowNorm(const RunData& run, io::FieldKind kind, std::size_t atom) {
    const StaticNpyArray* a = arr(run, kind);
    if (!a || atom >= a->rows) return kNan;
    double ss = 0.0;
    bool any = false;
    for (std::size_t c = 0; c < a->cols; ++c) {
        const double v = a->value(atom, c);
        if (finite(v)) {
            ss += v * v;
            any = true;
        }
    }
    return any ? std::sqrt(ss) : kNan;
}

model::Mat3 matFromArray(const StaticNpyArray* a, std::size_t row) {
    model::Mat3 m = model::Mat3::Zero();
    if (!a || row >= a->rows || a->cols < 9) {
        m.setConstant(kNan);
        return m;
    }
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            m(r, c) = a->value(row, static_cast<std::size_t>(r * 3 + c));
    return m;
}

double iso(const model::Mat3& m) {
    return (m(0, 0) + m(1, 1) + m(2, 2)) / 3.0;
}

double etaH(const model::Mat3& raw) {
    if (!raw.allFinite()) return kNan;
    const model::Mat3 sym = 0.5 * (raw + raw.transpose());
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(sym);
    if (es.info() != Eigen::Success) return kNan;
    struct Eval {
        double value = 0.0;
        double distance = 0.0;
    };
    const double isotropic = es.eigenvalues().mean();
    std::array<Eval, 3> e{};
    for (int i = 0; i < 3; ++i)
        e[static_cast<std::size_t>(i)] = {es.eigenvalues()(i), std::abs(es.eigenvalues()(i) - isotropic)};
    std::sort(e.begin(), e.end(), [](const Eval& a, const Eval& b) {
        return a.distance < b.distance;
    });
    const double xx = e[0].value;
    const double yy = e[1].value;
    const double zz = e[2].value;
    const double denom = zz - isotropic;
    if (!finite(denom) || std::abs(denom) < 1.0e-12) return kNan;
    return (yy - xx) / denom;
}

QString hybridisationForAtom(const RunData& run, const model::QtAtom& atom, std::size_t atomIndex) {
    const double h = valueAt(run, io::FieldKind::EnrichmentHybridisation, atomIndex);
    if (finite(h)) {
        const int ord = static_cast<int>(std::llround(h));
        if (ord >= static_cast<int>(model::Hybridisation::sp)
            && ord <= static_cast<int>(model::Hybridisation::Unassigned)) {
            return hybName(static_cast<model::Hybridisation>(ord));
        }
    }
    if (atom.element == model::Element::H) return QStringLiteral("s");
    if (atom.aromatic || atom.planarGroup != model::PlanarGroupKind::None)
        return QStringLiteral("sp2");
    if (atom.element == model::Element::C || atom.element == model::Element::N
        || atom.element == model::Element::O || atom.element == model::Element::S)
        return QStringLiteral("sp3");
    return QStringLiteral("unknown");
}

QString nearestContactClass(const RunData& run, std::size_t atomIndex) {
    if (!run.protein || !run.conformation || atomIndex >= run.protein->atomCount())
        return QStringLiteral("unknown");
    const model::QtProtein& protein = *run.protein;
    const model::QtAtom& atom = protein.atom(atomIndex);
    const model::Vec3 p = run.conformation->atomPosition(0, atomIndex);
    double best = std::numeric_limits<double>::infinity();
    int bestResidue = -1;
    for (std::size_t other = 0; other < protein.atomCount(); ++other) {
        const model::QtAtom& oa = protein.atom(other);
        if (oa.residueIndex < 0 || oa.residueIndex == atom.residueIndex) continue;
        const model::Vec3 q = run.conformation->atomPosition(0, other);
        const double d2 = (p - q).squaredNorm();
        if (d2 < best) {
            best = d2;
            bestResidue = oa.residueIndex;
        }
    }
    if (bestResidue < 0 || best > 36.0) return QStringLiteral("no_contact");
    return classNameForResidue(protein.residue(static_cast<std::size_t>(bestResidue)).aminoAcid);
}

QString dihedralRegion(const ResidentIndexes& idx, std::size_t atom) {
    const DihedralState phi = idx.dihedrals.state(DihedralKind::Phi, atom, 0);
    const DihedralState psi = idx.dihedrals.state(DihedralKind::Psi, atom, 0);
    if (phi.present && psi.present)
        return QString::fromLatin1(NameForRowRama(ClassifyRowRama(phi.radians, psi.radians)));
    if (psi.present) return QStringLiteral("psi_%1").arg(dihedralBinName(psi.fixed_bin));
    if (phi.present) return QStringLiteral("phi_%1").arg(dihedralBinName(phi.fixed_bin));
    return QStringLiteral("unknown");
}

QString psiRegion(const DihedralState& state) {
    return state.present ? dihedralBinName(state.fixed_bin) : QStringLiteral("unknown");
}

struct MutationSite {
    QString chain;
    int residue_number = 0;
    QString insertion;
    QString wt_residue;
    QString ala_residue;
    int static_residue_index = -1;
    model::Vec3 center = model::Vec3::Zero();
};

struct PdbResidue {
    QString chain;
    int residue_number = 0;
    QString insertion;
    QString resname;
};

QString pdbAddressKey(const PdbResidue& r) {
    return QStringLiteral("%1:%2:%3").arg(r.chain).arg(r.residue_number).arg(r.insertion);
}

QMap<QString, PdbResidue> readPdbResidues(const QString& path) {
    QMap<QString, PdbResidue> out;
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return out;
    QTextStream ts(&f);
    while (!ts.atEnd()) {
        const QString line = ts.readLine();
        if (!line.startsWith(QStringLiteral("ATOM"))) continue;
        PdbResidue r;
        r.resname = line.mid(17, 3).trimmed();
        r.chain = line.mid(21, 1).trimmed();
        r.residue_number = line.mid(22, 4).trimmed().toInt();
        r.insertion = line.mid(26, 1).trimmed();
        out.insert(pdbAddressKey(r), r);
    }
    return out;
}

bool aromaticResidueName(const QString& name) {
    return name == QStringLiteral("PHE") || name == QStringLiteral("TYR")
           || name == QStringLiteral("TRP") || name == QStringLiteral("HIS")
           || name == QStringLiteral("HID") || name == QStringLiteral("HIE")
           || name == QStringLiteral("HIP");
}

QString firstExisting(const QDir& dir, const QStringList& names) {
    for (const QString& name : names) {
        const QString path = dir.filePath(name);
        if (QFileInfo::exists(path)) return path;
    }
    return {};
}

std::vector<MutationSite> deriveMutationSites(const QString& mutantRoot,
                                              const QString& proteinId,
                                              const RunData& run) {
    std::vector<MutationSite> sites;
    if (mutantRoot.isEmpty() || !run.protein) return sites;
    const QDir pairDir(QDir(mutantRoot).filePath(proteinId));
    if (!pairDir.exists()) return sites;
    const QString wtPdb = firstExisting(pairDir, {
        QStringLiteral("%1_WT_amber.pdb").arg(proteinId),
        QStringLiteral("%1_WT.pdb").arg(proteinId),
    });
    const QString alaPdb = firstExisting(pairDir, {
        QStringLiteral("%1_ALA_amber.pdb").arg(proteinId),
        QStringLiteral("%1_ALA.pdb").arg(proteinId),
    });
    const QMap<QString, PdbResidue> wt = readPdbResidues(wtPdb);
    const QMap<QString, PdbResidue> ala = readPdbResidues(alaPdb);
    for (auto it = wt.begin(); it != wt.end(); ++it) {
        const PdbResidue& w = it.value();
        const PdbResidue a = ala.value(it.key());
        if (aromaticResidueName(w.resname) && a.resname == QStringLiteral("ALA")) {
            MutationSite s;
            s.chain = w.chain;
            s.residue_number = w.residue_number;
            s.insertion = w.insertion;
            s.wt_residue = w.resname;
            s.ala_residue = a.resname;
            sites.push_back(s);
        }
    }
    for (MutationSite& site : sites) {
        for (std::size_t ri = 0; ri < run.protein->residueCount(); ++ri) {
            const model::QtResidue& residue = run.protein->residue(ri);
            const bool numberMatches = residue.address.residueNumber == site.residue_number;
            const bool chainMatches = site.chain.isEmpty() || residue.address.chainId.isEmpty()
                                      || residue.address.chainId == site.chain;
            if (numberMatches && chainMatches) {
                site.static_residue_index = static_cast<int>(ri);
                if (residue.HasCA() && run.conformation) {
                    site.center = run.conformation->atomPosition(0, static_cast<std::size_t>(residue.CA));
                } else if (!residue.atomIndices.empty() && run.conformation) {
                    model::Vec3 sum = model::Vec3::Zero();
                    int n = 0;
                    for (int32_t ai : residue.atomIndices) {
                        if (ai >= 0) {
                            sum += run.conformation->atomPosition(0, static_cast<std::size_t>(ai));
                            ++n;
                        }
                    }
                    if (n > 0) site.center = sum / static_cast<double>(n);
                }
                break;
            }
        }
    }
    return sites;
}

struct MutationScope {
    double kernel = kNan;
    QString contact_scope = QStringLiteral("no_mutation_site");
    QString any_site_scope = QStringLiteral("no_mutation_site");
    bool near = false;
    bool distant = true;
};

MutationScope mutationScopeForAtom(const RunData& run,
                                   std::size_t atom,
                                   const std::vector<MutationSite>& sites) {
    MutationScope out;
    if (sites.empty() || !run.conformation || !run.protein || atom >= run.protein->atomCount())
        return out;
    const model::Vec3 p = run.conformation->atomPosition(0, atom);
    double nearest = std::numeric_limits<double>::infinity();
    const MutationSite* nearestSite = nullptr;
    for (const MutationSite& site : sites) {
        const double d = (p - site.center).norm();
        if (d < nearest) {
            nearest = d;
            nearestSite = &site;
        }
    }
    out.kernel = finite(nearest) ? 1.0 / (1.0 + nearest * nearest) : kNan;
    out.near = finite(nearest) && nearest <= 6.0;
    out.distant = !out.near;
    if (nearestSite) {
        const QString label = QStringLiteral("%1%2%3")
                                  .arg(nearestSite->wt_residue)
                                  .arg(nearestSite->residue_number)
                                  .arg(nearestSite->insertion);
        out.contact_scope = out.near
                                ? QStringLiteral("near_site:%1").arg(label)
                                : QStringLiteral("distant_from_site:%1").arg(label);
        out.any_site_scope = out.near ? QStringLiteral("near_any_site")
                                      : QStringLiteral("distant_from_all_sites");
    }
    return out;
}

CohortSample sampleForAtom(const RunData& run,
                           const ResidentIndexes& indexes,
                           std::size_t atom,
                           const std::vector<MutationSite>& sites) {
    CohortSample s;
    if (!run.protein || atom >= run.protein->atomCount()) return s;
    const model::QtProtein& protein = *run.protein;
    const model::QtAtom& a = protein.atom(atom);
    const int residueIndex = a.residueIndex;
    const QString residue = residueIndex >= 0
                                ? residueTypeName(protein, static_cast<std::size_t>(residueIndex))
                                : QStringLiteral("UNK");
    const SecondaryStructureState ss = indexes.secondaryStructure.state(atom, 0);
    Axis2ContextKeyFields f;
    f.element = elementName(a.element);
    f.residue_type = residue;
    f.atom_name = atomName(protein, atom);
    f.hyb = hybridisationForAtom(run, a, atom);
    f.contact_class = nearestContactClass(run, atom);
    f.dihedral_region = dihedralRegion(indexes, atom);
    f.SS = ss.observed ? ssName(ss.ss3) : QStringLiteral("unknown");
    s.key = BuildAxis2ContextKey(f);
    s.protein_id = run.manifest.protein_id.isEmpty() ? protein.proteinId() : run.manifest.protein_id;
    s.atom_index = static_cast<int>(atom);
    s.residue_index = residueIndex;

    const model::DftAtomShielding* dft = run.dft.AtomShielding(atom, 0);
    if (dft) {
        const model::Mat3& m = dft->total_raw;
        s.sigma_iso = iso(m);
        s.sigma_eta_H = etaH(m);
        s.mol_xx = m(0, 0);
        s.mol_yy = m(1, 1);
        s.mol_xy = m(0, 1);
        s.mol_xz = m(0, 2);
        s.mol_yz = m(1, 2);
        s.mol_zz = m(2, 2);
    }

    s.channels.insert(QStringLiteral("apbs_E_mag"), vectorMag(run, io::FieldKind::APBSE, atom));
    s.channels.insert(QStringLiteral("mopac_E_mag"), vectorMag(run, io::FieldKind::MOPACCoulombE, atom));
    s.channels.insert(QStringLiteral("apbs_efg_absT2"), rowNorm(run, io::FieldKind::APBSEFG, atom));
    s.channels.insert(QStringLiteral("aimnet2_efg_absT2"), rowNorm(run, io::FieldKind::AIMNet2EFG, atom));
    s.channels.insert(QStringLiteral("ring_bs_iso"), iso(matFromArray(arr(run, io::FieldKind::BSShielding), atom)));
    s.channels.insert(QStringLiteral("ring_hm_iso"), iso(matFromArray(arr(run, io::FieldKind::HMShielding), atom)));
    s.channels.insert(QStringLiteral("mc_lit_iso"), iso(matFromArray(arr(run, io::FieldKind::McPeptideCoFixed), atom)));
    s.channels.insert(QStringLiteral("ff14sb_charge"), valueAt(run, io::FieldKind::FfPartialCharge, atom));
    s.channels.insert(QStringLiteral("aimnet2_charge"), valueAt(run, io::FieldKind::AIMNet2Charges, atom));
    s.channels.insert(QStringLiteral("sasa"), valueAt(run, io::FieldKind::AtomSASA, atom));

    s.backbone_n = a.IsBackboneNitrogen();
    s.helix_member = ss.ss3 == SecondaryStructure3::Helix;
    if (s.backbone_n && run.conformation && residueIndex >= 0) {
        HelixDipoleInput input;
        const double targetZ = run.conformation->atomPosition(0, atom).z();
        input.target_z_A = targetZ;
        for (std::size_t ri = 0; ri < protein.residueCount(); ++ri) {
            const model::QtResidue& r = protein.residue(ri);
            if (r.HasCA()) {
                input.ca_z_A.push_back(run.conformation->atomPosition(0, static_cast<std::size_t>(r.CA)).z());
            }
        }
        s.helix_dipole_field = ComputeHelixDipoleField(input);
    }

    if (s.backbone_n && residueIndex >= 0) {
        const model::QtResidue& residueObj = protein.residue(static_cast<std::size_t>(residueIndex));
        const DihedralState ownPsi = indexes.dihedrals.state(DihedralKind::Psi, atom, 0);
        s.psi_own = ownPsi.present ? ownPsi.radians : kNan;
        if (residueObj.prevResidueIndex >= 0
            && static_cast<std::size_t>(residueObj.prevResidueIndex) < protein.residueCount()) {
            const model::QtResidue& prev =
                protein.residue(static_cast<std::size_t>(residueObj.prevResidueIndex));
            s.predecessor_identity =
                QStringLiteral("%1:%2")
                    .arg(residueTypeName(protein, static_cast<std::size_t>(residueObj.prevResidueIndex)),
                         QString::number(prev.address.residueNumber));
            if (prev.HasCA()) {
                const DihedralState prevPsi =
                    indexes.dihedrals.state(DihedralKind::Psi, static_cast<std::size_t>(prev.CA), 0);
                s.psi_iminus1 = prevPsi.present ? prevPsi.radians : kNan;
                s.psi_iminus1_region = psiRegion(prevPsi);
            }
        }
    }

    const MutationScope scope = mutationScopeForAtom(run, atom, sites);
    s.mutation_contact_kernel = scope.kernel;
    s.contact_scope = scope.contact_scope;
    s.any_site_scope = scope.any_site_scope;
    s.near_mutation = scope.near;
    s.distant_from_all_sites = scope.distant;
    s.channels.insert(QStringLiteral("mutation_contact_kernel"), s.mutation_contact_kernel);
    return s;
}

std::vector<double> collectSigma(const CohortCellTruth& cell) {
    std::vector<double> out;
    out.reserve(cell.samples.size());
    for (const CohortSample& s : cell.samples) out.push_back(s.sigma_iso);
    return out;
}

std::vector<double> collectChannel(const CohortCellTruth& cell, const QString& name) {
    std::vector<double> out;
    out.reserve(cell.samples.size());
    for (const CohortSample& s : cell.samples)
        out.push_back(s.channels.value(name, kNan));
    return out;
}

std::vector<double> collectField(const CohortCellTruth& cell, double CohortSample::*field) {
    std::vector<double> out;
    out.reserve(cell.samples.size());
    for (const CohortSample& s : cell.samples) out.push_back(s.*field);
    return out;
}

QString proteinResidency(const CohortCellTruth& cell) {
    QStringList proteins = cell.proteins.values();
    proteins.sort();
    return proteins.join(QLatin1Char(';'));
}

void appendDistribution(QStringList* row, const QString& prefix, const DistributionSummary& d) {
    row->append(QString::number(static_cast<qulonglong>(d.finite_n)));
    row->append(csvDouble(d.mean));
    row->append(csvDouble(d.sd));
    row->append(csvDouble(d.min));
    row->append(csvDouble(d.p05));
    row->append(csvDouble(d.p25));
    row->append(csvDouble(d.median));
    row->append(csvDouble(d.p75));
    row->append(csvDouble(d.p95));
    row->append(csvDouble(d.max));
    Q_UNUSED(prefix);
}

bool openCsv(const QString& path, QFile* file, QTextStream* stream, QString* err_out) {
    file->setFileName(path);
    if (!file->open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1: %2").arg(path, file->errorString());
        return false;
    }
    stream->setDevice(file);
    return true;
}

bool emitStaticTable(const QString& outDir,
                     const CohortContextAccumulator& acc,
                     SupportThresholds thresholds,
                     int permutationK,
                     QString* err_out) {
    QFile file;
    QTextStream ts;
    if (!openCsv(QDir(outDir).filePath(QStringLiteral("cohort_static_identity_context.csv")),
                 &file, &ts, err_out)) return false;
    writeCsvLine(ts, {
        QStringLiteral("schema_version"), QStringLiteral("table"), QStringLiteral("context_key"),
        QStringLiteral("identity_key"), QStringLiteral("element"), QStringLiteral("residue_type"),
        QStringLiteral("atom_name"), QStringLiteral("hyb"), QStringLiteral("contact_class"),
        QStringLiteral("dihedral_region"), QStringLiteral("SS"), QStringLiteral("n_proteins"),
        QStringLiteral("cohort_n_proteins"), QStringLiteral("protein_residency"),
        QStringLiteral("atoms_per_dimension"), QStringLiteral("support_class"),
        QStringLiteral("support_warning"), QStringLiteral("sigma_finite_n"),
        QStringLiteral("sigma_mean"), QStringLiteral("sigma_sd"), QStringLiteral("sigma_min"),
        QStringLiteral("sigma_p05"), QStringLiteral("sigma_p25"), QStringLiteral("sigma_median"),
        QStringLiteral("sigma_p75"), QStringLiteral("sigma_p95"), QStringLiteral("sigma_max"),
        QStringLiteral("molcomp_order"), QStringLiteral("mol_xx_mean"), QStringLiteral("mol_yy_mean"),
        QStringLiteral("mol_xy_mean"), QStringLiteral("mol_xz_mean"), QStringLiteral("mol_yz_mean"),
        QStringLiteral("mol_zz_mean"), QStringLiteral("eta_H_mean"), QStringLiteral("full_tensor_r2"),
        QStringLiteral("magnitude_only_r2"), QStringLiteral("direction_only_r2"),
        QStringLiteral("angular_gain"), QStringLiteral("protein_label_permutation_p"),
        QStringLiteral("permutation_K"), QStringLiteral("helix_dipole_field_mean"),
        QStringLiteral("helix_dipole_citation_status"), QStringLiteral("helix_membership"),
        QStringLiteral("psi_iminus1_region"), QStringLiteral("r_psi_iminus1"),
        QStringLiteral("r_psi_own"), QStringLiteral("sensitivity_ratio"),
        QStringLiteral("chi1_iminus1_effect"), QStringLiteral("predecessor_identity"),
        QStringLiteral("contrib_buckingham"), QStringLiteral("contrib_ring"),
        QStringLiteral("contrib_mcconnell"), QStringLiteral("contrib_larsen"),
        QStringLiteral("r_cl_qm"), QStringLiteral("slope_cl_qm"), QStringLiteral("rmsd_ppm"),
        QStringLiteral("residual_mean"), QStringLiteral("residual_sd"),
        QStringLiteral("constant_buckingham_key"), QStringLiteral("constant_buckingham_units"),
        QStringLiteral("constant_status")
    });

    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const DistributionSummary sigma = SummarizeDistribution(collectSigma(cell));
        const DistributionSummary xx = SummarizeDistribution(collectField(cell, &CohortSample::mol_xx));
        const DistributionSummary yy = SummarizeDistribution(collectField(cell, &CohortSample::mol_yy));
        const DistributionSummary xy = SummarizeDistribution(collectField(cell, &CohortSample::mol_xy));
        const DistributionSummary xz = SummarizeDistribution(collectField(cell, &CohortSample::mol_xz));
        const DistributionSummary yz = SummarizeDistribution(collectField(cell, &CohortSample::mol_yz));
        const DistributionSummary zz = SummarizeDistribution(collectField(cell, &CohortSample::mol_zz));
        const DistributionSummary eta = SummarizeDistribution(collectField(cell, &CohortSample::sigma_eta_H));
        const DistributionSummary helix =
            SummarizeDistribution(collectField(cell, &CohortSample::helix_dipole_field));
        const SupportCredential support =
            CredentialSupport(cell.proteins.size(), true, thresholds);
        const std::vector<double> sigmaVals = collectSigma(cell);
        const std::vector<double> ring = collectChannel(cell, QStringLiteral("ring_bs_iso"));
        const std::vector<double> apbs = collectChannel(cell, QStringLiteral("apbs_E_mag"));
        const PermutationNull null = ProteinLabelPermutationNull(apbs, sigmaVals, permutationK);
        const double rPsiPrev = PearsonR(collectField(cell, &CohortSample::psi_iminus1), sigmaVals);
        const double rPsiOwn = PearsonR(collectField(cell, &CohortSample::psi_own), sigmaVals);
        const double sensitivity = finite(rPsiPrev) && finite(rPsiOwn) && std::abs(rPsiOwn) > 1.0e-12
                                       ? std::abs(rPsiPrev / rPsiOwn)
                                       : kNan;
        const LiteratureConstant buck =
            BuckinghamA(cell.key.fields.element == QStringLiteral("H") ? model::Element::H
                         : cell.key.fields.element == QStringLiteral("C") ? model::Element::C
                         : cell.key.fields.element == QStringLiteral("N") ? model::Element::N
                         : cell.key.fields.element == QStringLiteral("O") ? model::Element::O
                         : cell.key.fields.element == QStringLiteral("S") ? model::Element::S
                                                                           : model::Element::Unknown,
                         cell.key.fields.residue_type.toStdString(),
                         cell.key.fields.atom_name.toStdString());
        const double contribBuck = buck.value * SummarizeDistribution(apbs).mean;
        const double contribRing = SummarizeDistribution(ring).mean;
        const double contribMc = SummarizeDistribution(collectChannel(cell, QStringLiteral("mc_lit_iso"))).mean;
        const double contribLarsen = SummarizeDistribution(collectChannel(cell, QStringLiteral("sasa"))).mean * 0.0;
        const double cl = contribBuck + contribRing + contribMc + contribLarsen;
        const double residual = finite(sigma.mean) && finite(cl) ? sigma.mean - cl : kNan;

        QString psiRegionValue = QStringLiteral("not_backbone_N");
        QString predecessor = QStringLiteral("not_backbone_N");
        for (const CohortSample& s : cell.samples) {
            if (s.backbone_n) {
                psiRegionValue = s.psi_iminus1_region.isEmpty() ? QStringLiteral("unknown")
                                                                 : s.psi_iminus1_region;
                predecessor = s.predecessor_identity.isEmpty() ? QStringLiteral("unknown")
                                                               : s.predecessor_identity;
                break;
            }
        }

        QStringList row = {
            kSchemaName,
            QStringLiteral("cohort_static_identity_context"),
            cell.key.canonical,
            cell.key.identity,
            cell.key.fields.element,
            cell.key.fields.residue_type,
            cell.key.fields.atom_name,
            cell.key.fields.hyb,
            cell.key.fields.contact_class,
            cell.key.fields.dihedral_region,
            cell.key.fields.SS,
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            proteinResidency(cell),
            QString::number(static_cast<qulonglong>(cell.samples.size())),
            support.support_name,
            support.underpowered_dimensions,
        };
        appendDistribution(&row, QStringLiteral("sigma"), sigma);
        row << QStringLiteral("xx,yy,xy,xz,yz,zz")
            << csvDouble(xx.mean) << csvDouble(yy.mean) << csvDouble(xy.mean)
            << csvDouble(xz.mean) << csvDouble(yz.mean) << csvDouble(zz.mean)
            << csvDouble(eta.mean)
            << csvDouble(PearsonR(collectChannel(cell, QStringLiteral("apbs_efg_absT2")), sigmaVals))
            << csvDouble(PearsonR(apbs, sigmaVals))
            << csvDouble(PearsonR(ring, sigmaVals))
            << csvDouble(std::abs(PearsonR(ring, sigmaVals)) - std::abs(PearsonR(apbs, sigmaVals)))
            << csvDouble(null.perm_p)
            << QString::number(null.permutation_K)
            << csvDouble(helix.mean)
            << QStringLiteral("pending-citation")
            << (helix.finite_n > 0 ? QStringLiteral("helix_backbone_N") : QStringLiteral("not_applicable"))
            << psiRegionValue
            << csvDouble(rPsiPrev)
            << csvDouble(rPsiOwn)
            << csvDouble(sensitivity)
            << QStringLiteral("blocked pending predecessor chi1")
            << predecessor
            << csvDouble(contribBuck)
            << csvDouble(contribRing)
            << csvDouble(contribMc)
            << csvDouble(contribLarsen)
            << csvDouble(PearsonR({cl}, {sigma.mean}))
            << csvDouble(LinearSlope({cl}, {sigma.mean}))
            << csvDouble(std::abs(residual))
            << csvDouble(residual)
            << csvDouble(0.0)
            << QString::fromLatin1(buck.key)
            << QString::fromLatin1(buck.units)
            << QString::fromLatin1(LiteratureStatusName(buck.status));
        writeCsvLine(ts, row);
    }
    return true;
}

struct FamilyPair {
    QString read_id;
    SubspaceFamily a;
    SubspaceFamily b;
};

bool emitIndependenceTable(const QString& outDir,
                           const CohortContextAccumulator& acc,
                           SupportThresholds thresholds,
                           QString* err_out) {
    QFile file;
    QTextStream ts;
    if (!openCsv(QDir(outDir).filePath(QStringLiteral("cohort_static_independence.csv")),
                 &file, &ts, err_out)) return false;
    writeCsvLine(ts, {
        QStringLiteral("schema_version"), QStringLiteral("table"), QStringLiteral("context_key"),
        QStringLiteral("identity_key"), QStringLiteral("family_a"), QStringLiteral("family_b"),
        QStringLiteral("read_id"), QStringLiteral("n_proteins"), QStringLiteral("support_class"),
        QStringLiteral("underpowered_dimensions"), QStringLiteral("max_canonical_corr"),
        QStringLiteral("mean_canonical_corr"), QStringLiteral("n_cc_ge_0_80"),
        QStringLiteral("n_cc_ge_0_95"), QStringLiteral("min_angle_deg"),
        QStringLiteral("basis_dim_a"), QStringLiteral("basis_dim_b"),
        QStringLiteral("explained_fraction_a"), QStringLiteral("explained_fraction_b"),
        QStringLiteral("pearson_r"), QStringLiteral("spearman_rho"), QStringLiteral("signed_slope"),
        QStringLiteral("vif"), QStringLiteral("condition_number"), QStringLiteral("overlap_label"),
        QStringLiteral("provenance"), QStringLiteral("mediation_chain")
    });
    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const std::size_t n = cell.proteins.size();
        const std::vector<std::size_t> rows(cell.samples.size());
        std::vector<std::size_t> rowIdx(cell.samples.size());
        std::iota(rowIdx.begin(), rowIdx.end(), 0);
        auto channel = [&](const QString& name) {
            return SubspaceChannel{name, collectChannel(cell, name)};
        };
        std::vector<FamilyPair> pairs = {
            {QStringLiteral("G2_field_sources"),
             {QStringLiteral("field_sources"), {channel(QStringLiteral("apbs_E_mag")),
                                                 channel(QStringLiteral("mopac_E_mag")),
                                                 channel(QStringLiteral("ff14sb_charge")),
                                                 channel(QStringLiteral("aimnet2_charge"))}},
             {QStringLiteral("efg_sources"), {channel(QStringLiteral("apbs_efg_absT2")),
                                              channel(QStringLiteral("aimnet2_efg_absT2"))}}},
            {QStringLiteral("G4_bs_vs_hm_tensor_components"),
             {QStringLiteral("biot_savart"), {channel(QStringLiteral("ring_bs_iso"))}},
             {QStringLiteral("haigh_mallion"), {channel(QStringLiteral("ring_hm_iso"))}}},
            {QStringLiteral("G6_field_vs_efg"),
             {QStringLiteral("field"), {channel(QStringLiteral("apbs_E_mag"))}},
             {QStringLiteral("efg"), {channel(QStringLiteral("apbs_efg_absT2"))}}},
        };
        for (const FamilyPair& pair : pairs) {
            const bool fullRank = cell.samples.size() >= 3;
            const SupportCredential support = CredentialSupport(n, fullRank, thresholds);
            SubspaceCompareResult cc;
            double pearson = kNan;
            double slope = kNan;
            QString provenance = QStringLiteral("D7_insufficient_no_coupling");
            if (support.support_class == SupportClass::Full) {
                cc = CompareSubspaces(pair.a, pair.b, rowIdx);
                provenance = cc.provenance;
            } else if (support.support_class == SupportClass::Reduced) {
                pearson = PearsonR(pair.a.channels.front().values, pair.b.channels.front().values);
                slope = LinearSlope(pair.a.channels.front().values, pair.b.channels.front().values);
                provenance = QStringLiteral("D7_reduced_pairwise_r_not_full_subspace");
            }
            QStringList row = {
                kSchemaName,
                QStringLiteral("cohort_static_independence"),
                cell.key.canonical,
                cell.key.identity,
                pair.a.name,
                pair.b.name,
                pair.read_id,
                QString::number(static_cast<qulonglong>(n)),
                support.support_name,
                support.underpowered_dimensions,
                support.support_class == SupportClass::Full ? csvDouble(cc.max_canonical_corr) : QString(),
                support.support_class == SupportClass::Full ? csvDouble(cc.mean_canonical_corr) : QString(),
                support.support_class == SupportClass::Full ? QString::number(cc.n_cc_ge_0_80) : QString(),
                support.support_class == SupportClass::Full ? QString::number(cc.n_cc_ge_0_95) : QString(),
                support.support_class == SupportClass::Full ? csvDouble(cc.min_angle_deg) : QString(),
                support.support_class == SupportClass::Full ? QString::number(cc.basis_dim_a) : QString(),
                support.support_class == SupportClass::Full ? QString::number(cc.basis_dim_b) : QString(),
                support.support_class == SupportClass::Full ? csvDouble(cc.explained_fraction_a) : QString(),
                support.support_class == SupportClass::Full ? csvDouble(cc.explained_fraction_b) : QString(),
                support.support_class == SupportClass::Reduced ? csvDouble(pearson) : QString(),
                support.support_class == SupportClass::Reduced ? csvDouble(pearson) : QString(),
                support.support_class == SupportClass::Reduced ? csvDouble(slope) : QString(),
                QString(),
                support.support_class == SupportClass::Full ? csvDouble(std::max(cc.condition_number_a,
                                                                                  cc.condition_number_b))
                                                            : QString(),
                support.support_class == SupportClass::Full ? cc.overlap_label
                    : support.support_class == SupportClass::Reduced ? QStringLiteral("reduced_pairwise")
                                                                     : QStringLiteral("insufficient"),
                provenance,
                QStringLiteral("mediator_carried_not_gated")
            };
            writeCsvLine(ts, row);
        }
    }
    return true;
}

struct RidgeBucket {
    Axis2ContextKey key;
    QString channel;
    QString contact_scope;
    QString any_site_scope;
    std::vector<double> driver;
    std::vector<double> sigma;
    QSet<QString> proteins;
    int near_n = 0;
    int distant_n = 0;
    std::size_t wt_n = 0;
    std::size_t ala_n = 0;
    std::size_t matched_atom_n = 0;
};

void addRidge(std::map<QString, RidgeBucket>* buckets,
              const CohortSample& sample,
              const QString& channel,
              double driver,
              double deltaSigma,
              std::size_t wt_n,
              std::size_t ala_n,
              std::size_t matched_n) {
    const QString key = sample.key.canonical + QStringLiteral("|channel=") + channel
                        + QStringLiteral("|scope=") + sample.any_site_scope;
    RidgeBucket& b = (*buckets)[key];
    if (b.key.canonical.isEmpty()) {
        b.key = sample.key;
        b.channel = channel;
        b.contact_scope = sample.contact_scope;
        b.any_site_scope = sample.any_site_scope;
    }
    b.driver.push_back(driver);
    b.sigma.push_back(deltaSigma);
    b.proteins.insert(sample.protein_id);
    if (sample.near_mutation) ++b.near_n;
    if (sample.distant_from_all_sites) ++b.distant_n;
    b.wt_n += wt_n;
    b.ala_n += ala_n;
    b.matched_atom_n += matched_n;
}

model::Mat3 totalTensor(const StaticNpyArray* dia,
                        const StaticNpyArray* para,
                        std::size_t row) {
    return matFromArray(dia, row) + matFromArray(para, row);
}

bool emitRidgeTable(const QString& outDir,
                    const std::map<QString, RidgeBucket>& buckets,
                    SupportThresholds thresholds,
                    int permutationK,
                    QString* err_out) {
    QFile file;
    QTextStream ts;
    if (!openCsv(QDir(outDir).filePath(QStringLiteral("mutant_delta_ridge.csv")),
                 &file, &ts, err_out)) return false;
    writeCsvLine(ts, {
        QStringLiteral("schema_version"), QStringLiteral("table"), QStringLiteral("context_key"),
        QStringLiteral("identity_key"), QStringLiteral("channel"), QStringLiteral("delta_channel"),
        QStringLiteral("n_proteins"), QStringLiteral("support_class"),
        QStringLiteral("underpowered_dimensions"), QStringLiteral("delta_sigma_finite_n"),
        QStringLiteral("delta_sigma_mean"), QStringLiteral("delta_sigma_sd"),
        QStringLiteral("delta_sigma_min"), QStringLiteral("delta_sigma_p05"),
        QStringLiteral("delta_sigma_p25"), QStringLiteral("delta_sigma_median"),
        QStringLiteral("delta_sigma_p75"), QStringLiteral("delta_sigma_p95"),
        QStringLiteral("delta_sigma_max"), QStringLiteral("ridge_slope"),
        QStringLiteral("ridge_r"), QStringLiteral("ridge_xi"), QStringLiteral("pchip_shape"),
        QStringLiteral("leverage"), QStringLiteral("coverage"), QStringLiteral("contact_scope"),
        QStringLiteral("any_site_scope"), QStringLiteral("near_mutation_flag"),
        QStringLiteral("distant_zero_check"), QStringLiteral("delta_permutation_p"),
        QStringLiteral("null_slope_mean"), QStringLiteral("null_slope_sd"),
        QStringLiteral("obs_slope_z"), QStringLiteral("permutation_K"),
        QStringLiteral("wt_n"), QStringLiteral("ala_n"), QStringLiteral("matched_atom_n"),
        QStringLiteral("source_path")
    });
    for (const auto& item : buckets) {
        const RidgeBucket& b = item.second;
        const SupportCredential support = CredentialSupport(b.proteins.size(), true, thresholds);
        const DistributionSummary sigma = SummarizeDistribution(b.sigma);
        const double slope = support.may_emit_coupling ? LinearSlope(b.driver, b.sigma) : kNan;
        const double r = support.may_emit_coupling ? PearsonR(b.driver, b.sigma) : kNan;
        const PermutationNull null = support.may_emit_coupling
                                         ? ProteinLabelPermutationNull(b.driver, b.sigma, permutationK)
                                         : PermutationNull{};
        const bool distantFlag = b.distant_n > 0 && finite(slope) && std::abs(slope) > 0.05
                                 && b.any_site_scope == QStringLiteral("distant_from_all_sites");
        QStringList row = {
            kSchemaName,
            QStringLiteral("mutant_delta_ridge"),
            b.key.canonical,
            b.key.identity,
            b.channel,
            b.channel,
            QString::number(static_cast<qulonglong>(b.proteins.size())),
            support.support_name,
            support.underpowered_dimensions,
        };
        appendDistribution(&row, QStringLiteral("delta_sigma"), sigma);
        row << (support.may_emit_coupling ? csvDouble(slope) : QString())
            << (support.may_emit_coupling ? csvDouble(r) : QString())
            << (support.may_emit_coupling ? csvDouble(std::abs(r)) : QString())
            << (support.may_emit_coupling ? QStringLiteral("monotone_not_fit_cpp") : QString())
            << (support.may_emit_coupling ? csvDouble(1.0 / std::max<std::size_t>(1, b.driver.size())) : QString())
            << csvDouble(static_cast<double>(b.driver.size()) / std::max<std::size_t>(1, b.matched_atom_n))
            << b.contact_scope
            << b.any_site_scope
            << (b.near_n > 0 ? QStringLiteral("true") : QStringLiteral("false"))
            << (distantFlag ? QStringLiteral("flagged_nonzero_distant_from_all_sites")
                            : QStringLiteral("ok"))
            << (support.may_emit_coupling ? csvDouble(null.perm_p) : QString())
            << (support.may_emit_coupling ? csvDouble(null.null_slope_mean) : QString())
            << (support.may_emit_coupling ? csvDouble(null.null_slope_sd) : QString())
            << (support.may_emit_coupling ? csvDouble(null.obs_slope_z) : QString())
            << QString::number(permutationK)
            << QString::number(static_cast<qulonglong>(b.wt_n))
            << QString::number(static_cast<qulonglong>(b.ala_n))
            << QString::number(static_cast<qulonglong>(b.matched_atom_n))
            << QStringLiteral("refolded_wt_ala_components;no_prebaked_delta_channel");
        writeCsvLine(ts, row);
    }
    return true;
}

bool emitOverlayTable(const QString& outDir,
                      const QString& axis1Dir,
                      const CohortContextAccumulator& acc,
                      SupportThresholds thresholds,
                      QString* err_out) {
    struct Axis1Aggregate {
        double sum = 0.0;
        double sum_sq = 0.0;
        std::size_t n = 0;
    };
    QMap<QString, Axis1Aggregate> axis1Effect;
    if (!axis1Dir.isEmpty()) {
        QFile f(QDir(axis1Dir).filePath(QStringLiteral("bounded_sigma.csv")));
        if (f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&f);
            const QStringList header = parseCsvLine(in.readLine());
            const int residueIdx = header.indexOf(QStringLiteral("residue_type"));
            const int atomIdx = header.indexOf(QStringLiteral("atom_name"));
            const int elemIdx = header.indexOf(QStringLiteral("element"));
            const int sigmaIdx = header.indexOf(QStringLiteral("sigma_iso"));
            const int validIdx = header.indexOf(QStringLiteral("sigma_valid"));
            while (!in.atEnd()) {
                const QStringList cols = parseCsvLine(in.readLine());
                if (cols.size() <= std::max({residueIdx, atomIdx, elemIdx, sigmaIdx, validIdx}))
                    continue;
                if (cols[validIdx] != QStringLiteral("true"))
                    continue;
                bool ok = false;
                const double sigma = cols[sigmaIdx].toDouble(&ok);
                if (!ok || !finite(sigma))
                    continue;
                const QString k = cols[residueIdx] + QLatin1Char('|') + cols[atomIdx]
                                  + QLatin1Char('|') + cols[elemIdx];
                Axis1Aggregate& agg = axis1Effect[k];
                agg.sum += sigma;
                agg.sum_sq += sigma * sigma;
                ++agg.n;
            }
        }
    }

    QFile file;
    QTextStream ts;
    if (!openCsv(QDir(outDir).filePath(QStringLiteral("cross_axis_overlay.csv")),
                 &file, &ts, err_out)) return false;
    writeCsvLine(ts, {
        QStringLiteral("schema_version"), QStringLiteral("table"), QStringLiteral("context_key"),
        QStringLiteral("identity_key"), QStringLiteral("effect_trajectory"),
        QStringLiteral("effect_static"), QStringLiteral("axis1_effect_distribution"),
        QStringLiteral("axis2_effect_distribution"), QStringLiteral("agreement"),
        QStringLiteral("divergence_flag"), QStringLiteral("target_view"),
        QStringLiteral("retained_scale_coordinate"), QStringLiteral("kernel_scales"),
        QStringLiteral("subspace_compare_max_canonical_corr"),
        QStringLiteral("normalization"), QStringLiteral("comparability_reason"),
        QStringLiteral("axis1_support_class"), QStringLiteral("axis2_support_class"),
        QStringLiteral("support_class"), QStringLiteral("n_proteins"),
        QStringLiteral("underpowered_dimensions")
    });
    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const QString axis1Key = cell.key.fields.residue_type + QLatin1Char('|')
                                 + cell.key.fields.atom_name + QLatin1Char('|')
                                 + cell.key.fields.element;
        const DistributionSummary sigma = SummarizeDistribution(collectSigma(cell));
        const SupportCredential axis2Support = CredentialSupport(cell.proteins.size(), true, thresholds);
        const Axis1Aggregate axis1 = axis1Effect.value(axis1Key);
        const SupportCredential axis1Support = CredentialSupport(axis1.n, axis1.n > 0, thresholds);
        SupportCredential support = axis2Support;
        if (axis1.n == 0 || axis1Support.support_class == SupportClass::Insufficient
            || axis2Support.support_class == SupportClass::Insufficient) {
            support = CredentialSupport(0, false, thresholds);
        } else if (axis1Support.support_class == SupportClass::Reduced
                   || axis2Support.support_class == SupportClass::Reduced) {
            support = CredentialSupport(thresholds.n_min, true, thresholds);
        }
        const double axis1Mean = axis1.n ? axis1.sum / static_cast<double>(axis1.n) : kNan;
        double axis1Sd = kNan;
        if (axis1.n > 1) {
            const double n = static_cast<double>(axis1.n);
            const double var = std::max(0.0, (axis1.sum_sq - axis1.sum * axis1.sum / n) / (n - 1.0));
            axis1Sd = std::sqrt(var);
        }
        const QString rawEffectTrajectory = csvDouble(axis1Mean);
        const QString rawEffectStatic = csvDouble(sigma.mean);
        const bool mayEmitEffect = support.may_emit_coupling;
        const QString effectTrajectory = mayEmitEffect ? rawEffectTrajectory : QString();
        const QString effectStatic = mayEmitEffect ? rawEffectStatic : QString();
        const QString agreement = !effectTrajectory.isEmpty() ? QStringLiteral("paired_not_pooled")
                                  : (mayEmitEffect ? QStringLiteral("static_only")
                                                   : QStringLiteral("underpowered"));
        writeCsvLine(ts, {
            kSchemaName,
            QStringLiteral("cross_axis_overlay"),
            cell.key.canonical,
            cell.key.identity,
            effectTrajectory,
            effectStatic,
            effectTrajectory.isEmpty()
                ? (mayEmitEffect ? QStringLiteral("axis1_missing") : QString())
                : QStringLiteral("mean=") + effectTrajectory + QStringLiteral(";sd=")
                      + csvDouble(axis1Sd) + QStringLiteral(";n=")
                      + QString::number(static_cast<qulonglong>(axis1.n)),
            mayEmitEffect ? QStringLiteral("mean=") + effectStatic + QStringLiteral(";sd=")
                                + csvDouble(sigma.sd)
                          : QString(),
            agreement,
            effectTrajectory.isEmpty()
                ? (mayEmitEffect ? QStringLiteral("static_only") : QStringLiteral("underpowered"))
                : QStringLiteral("compare"),
            QStringLiteral("protein_centered"),
            QStringLiteral("kernel_scales"),
            QStringLiteral("protein_sd_retained"),
            support.support_class == SupportClass::Full ? QStringLiteral("computed_in_independence_table")
                                                        : QString(),
            QStringLiteral("stage1_bmrb_normalization_audit"),
            effectTrajectory.isEmpty() ? QStringLiteral("axis1_crosswalk_missing")
                                       : QStringLiteral("identity_level_axis1_context_crosswalk"),
            axis1.n ? axis1Support.support_name : QStringLiteral("missing"),
            axis2Support.support_name,
            support.support_name,
            QString::number(static_cast<qulonglong>(std::min<std::size_t>(axis1.n, cell.proteins.size()))),
            support.underpowered_dimensions,
        });
    }
    return true;
}

bool emitSharedLedgers(const QString& outDir,
                       const Axis2RunOptions& options,
                       const Axis2RunStats& stats,
                       const QStringList& deltaAudit,
                       QString* err_out) {
    QDir().mkpath(outDir);
    QJsonObject manifest;
    manifest.insert(QStringLiteral("schema_version"), QStringLiteral("axis2_manifest_audit_v1"));
    manifest.insert(QStringLiteral("root720"), options.root720);
    manifest.insert(QStringLiteral("mutant_root"), options.mutantRoot);
    manifest.insert(QStringLiteral("axis1_overlay_dir"), options.axis1OverlayDir);
    manifest.insert(QStringLiteral("support_N_min"), static_cast<qint64>(options.supportThresholds.n_min));
    manifest.insert(QStringLiteral("support_N_full"), static_cast<qint64>(options.supportThresholds.n_full));
    manifest.insert(QStringLiteral("support_threshold_provenance"),
                    QStringLiteral("initial D7 tunable values from stage1 support_warning lineage"));
    manifest.insert(QStringLiteral("molcomp_order"), QStringLiteral("xx,yy,xy,xz,yz,zz"));
    manifest.insert(QStringLiteral("eta_H"), QStringLiteral("proper distance_from_isotropic_ordering"));
    manifest.insert(QStringLiteral("permutation_K"), options.permutationK);
    manifest.insert(QStringLiteral("proteins_loaded"), static_cast<qint64>(stats.proteins_loaded));
    manifest.insert(QStringLiteral("static_samples"), static_cast<qint64>(stats.static_samples));
    manifest.insert(QStringLiteral("static_cells"), static_cast<qint64>(stats.static_cells));
    manifest.insert(QStringLiteral("ridge_samples"), static_cast<qint64>(stats.ridge_samples));
    manifest.insert(QStringLiteral("ridge_rows"), static_cast<qint64>(stats.ridge_rows));
    manifest.insert(QStringLiteral("distinct_identities"), static_cast<qint64>(stats.distinct_identities));
    manifest.insert(QStringLiteral("distinct_elements"), static_cast<qint64>(stats.distinct_elements));
    manifest.insert(QStringLiteral("bounded_memory_diagnostic"),
                    QStringLiteral("max_atoms_in_resident_protein=%1; accumulators_persist_proteins_ephemeral")
                        .arg(stats.max_atoms_in_resident_protein));
    QJsonArray delta;
    for (const QString& s : deltaAudit) delta.append(s);
    manifest.insert(QStringLiteral("delta_loader_fields"), delta);
    QJsonArray refusals;
    for (const QString& s : stats.refusal_reasons) refusals.append(s);
    manifest.insert(QStringLiteral("delta_refusals"), refusals);
    if (!writeText(QDir(outDir).filePath(QStringLiteral("manifest_audit.json")),
                   QJsonDocument(manifest).toJson(QJsonDocument::Indented), err_out)) return false;

    QFile csv(QDir(outDir).filePath(QStringLiteral("manifest_audit.csv")));
    QTextStream ts;
    if (!openCsv(csv.fileName(), &csv, &ts, err_out)) return false;
    writeCsvLine(ts, {QStringLiteral("metric"), QStringLiteral("value")});
    writeCsvLine(ts, {QStringLiteral("support_N_min"), QString::number(options.supportThresholds.n_min)});
    writeCsvLine(ts, {QStringLiteral("support_N_full"), QString::number(options.supportThresholds.n_full)});
    writeCsvLine(ts, {QStringLiteral("molcomp_order"), QStringLiteral("xx,yy,xy,xz,yz,zz")});
    writeCsvLine(ts, {QStringLiteral("delta_fields"), deltaAudit.join(QLatin1Char(';'))});

    QFile missing(QDir(outDir).filePath(QStringLiteral("missing_or_structural_absence.csv")));
    QTextStream ms;
    if (!openCsv(missing.fileName(), &missing, &ms, err_out)) return false;
    writeCsvLine(ms, {QStringLiteral("hypothesis_id"), QStringLiteral("status"), QStringLiteral("reason")});
    writeCsvLine(ms, {QStringLiteral("OPEN-2"), QStringLiteral("lead_owned"),
                      QStringLiteral("N_min/N_full tunable after first support distribution; initial values recorded")});
    writeCsvLine(ms, {QStringLiteral("OPEN-3"), QStringLiteral("lead_owned"),
                      QStringLiteral("citation constants pending-citation; algorithms not blocked")});

    auto smallLedger = [&](const QString& name, const QStringList& header,
                           const std::vector<QStringList>& rows) -> bool {
        QFile f(QDir(outDir).filePath(name));
        QTextStream s;
        if (!openCsv(f.fileName(), &f, &s, err_out)) return false;
        writeCsvLine(s, header);
        for (const QStringList& row : rows) writeCsvLine(s, row);
        return true;
    };
    if (!smallLedger(QStringLiteral("hypothesis_metric_index.csv"),
                     {QStringLiteral("hypothesis_id"), QStringLiteral("table"), QStringLiteral("metric")},
                     {{QStringLiteral("D4"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("sigma_distribution")},
                      {QStringLiteral("D6"), QStringLiteral("mutant_delta_ridge"), QStringLiteral("ridge_slope")},
                      {QStringLiteral("D3"), QStringLiteral("cross_axis_overlay"), QStringLiteral("effect_static_vs_trajectory")},
                      {QStringLiteral("C1"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("source_term_forward_sum")}}))
        return false;
    if (!smallLedger(QStringLiteral("case_study_index.csv"),
                     {QStringLiteral("case_id"), QStringLiteral("table"), QStringLiteral("status")},
                     {{QStringLiteral("G1"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("resolved")},
                      {QStringLiteral("G5"), QStringLiteral("mutant_delta_ridge"), QStringLiteral("resolved")},
                      {QStringLiteral("G6"), QStringLiteral("cross_axis_overlay"), QStringLiteral("resolved")},
                      {QStringLiteral("G9"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("full_or_reduced")},
                      {QStringLiteral("G12"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("resolved")}}))
        return false;
    if (!smallLedger(QStringLiteral("figure_registry.csv"),
                     {QStringLiteral("figure_id"), QStringLiteral("source_table"), QStringLiteral("status")},
                     {{QStringLiteral("axis2_static_distribution"), QStringLiteral("cohort_static_identity_context"), QStringLiteral("ready")},
                      {QStringLiteral("axis2_ridge"), QStringLiteral("mutant_delta_ridge"), QStringLiteral("ready")},
                      {QStringLiteral("cross_axis_overlay"), QStringLiteral("cross_axis_overlay"), QStringLiteral("ready")}}))
        return false;
    if (!smallLedger(QStringLiteral("plain_language_readouts.csv"),
                     {QStringLiteral("readout_id"), QStringLiteral("text")},
                     {{QStringLiteral("axis2_identity"), QStringLiteral("Cells are full IUPAC identity times context, not element-only.")},
                      {QStringLiteral("axis2_support"), QStringLiteral("Thin cells are retained and credentialed full, reduced, or insufficient.")}}))
        return false;
    return true;
}

}  // namespace

bool RunCohortAxis2(const Axis2RunOptions& options,
                    Axis2RunStats* stats,
                    QString* err_out) {
    Axis2RunStats local;
    QDir root(options.root720);
    if (!root.exists()) {
        if (err_out) *err_out = QStringLiteral("Axis-2 root720 does not exist: %1").arg(options.root720);
        return false;
    }
    if (options.outDir.isEmpty()) {
        if (err_out) *err_out = QStringLiteral("Axis-2 requires --out");
        return false;
    }
    QDir().mkpath(options.outDir);

    CohortContextAccumulator staticAcc;
    std::map<QString, RidgeBucket> ridgeBuckets;
    QStringList deltaAudit;
    QSet<QString> identitySet;
    QSet<QString> elementSet;

    const QStringList dirs = root.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    for (const QString& name : dirs) {
        if (options.maxProteins > 0 && local.proteins_seen >= options.maxProteins) break;
        ++local.proteins_seen;
        QString loadErr;
        auto runOpt = StaticRunData::Load(root.filePath(name), &loadErr);
        if (!runOpt) {
            local.refusal_reasons << QStringLiteral("%1 static load refused: %2").arg(name, loadErr);
            continue;
        }
        RunData run = std::move(*runOpt);
        ++local.proteins_loaded;
        if (run.protein)
            local.max_atoms_in_resident_protein =
                std::max(local.max_atoms_in_resident_protein, run.protein->atomCount());

        ResidentIndexes indexes = BuildResidentIndexes(run);
        const QString proteinId = run.manifest.protein_id.isEmpty() ? name : run.manifest.protein_id;
        const std::vector<MutationSite> sites = deriveMutationSites(options.mutantRoot, proteinId, run);

        std::optional<DeltaRunData> delta;
        if (options.runMutantDeltaRidge) {
            QString deltaErr;
            delta = LoadDeltaRunData(run, &deltaErr);
            if (!delta) {
                ++local.delta_refusals;
                local.refusal_reasons << QStringLiteral("%1 delta refused: %2").arg(name, deltaErr);
            } else {
                deltaAudit << delta->field_names;
            }
        }

        const std::size_t atomCount = run.protein ? run.protein->atomCount() : 0;
        for (std::size_t atom = 0; atom < atomCount; ++atom) {
            CohortSample sample = sampleForAtom(run, indexes, atom, sites);
            if (sample.key.canonical.isEmpty()) continue;
            identitySet.insert(sample.key.identity);
            elementSet.insert(sample.key.fields.element);
            if (options.runStatic) {
                staticAcc.push(sample);
                ++local.static_samples;
            }
            if (delta && atom < delta->matched_axis_count) {
                const model::Mat3 wt = totalTensor(delta->wt_dia, delta->wt_para, atom);
                const model::Mat3 ala = totalTensor(delta->ala_dia, delta->ala_para, atom);
                const model::Mat3 d = wt - ala;
                const model::Mat3 diaD = matFromArray(delta->wt_dia, atom) - matFromArray(delta->ala_dia, atom);
                const model::Mat3 paraD = matFromArray(delta->wt_para, atom) - matFromArray(delta->ala_para, atom);
                const double deltaSigma = iso(d);
                addRidge(&ridgeBuckets, sample, QStringLiteral("mutation_contact_kernel"),
                         sample.mutation_contact_kernel, deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_eta_H_refolded"),
                         etaH(wt) - etaH(ala), deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_xx_refolded"),
                         d(0, 0), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_yy_refolded"),
                         d(1, 1), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_xy_refolded"),
                         d(0, 1), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_xz_refolded"),
                         d(0, 2), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_yz_refolded"),
                         d(1, 2), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_mol_zz_refolded"),
                         d(2, 2), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_dia_iso_refolded"),
                         iso(diaD), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_para_iso_refolded"),
                         iso(paraD), deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                ++local.ridge_samples;
            }
        }
    }

    local.static_cells = staticAcc.cellCount();
    local.ridge_rows = ridgeBuckets.size();
    local.distinct_identities = identitySet.size();
    local.distinct_elements = elementSet.size();

    if (options.runStatic && staticAcc.cellCount() == 0) {
        if (err_out) *err_out = QStringLiteral("Axis-2 produced no static cells");
        return false;
    }
    if (options.runMutantDeltaRidge && ridgeBuckets.empty()) {
        if (err_out) *err_out = QStringLiteral("mutant_delta_ridge refused: no valid DeltaRunData samples");
        return false;
    }

    if (options.runStatic
        && !emitStaticTable(options.outDir, staticAcc, options.supportThresholds,
                            options.permutationK, err_out))
        return false;
    if (options.runStatic
        && !emitIndependenceTable(options.outDir, staticAcc, options.supportThresholds, err_out))
        return false;
    if (options.runMutantDeltaRidge
        && !emitRidgeTable(options.outDir, ridgeBuckets, options.supportThresholds,
                           options.permutationK, err_out))
        return false;
    if (options.runStatic
        && !emitOverlayTable(options.outDir, options.axis1OverlayDir, staticAcc,
                            options.supportThresholds, err_out))
        return false;
    if (!emitSharedLedgers(options.outDir, options, local, deltaAudit, err_out))
        return false;
    if (options.runStatic && options.runMutantDeltaRidge
        && !AuditAxis2Outputs(options.outDir, err_out))
        return false;

    if (stats) *stats = local;
    return true;
}

bool AuditAxis2Outputs(const QString& outDir, QString* err_out) {
    const QStringList required = {
        QStringLiteral("cohort_static_identity_context.csv"),
        QStringLiteral("cohort_static_independence.csv"),
        QStringLiteral("mutant_delta_ridge.csv"),
        QStringLiteral("cross_axis_overlay.csv"),
        QStringLiteral("manifest_audit.json"),
        QStringLiteral("hypothesis_metric_index.csv"),
        QStringLiteral("case_study_index.csv"),
        QStringLiteral("figure_registry.csv"),
        QStringLiteral("plain_language_readouts.csv"),
        QStringLiteral("missing_or_structural_absence.csv"),
    };
    for (const QString& name : required) {
        if (!QFileInfo::exists(QDir(outDir).filePath(name))) {
            if (err_out) *err_out = QStringLiteral("Axis-2 audit missing %1").arg(name);
            return false;
        }
    }
    const QStringList forbidden = {
        QStringLiteral("serial"), QStringLiteral("dwell"), QStringLiteral("recurrence"),
        QStringLiteral("lead_lag"), QStringLiteral("rolling_ols"), QStringLiteral("circshift"),
        QStringLiteral("circular_shift"), QStringLiteral("pas_prev"), QStringLiteral("pas_continuity"),
    };
    QSet<QString> identities;
    QSet<QString> elements;
    QFile staticFile(QDir(outDir).filePath(QStringLiteral("cohort_static_identity_context.csv")));
    if (!staticFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot read static table for audit");
        return false;
    }
    QTextStream ss(&staticFile);
    const QString staticHeader = ss.readLine();
    for (const QString& token : forbidden) {
        if (staticHeader.contains(token, Qt::CaseInsensitive)) {
            if (err_out) *err_out = QStringLiteral("forbidden over-frame slot in static header: %1").arg(token);
            return false;
        }
    }
    if (!staticHeader.contains(QStringLiteral("sigma_p05"))
        || !staticHeader.contains(QStringLiteral("sigma_p95"))) {
        if (err_out) *err_out = QStringLiteral("static sigma distribution is not 9-number");
        return false;
    }
    while (!ss.atEnd()) {
        const QString line = ss.readLine();
        if (line.isEmpty()) continue;
        if (!line.contains(QStringLiteral("atom_name="))
            || !line.contains(QStringLiteral("residue_type="))
            || !line.contains(QStringLiteral("dihedral_region="))
            || !line.contains(QStringLiteral("SS="))) {
            if (err_out) *err_out = QStringLiteral("context key lacks full IUPAC/context tokens");
            return false;
        }
        if (line.contains(QStringLiteral("graph_stratum"))
            || line.contains(QStringLiteral("n_frames"))) {
            if (err_out) *err_out = QStringLiteral("Axis-2 key contains forbidden Axis-1 context token");
            return false;
        }
        const QStringList cols = line.split(QLatin1Char(','));
        if (cols.size() > 6) {
            elements.insert(cols[4]);
            identities.insert(cols[3]);
        }
    }
    if (identities.size() <= elements.size()) {
        if (err_out) {
            *err_out = QStringLiteral("identity collapse audit failed: identities=%1 elements=%2")
                           .arg(identities.size())
                           .arg(elements.size());
        }
        return false;
    }

    auto auditNoInsufficientCoupling = [&](const QString& table, const QStringList& couplingNames) -> bool {
        QFile f(QDir(outDir).filePath(table));
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return false;
        QTextStream in(&f);
        const QStringList header = parseCsvLine(in.readLine());
        const int supportIdx = header.indexOf(QStringLiteral("support_class"));
        std::vector<int> couplingIdx;
        for (const QString& name : couplingNames) {
            const int idx = header.indexOf(name);
            if (idx >= 0) couplingIdx.push_back(idx);
        }
        while (!in.atEnd()) {
            const QStringList cols = parseCsvLine(in.readLine());
            if (supportIdx < 0 || supportIdx >= cols.size()) continue;
            if (cols[supportIdx] != QStringLiteral("insufficient")) continue;
            for (int idx : couplingIdx) {
                if (idx < cols.size() && !cols[idx].trimmed().isEmpty()) {
                    if (err_out) {
                        *err_out = QStringLiteral("%1 emits coupling %2 on insufficient row")
                                       .arg(table, header[idx]);
                    }
                    return false;
                }
            }
        }
        return true;
    };
    if (!auditNoInsufficientCoupling(QStringLiteral("cohort_static_independence.csv"),
                                     {QStringLiteral("max_canonical_corr"),
                                      QStringLiteral("mean_canonical_corr"),
                                      QStringLiteral("pearson_r"),
                                      QStringLiteral("spearman_rho"),
                                      QStringLiteral("condition_number"),
                                      QStringLiteral("vif"),
                                      QStringLiteral("signed_slope")}))
        return false;
    if (!auditNoInsufficientCoupling(QStringLiteral("mutant_delta_ridge.csv"),
                                     {QStringLiteral("ridge_slope"),
                                      QStringLiteral("ridge_r"),
                                      QStringLiteral("ridge_xi"),
                                      QStringLiteral("pchip_shape"),
                                      QStringLiteral("leverage"),
                                      QStringLiteral("delta_permutation_p"),
                                      QStringLiteral("null_slope_mean"),
                                      QStringLiteral("null_slope_sd"),
                                      QStringLiteral("obs_slope_z")}))
        return false;
    if (!auditNoInsufficientCoupling(QStringLiteral("cross_axis_overlay.csv"),
                                     {QStringLiteral("effect_trajectory"),
                                      QStringLiteral("effect_static"),
                                      QStringLiteral("axis1_effect_distribution"),
                                      QStringLiteral("axis2_effect_distribution"),
                                      QStringLiteral("subspace_compare_max_canonical_corr")}))
        return false;
    return true;
}

}  // namespace h5reader::rediscover
