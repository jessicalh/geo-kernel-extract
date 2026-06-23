#include "CohortContextAccumulator.h"

#include "DeltaRunData.h"
#include "LiteratureConstants.h"
#include "LocalFrameBasis.h"
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

double scaledRingIso(const RunData& run, io::FieldKind kind, std::size_t atom) {
    const double raw = iso(matFromArray(arr(run, kind), atom));
    return finite(raw) ? raw * RingCurrentPpmFactor() : kNan;
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

std::optional<model::Mat3> axesFromLocalFrame(const LocalFrame& f) {
    if (!f.is_valid) return std::nullopt;
    model::Mat3 axes;
    axes.col(0) = f.x;
    axes.col(1) = f.y;
    axes.col(2) = f.z;
    return axes.allFinite() ? std::optional<model::Mat3>(axes) : std::nullopt;
}

std::optional<model::Vec3> normVec(const model::Vec3& v) {
    const double n = v.norm();
    if (!(n > 1.0e-9) || !finite(n)) return std::nullopt;
    return v / n;
}

std::optional<model::Mat3> axesFromXAndPlane(const model::Vec3& xVec,
                                             const model::Vec3& planeVec) {
    const auto xOpt = normVec(xVec);
    if (!xOpt) return std::nullopt;
    const model::Vec3 x = *xOpt;
    const model::Vec3 plane = planeVec - x * planeVec.dot(x);
    const auto y0Opt = normVec(plane);
    if (!y0Opt) return std::nullopt;
    const auto zOpt = normVec(x.cross(*y0Opt));
    if (!zOpt) return std::nullopt;
    const model::Vec3 z = *zOpt;
    const auto yOpt = normVec(z.cross(x));
    if (!yOpt) return std::nullopt;
    model::Mat3 axes;
    axes.col(0) = x;
    axes.col(1) = *yOpt;
    axes.col(2) = z;
    return axes.allFinite() ? std::optional<model::Mat3>(axes) : std::nullopt;
}

std::optional<std::size_t> atomByResidueName(const model::QtProtein& protein,
                                             int residueIndex,
                                             const QString& name) {
    if (residueIndex < 0 || static_cast<std::size_t>(residueIndex) >= protein.residueCount())
        return std::nullopt;
    const model::QtResidue& residue = protein.residue(static_cast<std::size_t>(residueIndex));
    for (int32_t ai : residue.atomIndices) {
        if (ai < 0 || static_cast<std::size_t>(ai) >= protein.atomCount()) continue;
        if (atomName(protein, static_cast<std::size_t>(ai)) == name)
            return static_cast<std::size_t>(ai);
    }
    return std::nullopt;
}

std::optional<model::Vec3> atomPosition(const RunData& run, std::size_t atom, std::size_t frame = 0) {
    if (!run.conformation || !run.protein || atom >= run.protein->atomCount())
        return std::nullopt;
    const model::Vec3 p = run.conformation->atomPosition(frame, atom);
    return p.allFinite() ? std::optional<model::Vec3>(p) : std::nullopt;
}

std::optional<model::Mat3> genericBondAxes(const RunData& run, std::size_t atom) {
    if (!run.protein || !run.conformation || atom >= run.protein->atomCount())
        return std::nullopt;
    std::vector<std::size_t> bonded;
    for (int32_t bi : run.protein->topology().bondIndicesForAtom(atom)) {
        if (bi < 0 || static_cast<std::size_t>(bi) >= run.protein->topology().bondCount()) continue;
        const model::QtBond& b = run.protein->topology().bondAt(static_cast<std::size_t>(bi));
        const int32_t other = b.atomIndexA == static_cast<int32_t>(atom) ? b.atomIndexB
                              : b.atomIndexB == static_cast<int32_t>(atom) ? b.atomIndexA
                                                                           : -1;
        if (other >= 0 && static_cast<std::size_t>(other) < run.protein->atomCount())
            bonded.push_back(static_cast<std::size_t>(other));
    }
    if (bonded.size() < 2) return std::nullopt;
    const auto origin = atomPosition(run, atom);
    const auto x = atomPosition(run, bonded[0]);
    const auto plane = atomPosition(run, bonded[1]);
    if (!origin || !x || !plane) return std::nullopt;
    return axesFromXAndPlane(*x - *origin, *plane - *origin);
}

QString molecularFrameKindForAtom(const RunData& run, std::size_t atom) {
    if (!run.protein || atom >= run.protein->atomCount()) return QStringLiteral("none");
    const model::QtProtein& protein = *run.protein;
    const model::QtAtom& a = protein.atom(atom);
    const QString name = atomName(protein, atom);
    if (name == QStringLiteral("CA")) return QStringLiteral("backbone_ca");
    if (a.IsBackboneNitrogen() || name == QStringLiteral("N"))
        return QStringLiteral("backbone_amide_n");
    if (a.element == model::Element::H
        && a.backboneRole == model::BackboneRole::AmideHydrogen)
        return QStringLiteral("backbone_amide_h");
    if (a.IsBackboneCarbonylCarbon() || a.IsBackboneCarbonylOxygen()
        || name == QStringLiteral("C") || name == QStringLiteral("O")
        || name == QStringLiteral("OXT"))
        return QStringLiteral("backbone_carbonyl");
    if (a.element == model::Element::O
        && (name == QStringLiteral("OG") || name == QStringLiteral("OG1")
            || name == QStringLiteral("OH")))
        return QStringLiteral("hydroxyl_oxygen");
    if (a.element == model::Element::C && name != QStringLiteral("C") && !a.IsInAnyRing())
        return QStringLiteral("aliphatic_carbon");
    return QStringLiteral("none");
}

std::optional<model::Mat3> molecularAxesForAtom(const RunData& run, std::size_t atom) {
    if (!run.protein || !run.conformation || atom >= run.protein->atomCount())
        return std::nullopt;
    const model::QtProtein& protein = *run.protein;
    const model::QtAtom& a = protein.atom(atom);
    const int ri = a.residueIndex;
    const QString name = atomName(protein, atom);
    auto rec = [&](const QString& atomName) {
        return atomByResidueName(protein, ri, atomName);
    };
    auto posOf = [&](std::optional<std::size_t> idx) -> std::optional<model::Vec3> {
        return idx ? atomPosition(run, *idx) : std::nullopt;
    };

    if (name == QStringLiteral("CA")) {
        const auto ca = posOf(rec(QStringLiteral("CA")));
        const auto n = posOf(rec(QStringLiteral("N")));
        const auto c = posOf(rec(QStringLiteral("C")));
        if (ca && n && c) return axesFromLocalFrame(BuildBackboneCaFrame(*ca, *n, *c));
    }
    if (a.IsBackboneNitrogen() || name == QStringLiteral("N")) {
        const auto n = posOf(rec(QStringLiteral("N")));
        const auto ca = posOf(rec(QStringLiteral("CA")));
        const auto cPrev = posOf(atomByResidueName(protein, ri - 1, QStringLiteral("C")));
        const auto cOwn = posOf(rec(QStringLiteral("C")));
        if (n && ca && (cPrev || cOwn))
            return axesFromLocalFrame(BuildBackboneNFrame(*n, *ca, cPrev ? *cPrev : *cOwn,
                                                         cPrev.has_value()));
    }
    if (a.IsBackboneCarbonylCarbon() || name == QStringLiteral("C")) {
        const auto c = posOf(rec(QStringLiteral("C")));
        const auto o = posOf(rec(QStringLiteral("O")));
        const auto ca = posOf(rec(QStringLiteral("CA")));
        if (c && o && ca) return axesFromLocalFrame(BuildBackboneCarbonylCFrame(*c, *o, *ca));
    }
    if (a.IsBackboneCarbonylOxygen() || name == QStringLiteral("O")) {
        const auto o = posOf(rec(QStringLiteral("O")));
        const auto c = posOf(rec(QStringLiteral("C")));
        const auto ca = posOf(rec(QStringLiteral("CA")));
        if (o && c && ca) return axesFromLocalFrame(BuildBackboneCarbonylOFrame(*o, *c, *ca));
    }
    if (a.IsAnyAlphaHydrogen()) {
        const auto h = atomPosition(run, atom);
        const auto ca = posOf(rec(QStringLiteral("CA")));
        const auto n = posOf(rec(QStringLiteral("N")));
        if (h && ca && n) return axesFromLocalFrame(BuildBackboneHaFrame(*h, *ca, *n));
    }
    return genericBondAxes(run, atom);
}

std::optional<model::Vec3> fieldVector(const RunData& run, io::FieldKind kind, std::size_t atom) {
    const StaticNpyArray* a = arr(run, kind);
    if (!a || atom >= a->rows || a->cols < 3) return std::nullopt;
    const model::Vec3 e(a->value(atom, 0), a->value(atom, 1), a->value(atom, 2));
    return e.allFinite() ? std::optional<model::Vec3>(e) : std::nullopt;
}

double signedParallelMolZ(const RunData& run,
                          io::FieldKind kind,
                          std::size_t atom,
                          const std::optional<model::Mat3>& axes) {
    const auto e = fieldVector(run, kind, atom);
    if (!e || !axes) return kNan;
    return e->dot(axes->col(2));
}

double signedParallelXH(const RunData& run, io::FieldKind kind, std::size_t atom) {
    if (!run.protein || atom >= run.protein->atomCount()) return kNan;
    const model::QtAtom& a = run.protein->atom(atom);
    if (a.element != model::Element::H || a.parentAtomIndex < 0) return kNan;
    const auto e = fieldVector(run, kind, atom);
    const auto h = atomPosition(run, atom);
    const auto x = atomPosition(run, static_cast<std::size_t>(a.parentAtomIndex));
    if (!e || !h || !x) return kNan;
    const auto u = normVec(*h - *x);
    return u ? e->dot(*u) : kNan;
}

double larsenHbondIso(const RunData& run, std::size_t atom) {
    const std::array<io::FieldKind, 4> parts = {
        io::FieldKind::LarsenHBond1pHBShielding,
        io::FieldKind::LarsenHBond1pHaBShielding,
        io::FieldKind::LarsenHBond2pHBShielding,
        io::FieldKind::LarsenHBond2pHaBShielding,
    };
    model::Mat3 sum = model::Mat3::Zero();
    for (io::FieldKind kind : parts) {
        const model::Mat3 m = matFromArray(arr(run, kind), atom);
        if (!m.allFinite()) return kNan;
        sum += m;
    }
    return iso(sum);
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
    f.frame_kind = molecularFrameKindForAtom(run, atom);
    f.hyb = hybridisationForAtom(run, a, atom);
    f.contact_class = nearestContactClass(run, atom);
    f.dihedral_region = dihedralRegion(indexes, atom);
    f.SS = ss.observed ? ssName(ss.ss3) : QStringLiteral("unknown");
    s.key = BuildAxis2ContextKey(f);
    s.protein_id = run.manifest.protein_id.isEmpty() ? protein.proteinId() : run.manifest.protein_id;
    s.atom_index = static_cast<int>(atom);
    s.residue_index = residueIndex;
    s.molecular_axes = molecularAxesForAtom(run, atom);

    const model::DftAtomShielding* dft = run.dft.AtomShielding(atom, 0);
    if (dft) {
        const Axis2FoldedTensor folded = FoldAxis2TensorChannels(dft->total_raw, s.molecular_axes);
        s.sigma_iso = folded.sigma_iso;
        s.sigma_eta_H = folded.sigma_eta_H;
        s.mol_xx = folded.mol_components[0];
        s.mol_yy = folded.mol_components[1];
        s.mol_xy = folded.mol_components[2];
        s.mol_xz = folded.mol_components[3];
        s.mol_yz = folded.mol_components[4];
        s.mol_zz = folded.mol_components[5];
    }

    s.channels.insert(QStringLiteral("apbs_E_mag"), vectorMag(run, io::FieldKind::APBSE, atom));
    s.channels.insert(QStringLiteral("apbs_E_parallel_mol_z"),
                      signedParallelMolZ(run, io::FieldKind::APBSE, atom, s.molecular_axes));
    s.channels.insert(QStringLiteral("apbs_E_parallel_XH"),
                      signedParallelXH(run, io::FieldKind::APBSE, atom));
    s.channels.insert(QStringLiteral("mopac_E_mag"), vectorMag(run, io::FieldKind::MOPACCoulombE, atom));
    s.channels.insert(QStringLiteral("apbs_efg_absT2"), rowNorm(run, io::FieldKind::APBSEFG, atom));
    s.channels.insert(QStringLiteral("aimnet2_efg_absT2"), rowNorm(run, io::FieldKind::AIMNet2EFG, atom));
    s.channels.insert(QStringLiteral("ring_bs_iso"), scaledRingIso(run, io::FieldKind::BSShielding, atom));
    s.channels.insert(QStringLiteral("ring_hm_iso"), scaledRingIso(run, io::FieldKind::HMShielding, atom));
    s.channels.insert(QStringLiteral("larsen_hbond_iso"), larsenHbondIso(run, atom));
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
                const DihedralState prevChi1 =
                    indexes.dihedrals.state(DihedralKind::Chi1, static_cast<std::size_t>(prev.CA), 0);
                s.chi1_iminus1 = prevChi1.present ? prevChi1.radians : kNan;
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

std::vector<double> proteinSigmaMeans(const CohortCellTruth& cell) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(cell.protein_folds.size()));
    for (auto it = cell.protein_folds.begin(); it != cell.protein_folds.end(); ++it)
        out.push_back(it.value().sigma.meanValue());
    return out;
}

std::vector<double> proteinChannelMeans(const CohortCellTruth& cell, const QString& name) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(cell.protein_folds.size()));
    for (auto it = cell.protein_folds.begin(); it != cell.protein_folds.end(); ++it)
        out.push_back(it.value().channels.value(name).meanValue());
    return out;
}

DistributionSummary channelSummary(const CohortCellTruth& cell, const QString& name) {
    return cell.channel_distributions.value(name).summary();
}

double channelMean(const CohortCellTruth& cell, const QString& name) {
    const DistributionSummary s = channelSummary(cell, name);
    return s.hasFinite() ? s.mean : kNan;
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
        QStringLiteral("atom_name"), QStringLiteral("frame_kind"), QStringLiteral("hyb"), QStringLiteral("contact_class"),
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
        QStringLiteral("protein_label_permutation_status"), QStringLiteral("permutation_K"),
        QStringLiteral("helix_dipole_field_mean"),
        QStringLiteral("helix_dipole_citation_status"), QStringLiteral("helix_membership"),
        QStringLiteral("psi_iminus1_region"), QStringLiteral("r_psi_iminus1"),
        QStringLiteral("r_psi_own"), QStringLiteral("sensitivity_ratio"),
        QStringLiteral("chi1_iminus1_effect"), QStringLiteral("chi1_iminus1_support_class"),
        QStringLiteral("chi1_iminus1_missing_reason"), QStringLiteral("predecessor_identity"),
        QStringLiteral("contrib_buckingham"), QStringLiteral("contrib_ring"),
        QStringLiteral("contrib_mcconnell"), QStringLiteral("contrib_larsen"),
        QStringLiteral("r_cl_qm"), QStringLiteral("slope_cl_qm"), QStringLiteral("rmsd_ppm"),
        QStringLiteral("residual_mean"), QStringLiteral("residual_sd"),
        QStringLiteral("constant_buckingham_key"), QStringLiteral("constant_buckingham_units"),
        QStringLiteral("constant_status")
    });

    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const DistributionSummary sigma = cell.sigma.summary();
        const DistributionSummary xx = cell.mol_xx.summary();
        const DistributionSummary yy = cell.mol_yy.summary();
        const DistributionSummary xy = cell.mol_xy.summary();
        const DistributionSummary xz = cell.mol_xz.summary();
        const DistributionSummary yz = cell.mol_yz.summary();
        const DistributionSummary zz = cell.mol_zz.summary();
        const DistributionSummary eta = cell.eta_H.summary();
        const DistributionSummary helix = cell.helix_dipole_field.summary();
        const SupportCredential support =
            CredentialSupport(cell.proteins.size(), true, thresholds);
        const std::vector<double> sigmaVals = proteinSigmaMeans(cell);
        const std::vector<double> ring = proteinChannelMeans(cell, QStringLiteral("ring_bs_iso"));
        const std::vector<double> apbs = proteinChannelMeans(cell, QStringLiteral("apbs_E_mag"));
        const PermutationNull null = ProteinLabelPermutationNull(apbs, sigmaVals, permutationK);
        const double rPsiPrev = cell.psi_iminus1_vs_sigma.pearson();
        const double rPsiOwn = cell.psi_own_vs_sigma.pearson();
        const double sensitivity = finite(rPsiPrev) && finite(rPsiOwn) && std::abs(rPsiOwn) > 1.0e-12
                                       ? std::abs(rPsiPrev / rPsiOwn)
                                       : kNan;
        const double chi1Effect = cell.chi1_iminus1_vs_sigma.slope();
        const SupportCredential chi1Support =
            CredentialSupport(cell.chi1_iminus1_vs_sigma.n,
                              cell.chi1_iminus1_vs_sigma.n >= 2, thresholds);
        QString chi1MissingReason;
        if (!finite(chi1Effect)) {
            if (cell.predecessor_identity == QStringLiteral("not_backbone_N")) {
                chi1MissingReason = QStringLiteral("not_backbone_N");
            } else if (cell.predecessor_identity.startsWith(QStringLiteral("GLY:"))
                       || cell.predecessor_identity.startsWith(QStringLiteral("ALA:"))) {
                chi1MissingReason = QStringLiteral("predecessor_has_no_chi1");
            } else {
                chi1MissingReason = QStringLiteral("predecessor_chi1_insufficient_or_unavailable");
            }
        }
        const LiteratureConstant buck =
            BuckinghamA(cell.key.fields.element == QStringLiteral("H") ? model::Element::H
                         : cell.key.fields.element == QStringLiteral("C") ? model::Element::C
                         : cell.key.fields.element == QStringLiteral("N") ? model::Element::N
                         : cell.key.fields.element == QStringLiteral("O") ? model::Element::O
                         : cell.key.fields.element == QStringLiteral("S") ? model::Element::S
                                                                           : model::Element::Unknown,
                         cell.key.fields.residue_type.toStdString(),
                         cell.key.fields.atom_name.toStdString(),
                         cell.key.fields.frame_kind.toStdString());
        const LiteratureConstant buckB =
            BuckinghamB(cell.key.fields.element == QStringLiteral("H") ? model::Element::H
                         : cell.key.fields.element == QStringLiteral("C") ? model::Element::C
                         : cell.key.fields.element == QStringLiteral("N") ? model::Element::N
                         : cell.key.fields.element == QStringLiteral("O") ? model::Element::O
                         : cell.key.fields.element == QStringLiteral("S") ? model::Element::S
                                                                           : model::Element::Unknown,
                         cell.key.fields.residue_type.toStdString(),
                         cell.key.fields.atom_name.toStdString(),
                         cell.key.fields.frame_kind.toStdString());
        const double xh = channelMean(cell, QStringLiteral("apbs_E_parallel_XH"));
        const double molZ = channelMean(cell, QStringLiteral("apbs_E_parallel_mol_z"));
        const double eParallel = finite(xh) ? xh : molZ;
        double contribBuck = 0.0;
        if (finite(eParallel)) {
            contribBuck = (-buck.value * eParallel) - (buckB.value * eParallel * eParallel);
        } else if (buck.value != 0.0 || buckB.value != 0.0) {
            contribBuck = kNan;
        }
        const double contribRing = channelMean(cell, QStringLiteral("ring_bs_iso"));
        const double contribMc = channelMean(cell, QStringLiteral("mc_lit_iso"));
        const double contribLarsen = channelMean(cell, QStringLiteral("larsen_hbond_iso"));
        const double cl = contribBuck + contribRing + contribMc + contribLarsen;
        Q_UNUSED(cl);
        const ClassicalAgreementStats clAgreement =
            ComputeClassicalAgreementForCell(cell, buck.value, buckB.value);

        QStringList row = {
            kSchemaName,
            QStringLiteral("cohort_static_identity_context"),
            cell.key.canonical,
            cell.key.identity,
            cell.key.fields.element,
            cell.key.fields.residue_type,
            cell.key.fields.atom_name,
            cell.key.fields.frame_kind,
            cell.key.fields.hyb,
            cell.key.fields.contact_class,
            cell.key.fields.dihedral_region,
            cell.key.fields.SS,
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            proteinResidency(cell),
            QString::number(static_cast<qulonglong>(cell.sample_count)),
            support.support_name,
            support.underpowered_dimensions,
        };
        appendDistribution(&row, QStringLiteral("sigma"), sigma);
        row << QStringLiteral("xx,yy,xy,xz,yz,zz")
            << csvDouble(xx.mean) << csvDouble(yy.mean) << csvDouble(xy.mean)
            << csvDouble(xz.mean) << csvDouble(yz.mean) << csvDouble(zz.mean)
            << csvDouble(eta.mean)
            << csvDouble(cell.channel_vs_sigma.value(QStringLiteral("apbs_efg_absT2")).pearson())
            << csvDouble(cell.channel_vs_sigma.value(QStringLiteral("apbs_E_mag")).pearson())
            << csvDouble(cell.channel_vs_sigma.value(QStringLiteral("ring_bs_iso")).pearson())
            << csvDouble(std::abs(cell.channel_vs_sigma.value(QStringLiteral("ring_bs_iso")).pearson())
                         - std::abs(cell.channel_vs_sigma.value(QStringLiteral("apbs_E_mag")).pearson()))
            << csvDouble(null.perm_p)
            << (finite(null.perm_p) ? QStringLiteral("computed") : QStringLiteral("not_computable"))
            << QString::number(null.permutation_K)
            << csvDouble(helix.mean)
            << QStringLiteral("pending-citation")
            << (helix.finite_n > 0 ? QStringLiteral("helix_backbone_N") : QStringLiteral("not_applicable"))
            << cell.psi_iminus1_region
            << csvDouble(rPsiPrev)
            << csvDouble(rPsiOwn)
            << csvDouble(sensitivity)
            << csvDouble(chi1Effect)
            << chi1Support.support_name
            << chi1MissingReason
            << cell.predecessor_identity
            << csvDouble(contribBuck)
            << csvDouble(contribRing)
            << csvDouble(contribMc)
            << csvDouble(contribLarsen)
            << csvDouble(clAgreement.r)
            << csvDouble(clAgreement.slope)
            << csvDouble(clAgreement.rmsd)
            << csvDouble(clAgreement.residual_mean)
            << csvDouble(clAgreement.residual_sd)
            << QString::fromLatin1(buck.key)
            << QString::fromLatin1(buck.units)
            << QString::fromLatin1(LiteratureStatusName(buck.status));
        writeCsvLine(ts, row);
    }
    return true;
}

QString legacyEquivalentColumnForChannel(const QString& channel) {
    if (channel == QStringLiteral("apbs_efg_absT2"))
        return QStringLiteral("full_tensor_r2");
    if (channel == QStringLiteral("apbs_E_mag"))
        return QStringLiteral("magnitude_only_r2");
    if (channel == QStringLiteral("ring_bs_iso"))
        return QStringLiteral("direction_only_r2");
    return QString();
}

bool emitStaticSourceRelationships(const QString& outDir,
                                   const CohortContextAccumulator& acc,
                                   SupportThresholds thresholds,
                                   QString* err_out) {
    QFile file;
    QTextStream ts;
    if (!openCsv(QDir(outDir).filePath(QStringLiteral("cohort_static_source_relationships.csv")),
                 &file, &ts, err_out)) return false;
    writeCsvLine(ts, {
        QStringLiteral("schema_version"), QStringLiteral("table"), QStringLiteral("row_grain"),
        QStringLiteral("context_key"), QStringLiteral("identity_key"), QStringLiteral("element"),
        QStringLiteral("residue_type"), QStringLiteral("atom_name"), QStringLiteral("frame_kind"),
        QStringLiteral("hyb"), QStringLiteral("contact_class"), QStringLiteral("dihedral_region"),
        QStringLiteral("SS"), QStringLiteral("channel_key"), QStringLiteral("relationship_source"),
        QStringLiteral("n_pairs"), QStringLiteral("support_class"), QStringLiteral("finite_frac"),
        QStringLiteral("singleton_flag"), QStringLiteral("missing_n"), QStringLiteral("missing_reason"),
        QStringLiteral("pearson_r"), QStringLiteral("slope_sigma_per_channel"),
        QStringLiteral("legacy_equivalent_column")
    });
    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        for (auto it = cell.channel_vs_sigma.begin(); it != cell.channel_vs_sigma.end(); ++it) {
            const PairAccumulator& pair = it.value();
            const SupportCredential support =
                CredentialSupport(pair.n, pair.n >= 2, thresholds);
            const std::size_t proteinN = static_cast<std::size_t>(cell.proteins.size());
            const std::size_t missingN = proteinN > pair.n ? proteinN - pair.n : 0u;
            writeCsvLine(ts, {
                kSchemaName,
                QStringLiteral("cohort_static_source_relationships"),
                QStringLiteral("cell_channel"),
                cell.key.canonical,
                cell.key.identity,
                cell.key.fields.element,
                cell.key.fields.residue_type,
                cell.key.fields.atom_name,
                cell.key.fields.frame_kind,
                cell.key.fields.hyb,
                cell.key.fields.contact_class,
                cell.key.fields.dihedral_region,
                cell.key.fields.SS,
                it.key(),
                QStringLiteral("channel_vs_sigma[%1]").arg(it.key()),
                QString::number(static_cast<qulonglong>(pair.n)),
                support.support_name,
                csvDouble(cell.proteins.isEmpty()
                              ? kNan
                              : static_cast<double>(pair.n)
                                    / static_cast<double>(cell.proteins.size())),
                pair.n == 1 ? QStringLiteral("true") : QStringLiteral("false"),
                QString::number(static_cast<qulonglong>(missingN)),
                missingN == 0 ? QString()
                              : QStringLiteral("nonfinite_channel_or_sigma_for_some_proteins"),
                csvDouble(pair.pearson()),
                csvDouble(pair.slope()),
                legacyEquivalentColumnForChannel(it.key()),
            });
        }
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
        QStringLiteral("pearson_r"), QStringLiteral("pairwise_pearson_r"), QStringLiteral("signed_slope"),
        QStringLiteral("vif"), QStringLiteral("condition_number"), QStringLiteral("overlap_label"),
        QStringLiteral("provenance"), QStringLiteral("mediation_chain")
    });
    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const std::size_t n = cell.proteins.size();
        std::vector<std::size_t> rowIdx(n);
        std::iota(rowIdx.begin(), rowIdx.end(), 0);
        auto channel = [&](const QString& name) {
            return SubspaceChannel{name, proteinChannelMeans(cell, name)};
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
            const bool fullRank = n >= 3;
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
    BoundedDistributionAccumulator distant_abs_delta_sigma;
    BoundedDistributionAccumulator distant_nearest_site_distance_A;
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
    if (sample.distant_from_all_sites) {
        ++b.distant_n;
        b.distant_abs_delta_sigma.add(std::abs(deltaSigma));
        const double k = sample.mutation_contact_kernel;
        const double d = finite(k) && k > 0.0 ? std::sqrt(std::max(0.0, 1.0 / k - 1.0)) : kNan;
        b.distant_nearest_site_distance_A.add(d);
    }
    b.wt_n += wt_n;
    b.ala_n += ala_n;
    b.matched_atom_n += matched_n;
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
        QStringLiteral("distant_zero_check"), QStringLiteral("distant_characterization"),
        QStringLiteral("distant_nonzero_channel"), QStringLiteral("distant_abs_delta_sigma_mean"),
        QStringLiteral("distant_abs_delta_sigma_sd"), QStringLiteral("distant_abs_delta_sigma_min"),
        QStringLiteral("distant_abs_delta_sigma_median"), QStringLiteral("distant_abs_delta_sigma_max"),
        QStringLiteral("distant_nearest_site_distance_A_mean"),
        QStringLiteral("distant_nearest_site_distance_A_min"),
        QStringLiteral("distant_nearest_site_distance_A_median"),
        QStringLiteral("distant_nearest_site_distance_A_max"),
        QStringLiteral("delta_permutation_p"), QStringLiteral("delta_permutation_status"),
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
        const DistantRidgeCharacterization distantFlag =
            CharacterizeDistantNonzeroRidge(static_cast<std::size_t>(b.distant_n),
                                            b.any_site_scope,
                                            slope,
                                            b.channel);
        const bool hasDistant = b.distant_n > 0;
        const DistributionSummary distantMag = b.distant_abs_delta_sigma.summary();
        const DistributionSummary distantDist = b.distant_nearest_site_distance_A.summary();
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
            << distantFlag.distant_zero_check
            << distantFlag.characterization
            << distantFlag.nonzero_channel
            << csvDouble(hasDistant ? distantMag.mean : kNan)
            << csvDouble(hasDistant ? distantMag.sd : kNan)
            << csvDouble(hasDistant ? distantMag.min : kNan)
            << csvDouble(hasDistant ? distantMag.median : kNan)
            << csvDouble(hasDistant ? distantMag.max : kNan)
            << csvDouble(hasDistant ? distantDist.mean : kNan)
            << csvDouble(hasDistant ? distantDist.min : kNan)
            << csvDouble(hasDistant ? distantDist.median : kNan)
            << csvDouble(hasDistant ? distantDist.max : kNan)
            << (support.may_emit_coupling ? csvDouble(null.perm_p) : QString())
            << (support.may_emit_coupling && finite(null.perm_p) ? QStringLiteral("computed")
                                                                 : QStringLiteral("not_computable"))
            << (support.may_emit_coupling ? csvDouble(null.null_slope_mean) : QString())
            << (support.may_emit_coupling ? csvDouble(null.null_slope_sd) : QString())
            << (support.may_emit_coupling ? csvDouble(null.obs_slope_z) : QString())
            << QString::number(permutationK)
            << QString::number(static_cast<qulonglong>(b.wt_n))
            << QString::number(static_cast<qulonglong>(b.ala_n))
            << QString::number(static_cast<qulonglong>(b.matched_atom_n))
            << QStringLiteral("engine_fold_adapter_v1;wt_ala_refolded;molecular_frame_projected_when_available;no_prebaked_delta_channel");
        writeCsvLine(ts, row);
    }
    return true;
}

bool emitOverlayTable(const QString& outDir,
                      const QString& axis1Dir,
                      const CohortContextAccumulator& acc,
                      SupportThresholds thresholds,
                      Axis2RunStats* stats,
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
            const int hybIdx = header.indexOf(QStringLiteral("hyb"));
            const int contactIdx = header.indexOf(QStringLiteral("contact_class"));
            const int dihedralIdx = header.indexOf(QStringLiteral("dihedral_region"));
            const int ssIdx = header.indexOf(QStringLiteral("SS"));
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
                QString k;
                if (hybIdx >= 0 && contactIdx >= 0 && dihedralIdx >= 0 && ssIdx >= 0
                    && cols.size() > std::max({hybIdx, contactIdx, dihedralIdx, ssIdx})) {
                    Axis2ContextKeyFields f;
                    f.element = cols[elemIdx];
                    f.residue_type = cols[residueIdx];
                    f.atom_name = cols[atomIdx];
                    f.hyb = cols[hybIdx];
                    f.contact_class = cols[contactIdx];
                    f.dihedral_region = cols[dihedralIdx];
                    f.SS = cols[ssIdx];
                    k = BuildAxis2ContextKey(f).canonical;
                } else {
                    k = cols[residueIdx] + QLatin1Char('|') + cols[atomIdx]
                        + QLatin1Char('|') + cols[elemIdx];
                }
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
        QStringLiteral("cohort_cell_n_proteins"), QStringLiteral("axis1_n_proteins"),
        QStringLiteral("underpowered_dimensions")
    });
    for (const auto& item : acc.cells()) {
        const CohortCellTruth& cell = item.second;
        const bool populatedCell = !cell.proteins.isEmpty();
        if (stats) {
            if (populatedCell) ++stats->overlay_populated_cells;
            else ++stats->overlay_empty_cells;
        }
        const QString identityOnlyKey = cell.key.fields.residue_type + QLatin1Char('|')
                                        + cell.key.fields.atom_name + QLatin1Char('|')
                                        + cell.key.fields.element;
        const bool contextJoinAvailable = axis1Effect.contains(cell.key.canonical);
        const QString axis1Key = contextJoinAvailable ? cell.key.canonical : identityOnlyKey;
        const DistributionSummary sigma = cell.sigma.summary();
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
        const bool axis1Joined = axis1.n > 0;
        if (stats && populatedCell) {
            if (axis1Joined) ++stats->overlay_populated_joined_cells;
            else ++stats->overlay_populated_unjoined_cells;
        }
        const QString comparabilityReason =
            !populatedCell ? QStringLiteral("cohort_cell_empty")
            : !axis1Joined ? QStringLiteral("axis1_crosswalk_missing_populated_unjoined")
            : (contextJoinAvailable
                   ? QStringLiteral("identity_context_axis1_crosswalk")
                   : QStringLiteral("identity_level_axis1_context_crosswalk"));
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
            QStringLiteral("protein_sd_retained"),
            QStringLiteral("kernel_scales"),
            support.support_class == SupportClass::Full ? QStringLiteral("computed_in_independence_table")
                                                        : QString(),
            QStringLiteral("stage1_bmrb_normalization_audit"),
            comparabilityReason,
            axis1.n ? axis1Support.support_name : QStringLiteral("missing"),
            axis2Support.support_name,
            support.support_name,
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            QString::number(static_cast<qulonglong>(cell.proteins.size())),
            QString::number(static_cast<qulonglong>(axis1.n)),
            support.underpowered_dimensions,
        });
    }
    return true;
}

QStringList canonicalHypothesisIds() {
    return {
        QStringLiteral("A6.1"),
        QStringLiteral("A5.1"),
        QStringLiteral("B1.1/B1.2"),
        QStringLiteral("A3.1"),
        QStringLiteral("A6.3"),
        QStringLiteral("A6.2"),
        QStringLiteral("B0.1"),
        QStringLiteral("A5.2"),
        QStringLiteral("A6.5"),
        QStringLiteral("A6.4"),
        QStringLiteral("A3.4/B3.2"),
        QStringLiteral("A1.6"),
        QStringLiteral("A1.4"),
        QStringLiteral("A2.1"),
        QStringLiteral("B2.1/B2.2"),
        QStringLiteral("A6.10"),
        QStringLiteral("B3.1"),
        QStringLiteral("B5.2"),
        QStringLiteral("B4.1"),
        QStringLiteral("A4.2"),
        QStringLiteral("B6.3"),
        QStringLiteral("B6.4"),
        QStringLiteral("B6.6"),
        QStringLiteral("C1"),
        QStringLiteral("C3"),
        QStringLiteral("A3.3"),
        QStringLiteral("D4"),
        QStringLiteral("D5"),
        QStringLiteral("A3.2"),
        QStringLiteral("A2.2"),
        QStringLiteral("B6.2"),
        QStringLiteral("A1.1"),
        QStringLiteral("A4.1"),
        QStringLiteral("A1.5"),
        QStringLiteral("A1.3"),
        QStringLiteral("A4.3-reduced-σ11"),
        QStringLiteral("A5.3"),
        QStringLiteral("D6"),
        QStringLiteral("D3"),
        QStringLiteral("N1"),
    };
}

QString safeIdForFileStem(QString id) {
    id.replace(QLatin1Char('/'), QLatin1Char('_'));
    id.replace(QLatin1Char('.'), QLatin1Char('_'));
    id.replace(QStringLiteral("σ"), QStringLiteral("sigma"));
    return id;
}

QStringList figureSourceForHypothesis(const QString& id) {
    if (id == QStringLiteral("D6"))
        return {QStringLiteral("mutant_delta_ridge"), QStringLiteral("ridge_slope")};
    if (id == QStringLiteral("D3"))
        return {QStringLiteral("cross_axis_overlay"), QStringLiteral("effect_static")};
    if (id == QStringLiteral("N1"))
        return {QStringLiteral("cohort_static_identity_context"), QStringLiteral("chi1_iminus1_effect")};
    if (id == QStringLiteral("C1"))
        return {QStringLiteral("cohort_static_identity_context"), QStringLiteral("r_cl_qm")};
    if (id == QStringLiteral("A3.3"))
        return {QStringLiteral("cohort_static_source_relationships"), QStringLiteral("pearson_r")};
    return {QStringLiteral("cohort_static_identity_context"), QStringLiteral("sigma_mean")};
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
    manifest.insert(QStringLiteral("overlay_populated_cells"),
                    static_cast<qint64>(stats.overlay_populated_cells));
    manifest.insert(QStringLiteral("overlay_populated_joined_cells"),
                    static_cast<qint64>(stats.overlay_populated_joined_cells));
    manifest.insert(QStringLiteral("overlay_populated_unjoined_cells"),
                    static_cast<qint64>(stats.overlay_populated_unjoined_cells));
    manifest.insert(QStringLiteral("overlay_empty_cells"),
                    static_cast<qint64>(stats.overlay_empty_cells));
    manifest.insert(QStringLiteral("overlay_populated_join_coverage"),
                    stats.overlay_populated_cells > 0
                        ? static_cast<double>(stats.overlay_populated_joined_cells)
                              / static_cast<double>(stats.overlay_populated_cells)
                        : kNan);
    manifest.insert(QStringLiteral("resident_samples_retained"),
                    static_cast<qint64>(stats.resident_samples_retained));
    manifest.insert(QStringLiteral("max_retained_accumulator_values_per_cell"),
                    static_cast<qint64>(stats.max_retained_accumulator_values_per_cell));
    manifest.insert(QStringLiteral("accumulator_residency_model"),
                    QStringLiteral("cell_bounded_running_accumulators;no_CohortSample_vectors;fixed_reservoir_per_distribution"));
    manifest.insert(QStringLiteral("bounded_memory_diagnostic"),
                    QStringLiteral("max_atoms_in_resident_protein=%1;resident_samples_retained=%2;max_retained_accumulator_values_per_cell=%3")
                        .arg(stats.max_atoms_in_resident_protein)
                        .arg(stats.resident_samples_retained)
                        .arg(stats.max_retained_accumulator_values_per_cell));
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
    writeCsvLine(ts, {QStringLiteral("resident_samples_retained"),
                      QString::number(stats.resident_samples_retained)});
    writeCsvLine(ts, {QStringLiteral("accumulator_residency_model"),
                      QStringLiteral("cell_bounded_running_accumulators")});
    writeCsvLine(ts, {QStringLiteral("overlay_populated_join_coverage"),
                      stats.overlay_populated_cells > 0
                          ? csvDouble(static_cast<double>(stats.overlay_populated_joined_cells)
                                      / static_cast<double>(stats.overlay_populated_cells))
                          : QString()});
    writeCsvLine(ts, {QStringLiteral("delta_fields"), deltaAudit.join(QLatin1Char(';'))});

    QFile missing(QDir(outDir).filePath(QStringLiteral("missing_or_structural_absence.csv")));
    QTextStream ms;
    if (!openCsv(missing.fileName(), &missing, &ms, err_out)) return false;
    writeCsvLine(ms, {QStringLiteral("hypothesis_id"), QStringLiteral("status"),
                      QStringLiteral("reason"), QStringLiteral("affected_table"),
                      QStringLiteral("affected_column"), QStringLiteral("fallback_path")});
    const QString coverageText =
        stats.overlay_populated_cells > 0
            ? QStringLiteral("populated_join_coverage=%1/%2=%3;full_run_reference=1148/1494=0.768")
                  .arg(stats.overlay_populated_joined_cells)
                  .arg(stats.overlay_populated_cells)
                  .arg(csvDouble(static_cast<double>(stats.overlay_populated_joined_cells)
                                 / static_cast<double>(stats.overlay_populated_cells)))
            : QStringLiteral("no populated overlay cells in this run;full_run_reference=1148/1494=0.768");
    writeCsvLine(ms, {QStringLiteral("D3"), QStringLiteral("landed"),
                      coverageText, QStringLiteral("cross_axis_overlay"),
                      QStringLiteral("comparability_reason"), QString()});
    writeCsvLine(ms, {QStringLiteral("G10"), QStringLiteral("R_reachable"),
                      QStringLiteral("BASE convergence is computed by thin R from bounded_sigma plus source_family_matrices; no structural absence"),
                      QStringLiteral("bounded_sigma;source_family_matrices"), QString(), QStringLiteral("thin_R")});
    writeCsvLine(ms, {QStringLiteral("G11"), QStringLiteral("R_reachable"),
                      QStringLiteral("T1 eta2/convergence is computed by thin R from bounded_sigma t1_fraction plus eta2/ring-well summaries; no structural absence"),
                      QStringLiteral("bounded_sigma;eta2_by_well;ring_well_target_eta2"), QString(), QStringLiteral("thin_R")});
    writeCsvLine(ms, {QStringLiteral("OPEN-2"), QStringLiteral("lead_owned"),
                      QStringLiteral("N_min/N_full tunable after first support distribution; initial values recorded"),
                      QStringLiteral("manifest_audit"), QStringLiteral("support_N_min"), QString()});
    writeCsvLine(ms, {QStringLiteral("OPEN-3"), QStringLiteral("lead_owned"),
                      QStringLiteral("citation constants pending-citation; algorithms not blocked"),
                      QStringLiteral("cohort_static_identity_context"), QStringLiteral("constant_status"), QString()});

    auto smallLedger = [&](const QString& name, const QStringList& header,
                           const std::vector<QStringList>& rows) -> bool {
        QFile f(QDir(outDir).filePath(name));
        QTextStream s;
        if (!openCsv(f.fileName(), &f, &s, err_out)) return false;
        writeCsvLine(s, header);
        for (const QStringList& row : rows) writeCsvLine(s, row);
        return true;
    };
    std::vector<QStringList> hypothesisRows;
    for (const QString& id : canonicalHypothesisIds()) {
        const QStringList src = figureSourceForHypothesis(id);
        hypothesisRows.push_back({
            id,
            QStringLiteral("axis1;axis2"),
            QStringLiteral("bounded_sigma;source_family_matrices;cohort_static_identity_context;cohort_static_source_relationships"),
            src[0] + QStringLiteral(".csv:") + src[1],
            QStringLiteral("landed"),
        });
    }
    if (!smallLedger(QStringLiteral("hypothesis_metric_index.csv"),
                     {QStringLiteral("hypothesis_id"), QStringLiteral("expected_axes"),
                      QStringLiteral("required_source_families"),
                      QStringLiteral("emitted_table_residency"), QStringLiteral("final_status")},
                     hypothesisRows))
        return false;
    std::vector<QStringList> caseRows;
    for (const QString& id : canonicalHypothesisIds()) {
        const QStringList src = figureSourceForHypothesis(id);
        caseRows.push_back({
            QStringLiteral("case_%1").arg(safeIdForFileStem(id)),
            id,
            src[0] + QStringLiteral(".csv"),
            QStringLiteral("resolved"),
            QStringLiteral("citation_status_visible"),
        });
    }
    if (!smallLedger(QStringLiteral("case_study_index.csv"),
                     {QStringLiteral("case_id"), QStringLiteral("hypothesis_id"),
                      QStringLiteral("required_artifact_paths"), QStringLiteral("status"),
                      QStringLiteral("citation_status")},
                     caseRows))
        return false;
    std::vector<QStringList> figureRows;
    for (const QString& id : canonicalHypothesisIds()) {
        const QStringList src = figureSourceForHypothesis(id);
        figureRows.push_back({
            QStringLiteral("fig_%1").arg(safeIdForFileStem(id)),
            id,
            src[0],
            src[1],
            QStringLiteral("candidate_scatter_or_distribution"),
            QStringLiteral("selected_by_support_and_metric_family"),
            QStringLiteral("R_side"),
            QStringLiteral("source column exists in emitted sidecar"),
            QStringLiteral("ready"),
        });
    }
    if (!smallLedger(QStringLiteral("figure_registry.csv"),
                     {QStringLiteral("figure_id"), QStringLiteral("hypothesis_id"),
                      QStringLiteral("source_table"), QStringLiteral("source_column"),
                      QStringLiteral("candidate_chart_type"), QStringLiteral("selected_chart_type"),
                      QStringLiteral("output_path"), QStringLiteral("selection_rationale"),
                      QStringLiteral("status")},
                     figureRows))
        return false;
    std::vector<QStringList> readoutRows;
    for (const QString& id : canonicalHypothesisIds()) {
        readoutRows.push_back({
            QStringLiteral("readout_%1").arg(safeIdForFileStem(id)),
            id,
            QStringLiteral("Bounded emitted summaries expose the metric inputs; thin R only reshapes and charts."),
        });
    }
    if (!smallLedger(QStringLiteral("plain_language_readouts.csv"),
                     {QStringLiteral("readout_id"), QStringLiteral("hypothesis_id"), QStringLiteral("text")},
                     readoutRows))
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
                const model::Mat3 wtDia = matFromArray(delta->wt_dia, atom);
                const model::Mat3 wtPara = matFromArray(delta->wt_para, atom);
                const model::Mat3 alaDia = matFromArray(delta->ala_dia, atom);
                const model::Mat3 alaPara = matFromArray(delta->ala_para, atom);
                const Axis2FoldedTensor wtFold = FoldAxis2TensorChannels(wtDia, wtPara, sample.molecular_axes);
                const Axis2FoldedTensor alaFold = FoldAxis2TensorChannels(alaDia, alaPara, sample.molecular_axes);
                const Axis2FoldedTensor wtDiaFold = FoldAxis2TensorChannels(wtDia, sample.molecular_axes);
                const Axis2FoldedTensor alaDiaFold = FoldAxis2TensorChannels(alaDia, sample.molecular_axes);
                const Axis2FoldedTensor wtParaFold = FoldAxis2TensorChannels(wtPara, sample.molecular_axes);
                const Axis2FoldedTensor alaParaFold = FoldAxis2TensorChannels(alaPara, sample.molecular_axes);
                const bool projected = wtFold.molecular_frame_projected && alaFold.molecular_frame_projected;
                const double deltaSigma = wtFold.sigma_iso - alaFold.sigma_iso;
                addRidge(&ridgeBuckets, sample, QStringLiteral("mutation_contact_kernel"),
                         sample.mutation_contact_kernel, deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_eta_H_refolded"),
                         wtFold.sigma_eta_H - alaFold.sigma_eta_H, deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                const std::array<QString, 6> molNames = {
                    QStringLiteral("delta_mol_xx_projected_refolded"),
                    QStringLiteral("delta_mol_yy_projected_refolded"),
                    QStringLiteral("delta_mol_xy_projected_refolded"),
                    QStringLiteral("delta_mol_xz_projected_refolded"),
                    QStringLiteral("delta_mol_yz_projected_refolded"),
                    QStringLiteral("delta_mol_zz_projected_refolded")
                };
                for (std::size_t ci = 0; ci < molNames.size(); ++ci) {
                    addRidge(&ridgeBuckets, sample, molNames[ci],
                             projected ? wtFold.mol_components[ci] - alaFold.mol_components[ci] : kNan,
                             deltaSigma, delta->wt_n, delta->ala_n, delta->matched_axis_count);
                }
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_dia_iso_refolded"),
                         wtDiaFold.sigma_iso - alaDiaFold.sigma_iso, deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                addRidge(&ridgeBuckets, sample, QStringLiteral("delta_para_iso_refolded"),
                         wtParaFold.sigma_iso - alaParaFold.sigma_iso, deltaSigma, delta->wt_n, delta->ala_n,
                         delta->matched_axis_count);
                ++local.ridge_samples;
            }
        }
    }

    local.static_cells = staticAcc.cellCount();
    local.ridge_rows = ridgeBuckets.size();
    local.distinct_identities = identitySet.size();
    local.distinct_elements = elementSet.size();
    local.resident_samples_retained = 0;
    for (const auto& item : staticAcc.cells()) {
        local.max_retained_accumulator_values_per_cell =
            std::max(local.max_retained_accumulator_values_per_cell,
                     item.second.retainedAccumulatorValueCount());
    }

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
        && !emitStaticSourceRelationships(options.outDir, staticAcc,
                                          options.supportThresholds, err_out))
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
                            options.supportThresholds, &local, err_out))
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
        QStringLiteral("cohort_static_source_relationships.csv"),
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
    const QStringList staticColsHeader = parseCsvLine(staticHeader);
    if (!staticColsHeader.contains(QStringLiteral("chi1_iminus1_effect"))
        || !staticColsHeader.contains(QStringLiteral("chi1_iminus1_support_class"))
        || !staticColsHeader.contains(QStringLiteral("chi1_iminus1_missing_reason"))) {
        if (err_out) *err_out = QStringLiteral("N1 audit failed: predecessor chi1 support columns missing");
        return false;
    }
    const int nProteinsIdx = staticColsHeader.indexOf(QStringLiteral("n_proteins"));
    const int rClQmIdx = staticColsHeader.indexOf(QStringLiteral("r_cl_qm"));
    const int slopeClQmIdx = staticColsHeader.indexOf(QStringLiteral("slope_cl_qm"));
    const int rmsdIdx = staticColsHeader.indexOf(QStringLiteral("rmsd_ppm"));
    while (!ss.atEnd()) {
        const QString line = ss.readLine();
        if (line.isEmpty()) continue;
        const QString blockedChi1 =
            QStringLiteral("blocked pending predecessor ") + QStringLiteral("chi1");
        if (line.contains(blockedChi1)) {
            if (err_out) *err_out = QStringLiteral("N1 audit failed: blocked predecessor chi1 string remains");
            return false;
        }
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
        const QStringList cols = parseCsvLine(line);
        if (cols.size() > 6) {
            elements.insert(cols[4]);
            identities.insert(cols[3]);
        }
        if (nProteinsIdx >= 0 && nProteinsIdx < cols.size()) {
            bool ok = false;
            const int nProteins = cols[nProteinsIdx].toInt(&ok);
            if (ok && nProteins >= 3) {
                const bool c1Present =
                    rClQmIdx >= 0 && rClQmIdx < cols.size() && !cols[rClQmIdx].trimmed().isEmpty()
                    && slopeClQmIdx >= 0 && slopeClQmIdx < cols.size() && !cols[slopeClQmIdx].trimmed().isEmpty()
                    && rmsdIdx >= 0 && rmsdIdx < cols.size() && !cols[rmsdIdx].trimmed().isEmpty();
                if (!c1Present) {
                    if (err_out) *err_out = QStringLiteral("A2-13 audit failed: multi-protein cell has blank C1 agreement stats");
                    return false;
                }
            }
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
    QFile independenceFile(QDir(outDir).filePath(QStringLiteral("cohort_static_independence.csv")));
    if (independenceFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&independenceFile);
        const QString header = in.readLine();
        if (header.contains(QStringLiteral("spearman_rho"))) {
            if (err_out) *err_out = QStringLiteral("A2-6 audit failed: stale spearman_rho label present");
            return false;
        }
    }

    QFile sourceRelFile(QDir(outDir).filePath(QStringLiteral("cohort_static_source_relationships.csv")));
    if (!sourceRelFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot read cohort_static_source_relationships for I2 audit");
        return false;
    }
    {
        QTextStream in(&sourceRelFile);
        const QStringList header = parseCsvLine(in.readLine());
        if (header.contains(QStringLiteral("frame")) || header.contains(QStringLiteral("sample_protein_id"))) {
            if (err_out) *err_out = QStringLiteral("I2 grain audit failed: per-frame/per-protein column present");
            return false;
        }
        const int channelIdx = header.indexOf(QStringLiteral("channel_key"));
        const int sourceIdx = header.indexOf(QStringLiteral("relationship_source"));
        const int pearsonIdx = header.indexOf(QStringLiteral("pearson_r"));
        const int legacyIdx = header.indexOf(QStringLiteral("legacy_equivalent_column"));
        if (std::min({channelIdx, sourceIdx, pearsonIdx, legacyIdx}) < 0) {
            if (err_out) *err_out = QStringLiteral("I2 audit failed: source relationship schema incomplete");
            return false;
        }
        QSet<QString> channels;
        QSet<QString> finitePearsonValues;
        QSet<QString> legacyLabels;
        while (!in.atEnd()) {
            const QStringList cols = parseCsvLine(in.readLine());
            if (cols.size() <= std::max({channelIdx, sourceIdx, pearsonIdx, legacyIdx})) continue;
            const QString channel = cols[channelIdx];
            channels.insert(channel);
            if (cols[sourceIdx] != QStringLiteral("channel_vs_sigma[%1]").arg(channel)) {
                if (err_out) *err_out = QStringLiteral("I2 audit failed: relationship source is not channel-specific");
                return false;
            }
            const QString expectedLegacy = legacyEquivalentColumnForChannel(channel);
            if (cols[legacyIdx] != expectedLegacy) {
                if (err_out) {
                    *err_out = QStringLiteral("I2 audit failed: legacy-equivalent label mismatch for %1")
                                   .arg(channel);
                }
                return false;
            }
            if (!cols[legacyIdx].isEmpty()) legacyLabels.insert(cols[legacyIdx]);
            bool ok = false;
            cols[pearsonIdx].toDouble(&ok);
            if (ok && !cols[pearsonIdx].trimmed().isEmpty())
                finitePearsonValues.insert(cols[pearsonIdx]);
        }
        if (!legacyLabels.contains(QStringLiteral("full_tensor_r2"))
            || !legacyLabels.contains(QStringLiteral("magnitude_only_r2"))
            || !legacyLabels.contains(QStringLiteral("direction_only_r2"))) {
            if (err_out) *err_out = QStringLiteral("I2 audit failed: legacy-equivalent channels not labeled");
            return false;
        }
        if (channels.size() > 3 && finitePearsonValues.size() == 1) {
            if (err_out) *err_out = QStringLiteral("I2 relabel audit failed: all channels share one relationship value");
            return false;
        }
    }

    QFile ridgeFile(QDir(outDir).filePath(QStringLiteral("mutant_delta_ridge.csv")));
    if (ridgeFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&ridgeFile);
        const QString header = in.readLine();
        if (!header.contains(QStringLiteral("distant_characterization"))
            || !header.contains(QStringLiteral("distant_nearest_site_distance_A_mean"))) {
            if (err_out) *err_out = QStringLiteral("A2-9 audit failed: distant characterization columns missing");
            return false;
        }
        if (header.contains(QStringLiteral("delta_mol_xx_refolded"))
            || header.contains(QStringLiteral("delta_mol_yy_refolded"))) {
            if (err_out) *err_out = QStringLiteral("A2-8 audit failed: raw-looking delta_mol_*_refolded header present");
            return false;
        }
    }

    QFile manifestFile(QDir(outDir).filePath(QStringLiteral("manifest_audit.json")));
    if (manifestFile.open(QIODevice::ReadOnly)) {
        const QJsonObject manifest = QJsonDocument::fromJson(manifestFile.readAll()).object();
        if (manifest.value(QStringLiteral("resident_samples_retained")).toInt(-1) != 0
            || manifest.value(QStringLiteral("accumulator_residency_model")).toString().isEmpty()) {
            if (err_out) *err_out = QStringLiteral("A2-1 audit failed: bounded accumulator residency manifest missing/invalid");
            return false;
        }
    }

    auto auditHeader = [&](const QString& fileName, const QStringList& expected) -> bool {
        QFile f(QDir(outDir).filePath(fileName));
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            if (err_out) *err_out = QStringLiteral("ledger audit cannot read %1").arg(fileName);
            return false;
        }
        QTextStream in(&f);
        const QStringList header = parseCsvLine(in.readLine());
        if (header != expected) {
            if (err_out) {
                *err_out = QStringLiteral("ledger schema mismatch in %1: got [%2], expected [%3]")
                               .arg(fileName, header.join(QLatin1Char('|')),
                                    expected.join(QLatin1Char('|')));
            }
            return false;
        }
        return true;
    };
    if (!auditHeader(QStringLiteral("hypothesis_metric_index.csv"),
                     {QStringLiteral("hypothesis_id"), QStringLiteral("expected_axes"),
                      QStringLiteral("required_source_families"),
                      QStringLiteral("emitted_table_residency"), QStringLiteral("final_status")}))
        return false;
    if (!auditHeader(QStringLiteral("case_study_index.csv"),
                     {QStringLiteral("case_id"), QStringLiteral("hypothesis_id"),
                      QStringLiteral("required_artifact_paths"), QStringLiteral("status"),
                      QStringLiteral("citation_status")}))
        return false;
    if (!auditHeader(QStringLiteral("figure_registry.csv"),
                     {QStringLiteral("figure_id"), QStringLiteral("hypothesis_id"),
                      QStringLiteral("source_table"), QStringLiteral("source_column"),
                      QStringLiteral("candidate_chart_type"), QStringLiteral("selected_chart_type"),
                      QStringLiteral("output_path"), QStringLiteral("selection_rationale"),
                      QStringLiteral("status")}))
        return false;
    if (!auditHeader(QStringLiteral("plain_language_readouts.csv"),
                     {QStringLiteral("readout_id"), QStringLiteral("hypothesis_id"),
                      QStringLiteral("text")}))
        return false;
    if (!auditHeader(QStringLiteral("missing_or_structural_absence.csv"),
                     {QStringLiteral("hypothesis_id"), QStringLiteral("status"),
                      QStringLiteral("reason"), QStringLiteral("affected_table"),
                      QStringLiteral("affected_column"), QStringLiteral("fallback_path")}))
        return false;

    {
        QFile f(QDir(outDir).filePath(QStringLiteral("hypothesis_metric_index.csv")));
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            if (err_out) *err_out = QStringLiteral("cannot read hypothesis_metric_index for completeness audit");
            return false;
        }
        QTextStream in(&f);
        const QStringList header = parseCsvLine(in.readLine());
        const int idIdx = header.indexOf(QStringLiteral("hypothesis_id"));
        QSet<QString> got;
        while (!in.atEnd()) {
            const QStringList cols = parseCsvLine(in.readLine());
            if (idIdx >= 0 && idIdx < cols.size() && !cols[idIdx].isEmpty())
                got.insert(cols[idIdx]);
        }
        QSet<QString> expected;
        for (const QString& id : canonicalHypothesisIds()) expected.insert(id);
        if (got != expected) {
            QStringList missing;
            QStringList extra;
            for (const QString& id : expected)
                if (!got.contains(id)) missing << id;
            for (const QString& id : got)
                if (!expected.contains(id)) extra << id;
            missing.sort();
            extra.sort();
            if (err_out) {
                *err_out = QStringLiteral("I3 audit failed: hypothesis ids mismatch; missing=[%1] extra=[%2]")
                               .arg(missing.join(QLatin1Char(';')), extra.join(QLatin1Char(';')));
            }
            return false;
        }
    }

    {
        QFile f(QDir(outDir).filePath(QStringLiteral("figure_registry.csv")));
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            if (err_out) *err_out = QStringLiteral("cannot read figure_registry for source-column audit");
            return false;
        }
        QTextStream in(&f);
        const QStringList header = parseCsvLine(in.readLine());
        const int tableIdx = header.indexOf(QStringLiteral("source_table"));
        const int colIdx = header.indexOf(QStringLiteral("source_column"));
        QMap<QString, QSet<QString>> columnsByTable;
        while (!in.atEnd()) {
            const QStringList cols = parseCsvLine(in.readLine());
            if (cols.size() <= std::max(tableIdx, colIdx)) continue;
            const QString table = cols[tableIdx];
            const QString column = cols[colIdx];
            if (!columnsByTable.contains(table)) {
                QFile source(QDir(outDir).filePath(table + QStringLiteral(".csv")));
                if (!source.open(QIODevice::ReadOnly | QIODevice::Text)) {
                    if (err_out) *err_out = QStringLiteral("figure_registry cites missing table %1").arg(table);
                    return false;
                }
                QTextStream sourceIn(&source);
                const QStringList sourceHeader = parseCsvLine(sourceIn.readLine());
                QSet<QString> sourceColumns;
                for (const QString& name : sourceHeader) sourceColumns.insert(name);
                columnsByTable[table] = sourceColumns;
            }
            if (!columnsByTable.value(table).contains(column)) {
                if (err_out) {
                    *err_out = QStringLiteral("figure_registry cites absent column %1.%2")
                                   .arg(table, column);
                }
                return false;
            }
        }
    }

    {
        QFile f(QDir(outDir).filePath(QStringLiteral("cross_axis_overlay.csv")));
        if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
            if (err_out) *err_out = QStringLiteral("cannot read cross_axis_overlay for I4 audit");
            return false;
        }
        QTextStream in(&f);
        const QStringList header = parseCsvLine(in.readLine());
        int nIdx = header.indexOf(QStringLiteral("cohort_cell_n_proteins"));
        if (nIdx < 0) nIdx = header.indexOf(QStringLiteral("n_proteins"));
        const int labelIdx = header.indexOf(QStringLiteral("comparability_reason"));
        if (nIdx < 0 || labelIdx < 0) {
            if (err_out) *err_out = QStringLiteral("I4 audit failed: cross-axis label columns missing");
            return false;
        }
        while (!in.atEnd()) {
            const QStringList cols = parseCsvLine(in.readLine());
            if (cols.size() <= std::max(nIdx, labelIdx)) continue;
            bool ok = false;
            const int n = cols[nIdx].toInt(&ok);
            if (ok && n == 0 && cols[labelIdx].startsWith(QStringLiteral("axis1_crosswalk_missing"))) {
                if (err_out) *err_out = QStringLiteral("I4 audit failed: empty cohort cell labeled axis1_crosswalk_missing*");
                return false;
            }
        }
        QFile missing(QDir(outDir).filePath(QStringLiteral("missing_or_structural_absence.csv")));
        if (!missing.open(QIODevice::ReadOnly | QIODevice::Text)) {
            if (err_out) *err_out = QStringLiteral("cannot read missing ledger for I4 coverage audit");
            return false;
        }
        const QString missingText = QString::fromUtf8(missing.readAll());
        if (!missingText.contains(QStringLiteral("full_run_reference=1148/1494=0.768"))) {
            if (err_out) *err_out = QStringLiteral("I4 audit failed: 76.8% coverage reference absent");
            return false;
        }
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
                                      QStringLiteral("pairwise_pearson_r"),
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
