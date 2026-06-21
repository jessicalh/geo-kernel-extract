#include "NewmanProjection.h"

#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"   // DihedralDegrees
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtResidueNames.h"          // IupacResidue3LetterFor
#include "../model/Types.h"                   // Vec3 (Eigen::Vector3d)

#include <Eigen/Dense>

#include <cmath>

namespace h5reader::app {

using model::Vec3;

namespace {

// Angle of a substituent atom (bonded to `hub`) in the projection basis (u, v),
// i.e. its direction projected onto the plane perpendicular to the sight axis.
double SpokeAngle(const Vec3& hub, const Vec3& sub,
                  const Vec3& sight, const Vec3& u, const Vec3& v) {
    Vec3 d = sub - hub;
    d = d - d.dot(sight) * sight;
    return std::atan2(d.dot(v), d.dot(u));
}

constexpr int NONE = model::QtResidue::NONE;

}  // namespace

NewmanProjection ComputeNewmanProjection(const model::QtProtein& protein,
                                         const model::Conformation& conf,
                                         std::size_t residueIndex, std::size_t frame,
                                         NewmanKind kind) {
    NewmanProjection out;
    out.kind = kind;
    out.residueIndex = residueIndex;

    if (residueIndex >= protein.residueCount()) {
        out.invalidReason = QStringLiteral("residue index out of range");
        return out;
    }
    const model::QtResidue& res = protein.residue(residueIndex);
    // Match the Inspector's residue identity ("LEU #28"): IUPAC 3-letter +
    // the residue NUMBER from the chain address, not the 0-based index.
    out.residueLabel = QStringLiteral("%1 #%2")
        .arg(QString::fromLatin1(model::IupacResidue3LetterFor(res.aminoAcid)))
        .arg(res.address.residueNumber);

    // Resolve the four dihedral atoms a-b-c-d (sight down b->c) and the
    // substituent spokes of the front (b) and back (c) atoms.
    int a = NONE, b = NONE, c = NONE, d = NONE;
    struct Sub { int idx; const char* label; bool front; bool reference; };
    std::vector<Sub> subs;

    if (!res.HasN() || !res.HasCA() || !res.HasC()) {
        out.invalidReason = QStringLiteral("missing backbone N/CA/C");
        return out;
    }

    if (kind == NewmanKind::Phi) {
        const int prev = res.prevResidueIndex;
        if (prev < 0 || static_cast<std::size_t>(prev) >= protein.residueCount()
            || !protein.residue(static_cast<std::size_t>(prev)).HasC()) {
            out.invalidReason = QStringLiteral("N-terminus / no preceding carbonyl C");
            return out;
        }
        a = protein.residue(static_cast<std::size_t>(prev)).C;
        b = res.N; c = res.CA; d = res.C;
        out.frontLabel = QStringLiteral("N");
        out.backLabel  = QStringLiteral("CA");
        // front (N): C(-1) [ref], H
        subs.push_back({a, "C(-1)", true, true});
        if (res.H != NONE)  subs.push_back({res.H, "H", true, false});
        // back (CA): C [ref], CB, HA
        subs.push_back({d, "C", false, true});
        if (res.CB != NONE) subs.push_back({res.CB, "CB", false, false});
        if (res.HA != NONE) subs.push_back({res.HA, "HA", false, false});
    } else {
        const int next = res.nextResidueIndex;
        if (next < 0 || static_cast<std::size_t>(next) >= protein.residueCount()
            || !protein.residue(static_cast<std::size_t>(next)).HasN()) {
            out.invalidReason = QStringLiteral("C-terminus / no following amide N");
            return out;
        }
        a = res.N; b = res.CA; c = res.C;
        d = protein.residue(static_cast<std::size_t>(next)).N;
        out.frontLabel = QStringLiteral("CA");
        out.backLabel  = QStringLiteral("C");
        // front (CA): N [ref], CB, HA
        subs.push_back({a, "N", true, true});
        if (res.CB != NONE) subs.push_back({res.CB, "CB", true, false});
        if (res.HA != NONE) subs.push_back({res.HA, "HA", true, false});
        // back (C): N(+1) [ref], O
        subs.push_back({d, "N(+1)", false, true});
        if (res.O != NONE)  subs.push_back({res.O, "O", false, false});
    }

    auto pos = [&](int idx) {
        return conf.atomPosition(frame, static_cast<std::size_t>(idx));
    };
    const Vec3 pa = pos(a), pb = pos(b), pc = pos(c), pd = pos(d);

    out.torsionDeg = model::DihedralDegrees(pa, pb, pc, pd);
    if (std::isnan(out.torsionDeg)) {
        out.invalidReason = QStringLiteral("degenerate geometry");
        return out;
    }

    // Sight axis b->c; anchor u to the front reference bond so the front-ref
    // spoke sits at angle 0 and the back-ref spoke reads ~torsion.
    Vec3 sight = pc - pb;
    if (sight.norm() < 1e-9) {
        out.invalidReason = QStringLiteral("zero sight axis");
        return out;
    }
    sight.normalize();
    Vec3 refDir = pa - pb;
    refDir = refDir - refDir.dot(sight) * sight;
    Vec3 u;
    if (refDir.norm() > 1e-9) {
        u = refDir.normalized();
    } else {
        const Vec3 t = (std::abs(sight.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        u = (t - t.dot(sight) * sight).normalized();
    }
    // Orient v so the back-reference spoke reads the signed torsion angle
    // directly (matches the Blondel-Karplus convention DihedralDegrees uses),
    // giving the panel a clean invariant: back-ref spoke angle == torsionDeg.
    const Vec3 v = u.cross(sight);

    out.spokes.reserve(subs.size());
    for (const Sub& s : subs) {
        NewmanSpoke spoke;
        spoke.angleRad  = SpokeAngle(s.front ? pb : pc, pos(s.idx), sight, u, v);
        spoke.label     = QString::fromLatin1(s.label);
        spoke.front     = s.front;
        spoke.reference = s.reference;
        out.spokes.push_back(spoke);
    }

    out.valid = true;
    return out;
}

}  // namespace h5reader::app
