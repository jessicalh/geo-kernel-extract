// csa_math_tests -- headless QtTest for the viewer-side CSA / molecular-frame
// math (model/CsaShape.h + model/MolecularFrame.h). Pure functions over
// model::Vec3/Mat3; links Qt6::Test + Eigen3 only, like scene_math_tests.
//
// Verifies the decomposition/projection grabbed from AnalysisAtom actually
// compiles on MSVC and is numerically correct: a known tensor -> known PAS /
// Haeberlen shape; the molecular-frame builders -> orthonormal right-handed
// axes; symmetricComponents in the engine's {xx,yy,xy,xz,yz,zz} order.

#include "model/CsaShape.h"
#include "model/MolecularFrame.h"

#include <QObject>
#include <QtTest>

#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

using namespace h5reader;

namespace {
bool nearly(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }

bool orthonormalRH(const model::Mat3& m) {
    const model::Mat3 g = m.transpose() * m;
    return (g - model::Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9
           && std::abs(m.determinant() - 1.0) < 1e-9;
}

// Build per-frame axes from a fixed position table; the single ring is centred
// at the origin with a +z normal (for the aromatic-ring-local kind).
std::optional<model::Mat3> buildAxes(const model::MolFrameSpec& spec,
                                     const std::vector<model::Vec3>& pos) {
    const model::PositionLookup posOf = [&](std::int32_t i) -> std::optional<model::Vec3> {
        if (i < 0 || static_cast<std::size_t>(i) >= pos.size()) return std::nullopt;
        return pos[static_cast<std::size_t>(i)];
    };
    const model::RingFrameLookup ringOf =
        [](std::int32_t) -> std::optional<std::pair<model::Vec3, model::Vec3>> {
        return std::make_pair(model::Vec3(0, 0, 0), model::Vec3(0, 0, 1));
    };
    model::MolFrameContinuity cont;
    return model::BuildMolecularFrameAxes(spec, posOf, ringOf, cont);
}
}  // namespace

class CsaMathTests : public QObject {
    Q_OBJECT

private slots:
    void csaDiagonalKnownShape();
    void csaIsotropicIsInvalid();
    void symmetricComponentsOrder();
    void projectIdentityAndRotationTraceInvariant();
    void frameFromXAndPlaneIsOrthonormal();
    void buildCarbonylFrameOrthonormal();
    // Every molecular-frame kind builds an orthonormal right-handed frame.
    void buildAmideNFrame();
    void buildAromaticRingFrame();
    void buildMetSdFrame();
    void buildCarboxylateFrame();
    void buildGuanidiniumFrame();
    void buildCarboxamideFrame();
};

void CsaMathTests::csaDiagonalKnownShape() {
    // diag(10, 20, 60) ppm. iso=30, span=50; zz=60 (farthest), xx=10, yy=20;
    // eta=(20-10)/(60-30)=1/3; skew=3(20-30)/50=-0.6.
    model::Mat3 t = model::Mat3::Zero();
    t(0, 0) = 10.0;
    t(1, 1) = 20.0;
    t(2, 2) = 60.0;
    const model::CsaShape s = model::ComputeCsaShape(t);
    QVERIFY(s.valid);
    QVERIFY(nearly(s.principal_values[0], 10.0));
    QVERIFY(nearly(s.principal_values[1], 20.0));
    QVERIFY(nearly(s.principal_values[2], 60.0));
    QVERIFY(nearly(s.sigma_iso, 30.0));
    QVERIFY(nearly(s.span, 50.0));
    QVERIFY(nearly(s.eta, 1.0 / 3.0));
    QVERIFY(nearly(s.skew, -0.6));
    QVERIFY(nearly(s.haeberlen_values[0], 10.0));
    QVERIFY(nearly(s.haeberlen_values[1], 20.0));
    QVERIFY(nearly(s.haeberlen_values[2], 60.0));
    const model::Mat3 gram = s.pas_axes.transpose() * s.pas_axes;
    QVERIFY((gram - model::Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9);
}

void CsaMathTests::csaIsotropicIsInvalid() {
    const model::Mat3 t = 30.0 * model::Mat3::Identity();
    QVERIFY(!model::ComputeCsaShape(t).valid);
}

void CsaMathTests::symmetricComponentsOrder() {
    model::Mat3 m;
    m << 1, 2, 3,
         2, 4, 5,
         3, 5, 6;
    const std::array<double, 6> c = model::symmetricComponents(m);
    QVERIFY(nearly(c[0], 1.0));  // xx
    QVERIFY(nearly(c[1], 4.0));  // yy
    QVERIFY(nearly(c[2], 2.0));  // xy
    QVERIFY(nearly(c[3], 3.0));  // xz
    QVERIFY(nearly(c[4], 5.0));  // yz
    QVERIFY(nearly(c[5], 6.0));  // zz
    QVERIFY(nearly(model::component(m, model::MolComp::xz), 3.0));
}

void CsaMathTests::projectIdentityAndRotationTraceInvariant() {
    model::Mat3 t;
    t << 1, 2, 3,
         2, 4, 5,
         3, 5, 6;
    const model::Mat3 p = model::ProjectToMolecularFrame(model::Mat3::Identity(), t);
    QVERIFY((p - t).cwiseAbs().maxCoeff() < 1e-12);
    // +90 deg rotation about z: a similarity transform preserves the trace.
    model::Mat3 rz;
    rz << 0, -1, 0,
          1,  0, 0,
          0,  0, 1;
    const model::Mat3 pr = model::ProjectToMolecularFrame(rz, t);
    QVERIFY(nearly(pr.trace(), t.trace()));
}

void CsaMathTests::frameFromXAndPlaneIsOrthonormal() {
    const auto f = model::frameFromXAndPlane(model::Vec3(1, 0, 0), model::Vec3(0, 1, 0));
    QVERIFY(f.has_value());
    const model::Mat3 gram = f->transpose() * (*f);
    QVERIFY((gram - model::Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9);
    QVERIFY(nearly(f->determinant(), 1.0));
    QVERIFY(nearly(f->col(0).x(), 1.0));
}

void CsaMathTests::buildCarbonylFrameOrthonormal() {
    // C at origin, O along +x, next-N along +y -> a peptide-plane frame.
    model::MolFrameSpec spec;
    spec.kind = model::MolecularFrameKind::BackboneCarbonyl;
    spec.origin = 0;
    spec.xAnchor = 1;
    spec.planeAnchor = 2;
    const model::PositionLookup posOf = [](std::int32_t i) -> std::optional<model::Vec3> {
        switch (i) {
        case 0: return model::Vec3(0, 0, 0);
        case 1: return model::Vec3(1, 0, 0);
        case 2: return model::Vec3(0, 1, 0);
        }
        return std::nullopt;
    };
    const model::RingFrameLookup ringOf =
        [](std::int32_t) -> std::optional<std::pair<model::Vec3, model::Vec3>> {
        return std::nullopt;
    };
    model::MolFrameContinuity cont;
    const auto axes = model::BuildMolecularFrameAxes(spec, posOf, ringOf, cont);
    QVERIFY(axes.has_value());
    const model::Mat3 gram = axes->transpose() * (*axes);
    QVERIFY((gram - model::Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9);
    QVERIFY(nearly(axes->determinant(), 1.0));
    QVERIFY(nearly(axes->col(0).x(), 1.0));  // x director along O - C = +x
}

void CsaMathTests::buildAmideNFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::BackboneAmideN;
    s.origin = 0; s.xAnchor = 1; s.planeAnchor = 2;
    const auto a = buildAxes(s, {model::Vec3(0, 0, 0), model::Vec3(1, 0, 0), model::Vec3(0, 1, 0)});
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x along N -> H
}

void CsaMathTests::buildAromaticRingFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::AromaticRingLocal;
    s.ring = 0; s.heavy = 0;
    const auto a = buildAxes(s, {model::Vec3(1, 0, 0)});  // heavy atom; ring at origin, +z normal
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(2).z(), 1.0));  // z along the ring normal
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x radial to the heavy atom
}

void CsaMathTests::buildMetSdFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::MetSd;
    s.origin = 0; s.xAnchor = 1; s.planeAnchor = 2;
    const auto a = buildAxes(s, {model::Vec3(0, 0, 0), model::Vec3(1, 0, 0), model::Vec3(0, 1, 0)});
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x along SD -> CE
}

void CsaMathTests::buildCarboxylateFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::SidechainCarboxylate;
    s.origin = 0; s.xAnchor = 1; s.secondAnchor = 2;
    // C at origin; O1 (1,1,0), O2 (1,-1,0): the O-C-O bisector is +x.
    const auto a = buildAxes(s, {model::Vec3(0, 0, 0), model::Vec3(1, 1, 0), model::Vec3(1, -1, 0)});
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x along the carboxylate bisector
}

void CsaMathTests::buildGuanidiniumFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::SidechainGuanidinium;
    s.origin = 0; s.xAnchor = 1; s.planeAnchor = 2;
    const auto a = buildAxes(s, {model::Vec3(0, 0, 0), model::Vec3(1, 0, 0), model::Vec3(0, 1, 0)});
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x along CZ -> N
}

void CsaMathTests::buildCarboxamideFrame() {
    model::MolFrameSpec s;
    s.kind = model::MolecularFrameKind::SidechainCarboxamide;
    s.origin = 0; s.xAnchor = 1; s.planeAnchor = 2;
    const auto a = buildAxes(s, {model::Vec3(0, 0, 0), model::Vec3(1, 0, 0), model::Vec3(0, 1, 0)});
    QVERIFY(a.has_value());
    QVERIFY(orthonormalRH(*a));
    QVERIFY(nearly(a->col(0).x(), 1.0));  // x along C -> O
}

QTEST_APPLESS_MAIN(CsaMathTests)
#include "csa_math_tests.moc"
