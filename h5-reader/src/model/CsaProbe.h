// CsaProbe -- compute one atom's DFT shielding-tensor CSA result for a given
// frame, aligned onto the displayed (stabilized) molecule. This is the SINGLE
// source of truth shared by the live glyph driver
// (ReaderMainWindow::updateCsaGlyph) and the REST probe (GET /csa), so the
// picture and the queried numbers can never disagree.
//
// Pulls the raw DFT tensor, builds the molecular frame in raw + display coords,
// lifts sigma into display coords (the alignment bridge: sigma_mol is rotation-
// invariant, so M_raw^T sigma M_raw then M_disp ... M_disp^T), and runs the PAS
// decomposition. Model-level orchestration; no Qt widgets / VTK.

#pragma once

#include "Conformation.h"
#include "ConformationGeometry.h"
#include "CsaShape.h"
#include "DftShieldingStore.h"
#include "MolecularFrame.h"
#include "MolecularFrameSelect.h"
#include "QtProtein.h"
#include "TransformedConformation.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>

namespace h5reader::model {

struct AtomCsaResult {
    bool dftPresent = false;                 // a DFT job exists for this frame
    bool valid = false;                      // a CSA shape was computed
    bool framed = false;                     // a molecular frame was built
    MolecularFrameKind frameKind = MolecularFrameKind::None;
    CsaShape shape;                          // PAS in display coords (valid only)
    std::optional<Mat3> molecularAxes;       // display-frame axes (framed only)
    Vec3 atomPos = Vec3::Zero();             // the atom's display position
};

// rawConf = the loader's untransformed trajectory; transformed = the
// Kabsch-stabilized display conformation the scene renders. requestDftFrame
// keeps REST/heroshot probes deterministic; live glyph refresh can pass false
// so frame advance never waits on ORCA output parsing.
inline AtomCsaResult ComputeAtomCsa(const QtProtein& protein,
                                    Conformation& rawConf,
                                    TransformedConformation& transformed,
                                    DftShieldingStore& dftStore,
                                    std::size_t atom,
                                    std::size_t frame,
                                    bool requestDftFrame = true) {
    AtomCsaResult out;
    if (atom >= protein.atomCount()) return out;
    out.atomPos = transformed.atomPosition(frame, atom);

    const std::size_t original = rawConf.originalFrameIndex(frame);
    if (!dftStore.hasJob(original)) return out;  // honest gap; dftPresent stays false
    if (requestDftFrame)
        dftStore.requestFrame(original);
    const DftShieldingFrame* dft = dftStore.frame(original);
    if (!dft || !dft->valid || atom >= dft->atoms.size()) return out;
    out.dftPresent = true;
    const Mat3 sigmaRaw = dft->atoms[atom].total_raw;

    if (const auto spec = SelectMolecularFrameSpec(protein, atom)) {
        const auto posRaw = [&](std::int32_t i) -> std::optional<Vec3> {
            if (i < 0 || static_cast<std::size_t>(i) >= protein.atomCount()) return std::nullopt;
            return rawConf.atomPosition(frame, static_cast<std::size_t>(i));
        };
        const auto ringRaw = [&](std::int32_t r) -> std::optional<std::pair<Vec3, Vec3>> {
            if (r < 0) return std::nullopt;
            const RingGeometry g = RingGeometryAt(rawConf, static_cast<std::size_t>(r), frame);
            return std::make_pair(g.center, g.normal);
        };
        const auto posDisp = [&](std::int32_t i) -> std::optional<Vec3> {
            if (i < 0 || static_cast<std::size_t>(i) >= protein.atomCount()) return std::nullopt;
            return transformed.atomPosition(frame, static_cast<std::size_t>(i));
        };
        const auto ringDisp = [&](std::int32_t r) -> std::optional<std::pair<Vec3, Vec3>> {
            if (r < 0) return std::nullopt;
            const RingGeometry g = RingGeometryAt(transformed, static_cast<std::size_t>(r), frame);
            return std::make_pair(g.center, g.normal);
        };
        MolFrameContinuity contRaw;
        MolFrameContinuity contDisp;
        const auto mRaw = BuildMolecularFrameAxes(*spec, posRaw, ringRaw, contRaw);
        const auto mDisp = BuildMolecularFrameAxes(*spec, posDisp, ringDisp, contDisp);
        if (mRaw && mDisp) {
            const Mat3 sigmaMol = mRaw->transpose() * sigmaRaw * (*mRaw);
            const Mat3 sigmaDisp = (*mDisp) * sigmaMol * mDisp->transpose();
            out.shape = ComputeCsaShape(sigmaDisp);
            if (out.shape.valid) {
                out.molecularAxes = *mDisp;
                out.framed = true;
                out.frameKind = spec->kind;
            }
        }
    }
    if (!out.shape.valid)
        out.shape = ComputeCsaShape(sigmaRaw);  // unframed fallback (raw-frame PAS)
    out.valid = out.shape.valid;
    return out;
}

}  // namespace h5reader::model
