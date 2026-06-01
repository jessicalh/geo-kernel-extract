#include "BuckinghamEfield.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "Relationship.h"
#include "RunData.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtTimeSeriesBuffers.h"

#include <QLoggingCategory>

#include <algorithm>
#include <cmath>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cBuckingham, "h5reader.rediscover.buckingham")

bool validResidue(const model::QtProtein& p, int32_t r) {
    return r >= 0 && static_cast<std::size_t>(r) < p.residueCount();
}

bool isBackboneFeatureAtom(const model::QtAtom& a) {
    return a.IsBackbone() || a.IsAnyAlphaHydrogen();
}

FrameResult backboneFrame(const Body& body, std::size_t atom, std::size_t frame) {
    FrameResult fr;
    const model::QtProtein& p = *body.run.protein;
    const model::QtAtom& a = p.atom(atom);
    if (!validResidue(p, a.residueIndex)) return fr;
    const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));

    auto posOf = [&](int32_t ai) {
        return verbs::pos(body, static_cast<std::size_t>(ai), frame);
    };

    if (a.IsBackboneNitrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        bool cPrevValid = false;
        Vec3 cRef = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cRef = posOf(prev.C);
                cPrevValid = true;
            }
        }
        if (!cPrevValid && r.C != model::QtResidue::NONE) cRef = posOf(r.C);
        fr.frame = BuildBackboneNFrame(posOf(r.N), posOf(r.CA), cRef, cPrevValid);
        fr.anchor_atom_index = cPrevValid && validResidue(p, r.prevResidueIndex)
                                   ? p.residue(static_cast<std::size_t>(r.prevResidueIndex)).C
                                   : r.C;
        return fr;
    }

    if (a.IsBackboneAlphaCarbon()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE
            || r.C == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCaFrame(posOf(r.CA), posOf(r.N), posOf(r.C));
        fr.anchor_atom_index = r.N;
        return fr;
    }

    if (a.IsBackboneCarbonylCarbon()) {
        if (r.C == model::QtResidue::NONE || r.O == model::QtResidue::NONE
            || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylCFrame(posOf(r.C), posOf(r.O), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }

    if (a.IsBackboneCarbonylOxygen()) {
        if (r.O == model::QtResidue::NONE || r.C == model::QtResidue::NONE
            || r.CA == model::QtResidue::NONE)
            return fr;
        fr.frame = BuildBackboneCarbonylOFrame(posOf(r.O), posOf(r.C), posOf(r.CA));
        fr.anchor_atom_index = r.CA;
        return fr;
    }

    if (a.IsBackboneAmideHydrogen()) {
        if (r.N == model::QtResidue::NONE || r.CA == model::QtResidue::NONE) return fr;
        bool cPrevValid = false;
        Vec3 cPrev = Vec3::Zero();
        if (validResidue(p, r.prevResidueIndex)) {
            const model::QtResidue& prev = p.residue(static_cast<std::size_t>(r.prevResidueIndex));
            if (prev.C != model::QtResidue::NONE) {
                cPrev = posOf(prev.C);
                cPrevValid = true;
            }
        }
        fr.frame = BuildHNFrame(posOf(r.N), verbs::pos(body, atom, frame), posOf(r.CA),
                                cPrev, cPrevValid);
        fr.anchor_atom_index = -1;
        return fr;
    }

    if (a.IsAnyAlphaHydrogen()) {
        if (r.CA == model::QtResidue::NONE || r.N == model::QtResidue::NONE) return fr;
        fr.frame = BuildBackboneHaFrame(verbs::pos(body, atom, frame), posOf(r.CA), posOf(r.N));
        fr.anchor_atom_index = r.N;
        return fr;
    }

    return fr;
}

bool apbsEfieldPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtVec3TimeSeries* efield = h5 ? h5->apbsEfield() : nullptr;
    return efield && body.catalog.present(body, ArrayId::ApbsEfield, atom, row)
           && efield->sourceAttachedAt(row);
}

void updateMagnitudeStats(BuckinghamEfieldStats& stats, double mag) {
    if (!std::isfinite(mag)) return;
    if (stats.finite_efield == 0) {
        stats.min_efield_magnitude = mag;
        stats.max_efield_magnitude = mag;
    } else {
        stats.min_efield_magnitude = std::min(stats.min_efield_magnitude, mag);
        stats.max_efield_magnitude = std::max(stats.max_efield_magnitude, mag);
    }
    ++stats.finite_efield;
}

}  // namespace

BuckinghamEfieldStats RunBuckinghamEfieldPerAtomFeature(const Body& body,
                                                        BuckinghamEfieldSink& sink) {
    BuckinghamEfieldStats stats;
    if (!body.run.protein || !body.run.trajectory()) return stats;

    const model::QtProtein& p = *body.run.protein;
    const QString units = body.catalog.has(ArrayId::ApbsEfield)
                              ? body.catalog.spec(ArrayId::ApbsEfield).unit
                              : QStringLiteral("V/Angstrom");

    qCInfo(cBuckingham).noquote()
        << "buckingham_efield per_atom_feature | atoms=" << p.atomCount()
        << "| dft rows=" << body.run.frameMap.dftRows().size()
        << "| apbs_efield=" << body.catalog.has(ArrayId::ApbsEfield);

    for (std::size_t row : body.run.frameMap.dftRows()) {
        const std::size_t orig = body.run.frameMap.originalIndex(row);
        for (std::size_t atom = 0; atom < p.atomCount(); ++atom) {
            const model::QtAtom& a = p.atom(atom);
            if (!isBackboneFeatureAtom(a)) continue;

            const FrameResult fr = backboneFrame(body, atom, row);
            const DftTarget target = BuildTarget(body.run, atom, orig, fr.frame);
            if (!target.present) continue;

            BuckinghamEfieldRow out;
            out.atom_index = a.atomIndex;
            out.residue_index = a.residueIndex;
            out.element = static_cast<int>(a.element);
            out.atom_name = p.atomLabel(atom, model::NamingConvention::Bmrb);
            if (validResidue(p, a.residueIndex)) {
                const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
                out.residue_number = r.address.residueNumber;
                out.amino_acid = static_cast<int>(r.aminoAcid);
            }
            out.frame_variant = static_cast<int>(fr.frame.variant);
            out.frame_valid = fr.frame.is_valid;
            out.frame_anchor_atom_index = fr.anchor_atom_index;
            out.frame_z = fr.frame.z;
            out.frame_x = fr.frame.x;
            out.frame_y = fr.frame.y;
            out.h5_row = static_cast<int32_t>(row);
            out.original_index = static_cast<int32_t>(orig);
            out.time_ps = body.run.trajectory()->timePicoseconds(row);
            out.dft_present = true;
            out.dft_total_raw = target.total_raw;
            out.dft_total_decomp = target.total_decomp;
            out.efield_units = units;

            if (out.frame_valid) ++stats.frame_valid;
            out.apbs_efield_present = apbsEfieldPresent(body, atom, row);
            if (out.apbs_efield_present) {
                const Vec3 eLab = body.catalog.valueVec3(body, ArrayId::ApbsEfield, atom, row);
                out.e_mag = eLab.norm();
                updateMagnitudeStats(stats, out.e_mag);
                if (out.frame_valid && FiniteVec3(eLab)) {
                    out.efield_local = fr.frame.ToLocal(eLab);
                    out.e_proj = out.efield_local.z();
                }
                ++stats.apbs_efield_present;
            }

            sink.Write(out);
            ++stats.rows;
            ++stats.dft_present;
        }
    }

    qCInfo(cBuckingham).noquote()
        << "buckingham_efield per_atom_feature rows=" << stats.rows
        << "| dft_present=" << stats.dft_present
        << "| frame_valid=" << stats.frame_valid
        << "| apbs_efield_present=" << stats.apbs_efield_present
        << "| finite_efield=" << stats.finite_efield
        << "| |E| min=" << stats.min_efield_magnitude
        << "max=" << stats.max_efield_magnitude;
    return stats;
}

}  // namespace h5reader::rediscover
