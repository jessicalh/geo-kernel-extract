#include "EfgFeature.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "RunData.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtTimeSeriesBuffers.h"

#include <QLoggingCategory>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cEfg, "h5reader.rediscover.efg")

bool validResidue(const model::QtProtein& p, int32_t r) {
    return r >= 0 && static_cast<std::size_t>(r) < p.residueCount();
}

FrameVariant classifyFrameVariant(const model::QtProtein& p, std::size_t atomIdx) {
    const model::QtAtom& a = p.atom(atomIdx);
    const bool hasResidue = validResidue(p, a.residueIndex);
    const model::QtResidue* r =
        hasResidue ? &p.residue(static_cast<std::size_t>(a.residueIndex)) : nullptr;
    const bool hasPrev = r && validResidue(p, r->prevResidueIndex);

    if (a.IsBackboneNitrogen())
        return hasPrev ? FrameVariant::BackboneN : FrameVariant::BackboneN_NTerminus;
    if (a.IsBackboneAlphaCarbon()) return FrameVariant::BackboneCA;
    if (a.IsBackboneCarbonylCarbon()) return FrameVariant::BackboneCarbonylC;
    if (a.IsBackboneCarbonylOxygen()) return FrameVariant::BackboneCarbonylO;
    if (a.IsBackboneAmideHydrogen())
        return hasPrev ? FrameVariant::HN_Standard : FrameVariant::HN_NTerminus;
    if (a.IsAnyAlphaHydrogen()) return FrameVariant::BackboneHA;
    return FrameVariant::Invalid;
}

bool apbsEfgPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtT2TimeSeries* efg = h5 ? h5->apbsEfg() : nullptr;
    return efg && body.catalog.present(body, ArrayId::ApbsEfg, atom, row)
           && efg->sourceAttachedAt(row);
}

void updateMagnitudeStats(EfgFeatureStats& stats, const std::array<double, 5>& t2) {
    if (!FiniteT2(t2)) return;
    const double mag = T2Magnitude(t2);
    if (!std::isfinite(mag)) return;
    if (stats.finite_efg == 0) {
        stats.min_efg_magnitude = mag;
        stats.max_efg_magnitude = mag;
    } else {
        stats.min_efg_magnitude = std::min(stats.min_efg_magnitude, mag);
        stats.max_efg_magnitude = std::max(stats.max_efg_magnitude, mag);
    }
    ++stats.finite_efg;
}

}  // namespace

EfgFeatureStats RunEfgPerAtomFeature(const Body& body, EfgFeatureSink& sink) {
    EfgFeatureStats stats;
    if (!body.run.protein || !body.run.trajectory()) return stats;

    const model::QtProtein& p = *body.run.protein;
    const QString units = body.catalog.has(ArrayId::ApbsEfg)
                              ? body.catalog.spec(ArrayId::ApbsEfg).unit
                              : QStringLiteral("V/Angstrom^2");

    qCInfo(cEfg).noquote() << "efg per_atom_feature | atoms=" << p.atomCount()
                           << "| dft rows=" << body.run.frameMap.dftRows().size()
                           << "| apbs_efg=" << body.catalog.has(ArrayId::ApbsEfg);

    const LocalFrame labFrame;  // invalid: BuildTarget emits lab-frame T2 only.
    for (std::size_t row : body.run.frameMap.dftRows()) {
        const std::size_t orig = body.run.frameMap.originalIndex(row);
        for (std::size_t atom = 0; atom < p.atomCount(); ++atom) {
            const DftTarget target = BuildTarget(body.run, atom, orig, labFrame);
            if (!target.present) continue;

            EfgFeatureRow out;
            const model::QtAtom& a = p.atom(atom);
            out.atom_index = a.atomIndex;
            out.residue_index = a.residueIndex;
            out.element = static_cast<int>(a.element);
            out.atom_name = p.atomLabel(atom, model::NamingConvention::Bmrb);
            if (validResidue(p, a.residueIndex)) {
                const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
                out.residue_number = r.address.residueNumber;
                out.amino_acid = static_cast<int>(r.aminoAcid);
            }
            out.frame_variant = static_cast<int>(classifyFrameVariant(p, atom));
            out.h5_row = static_cast<int32_t>(row);
            out.original_index = static_cast<int32_t>(orig);
            out.time_ps = body.run.trajectory()->timePicoseconds(row);
            out.dft_present = true;
            out.dft_target_T2 = target.total_decomp.T2;
            out.efg_units = units;

            out.apbs_efg_present = apbsEfgPresent(body, atom, row);
            if (out.apbs_efg_present) {
                // APBS EFG in H5 was written by ApbsFieldResult as
                // SphericalTensor::Decompose(EFG).PackT2, whose T2 order and
                // isometric normalization are exactly DecomposeLibrary's
                // [xy, yz, zz, xz, xx-yy] basis. Do not re-project in Python.
                out.efg_feature_T2 = body.catalog.valueT2(body, ArrayId::ApbsEfg, atom, row);
                ++stats.apbs_efg_present;
                updateMagnitudeStats(stats, out.efg_feature_T2);
            }

            sink.Write(out);
            ++stats.rows;
            ++stats.dft_present;
        }
    }

    qCInfo(cEfg).noquote() << "efg per_atom_feature rows=" << stats.rows
                           << "| dft_present=" << stats.dft_present
                           << "| apbs_efg_present=" << stats.apbs_efg_present
                           << "| finite_efg=" << stats.finite_efg
                           << "| |EFG| min=" << stats.min_efg_magnitude
                           << "max=" << stats.max_efg_magnitude;
    return stats;
}

}  // namespace h5reader::rediscover
