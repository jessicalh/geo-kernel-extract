#include "Aimnet2Feature.h"

#include "Catalog.h"
#include "ExtractionSupport.h"
#include "LocalFrameBasis.h"
#include "Relationship.h"
#include "RunData.h"
#include "SphericalBasis.h"
#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/QtTimeSeriesBuffers.h"

#include <QLoggingCategory>

#include <cmath>
#include <stdexcept>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cAimnet2, "h5reader.rediscover.aimnet2")

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

bool aimnetChargePresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtScalarTimeSeries* charge = h5 ? h5->aimnet2Charge() : nullptr;
    return charge && body.catalog.present(body, ArrayId::Aimnet2Charge, atom, row)
           && charge->sourceAttachedAt(row);
}

bool aimnetCrgPresent(const Body& body, std::size_t atom, std::size_t row) {
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtAimnet2ChargeResponseGradientTimeSeries* crg =
        h5 ? h5->aimnet2ChargeResponseGradient() : nullptr;
    return crg && body.catalog.present(body, ArrayId::Aimnet2ChargeRespScalar, atom, row)
           && body.catalog.present(body, ArrayId::Aimnet2ChargeRespVector, atom, row)
           && crg->meta.sourceAttachedAt(row);
}

bool aimnetEmbeddingPresent(const Body& body, std::size_t atom, std::size_t row,
                            std::size_t& dims, const float*& ptr) {
    dims = 0;
    ptr = body.catalog.valueEmbedding(body, ArrayId::Aimnet2Embedding, atom, row, dims);
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtEmbeddingTimeSeries* emb = h5 ? h5->aimnet2Embedding() : nullptr;
    return ptr && emb && body.catalog.present(body, ArrayId::Aimnet2Embedding, atom, row)
           && emb->meta.sourceAttachedAt(row);
}

}  // namespace

Aimnet2FeatureStats RunAimnet2PerAtomFeature(const Body& body, Aimnet2FeatureSink& sink) {
    Aimnet2FeatureStats stats;
    if (!body.run.protein || !body.run.trajectory()) return stats;

    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtEmbeddingTimeSeries* emb = h5 ? h5->aimnet2Embedding() : nullptr;
    if (emb && emb->n_dims != sink.embeddingDims()) {
        throw std::runtime_error(QStringLiteral("aimnet2_embedding dimension %1 != expected %2")
                                     .arg(static_cast<qulonglong>(emb->n_dims))
                                     .arg(static_cast<qulonglong>(sink.embeddingDims()))
                                     .toStdString());
    }

    const model::QtProtein& p = *body.run.protein;
    qCInfo(cAimnet2).noquote()
        << "aimnet2_features per_atom_feature | atoms=" << p.atomCount()
        << "| dft rows=" << body.run.frameMap.dftRows().size()
        << "| charge=" << body.catalog.has(ArrayId::Aimnet2Charge)
        << "| crg_scalar=" << body.catalog.has(ArrayId::Aimnet2ChargeRespScalar)
        << "| crg_vector=" << body.catalog.has(ArrayId::Aimnet2ChargeRespVector)
        << "| embedding=" << body.catalog.has(ArrayId::Aimnet2Embedding)
        << "| embedding_dims=" << (emb ? emb->n_dims : 0);

    for (std::size_t row : body.run.frameMap.dftRows()) {
        const std::size_t orig = body.run.frameMap.originalIndex(row);
        for (std::size_t atom = 0; atom < p.atomCount(); ++atom) {
            const model::QtAtom& a = p.atom(atom);
            if (!isBackboneFeatureAtom(a)) continue;

            const FrameResult fr = backboneFrame(body, atom, row);
            const DftTarget target = BuildTarget(body.run, atom, orig, fr.frame);
            if (!target.present) continue;

            Aimnet2FeatureRow out;
            out.atom_index = a.atomIndex;
            out.residue_index = a.residueIndex;
            out.element = static_cast<int>(a.element);
            out.atom_name = p.atomLabel(atom, model::NamingConvention::Bmrb);
            if (validResidue(p, a.residueIndex)) {
                const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
                out.residue_number = r.address.residueNumber;
                out.amino_acid = static_cast<int>(r.aminoAcid);
            }
            out.frame_variant =
                fr.frame.is_valid ? static_cast<int>(fr.frame.variant)
                                  : static_cast<int>(classifyFrameVariant(p, atom));
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
            out.dft_local_frame_valid = target.local_frame_valid;
            if (target.local_frame_valid) out.dft_target_local_T2 = DecomposeLibrary(target.total_local).T2;

            if (out.frame_valid) ++stats.frame_valid;

            out.aimnet2_charge_present = aimnetChargePresent(body, atom, row);
            if (out.aimnet2_charge_present) {
                out.aimnet2_charge = body.catalog.value(body, ArrayId::Aimnet2Charge, atom, row);
                if (std::isfinite(out.aimnet2_charge)) ++stats.charge_present;
            }

            out.crg_present = aimnetCrgPresent(body, atom, row);
            if (out.crg_present) {
                out.crg_scalar = body.catalog.value(body, ArrayId::Aimnet2ChargeRespScalar, atom, row);
                out.crg_vector_lab = body.catalog.valueVec3(body, ArrayId::Aimnet2ChargeRespVector, atom, row);
                if (out.frame_valid) out.crg_vector_local = fr.frame.ToLocal(out.crg_vector_lab);
                if (std::isfinite(out.crg_scalar) && FiniteAimnetVec3(out.crg_vector_lab))
                    ++stats.crg_present;
            }

            const float* embeddingPtr = nullptr;
            std::size_t embeddingDims = 0;
            out.embedding_present = aimnetEmbeddingPresent(body, atom, row, embeddingDims, embeddingPtr);
            out.embedding = embeddingPtr;
            out.embedding_dims = embeddingDims;
            if (out.embedding_present && embeddingDims == sink.embeddingDims()) ++stats.embedding_present;

            sink.Write(out);
            ++stats.rows;
            ++stats.dft_present;
        }
    }

    stats.embedding_dims = sink.embeddingDims();
    qCInfo(cAimnet2).noquote()
        << "aimnet2_features per_atom_feature rows=" << stats.rows
        << "| dft_present=" << stats.dft_present
        << "| frame_valid=" << stats.frame_valid
        << "| charge_present=" << stats.charge_present
        << "| crg_present=" << stats.crg_present
        << "| embedding_present=" << stats.embedding_present;
    return stats;
}

}  // namespace h5reader::rediscover
