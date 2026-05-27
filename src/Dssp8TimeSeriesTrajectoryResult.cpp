#include "Dssp8TimeSeriesTrajectoryResult.h"

#include "Atom.h"
#include "DsspResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Translate DsspResidue::HBondPartner.residue_index (size_t with
// SIZE_MAX sentinel) to int32 with -1 sentinel. Out-of-range residue
// indices are also treated as "no partner." The INT32_MAX cap prevents
// an unguarded narrowing cast.
std::int32_t PartnerToInt32(std::size_t partner_idx, std::size_t n_res) {
    if (partner_idx >= n_res) return Dssp8TimeSeriesTrajectoryResult::kNoPartner;
    if (partner_idx > static_cast<std::size_t>(
            std::numeric_limits<std::int32_t>::max()))
        return Dssp8TimeSeriesTrajectoryResult::kNoPartner;
    return static_cast<std::int32_t>(partner_idx);
}

}  // anonymous namespace


std::unique_ptr<Dssp8TimeSeriesTrajectoryResult>
Dssp8TimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<Dssp8TimeSeriesTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = tp.AtomCount();

    r->ss8_code_.assign(R, {});
    r->hbond_acceptor_partner_.assign(R, {});
    r->hbond_acceptor_energy_.assign(R, {});
    r->hbond_donor_partner_.assign(R, {});
    r->hbond_donor_energy_.assign(R, {});

    r->residue_index_per_atom_.assign(N, -1);
    for (std::size_t ai = 0; ai < N; ++ai) {
        r->residue_index_per_atom_[ai] =
            static_cast<std::int32_t>(protein.AtomAt(ai).residue_index);
    }
    return r;
}


// ── Compute ──────────────────────────────────────────────────────────
//
// Per frame: if DsspResult attached, copy per-residue SS code + H-bond
// partners/energies into growing buffers. If not attached, push absent
// sentinels (ss8_code=kSSUnassigned, partner=kNoPartner, energy=NaN).

void Dssp8TimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<DsspResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    const Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    if (source_attached) {
        const auto& dssp = conf.Result<DsspResult>();
        const auto& dssp_residues = dssp.AllResidues();
        const std::size_t dssp_R = dssp_residues.size();
        for (std::size_t ri = 0; ri < R; ++ri) {
            // Codex review 2026-05-19: index-bounds check is necessary
            // but NOT sufficient. DsspResult resizes residues_ to
            // ResidueCount() and leaves unmapped entries at default
            // (secondary_structure='C', observed=false). Treating
            // `ri < dssp_R` as observed would silently bias unmapped
            // residues (caps, insertion-coded mismatches, DSSP skips)
            // into real coil. Use the observed flag instead.
            const bool residue_observed =
                ri < dssp_R && dssp_residues[ri].observed;
            if (residue_observed) {
                const auto& dr = dssp_residues[ri];
                ss8_code_[ri].push_back(Ss8Code(dr.secondary_structure));
                hbond_acceptor_partner_[ri].push_back({
                    PartnerToInt32(dr.acceptors[0].residue_index, R),
                    PartnerToInt32(dr.acceptors[1].residue_index, R)});
                hbond_acceptor_energy_[ri].push_back({
                    dr.acceptors[0].residue_index < R ? dr.acceptors[0].energy : kNaN,
                    dr.acceptors[1].residue_index < R ? dr.acceptors[1].energy : kNaN});
                hbond_donor_partner_[ri].push_back({
                    PartnerToInt32(dr.donors[0].residue_index, R),
                    PartnerToInt32(dr.donors[1].residue_index, R)});
                hbond_donor_energy_[ri].push_back({
                    dr.donors[0].residue_index < R ? dr.donors[0].energy : kNaN,
                    dr.donors[1].residue_index < R ? dr.donors[1].energy : kNaN});
            } else {
                // DSSP didn't actually map this residue (size mismatch
                // OR observed=false). Treat as absent for this residue
                // this frame so consumers can distinguish "DSSP said
                // coil" from "DSSP never wrote here."
                ss8_code_[ri].push_back(kSSUnassigned);
                hbond_acceptor_partner_[ri].push_back({kNoPartner, kNoPartner});
                hbond_acceptor_energy_[ri].push_back({kNaN, kNaN});
                hbond_donor_partner_[ri].push_back({kNoPartner, kNoPartner});
                hbond_donor_energy_[ri].push_back({kNaN, kNaN});
            }
        }
    } else {
        // Source not attached this frame: push absent sentinels across
        // all residues. Source_attached_per_frame mask captures the
        // provenance.
        for (std::size_t ri = 0; ri < R; ++ri) {
            ss8_code_[ri].push_back(kSSUnassigned);
            hbond_acceptor_partner_[ri].push_back({kNoPartner, kNoPartner});
            hbond_acceptor_energy_[ri].push_back({kNaN, kNaN});
            hbond_donor_partner_[ri].push_back({kNoPartner, kNoPartner});
            hbond_donor_energy_[ri].push_back({kNaN, kNaN});
        }
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void Dssp8TimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "Dssp8TimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(ss8_code_.size()) + " residues.");
}


// ── WriteH5Group ─────────────────────────────────────────────────────
//
// Skip-emission when source never attached (absent, not faked). Float
// data NaN-filled on absent frames; int sentinels (-1 partner,
// 255 ss8_code) for absent. residue_index_per_atom emitted alongside
// for SDK/viewer atom-axis broadcast.

void Dssp8TimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t R = ss8_code_.size();
    const std::size_t T = n_frames_;
    const std::size_t N = tp.AtomCount();

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "Dssp8TimeSeriesTrajectoryResult::WriteH5Group",
            "DsspResult was not attached in any of " +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/dssp8_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    auto grp = file.createGroup("/trajectory/dssp8_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_residues",  R);
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("source_attached_count", source_attached_count);

    grp.createAttribute("ss8_legend", std::string(
        "H=0 (alpha helix), G=1 (3_10 helix), I=2 (pi helix), "
        "E=3 (extended strand), B=4 (beta bridge), T=5 (turn), "
        "S=6 (bend), C=7 (coil); 255 = no observation."));
    grp.createAttribute("ss8_unassigned_sentinel", std::uint8_t{kSSUnassigned});
    grp.createAttribute("hbond_partner_sentinel", std::string(
        "-1 (kNoPartner) = no partner at this slot. DSSP records up to "
        "2 acceptors and 2 donors per residue; slot [0] is the lower-"
        "energy partner."));
    grp.createAttribute("hbond_energy_units",            std::string("kcal/mol"));
    grp.createAttribute("hbond_energy_absent_sentinel",  std::string("NaN"));
    grp.createAttribute("hbond_threshold", std::string(
        "Kabsch-Sander 1983 (J. Biomol. Struct. Dyn. 1: 879) defines "
        "E < -0.5 kcal/mol as the threshold for COUNTING an H-bond "
        "toward SS classification (helix detection etc.). libdssp "
        "writes the two best (lowest-energy) acceptor and donor "
        "partners to the slot REGARDLESS of strict threshold, so "
        "consumers may see observed slot energies in [-0.5, 0] -- "
        "real (attractive) but not strict-K-S-counted. For consumers "
        "wanting strict K-S H-bond stats: filter on `hbond_*_energy "
        "< -0.5`. For consumers wanting all attractive interactions: "
        "filter on `hbond_*_energy < 0`."));
    grp.createAttribute("source", std::string(
        "DsspResult (libdssp via Joosten 2011 / Kabsch-Sander 1983). "
        "Phi/Psi/SASA from DsspResult are NOT mirrored here -- phi/psi "
        "live in DihedralTimeSeriesTrajectoryResult (which negates DSSP's "
        "libdssp values to recover IUPAC convention), and SASA lives in "
        "SasaResult."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- DsspResult attaches only when "
        "PerFrameRunOptions::skip_dssp == false. Reader contract: "
        "consult source_attached_per_frame mask before downstream stats; "
        "absent frames have ss8_code=255, partner=-1, energy=NaN."));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));
    grp.createAttribute("atom_axis",    std::string("protein_atom_index"));

    // ── Per-frame metadata (T,) ──────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));

    // ── ss8_code (R, T) uint8 ────────────────────────────────────────
    {
        std::vector<std::uint8_t> flat(R * T, kSSUnassigned);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = ss8_code_[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[ri * T + f] = row[f];
            }
        }
        const std::vector<std::size_t> dims = {R, T};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint8_t>("ss8_code", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("category"));
    }

    // ── Helpers for the four 2-slot (R, T, 2) datasets ───────────────
    auto emit_partner = [&](const std::string& name,
                             const std::vector<std::vector<
                                 std::array<std::int32_t, 2>>>& src) {
        std::vector<std::int32_t> flat(R * T * 2, kNoPartner);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = src[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[(ri * T + f) * 2 + 0] = row[f][0];
                flat[(ri * T + f) * 2 + 1] = row[f][1];
            }
        }
        const std::vector<std::size_t> dims = {R, T, std::size_t(2)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::int32_t>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("residue_index"));
    };
    auto emit_energy = [&](const std::string& name,
                            const std::vector<std::vector<
                                std::array<double, 2>>>& src) {
        std::vector<double> flat(R * T * 2, kNaN);
        for (std::size_t ri = 0; ri < R; ++ri) {
            const auto& row = src[ri];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[(ri * T + f) * 2 + 0] = row[f][0];
                flat[(ri * T + f) * 2 + 1] = row[f][1];
            }
        }
        const std::vector<std::size_t> dims = {R, T, std::size_t(2)};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("kcal/mol"));
    };

    emit_partner("hbond_acceptor_partner", hbond_acceptor_partner_);
    emit_energy ("hbond_acceptor_energy",  hbond_acceptor_energy_);
    emit_partner("hbond_donor_partner",    hbond_donor_partner_);
    emit_energy ("hbond_donor_energy",     hbond_donor_energy_);

    // ── Per-atom lookup (N,) ─────────────────────────────────────────
    grp.createDataSet("residue_index_per_atom", residue_index_per_atom_)
       .createAttribute("units", std::string("residue_index"));

    OperationLog::Info(LogCalcOther,
        "Dssp8TimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/dssp8_time_series with " +
        std::to_string(R) + " residues x " + std::to_string(T) +
        " frames (source attached in " +
        std::to_string(source_attached_count) + " / " +
        std::to_string(source_attached_per_frame_.size()) + " frames)");
}


}  // namespace nmr
