#include "HBondResult.h"
#include "Protein.h"
#include "DsspResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <optional>
#include <queue>
#include <set>

namespace nmr {

std::vector<std::type_index> HBondResult::Dependencies() const {
    return {
        std::type_index(typeid(DsspResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// An identified H-bond from DSSP, resolved to atom positions.
//
// DSSP identifies backbone H-bonds by the Kabsch-Sander energy criterion
// between residue pairs. Each H-bond has:
//   - donor residue (N-H donates)
//   - acceptor residue (C=O accepts)
//
// We resolve this to atoms:
//   - donor_N: the backbone N metadata endpoint of the donor residue
//   - donor_H: the explicit backbone H used as the physical source point
//   - acceptor_O: the backbone O of the acceptor residue
//   - h_hat: unit direction from donor H to acceptor O
//   - distance: |H...O| source extent
// ============================================================================

struct ResolvedHBond {
    size_t donor_N = SIZE_MAX;
    size_t donor_H = SIZE_MAX;
    size_t acceptor_O = SIZE_MAX;
    size_t donor_residue = SIZE_MAX;
    size_t acceptor_residue = SIZE_MAX;
    Vec3 source_point = Vec3::Zero();
    Vec3 h_hat = Vec3::Zero();        // donor H → acceptor O direction
    double distance = 0.0;            // H...O distance
    int sequence_separation = 0;
};


// ============================================================================
// The clean H-bond angular scalar from one H-bond at one atom.
// The former McConnell-form rank-2 kernel is intentionally no longer
// computed or emitted.
// ============================================================================

namespace hbond_result_detail {

int BackboneResidueSeparation(const Protein& protein,
                              size_t residue_a,
                              size_t residue_b) {
    const size_t residue_count = protein.ResidueCount();
    if (residue_a >= residue_count || residue_b >= residue_count) return -1;
    if (residue_a == residue_b) return 0;

    std::vector<int> distance(residue_count, -1);
    std::queue<size_t> frontier;
    distance[residue_a] = 0;
    frontier.push(residue_a);

    while (!frontier.empty()) {
        const size_t current = frontier.front();
        frontier.pop();

        const std::optional<size_t> neighbours[2] = {
            protein.BackbonePredecessor(current),
            protein.BackboneSuccessor(current),
        };
        for (const auto& neighbour : neighbours) {
            if (!neighbour || *neighbour >= residue_count) continue;
            if (distance[*neighbour] >= 0) continue;
            distance[*neighbour] = distance[current] + 1;
            if (*neighbour == residue_b) return distance[*neighbour];
            frontier.push(*neighbour);
        }
    }
    return -1;
}

namespace {

int MinimumKnownSeparation(int a, int b) {
    if (a < 0) return b;
    if (b < 0) return a;
    return std::min(a, b);
}

}  // namespace

KernelResult ComputeKernel(
        const Vec3& atom_pos,
        const Vec3& donor_h_pos,
        const Vec3& h_hat) {

    KernelResult result;

    // Field direction from the explicit donor H source point to the atom.
    Vec3 d = atom_pos - donor_h_pos;
    double r = d.norm();

    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;

    result.distance = r;
    double r3 = r * r * r;
    Vec3 d_hat = d / r;
    double cos_theta = d_hat.dot(h_hat);

    // angular factor
    result.f = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    return result;
}

}  // namespace hbond_result_detail


// ============================================================================
// HBondResult::Compute
// ============================================================================

std::unique_ptr<HBondResult> HBondResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("HBondResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& dssp = conf.Result<DsspResult>();
    // SpatialIndexResult is a declared dependency (ensures it's computed
    // before us) but we iterate H-bonds directly, not via spatial search.
    (void)conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_residues = protein.ResidueCount();

    auto result_ptr = std::make_unique<HBondResult>();
    result_ptr->conf_ = &conf;

    GeometryChoiceBuilder choices(conf);

    // ------------------------------------------------------------------
    // Step 1: Resolve DSSP H-bond partners to atom positions.
    //
    // DSSP provides up to 2 acceptor partners and 2 donor partners per
    // residue. Each partner is a residue index. We resolve to backbone
    // N/H (donor) and O (acceptor) atoms.
    //
    // Skip H-bonds where:
    //   - partner residue index is SIZE_MAX (no partner)
    //   - backbone N, donor H, or acceptor O atoms are not present
    //   - sequence separation below hbond_sequential_exclusion_residues
    //   - H...O distance outside the configured H-bond source range
    // ------------------------------------------------------------------

    std::vector<ResolvedHBond> hbonds;

    // Use a set to deduplicate: (donor_N, acceptor_O) pairs
    std::set<std::pair<size_t, size_t>> seen;
    size_t choice_record_seq = 0;

    for (size_t ri = 0; ri < n_residues; ++ri) {
        const auto& dr = dssp.AllResidues()[ri];
        const Residue& res = protein.ResidueAt(ri);

        // donor side: this residue's backbone N -> partner's acceptor O
        for (int partner_slot = 0; partner_slot < 2; ++partner_slot) {
            size_t acceptor_residue_idx = dr.acceptors[partner_slot].residue_index;
            if (acceptor_residue_idx == SIZE_MAX || acceptor_residue_idx >= n_residues) continue;

            const Residue& acc_res = protein.ResidueAt(acceptor_residue_idx);
            if (res.N == Residue::NONE || acc_res.O == Residue::NONE) continue;
            if (res.H == Residue::NONE) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, acceptor_residue_idx](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(acceptor_residue_idx), "index");
                        AddNumber(gc, "rejection", 1.0, "missing_explicit_donor_H");
                    });
                continue;
            }

            // sequence exclusion
            const int seq_sep =
                hbond_result_detail::BackboneResidueSeparation(
                    protein, ri, acceptor_residue_idx);
            if (seq_sep >= 0 &&
                seq_sep < static_cast<int>(CalculatorConfig::Get(
                    "hbond_sequential_exclusion_residues"))) {
                // record: seq-sep reject
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, acceptor_residue_idx, seq_sep](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(acceptor_residue_idx), "index");
                        AddNumber(gc, "sequence_separation", static_cast<double>(seq_sep), "residues");
                        AddNumber(gc, "rejection", 1.0, "seq_sep_too_small");
                    });
                continue;
            }

            auto key = std::make_pair(res.N, acc_res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 H_pos = conf.PositionAt(res.H);
            Vec3 O_pos = conf.PositionAt(acc_res.O);
            Vec3 d = O_pos - H_pos;
            double dist = d.norm();

            // source distance gate (singular or beyond max H...O)
            if (dist < CalculatorConfig::Get("singularity_guard_distance") || dist > CalculatorConfig::Get("hbond_dipolar_max_distance")) {
                // record: distance reject
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, acceptor_residue_idx, dist](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(acceptor_residue_idx), "index");
                        AddNumber(gc, "distance", dist, "A");
                        AddNumber(gc, "rejection", 1.0, "distance_out_of_range");
                    });
                continue;
            }

            ResolvedHBond hb;
            hb.donor_N = res.N;
            hb.donor_H = res.H;
            hb.acceptor_O = acc_res.O;
            hb.donor_residue = ri;
            hb.acceptor_residue = acceptor_residue_idx;
            hb.source_point = H_pos;
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }

        // acceptor side: partner's donor N -> this residue's backbone O
        for (int partner_slot = 0; partner_slot < 2; ++partner_slot) {
            size_t donor_residue_idx = dr.donors[partner_slot].residue_index;
            if (donor_residue_idx == SIZE_MAX || donor_residue_idx >= n_residues) continue;

            const Residue& don_res = protein.ResidueAt(donor_residue_idx);
            if (don_res.N == Residue::NONE || res.O == Residue::NONE) continue;
            if (don_res.H == Residue::NONE) {
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [donor_residue_idx, ri](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(donor_residue_idx), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "rejection", 1.0, "missing_explicit_donor_H");
                    });
                continue;
            }

            // sequence exclusion
            const int seq_sep =
                hbond_result_detail::BackboneResidueSeparation(
                    protein, donor_residue_idx, ri);
            if (seq_sep >= 0 &&
                seq_sep < static_cast<int>(CalculatorConfig::Get(
                    "hbond_sequential_exclusion_residues"))) {
                // record: seq-sep reject
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [ri, donor_residue_idx, seq_sep](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(donor_residue_idx), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "sequence_separation", static_cast<double>(seq_sep), "residues");
                        AddNumber(gc, "rejection", 1.0, "seq_sep_too_small");
                    });
                continue;
            }

            auto key = std::make_pair(don_res.N, res.O);
            if (seen.count(key)) continue;
            seen.insert(key);

            Vec3 H_pos = conf.PositionAt(don_res.H);
            Vec3 O_pos = conf.PositionAt(res.O);
            Vec3 d = O_pos - H_pos;
            double dist = d.norm();

            // source distance gate (singular or beyond max H...O)
            if (dist < CalculatorConfig::Get("singularity_guard_distance") || dist > CalculatorConfig::Get("hbond_dipolar_max_distance")) {
                // record: distance reject
                choices.Record(CalculatorId::HBond, choice_record_seq++, "hbond resolution",
                    [donor_residue_idx, ri, dist](GeometryChoice& gc) {
                        AddNumber(gc, "donor_residue", static_cast<double>(donor_residue_idx), "index");
                        AddNumber(gc, "acceptor_residue", static_cast<double>(ri), "index");
                        AddNumber(gc, "distance", dist, "A");
                        AddNumber(gc, "rejection", 1.0, "distance_out_of_range");
                    });
                continue;
            }

            ResolvedHBond hb;
            hb.donor_N = don_res.N;
            hb.donor_H = don_res.H;
            hb.acceptor_O = res.O;
            hb.donor_residue = donor_residue_idx;
            hb.acceptor_residue = ri;
            hb.source_point = H_pos;
            hb.h_hat = d / dist;
            hb.distance = dist;
            hb.sequence_separation = seq_sep;
            hbonds.push_back(hb);
        }
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "resolved " + std::to_string(hbonds.size()) +
        " unique backbone H-bonds from DSSP");

    result_ptr->hbond_count_ = hbonds.size();

    if (hbonds.empty()) {
        return result_ptr;
    }

    // ------------------------------------------------------------------
    // Step 2: Build the filter set for H-bond kernel evaluations.
    //
    // MinDistanceFilter: singularity guard near the donor H source point.
// SelfSourceFilter: atom cannot be a field point for an H-bond
    //   where it is the donor H or acceptor O.
// DipolarNearFieldFilter: point-source model invalid when field
    //   point is inside the H...O source distribution.
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    filters.Add(std::make_unique<SequentialExclusionFilter>());
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "filter set: " + filters.Describe());

    // ------------------------------------------------------------------
    // Step 3: For each atom, compute the dipolar tensor from all H-bonds
    // that pass the filter set. 1/r³ decay handles range naturally.
    // ------------------------------------------------------------------

    int total_pairs = 0;
    int filtered_out = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        // Residue index of this atom (for sequence separation)
        size_t ai_res = protein.AtomAt(ai).residue_index;

        double f_sum = 0.0;  // McConnell angular scalar, summed over contributing H-bonds
        double nearest_dist = NO_DATA_SENTINEL;
        size_t nearest_hb_idx = SIZE_MAX;
        int nearby_hbond_count = 0;

        for (size_t hi = 0; hi < hbonds.size(); ++hi) {
            const auto& hb = hbonds[hi];

            // kernel geometry (computed first so filters can read distance/extent)
            hbond_result_detail::KernelResult kernel =
                hbond_result_detail::ComputeKernel(
                    atom_pos, hb.source_point, hb.h_hat);

            // Build evaluation context from already-computed geometry
            KernelEvaluationContext ctx;
            ctx.distance = kernel.distance;
            ctx.source_extent = hb.distance;  // H...O distance
            ctx.atom_index = ai;
            ctx.source_atom_a = hb.donor_H;
            ctx.source_atom_b = hb.acceptor_O;

            // endpoint sequence gap (min to either donor/acceptor residue)
            const int sep_don =
                hbond_result_detail::BackboneResidueSeparation(
                    protein, ai_res, hb.donor_residue);
            const int sep_acc =
                hbond_result_detail::BackboneResidueSeparation(
                    protein, ai_res, hb.acceptor_residue);
            ctx.sequence_separation =
                hbond_result_detail::MinimumKnownSeparation(
                    sep_don, sep_acc);

            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: filter exclusion ----
                choices.Record(CalculatorId::HBond, hi, "filter exclusion",
                    [&](GeometryChoice& gc) {
                        AddAtom(gc, &conf.AtomAt(hb.donor_N), hb.donor_N,
                                EntityRole::Context, EntityOutcome::Included);
                        AddAtom(gc, &conf.AtomAt(hb.donor_H), hb.donor_H,
                                EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &conf.AtomAt(hb.acceptor_O), hb.acceptor_O,
                                EntityRole::Context, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", kernel.distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                filtered_out++;
                continue;
            }

            // count H-bonds within hbond_counting_radius (default 3.5 A) of this atom
            if (kernel.distance < CalculatorConfig::Get("hbond_counting_radius")) nearby_hbond_count++;

            // Track nearest (among accepted evaluations only)
            if (kernel.distance < nearest_dist) {
                nearest_dist = kernel.distance;
                nearest_hb_idx = hi;
            }

            f_sum += kernel.f;  // same accepted H-bond set as the tensor
            total_pairs++;
        }

        ca.hbond_count_within_3_5A = nearby_hbond_count;

        // ---- GeometryChoice: hbond neighbourhood ----
        choices.Record(CalculatorId::HBond, ai, "hbond neighbourhood",
            [&ca, ai, nearby_hbond_count](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "hbond_count_within_3_5A", static_cast<double>(nearby_hbond_count), "count");
            });

        if (nearest_hb_idx != SIZE_MAX) {
            const auto& nearest_hb = hbonds[nearest_hb_idx];
            ca.hbond_nearest_dist = nearest_dist;
            ca.hbond_nearest_dir =
                (atom_pos - nearest_hb.source_point).normalized();
            ca.hbond_is_backbone = true;  // all DSSP H-bonds are backbone

            // 1/r3 to the nearest donor-H source point (distinct from the
            // H...O source extent).
            ca.hbond_inv_d3 = 1.0 / (nearest_dist * nearest_dist * nearest_dist);
        }

        // endpoint flags: is this atom a donor N / acceptor O of any resolved bond
        for (const auto& hb : hbonds) {
            if (ai == hb.donor_N) ca.hbond_is_donor = true;
            if (ai == hb.acceptor_O) ca.hbond_is_acceptor = true;
        }

        ca.hbond_mcconnell_scalar = f_sum;
    }

    OperationLog::Info(LogCalcOther, "HBondResult::Compute",
        "atom_hbond_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " hbonds=" + std::to_string(hbonds.size()) +
        " atoms=" + std::to_string(n_atoms));

    return result_ptr;
}


int HBondResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> scalars(N * 4);
    std::vector<double> nearest_dir(N * 3);
    std::vector<int8_t> flags(N * 3, 0);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        scalars[i*4+0] = ca.hbond_nearest_dist;
        scalars[i*4+1] = ca.hbond_inv_d3;
        scalars[i*4+2] = static_cast<double>(ca.hbond_count_within_3_5A);
        scalars[i*4+3] = ca.hbond_mcconnell_scalar;
        nearest_dir[i*3+0] = ca.hbond_nearest_dir.x();
        nearest_dir[i*3+1] = ca.hbond_nearest_dir.y();
        nearest_dir[i*3+2] = ca.hbond_nearest_dir.z();
        flags[i*3 + 0] = ca.hbond_is_backbone ? 1 : 0;
        flags[i*3 + 1] = ca.hbond_is_donor ? 1 : 0;
        flags[i*3 + 2] = ca.hbond_is_acceptor ? 1 : 0;
    }

    NpyWriter::WriteFloat64(output_dir + "/hbond_scalars.npy", scalars.data(), N, 4);
    NpyWriter::WriteFloat64(output_dir + "/hbond_nearest_dir.npy",
                            nearest_dir.data(), N, 3);
    NpyWriter::WriteInt8(output_dir + "/hbond_flags.npy",
                         flags.data(), N, 3);
    return 3;
}

}  // namespace nmr
