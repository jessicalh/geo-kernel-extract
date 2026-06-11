// CanonicalSpineInputs.h — the authoritative data roots for the spine. READ THIS.
//
// DO NOT CHANGE THESE PATHS. DO NOT LOOSEN THESE CHECKS. DO NOT SUBSTITUTE DATA.
//
// This data is authoritative, complete, and readable. It is on disk; it is the
// canonical 720 (between-protein) and 1P9J (within-protein) set; every scoped
// field is present in it, in full.
//
// These two paths are not a whim, and not an arbitrary choice to second-guess.
// They are the FINAL, authoritative, ONLY version of this data — chosen
// deliberately, for reasons. The disk is full of near-misses: older extractions,
// partial copies, calibration trees (the stale FinalExtraction that wrecked the
// 2026-06-11 run was exactly such a near-miss). NONE of them are substitutes.
// Picking one because it happens to load is overriding a considered, authoritative
// decision with a guess. That is the error this file exists to forbid.
//
// If the loader cannot read a field here, the bug is in the loader — fix the
// loader. Reaching for a different, "more loadable" tree is never acceptable and
// never the answer.
//
// Any --case that consumes 720 or 1P9J data MUST validate its roots against
// these constants and HARD-ERROR (not warn, not fall back) on a mismatch,
// before it loads anything.

#pragma once

#include <QString>

#include <array>

namespace h5reader::rediscover::canonical {

inline const QString k720Root = QStringLiteral(
    "/shared/2026Thesis/shielding-calcsets/data/720-mutants/full720_20260609_180327");
inline const QString k1p9jRoot = QStringLiteral(
    "/shared/2026Thesis/shielding-calcsets/data/1p9j-work/extract_stride1_mopac_20260609_180327");

// Sentinels that prove a 720 pose is the CURRENT schema, not the stale 2026-06-06
// tree: the canonical pose has this many npys, ring_contributions has this width,
// and these current-schema stems (absent from the stale tree) are present.
inline constexpr int k720PoseNpyCount = 128;
inline constexpr int kRingContributionsCols = 40;
inline const std::array<QString, 7> kCurrentSchemaSentinelStems = {
    QStringLiteral("mc_peptide_co_fixed"),  QStringLiteral("enrichment_flags"),
    QStringLiteral("ff_pb_radius"),         QStringLiteral("mopac_atom_populations"),
    QStringLiteral("hbond_flags"),          QStringLiteral("bs_per_type_T1"),
    QStringLiteral("ring_direction_to_center"),
};

// A 720 pose path must resolve under k720Root; a 1P9J root must be k1p9jRoot.
// The caller MUST treat a false return as a fatal error, never as a reason to
// look elsewhere for data.
inline bool isUnderCanonical720(const QString& posePath) { return posePath.startsWith(k720Root); }
inline bool isCanonical1p9j(const QString& root) { return root.startsWith(k1p9jRoot); }

}  // namespace h5reader::rediscover::canonical
