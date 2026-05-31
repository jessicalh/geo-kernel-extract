// CalcsetManifest — the typed wrapper over a calcset's `.LGS` (Lowly
// Graduate Student) file. One per calcset directory; the consumer's
// entry point into the artifact tree.
//
// Authoritative spec: spec/CALCSET_MANIFEST.md (schema v1).
//
// The `.LGS` carries the top-level identity, the `kind` dispatch
// (trajectory / single_pose / mutant_pair), the artifact-pointer
// table, and — when DFT exists for this calcset — the typed
// `frame_index → meta_json` map. It is the index, not a parser of the
// artifacts it points at.
//
// Loader semantics — see CALCSET_MANIFEST.md §"Loader behaviour":
//   1. The argument may be a directory (look for the single `*.LGS`
//      inside; zero or > 1 matches → hard error) or a `.LGS` file path.
//   2. `schema_version` must equal `kSupportedSchemaVersion` (== 1);
//      mismatches reject with a clear error.
//   3. Required keys for the active `kind` must be present; every
//      declared file/dir path must exist on disk with the right kind.
//   4. Optional sub-blocks absent → silently absent; optional paths
//      declared but missing on disk → hard error (manifest lies).
//
// No exceptions cross the loader boundary: failures return
// `std::nullopt` and write a human-readable message to `err_out`.
// CalcsetManifest is a plain struct, NOT a QObject — no CENSUS, no
// signals/slots; per-frame meta lookups are pure functions on
// resolved paths.

#pragma once

#include <QString>

#include <cstdint>
#include <optional>
#include <vector>

namespace h5reader::io {

/**
 * One DFT job entry from the `.LGS`'s `dft.frames[]` array.
 *
 * Carries the typed frame index + the resolved abs path to the job's
 * `meta_json`. The meta.json is the per-frame truth — `framePs()`,
 * `orcaOutAbspath()`, `orcaExitCode()` lazy-load it on demand so the
 * common case (loader iterates frames, callers map to a hash) doesn't
 * touch every per-frame file.
 */
struct DftFrame {
    std::int32_t frame_index = 0;
    QString      meta_json_relpath;   // as written in the .LGS
    QString      meta_json_abspath;   // resolved against calcset_root

    // Read meta.json once on the calling thread and cache the parsed
    // fields needed for downstream wiring. Returns true on success;
    // writes a message to `err_out` (if non-null) on parse failure.
    bool LoadMeta(QString* err_out = nullptr) const;

    // Lazy accessors — call LoadMeta() first; otherwise they return
    // sentinel values (frame_ps_=NAN, etc.). Returning sentinels keeps
    // callers Qt-citizen (no exceptions across the boundary); strict
    // callers check LoadMeta() ok=true first.
    double  framePs() const;
    QString orcaOutAbspath() const;
    int     orcaExitCode() const;

private:
    // Mutable so const accessors can fill on demand. Honest about cache.
    mutable bool    meta_loaded_     = false;
    mutable double  frame_ps_        = 0.0;
    mutable QString orca_out_abspath_;
    mutable int     orca_exit_code_  = -1;
};

/**
 * Typed view of a `.LGS` calcset manifest (schema v1). Each `kind`
 * has its own sub-struct; only the one matching `kind` is populated.
 *
 * Resolved paths are all absolute (calcset_root_abspath joined with
 * the relative path written in the .LGS).
 */
struct CalcsetManifest {
    enum class Kind     { Trajectory, SinglePose, MutantPair };
    enum class PoseKind { Pdb, ProtonatedPdb, Orca };

    int     schema_version = 0;
    Kind    kind = Kind::Trajectory;
    QString dataset_id;
    QString protein_id;
    QString human_name;

    struct Trajectory {
        QString md_dir_abspath;
        QString topology_top_abspath;
        QString extraction_dir_abspath;
        QString trajectory_h5_abspath;
        QString extraction_manifest_abspath;
        double  frame_dt_ps = 0.0;
        QString frame_index_basis;
        QString reference_pdb_abspath;  // optional; empty if absent
    };
    struct SinglePose {
        PoseKind pose_kind = PoseKind::Pdb;
        QString  pose_dir_abspath;
        QString  extraction_manifest_abspath;
    };
    struct MutantPair {
        QString wt_lgs_abspath;
        QString ala_lgs_abspath;
    };
    struct Dft {
        QString method;
        int     campaign_target_frames = 0;
        struct Stride {
            int first = 0;
            int last  = 0;
            int step  = 1;
        };
        Stride                 frame_stride;
        std::vector<DftFrame>  frames;
    };

    std::optional<Trajectory> trajectory;
    std::optional<SinglePose> single_pose;
    std::optional<MutantPair> mutant_pair;
    std::optional<Dft>        dft;

    QString calcset_root_abspath;  // resolved at load time
    QString lgs_path_abspath;      // for diagnostics

    QString generated_at_utc;
    QString lgs_writer;
    QString producer_extractor_version;

    static constexpr int        kSupportedSchemaVersion = 1;
    static constexpr const char* kExtension = ".LGS";

    /**
     * Load + validate a `.LGS` from a calcset directory or directly
     * from a `.LGS` path.
     *
     * @param root_or_lgs_path
     *   A directory holding exactly one `*.LGS` file, or the `.LGS`
     *   file itself. Zero or > 1 `.LGS` files in a directory is a hard
     *   error — no glob-and-pick.
     * @param err_out
     *   On `nullopt` return, receives a human-readable error message
     *   (kept short enough for an error dialog). Pass `nullptr` to
     *   discard.
     * @return
     *   The populated manifest on success; `std::nullopt` on any
     *   failure. ErrorBus is also notified.
     */
    static std::optional<CalcsetManifest>
    Load(const QString& root_or_lgs_path, QString* err_out = nullptr);

    // Stringify Kind for log/error messages.
    static const char* NameForKind(Kind k);
};

}  // namespace h5reader::io
