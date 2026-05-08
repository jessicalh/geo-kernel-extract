#include "FramePdbEmitter.h"

#include "Atom.h"
#include "AminoAcidType.h"
#include "Bond.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Types.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace nmr {

// ============================================================================
// Singleton state -- file-local. Configure / OnFrame / Reset are the
// only access surface; nothing else needs to see this struct.
// ============================================================================

namespace {

struct EmitterState {
    bool                       active  = false;
    const Protein*             protein = nullptr;
    FramePdbEmitter::Config    config{};
};

EmitterState& State() {
    static EmitterState s;
    return s;
}


// Format an atom name into 4 chars, AMBER / PDB-v3 convention for the
// 1-char elements present in proteins (H, C, N, O, S):
//   - 4-char names (e.g. "HD11", "HD22") fill cols 13-16 verbatim.
//   - 1-3 char names get a leading space (col 13 = ' '), name in
//     cols 14-16 left-justified.
//
// Pre-condition: name is the canonical AMBER ff14SB atom name on
// `Atom.pdb_atom_name` (post-NamingApplicator). 2-char elements would
// occupy cols 13-14; we don't have any in the standard 20.
std::string FormatAtomName4(const std::string& name) {
    std::string out;
    out.reserve(4);
    if (name.size() >= 4) {
        out.assign(name.begin(), name.begin() + 4);
    } else {
        out.push_back(' ');
        out.append(name);
        while (out.size() < 4) out.push_back(' ');
    }
    return out;
}


// Element symbol right-justified in 2 chars, for cols 77-78.
const char* ElementSymbolRight2(Element e) {
    switch (e) {
        case Element::H:       return " H";
        case Element::C:       return " C";
        case Element::N:       return " N";
        case Element::O:       return " O";
        case Element::S:       return " S";
        case Element::Unknown: return "  ";
    }
    return "  ";
}


// Compute a/b/c/alpha/beta/gamma from a box matrix whose columns are
// the lattice vectors a, b, c (Å). Returns false on degenerate boxes.
bool BoxToCellParameters(const Eigen::Matrix3d& box,
                         double& a_len, double& b_len, double& c_len,
                         double& alpha_deg, double& beta_deg,
                         double& gamma_deg) {
    const Eigen::Vector3d a = box.col(0);
    const Eigen::Vector3d b = box.col(1);
    const Eigen::Vector3d c = box.col(2);
    a_len = a.norm();
    b_len = b.norm();
    c_len = c.norm();
    if (a_len < 1e-6 || b_len < 1e-6 || c_len < 1e-6) return false;
    constexpr double kRadToDeg = 180.0 / M_PI;
    // Clamp to [-1, 1] before acos for numerical safety.
    auto angle = [](double cosv) {
        return std::acos(std::min(1.0, std::max(-1.0, cosv))) * kRadToDeg;
    };
    alpha_deg = angle(b.dot(c) / (b_len * c_len));
    beta_deg  = angle(a.dot(c) / (a_len * c_len));
    gamma_deg = angle(a.dot(b) / (a_len * b_len));
    return true;
}


std::string MakeFilename(const FramePdbEmitter::Config& cfg,
                         std::size_t frame_idx, double time_ps) {
    char buf[512];
    if (cfg.decorator.empty()) {
        std::snprintf(buf, sizeof(buf), "%s_f%06zu_t%.1f.pdb",
                      cfg.stem.c_str(),
                      frame_idx, time_ps);
    } else {
        std::snprintf(buf, sizeof(buf), "%s_%s_f%06zu_t%.1f.pdb",
                      cfg.stem.c_str(),
                      cfg.decorator.c_str(),
                      frame_idx, time_ps);
    }
    return std::string(buf);
}


void WriteCryst1IfPresent(std::ofstream& out,
                          const Eigen::Matrix3d* box_matrix) {
    if (!box_matrix) return;
    if (box_matrix->isZero(1e-6)) return;
    double a, b, c, alpha, beta, gamma;
    if (!BoxToCellParameters(*box_matrix, a, b, c, alpha, beta, gamma)) return;
    char buf[96];
    std::snprintf(buf, sizeof(buf),
        "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
        a, b, c, alpha, beta, gamma);
    out << buf;
}


// Emit one ATOM record. `serial` is 1-indexed.
void WriteAtomLine(std::ofstream& out,
                   int serial,
                   const Atom& atom,
                   const Residue& res,
                   const Vec3& pos) {
    const std::string atom_name4 = FormatAtomName4(atom.pdb_atom_name);
    const char* res3 = res.AminoAcidInfo().three_letter_code;
    const char chain_char =
        res.chain_id.empty() ? 'A' : res.chain_id.front();
    const char ins_char =
        res.insertion_code.empty() ? ' ' : res.insertion_code.front();
    const char* elem2 = ElementSymbolRight2(atom.element);

    char line[96];
    std::snprintf(line, sizeof(line),
        "ATOM  %5d %4s %-3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
        serial,
        atom_name4.c_str(),
        res3 ? res3 : "UNK",
        chain_char,
        res.sequence_number,
        ins_char,
        pos.x(), pos.y(), pos.z(),
        1.00,  // occupancy
        0.00,  // B-factor
        elem2);
    out << line;
}


// Emit a TER record for the residue that ended the previous chain.
void WriteTer(std::ofstream& out, int serial, const Residue& res) {
    const char* res3 = res.AminoAcidInfo().three_letter_code;
    const char chain_char =
        res.chain_id.empty() ? 'A' : res.chain_id.front();
    const char ins_char =
        res.insertion_code.empty() ? ' ' : res.insertion_code.front();
    char line[96];
    std::snprintf(line, sizeof(line),
        "TER   %5d      %-3s %c%4d%c\n",
        serial,
        res3 ? res3 : "UNK",
        chain_char,
        res.sequence_number,
        ins_char);
    out << line;
}


// Body of one PDB emission. Called only when active and gates passed.
void EmitOnePdb(const Protein& protein,
                const ProteinConformation& conf,
                std::size_t frame_idx,
                double time_ps,
                const Eigen::Matrix3d* box_matrix) {
    const auto& cfg = State().config;

    std::error_code ec;
    std::filesystem::create_directories(cfg.output_dir, ec);
    if (ec) {
        OperationLog::Error("FramePdbEmitter",
            "create_directories failed for " + cfg.output_dir.string() +
            ": " + ec.message());
        return;
    }

    const std::filesystem::path path =
        cfg.output_dir / MakeFilename(cfg, frame_idx, time_ps);

    std::ofstream out(path);
    if (!out) {
        OperationLog::Error("FramePdbEmitter",
            "failed to open " + path.string() + " for write");
        return;
    }

    // ── HEADER + REMARK ─────────────────────────────────────────────
    {
        char buf[128];
        std::snprintf(buf, sizeof(buf),
            "HEADER    FRAME PDB FROM nmr_extract                       %s\n",
            cfg.stem.c_str());
        out << buf;
        std::snprintf(buf, sizeof(buf),
            "REMARK   1 frame_idx=%zu time_ps=%.3f\n",
            frame_idx, time_ps);
        out << buf;
        std::snprintf(buf, sizeof(buf),
            "REMARK   1 stem=%s\n", cfg.stem.c_str());
        out << buf;
        if (!cfg.decorator.empty()) {
            std::snprintf(buf, sizeof(buf),
                "REMARK   1 decorator=%s\n", cfg.decorator.c_str());
            out << buf;
        }
    }

    // ── CRYST1 if box matrix present ───────────────────────────────
    WriteCryst1IfPresent(out, box_matrix);

    // ── ATOM records + TER between biological chains ──────────────
    //
    // Chain break = residue.chain_id changes between consecutive atoms.
    // Emit TER for the LAST residue of the previous chain when a break
    // is detected. Emit a final TER after the last atom.
    const std::size_t natoms = protein.AtomCount();
    int serial = 0;
    std::string current_chain;
    std::size_t last_residue_index = SIZE_MAX;

    for (std::size_t ai = 0; ai < natoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(atom.residue_index);
        const std::string& this_chain =
            res.chain_id.empty() ? std::string("A") : res.chain_id;

        // Chain break: TER for the previous residue.
        if (last_residue_index != SIZE_MAX && this_chain != current_chain) {
            ++serial;
            const Residue& last_res = protein.ResidueAt(last_residue_index);
            WriteTer(out, serial, last_res);
        }
        current_chain = this_chain;
        last_residue_index = atom.residue_index;

        ++serial;
        WriteAtomLine(out, serial, atom, res, conf.AtomAt(ai).Position());
    }

    // Final TER after the last atom.
    if (last_residue_index != SIZE_MAX) {
        ++serial;
        const Residue& last_res = protein.ResidueAt(last_residue_index);
        WriteTer(out, serial, last_res);
    }

    // ── CONECT records: disulfides only ────────────────────────────
    //
    // MolProbity prefers to detect routine bonds itself; explicit CONECT
    // for disulfides is the convention because they cross sequence
    // distance and aren't inferable from proximity alone in every case.
    const auto& bonds = protein.LegacyAmber().BondList();
    for (const Bond& b : bonds) {
        if (b.category != BondCategory::Disulfide) continue;
        char line[64];
        std::snprintf(line, sizeof(line),
            "CONECT%5zu%5zu\n",
            b.atom_index_a + 1, b.atom_index_b + 1);
        out << line;
    }

    out << "END\n";
}

}  // namespace


// ============================================================================
// FramePdbEmitter -- public surface
// ============================================================================

void FramePdbEmitter::Configure(const Protein& protein, Config config) {
    auto& s = State();
    s.protein = &protein;
    s.config  = std::move(config);
    s.active  = !s.config.output_dir.empty();

    if (s.active) {
        std::string msg = "configured: dir=" + s.config.output_dir.string() +
            " stem=" + s.config.stem +
            " stride=" + std::to_string(s.config.stride);
        if (!s.config.decorator.empty()) {
            msg += " decorator=" + s.config.decorator;
        }
        OperationLog::Info(LogCalcOther, "FramePdbEmitter", msg);
    } else {
        OperationLog::Warn("FramePdbEmitter",
            "Configure called with empty output_dir; emitter remains inactive");
    }
}


void FramePdbEmitter::OnFrame(const ProteinConformation& conf,
                              std::size_t frame_idx,
                              double time_ps,
                              const Eigen::Matrix3d* box_matrix) {
    auto& s = State();
    if (!s.active) return;
    if (!s.protein) return;
    if (time_ps < s.config.from_ps || time_ps >= s.config.to_ps) return;
    if (s.config.stride > 1 && (frame_idx % s.config.stride) != 0) return;
    EmitOnePdb(*s.protein, conf, frame_idx, time_ps, box_matrix);
}


void FramePdbEmitter::Reset() {
    auto& s = State();
    s.active  = false;
    s.protein = nullptr;
    s.config  = Config{};
}


bool FramePdbEmitter::IsActive() {
    return State().active;
}

}  // namespace nmr
