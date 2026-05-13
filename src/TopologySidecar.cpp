#include "TopologySidecar.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "Bond.h"
#include "CovalentTopology.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "Protein.h"
#include "Residue.h"
#include "Ring.h"
#include "RingTopology.h"
#include "Types.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// Internal helpers (anonymous namespace)
// ============================================================================

namespace {

// ── NPY 1.0 writer ────────────────────────────────────────────────
//
// Same layout as CategoryInfoProjection.cpp::WriteStructuredNpy. Kept
// local rather than refactored into a shared utility (PATTERNS.md
// "Duplication is preferred over chaining" + "The utility namespace"
// anti-pattern) — three independent emitters with stable shapes.
//
bool WriteStructuredNpy(const fs::path& path,
                          const char* dtype_descr,
                          size_t n_records,
                          size_t bytes_per_record,
                          const std::vector<unsigned char>& bytes) {
    if (bytes.size() != n_records * bytes_per_record) {
        std::fprintf(stderr,
            "FATAL: TopologySidecar -- buffer size mismatch for %s: "
            "%zu bytes for %zu records (expected %zu).\n",
            path.string().c_str(), bytes.size(), n_records,
            n_records * bytes_per_record);
        std::abort();
    }

    std::ostringstream hdr;
    hdr << "{'descr': " << dtype_descr
        << ", 'fortran_order': False, 'shape': (" << n_records << ",), }";
    std::string header = hdr.str();

    constexpr size_t kPreHeader = 10;
    size_t total = kPreHeader + header.size() + 1;
    size_t pad = (64 - (total % 64)) % 64;
    header.append(pad, ' ');
    header.push_back('\n');

    if (header.size() > 0xFFFF) {
        std::fprintf(stderr,
            "FATAL: TopologySidecar -- NPY header too large for v1.0 "
            "format (%zu > 65535) for %s.\n",
            header.size(), path.string().c_str());
        std::abort();
    }
    const uint16_t header_len = static_cast<uint16_t>(header.size());

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        OperationLog::Error("TopologySidecar",
            "could not open " + path.string() + " for write");
        return false;
    }
    const char magic[] = {'\x93','N','U','M','P','Y'};
    out.write(magic, 6);
    const char ver[] = {'\x01','\x00'};
    out.write(ver, 2);
    out.write(reinterpret_cast<const char*>(&header_len), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char*>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
    return out.good();
}

// ── Residues emission ─────────────────────────────────────────────
//
// One row per Residue. Codex contract item: "Residue Table Minimum"
// plus prev/next links and Markley-style 1-letter / 3-letter
// renderings.
//
// Chain boundary detection: row i is at chain start if i == 0 OR
// chain_id changes from i-1 to i. Same for chain end with i+1.
//
constexpr const char* kResiduesDtypeDescr =
    "[('residue_index', '<i4'),"
    " ('chain_id', '|S2'),"
    " ('residue_number', '<i4'),"
    " ('insertion_code', '|S1'),"
    " ('residue_type', 'i1'),"
    " ('amber_residue_3letter', '|S4'),"
    " ('iupac_residue_3letter', '|S4'),"
    " ('one_letter', '|S1'),"
    " ('protonation_variant_index', 'i1'),"
    " ('terminal_state', 'i1'),"
    " ('prev_residue_index', '<i4'),"
    " ('next_residue_index', '<i4'),"
    " ('prev_residue_type', 'i1'),"
    " ('next_residue_type', 'i1'),"
    " ('atom_count', '<i4'),"
    " ('is_proline', 'i1'),"
    " ('is_aromatic', 'i1'),"
    " ('is_titratable', 'i1'),"
    " ('has_amide_h', 'i1'),"
    " ('is_xpro_context', 'i1')]";

constexpr size_t kResidueRecordSize =
    4 + 2 + 4 + 1 + 1 + 4 + 4 + 1 + 1 + 1 + 4 + 4 + 1 + 1 + 4 + 1 + 1 + 1 + 1 + 1;

void PackFixedString(unsigned char* dst, size_t dst_size,
                       const std::string& src,
                       const char* field_name, size_t row_index) {
    if (src.size() > dst_size) {
        std::fprintf(stderr,
            "FATAL: TopologySidecar -- field \"%s\" length %zu exceeds "
            "column width %zu at row %zu (value=\"%s\").\n",
            field_name, src.size(), dst_size, row_index, src.c_str());
        std::abort();
    }
    std::memset(dst, 0, dst_size);
    std::memcpy(dst, src.data(), src.size());
}

bool WriteResidues(const Protein& protein, const fs::path& out_dir) {
    const size_t N = protein.ResidueCount();
    std::vector<unsigned char> buf(N * kResidueRecordSize, 0);

    for (size_t ri = 0; ri < N; ++ri) {
        const Residue& r = protein.ResidueAt(ri);
        unsigned char* row = buf.data() + ri * kResidueRecordSize;
        size_t off = 0;

        const int32_t r_idx = static_cast<int32_t>(ri);
        std::memcpy(row + off, &r_idx, 4); off += 4;

        PackFixedString(row + off, 2, r.chain_id, "chain_id", ri);
        off += 2;

        const int32_t seq_num = static_cast<int32_t>(r.sequence_number);
        std::memcpy(row + off, &seq_num, 4); off += 4;

        PackFixedString(row + off, 1, r.insertion_code, "insertion_code", ri);
        off += 1;

        row[off++] = static_cast<int8_t>(r.type);

        const AminoAcidType& aat = (r.type == AminoAcid::Unknown)
            ? GetAminoAcidType(AminoAcid::ALA)  // dummy to satisfy reference
            : GetAminoAcidType(r.type);

        std::string amber_3 = (r.type == AminoAcid::Unknown) ? "UNK"
            : (r.protonation_variant_index >= 0
                  && static_cast<size_t>(r.protonation_variant_index) < aat.variants.size()
                  ? std::string(aat.variants[r.protonation_variant_index].name)
                  : std::string(aat.three_letter_code));
        std::string iupac_3 = (r.type == AminoAcid::Unknown) ? "UNK"
            : std::string(aat.three_letter_code);
        std::string one_letter(1,
            r.type == AminoAcid::Unknown ? 'X' : aat.one_letter_code);

        PackFixedString(row + off, 4, amber_3, "amber_residue_3letter", ri);
        off += 4;
        PackFixedString(row + off, 4, iupac_3, "iupac_residue_3letter", ri);
        off += 4;
        PackFixedString(row + off, 1, one_letter, "one_letter", ri);
        off += 1;

        row[off++] = static_cast<int8_t>(r.protonation_variant_index);
        row[off++] = static_cast<int8_t>(r.terminal_state);

        // Chain-aware prev/next links. -1 at chain boundaries.
        int32_t prev_idx = -1;
        int8_t  prev_type = static_cast<int8_t>(AminoAcid::Unknown);
        if (ri > 0) {
            const Residue& prev = protein.ResidueAt(ri - 1);
            if (prev.chain_id == r.chain_id) {
                prev_idx = static_cast<int32_t>(ri - 1);
                prev_type = static_cast<int8_t>(prev.type);
            }
        }
        std::memcpy(row + off, &prev_idx, 4); off += 4;

        int32_t next_idx = -1;
        int8_t  next_type = static_cast<int8_t>(AminoAcid::Unknown);
        if (ri + 1 < N) {
            const Residue& next = protein.ResidueAt(ri + 1);
            if (next.chain_id == r.chain_id) {
                next_idx = static_cast<int32_t>(ri + 1);
                next_type = static_cast<int8_t>(next.type);
            }
        }
        std::memcpy(row + off, &next_idx, 4); off += 4;
        row[off++] = prev_type;
        row[off++] = next_type;

        const int32_t atom_count = static_cast<int32_t>(r.atom_indices.size());
        std::memcpy(row + off, &atom_count, 4); off += 4;

        row[off++] = (r.type == AminoAcid::PRO) ? 1 : 0;
        row[off++] = (r.type != AminoAcid::Unknown && r.IsAromatic()) ? 1 : 0;
        row[off++] = (r.type != AminoAcid::Unknown && r.IsTitratable()) ? 1 : 0;
        row[off++] = (r.type != AminoAcid::Unknown && r.HasAmideH()) ? 1 : 0;
        // is_xpro_context: this residue's i+1 is PRO. Relevant for the
        // (i, i+1) peptide-bond omega — X→Pro permits cis isomerism that
        // standard non-Pro context does not. -1 next means chain end.
        const bool x_pro =
            (next_idx >= 0 && next_type == static_cast<int8_t>(AminoAcid::PRO));
        row[off++] = x_pro ? 1 : 0;

        if (off != kResidueRecordSize) {
            std::fprintf(stderr,
                "FATAL: TopologySidecar::WriteResidues record-size mismatch "
                "at row %zu: wrote %zu, expected %zu.\n",
                ri, off, kResidueRecordSize);
            std::abort();
        }
    }
    return WriteStructuredNpy(out_dir / "residues.npy", kResiduesDtypeDescr,
                                N, kResidueRecordSize, buf);
}

// ── Bonds emission ────────────────────────────────────────────────
//
// One row per Bond from LegacyAmberTopology::BondList(). Row i has
// bond_index == i (bond rows are listed in topology order; this is
// the canonical bond axis).
//
constexpr const char* kBondsDtypeDescr =
    "[('bond_index', '<i4'),"
    " ('atom_index_a', '<i4'),"
    " ('atom_index_b', '<i4'),"
    " ('bond_order', 'i1'),"
    " ('bond_category', 'i1'),"
    " ('is_rotatable', 'i1'),"
    " ('is_aromatic', 'i1'),"
    " ('is_peptide', 'i1'),"
    " ('is_backbone', 'i1')]";

constexpr size_t kBondRecordSize = 4 + 4 + 4 + 1 + 1 + 1 + 1 + 1 + 1;

bool WriteBonds(const Protein& protein, const fs::path& out_dir) {
    const auto& bonds = protein.LegacyAmber().BondList();
    const size_t N = bonds.size();
    std::vector<unsigned char> buf(N * kBondRecordSize, 0);

    for (size_t bi = 0; bi < N; ++bi) {
        const Bond& b = bonds[bi];
        unsigned char* row = buf.data() + bi * kBondRecordSize;
        size_t off = 0;
        const int32_t bi32 = static_cast<int32_t>(bi);
        std::memcpy(row + off, &bi32, 4); off += 4;
        const int32_t aa = static_cast<int32_t>(b.atom_index_a);
        std::memcpy(row + off, &aa, 4); off += 4;
        const int32_t bb = static_cast<int32_t>(b.atom_index_b);
        std::memcpy(row + off, &bb, 4); off += 4;
        row[off++] = static_cast<int8_t>(b.order);
        row[off++] = static_cast<int8_t>(b.category);
        row[off++] = b.is_rotatable ? 1 : 0;
        row[off++] = b.IsAromatic() ? 1 : 0;
        row[off++] = b.IsPeptideBond() ? 1 : 0;
        row[off++] = b.IsBackbone() ? 1 : 0;
    }
    return WriteStructuredNpy(out_dir / "bonds.npy", kBondsDtypeDescr,
                                N, kBondRecordSize, buf);
}

// ── Rings emission ────────────────────────────────────────────────
//
// One row per Ring. Aromatic rings come first, in aromatic-axis order
// (matching the existing ring_geometry.npy row order); saturated rings
// follow. ring_id is the absolute row index in rings.npy.
// native_axis_index is the index within the aromatic-only or
// saturated-only axis -- aromatic rows have native_axis_index ==
// ring_id, saturated rows have native_axis_index == ring_id -
// aromatic_count.
//
// fused_partner_ring_id is -1 when not fused, otherwise the absolute
// ring_id of the partner.
//
constexpr const char* kRingsDtypeDescr =
    "[('ring_id', '<i4'),"
    " ('ring_kind', 'i1'),"
    " ('ring_type_index', 'i1'),"
    " ('atom_count', 'i1'),"
    " ('_pad0', 'i1'),"
    " ('native_axis_index', '<i4'),"
    " ('parent_residue_index', '<i4'),"
    " ('parent_residue_number', '<i4'),"
    " ('fused_partner_ring_id', '<i4')]";

constexpr size_t kRingRecordSize = 4 + 1 + 1 + 1 + 1 + 4 + 4 + 4 + 4;

bool WriteRings(const Protein& protein, const fs::path& out_dir,
                size_t* out_total_membership_rows) {
    const RingTopology& rt = protein.LegacyAmber().Rings();
    const size_t n_arom = rt.AromaticCount();
    const size_t n_satd = rt.SaturatedCount();
    const size_t N = n_arom + n_satd;

    std::vector<unsigned char> buf(N * kRingRecordSize, 0);

    // Build aromatic-native to absolute ring_id map for fused-partner
    // resolution. The Ring struct's fused_partner_index is the
    // aromatic-axis index (SIZE_MAX if not fused) per Ring.h:57.
    // We translate it to absolute ring_id (== aromatic-axis index here).
    size_t membership_rows = 0;

    for (size_t ai = 0; ai < n_arom; ++ai) {
        const Ring& r = rt.AromaticAt(ai);
        unsigned char* row = buf.data() + ai * kRingRecordSize;
        size_t off = 0;
        const int32_t ring_id = static_cast<int32_t>(ai);
        std::memcpy(row + off, &ring_id, 4); off += 4;
        row[off++] = 0;  // ring_kind: 0 = aromatic
        row[off++] = static_cast<int8_t>(r.type_index);
        row[off++] = static_cast<int8_t>(r.atom_indices.size());
        row[off++] = 0;  // _pad0
        const int32_t native_idx = static_cast<int32_t>(ai);
        std::memcpy(row + off, &native_idx, 4); off += 4;
        const int32_t parent_idx = static_cast<int32_t>(r.parent_residue_index);
        std::memcpy(row + off, &parent_idx, 4); off += 4;
        const int32_t parent_num = static_cast<int32_t>(r.parent_residue_number);
        std::memcpy(row + off, &parent_num, 4); off += 4;
        const int32_t fp = (r.fused_partner_index == SIZE_MAX)
            ? -1
            : static_cast<int32_t>(r.fused_partner_index);
        std::memcpy(row + off, &fp, 4); off += 4;

        membership_rows += r.atom_indices.size();
    }

    for (size_t si = 0; si < n_satd; ++si) {
        const Ring& r = rt.SaturatedAt(si);
        const size_t absolute_id = n_arom + si;
        unsigned char* row = buf.data() + absolute_id * kRingRecordSize;
        size_t off = 0;
        const int32_t ring_id = static_cast<int32_t>(absolute_id);
        std::memcpy(row + off, &ring_id, 4); off += 4;
        row[off++] = 1;  // ring_kind: 1 = saturated
        row[off++] = static_cast<int8_t>(r.type_index);
        row[off++] = static_cast<int8_t>(r.atom_indices.size());
        row[off++] = 0;  // _pad0
        const int32_t native_idx = static_cast<int32_t>(si);
        std::memcpy(row + off, &native_idx, 4); off += 4;
        const int32_t parent_idx = static_cast<int32_t>(r.parent_residue_index);
        std::memcpy(row + off, &parent_idx, 4); off += 4;
        const int32_t parent_num = static_cast<int32_t>(r.parent_residue_number);
        std::memcpy(row + off, &parent_num, 4); off += 4;
        // Saturated rings don't carry fused-partner info (Ring.h documents
        // fused_partner_index as aromatic-only; Pro pyrrolidine is never
        // fused in our chemistry).
        const int32_t fp = -1;
        std::memcpy(row + off, &fp, 4); off += 4;

        membership_rows += r.atom_indices.size();
    }

    *out_total_membership_rows = membership_rows;
    return WriteStructuredNpy(out_dir / "rings.npy", kRingsDtypeDescr,
                                N, kRingRecordSize, buf);
}

// ── Ring membership emission ──────────────────────────────────────
//
// One row per (ring, atom) pair, in ring-walk order. All atoms in our
// Ring model are vertices (Ring::atom_indices stores the cyclic walk);
// is_substituent is reserved for future extension and currently 0.
//
constexpr const char* kRingMembershipDtypeDescr =
    "[('ring_id', '<i4'),"
    " ('atom_index', '<i4'),"
    " ('ring_atom_order', 'i1'),"
    " ('is_vertex', 'i1'),"
    " ('is_substituent', 'i1'),"
    " ('_pad0', 'i1')]";

constexpr size_t kRingMembershipRecordSize = 4 + 4 + 1 + 1 + 1 + 1;

bool WriteRingMembership(const Protein& protein, const fs::path& out_dir,
                          size_t total_rows) {
    const RingTopology& rt = protein.LegacyAmber().Rings();
    const size_t n_arom = rt.AromaticCount();
    const size_t n_satd = rt.SaturatedCount();

    std::vector<unsigned char> buf(total_rows * kRingMembershipRecordSize, 0);
    size_t cursor = 0;

    auto write_ring = [&](const Ring& r, size_t absolute_ring_id) {
        for (size_t k = 0; k < r.atom_indices.size(); ++k) {
            unsigned char* row = buf.data() + cursor * kRingMembershipRecordSize;
            size_t off = 0;
            const int32_t rid = static_cast<int32_t>(absolute_ring_id);
            std::memcpy(row + off, &rid, 4); off += 4;
            const int32_t aix = static_cast<int32_t>(r.atom_indices[k]);
            std::memcpy(row + off, &aix, 4); off += 4;
            row[off++] = static_cast<int8_t>(k);
            row[off++] = 1;  // is_vertex
            row[off++] = 0;  // is_substituent (reserved)
            row[off++] = 0;  // _pad0
            ++cursor;
        }
    };

    for (size_t ai = 0; ai < n_arom; ++ai) {
        write_ring(rt.AromaticAt(ai), ai);
    }
    for (size_t si = 0; si < n_satd; ++si) {
        write_ring(rt.SaturatedAt(si), n_arom + si);
    }

    if (cursor != total_rows) {
        std::fprintf(stderr,
            "FATAL: TopologySidecar::WriteRingMembership -- wrote %zu "
            "rows, expected %zu.\n", cursor, total_rows);
        std::abort();
    }

    return WriteStructuredNpy(out_dir / "ring_membership.npy",
                                kRingMembershipDtypeDescr,
                                total_rows, kRingMembershipRecordSize, buf);
}

// ── Manifest emission ─────────────────────────────────────────────
//
// Hand-written JSON; format-style direct to match the project's "no
// utility namespace, no abstractions" rule. Topology summary + axis
// sizes + axis-alignment statements.
//
// The axis-alignment block makes explicit what is implicit in the
// extraction layout: all atom-axis NPYs share row order, ring_geometry
// matches the aromatic prefix of rings, ring_contributions references
// the aromatic-ring axis.
//
std::string Iso8601UtcNow() {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const auto t = system_clock::to_time_t(now);
    std::tm tm{};
    gmtime_r(&t, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf);
}

bool WriteManifest(const Protein& protein, const fs::path& out_dir,
                    const std::string& protein_id,
                    size_t bond_count, size_t aromatic_ring_count,
                    size_t saturated_ring_count, size_t ring_membership_count) {
    const LegacyAmberTopology& lat = protein.LegacyAmber();
    const bool has_substrate = lat.HasAtomSemantic();
    const bool has_ff_types  = !lat.AtomtypeString().empty();
    const bool has_mass      = !lat.Mass().empty();

    std::ostringstream j;
    j << "{\n";
    j << "  \"schema_version\": \"1.0\",\n";
    j << "  \"extractor\": \"nmr_extract\",\n";
    j << "  \"extractor_version\": \"0.2.0\",\n";
    j << "  \"generated_at_utc\": \"" << Iso8601UtcNow() << "\",\n";
    j << "  \"protein_id\": \"" << protein_id << "\",\n";
    j << "  \"topology\": {\n";
    j << "    \"source\": \"amber-ff14SB+cifpp\",\n";
    j << "    \"has_atom_semantic\": " << (has_substrate ? "true" : "false") << ",\n";
    j << "    \"has_ff_atom_types\": " << (has_ff_types ? "true" : "false") << ",\n";
    j << "    \"has_ff_mass\": " << (has_mass ? "true" : "false") << "\n";
    j << "  },\n";
    // enum_vocab — integer-to-name mapping for the typed enum columns in
    // residues.npy / bonds.npy / rings.npy so a consumer reading
    // `bond_order == 4` knows it's `Peptide` without grepping the source.
    // Per codex contract `enum_vocab_refs`. Mirrors:
    //   ResidueTerminalState (Residue.h), BondOrder / BondCategory /
    //   RingTypeIndex / AminoAcid (Types.h).
    j << "  \"enum_vocab\": {\n";
    j << "    \"terminal_state\": {\"0\":\"Internal\",\"1\":\"NTerminus\","
      <<                          "\"2\":\"CTerminus\",\"3\":\"NAndCTerminus\","
      <<                          "\"4\":\"Unknown\"},\n";
    j << "    \"bond_order\": {\"0\":\"Single\",\"1\":\"Double\",\"2\":\"Triple\","
      <<                       "\"3\":\"Aromatic\",\"4\":\"Peptide\",\"5\":\"Unknown\"},\n";
    j << "    \"bond_category\": {\"0\":\"PeptideCO\",\"1\":\"PeptideCN\","
      <<                          "\"2\":\"BackboneOther\",\"3\":\"SidechainCO\","
      <<                          "\"4\":\"Aromatic\",\"5\":\"Disulfide\","
      <<                          "\"6\":\"SidechainOther\",\"7\":\"Unknown\"},\n";
    j << "    \"ring_kind\": {\"0\":\"aromatic\",\"1\":\"saturated\"},\n";
    j << "    \"ring_type_index\": {\"0\":\"PheBenzene\",\"1\":\"TyrPhenol\","
      <<                            "\"2\":\"TrpBenzene\",\"3\":\"TrpPyrrole\","
      <<                            "\"4\":\"TrpPerimeter\",\"5\":\"HisImidazole\","
      <<                            "\"6\":\"HidImidazole\",\"7\":\"HieImidazole\","
      <<                            "\"8\":\"ProPyrrolidine\"},\n";
    j << "    \"residue_type\": {\"0\":\"ALA\",\"1\":\"ARG\",\"2\":\"ASN\","
      <<                         "\"3\":\"ASP\",\"4\":\"CYS\",\"5\":\"GLN\","
      <<                         "\"6\":\"GLU\",\"7\":\"GLY\",\"8\":\"HIS\","
      <<                         "\"9\":\"ILE\",\"10\":\"LEU\",\"11\":\"LYS\","
      <<                         "\"12\":\"MET\",\"13\":\"PHE\",\"14\":\"PRO\","
      <<                         "\"15\":\"SER\",\"16\":\"THR\",\"17\":\"TRP\","
      <<                         "\"18\":\"TYR\",\"19\":\"VAL\",\"20\":\"Unknown\"}\n";
    j << "  },\n";
    j << "  \"axis_sizes\": {\n";
    j << "    \"atom\": " << protein.AtomCount() << ",\n";
    j << "    \"residue\": " << protein.ResidueCount() << ",\n";
    j << "    \"bond\": " << bond_count << ",\n";
    j << "    \"aromatic_ring\": " << aromatic_ring_count << ",\n";
    j << "    \"saturated_ring\": " << saturated_ring_count << ",\n";
    j << "    \"ring\": " << (aromatic_ring_count + saturated_ring_count) << ",\n";
    j << "    \"ring_membership\": " << ring_membership_count << "\n";
    j << "  },\n";
    j << "  \"axis_alignment\": {\n";
    j << "    \"atom\": \"All atom-axis NPYs share row order: "
      <<                "atoms_category_info.row[i] == pos.row[i] == "
      <<                "element.row[i] == residue_index.row[i] == atom_index i. "
      <<                "Calculator atom-axis NPYs (bs_shielding, hm_shielding, "
      <<                "mc_shielding, coulomb_shielding, hbond_shielding, "
      <<                "larsen_hbond_*, tripeptide_*, etc.) follow the same convention.\",\n";
    j << "    \"residue\": \"residues.npy is the canonical residue axis. "
      <<                  "residue_type.npy / residue_index.npy in the identity block are atom-axis "
      <<                  "(N atom rows, each carrying that atom's residue's type / index). "
      <<                  "atoms_category_info.residue_index references the residue axis.\",\n";
    j << "    \"bond\": \"bonds.npy is the canonical bond axis. "
      <<               "bonds.bond_index is the row index.\",\n";
    j << "    \"ring\": \"rings.npy lists aromatic rings first (rows 0..aromatic_ring-1) "
      <<               "then saturated rings (rows aromatic_ring..ring-1). ring_id is the absolute row index.\",\n";
    j << "    \"aromatic_ring\": \"ring_geometry.npy is aromatic-only. "
      <<                          "Its rows correspond to rings.npy entries where ring_kind == 0 (aromatic), in order. "
      <<                          "rings.native_axis_index gives that aromatic-axis index.\",\n";
    j << "    \"ring_contribution_pair\": \"ring_contributions.npy is per (atom, aromatic_ring) "
      <<                                  "pair. The ring_index column references the aromatic-ring axis.\",\n";
    j << "    \"ring_membership\": \"ring_membership.npy is per (ring, ring-vertex-atom) pair. "
      <<                            "ring_id references rings.npy; atom_index references the atom axis.\"\n";
    j << "  }\n";
    j << "}\n";

    const fs::path path = out_dir / "extraction_manifest.json";
    std::ofstream out(path);
    if (!out) {
        OperationLog::Error("TopologySidecar::WriteManifest",
            "could not open " + path.string() + " for write");
        return false;
    }
    out << j.str();
    return out.good();
}

}  // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

int TopologySidecar::WriteFeatures(const Protein& protein,
                                     const std::string& output_dir,
                                     const std::string& protein_id) {
    OperationLog::Scope scope("TopologySidecar::WriteFeatures",
        "atoms=" + std::to_string(protein.AtomCount()) +
        " bonds=" + std::to_string(protein.LegacyAmber().BondCount()) +
        " arom_rings=" + std::to_string(protein.LegacyAmber().AromaticRingCount()) +
        " satd_rings=" + std::to_string(protein.LegacyAmber().SaturatedRingCount()) +
        " dir=" + output_dir);

    fs::create_directories(output_dir);
    const fs::path out_dir(output_dir);

    // Effective protein_id: caller's value or fall back to output_dir's
    // basename. Either way it is just an identifier for joins;
    // misalignment is the consumer's problem.
    std::string effective_id = protein_id;
    if (effective_id.empty()) {
        effective_id = out_dir.filename().string();
    }

    int written = 0;
    if (WriteResidues(protein, out_dir)) ++written;
    if (WriteBonds(protein, out_dir)) ++written;
    size_t membership_rows = 0;
    if (WriteRings(protein, out_dir, &membership_rows)) ++written;
    if (WriteRingMembership(protein, out_dir, membership_rows)) ++written;
    if (WriteManifest(protein, out_dir, effective_id,
                       protein.LegacyAmber().BondCount(),
                       protein.LegacyAmber().AromaticRingCount(),
                       protein.LegacyAmber().SaturatedRingCount(),
                       membership_rows)) ++written;

    OperationLog::Info("TopologySidecar::WriteFeatures",
        "wrote " + std::to_string(written) + " files in " + output_dir);
    return written;
}

}  // namespace nmr
