#pragma once

#include "DistributionSummary.h"
#include "RunData.h"

#include <QMap>
#include <QSet>
#include <QString>
#include <QStringList>
#include <QtGlobal>

#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace h5reader::rediscover {

struct Axis2ContextKeyFields {
    QString element;
    QString residue_type;
    QString atom_name;
    QString hyb;
    QString contact_class;
    QString dihedral_region;
    QString SS;
    int schema_version = 1;
};

struct Axis2ContextKey {
    Axis2ContextKeyFields fields;
    QString identity;
    QString context;
    QString canonical;
};

Axis2ContextKey BuildAxis2ContextKey(const Axis2ContextKeyFields& fields);

struct SupportThresholds {
    std::size_t n_min = 3;
    std::size_t n_full = 12;
};

enum class SupportClass {
    Full,
    Reduced,
    Insufficient,
};

struct SupportCredential {
    SupportClass support_class = SupportClass::Insufficient;
    QString support_name = QStringLiteral("insufficient");
    bool may_emit_coupling = false;
    bool may_emit_full_subspace = false;
    QString underpowered_dimensions;
};

SupportCredential CredentialSupport(std::size_t n_proteins,
                                    bool full_rank,
                                    SupportThresholds thresholds = {});

struct PermutationNull {
    int permutation_K = 0;
    double null_slope_mean = 0.0;
    double null_slope_sd = 0.0;
    double obs_slope_z = 0.0;
    double perm_p = 1.0;
};

double LinearSlope(const std::vector<double>& x, const std::vector<double>& y);
double PearsonR(const std::vector<double>& x, const std::vector<double>& y);
PermutationNull ProteinLabelPermutationNull(const std::vector<double>& driver,
                                            const std::vector<double>& sigma,
                                            int k,
                                            quint32 seed = 0xA22026u);

struct HelixDipoleInput {
    std::vector<double> ca_z_A;
    double target_z_A = 0.0;
    double charge_per_residue = 0.35;
};

double ComputeHelixDipoleField(const HelixDipoleInput& input);

struct CohortSample {
    Axis2ContextKey key;
    QString protein_id;
    int atom_index = -1;
    int residue_index = -1;
    QString predecessor_identity;
    QString psi_iminus1_region;
    double psi_iminus1 = std::numeric_limits<double>::quiet_NaN();
    double psi_own = std::numeric_limits<double>::quiet_NaN();
    double sigma_iso = std::numeric_limits<double>::quiet_NaN();
    double sigma_eta_H = std::numeric_limits<double>::quiet_NaN();
    double mol_xx = std::numeric_limits<double>::quiet_NaN();
    double mol_yy = std::numeric_limits<double>::quiet_NaN();
    double mol_xy = std::numeric_limits<double>::quiet_NaN();
    double mol_xz = std::numeric_limits<double>::quiet_NaN();
    double mol_yz = std::numeric_limits<double>::quiet_NaN();
    double mol_zz = std::numeric_limits<double>::quiet_NaN();
    bool backbone_n = false;
    bool helix_member = false;
    double helix_dipole_field = std::numeric_limits<double>::quiet_NaN();
    double mutation_contact_kernel = std::numeric_limits<double>::quiet_NaN();
    QString contact_scope;
    QString any_site_scope;
    bool near_mutation = false;
    bool distant_from_all_sites = true;
    QMap<QString, double> channels;
};

struct CohortCellTruth {
    Axis2ContextKey key;
    QSet<QString> proteins;
    std::vector<CohortSample> samples;
};

class CohortContextAccumulator {
public:
    void push(const CohortSample& sample);
    const std::map<QString, CohortCellTruth>& cells() const { return cells_; }
    std::size_t sampleCount() const { return sample_count_; }
    std::size_t cellCount() const { return cells_.size(); }

private:
    std::map<QString, CohortCellTruth> cells_;
    std::size_t sample_count_ = 0;
};

struct Axis2RunOptions {
    QString root720;
    QString mutantRoot;
    QString axis1OverlayDir;
    QString outDir;
    std::size_t maxProteins = 0;
    int permutationK = 99;
    SupportThresholds supportThresholds;
    bool runStatic = true;
    bool runMutantDeltaRidge = true;
};

struct Axis2RunStats {
    std::size_t proteins_seen = 0;
    std::size_t proteins_loaded = 0;
    std::size_t static_samples = 0;
    std::size_t static_cells = 0;
    std::size_t ridge_samples = 0;
    std::size_t ridge_rows = 0;
    std::size_t distinct_identities = 0;
    std::size_t distinct_elements = 0;
    std::size_t max_atoms_in_resident_protein = 0;
    std::size_t delta_refusals = 0;
    QStringList refusal_reasons;
};

bool RunCohortAxis2(const Axis2RunOptions& options,
                    Axis2RunStats* stats,
                    QString* err_out = nullptr);

bool AuditAxis2Outputs(const QString& outDir, QString* err_out = nullptr);

}  // namespace h5reader::rediscover
