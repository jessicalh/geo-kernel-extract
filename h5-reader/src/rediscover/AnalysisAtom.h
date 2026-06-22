// AnalysisAtom -- atom branch of the analysis-object model.

#pragma once

#include "AnalysisElement.h"

#include <QJsonObject>
#include <QString>
#include <QTextStream>

#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

struct ClassicalSourceTermRecord {
    QString schema_version = QStringLiteral("classical_source_terms_v2");
    QString granularity = QStringLiteral("atom_leaf");
    QString dataset_id;
    QString protein_id;
    QString atom_uid;
    std::size_t atom_index = 0;
    int residue_index = -1;
    int residue_number = 0;
    QString residue_type;
    QString atom_name;
    QString element;
    QString backbone_role;
    QString frame_kind;
    QString frame_kind_ord;

    std::vector<double> sigma_qm;
    std::vector<double> sigma0;
    std::vector<double> buckingham_linear;
    std::vector<double> buckingham_quadratic;
    std::vector<double> ring;
    std::vector<double> mcconnell;
    std::vector<double> larsen;
    std::vector<double> sigma_cl;
    std::vector<double> residual;

    QString sigma0_status;
    QString buckingham_linear_status;
    QString buckingham_linear_key;
    double buckingham_linear_constant_value = std::numeric_limits<double>::quiet_NaN();
    QString buckingham_linear_units;
    QString buckingham_quadratic_status;
    QString buckingham_quadratic_key;
    double buckingham_quadratic_constant_value = std::numeric_limits<double>::quiet_NaN();
    QString buckingham_quadratic_units;
    QString ring_status = QStringLiteral("good_enough");
    QString mcconnell_status = QStringLiteral("good_enough");
    QString larsen_status = QStringLiteral("cited");
    int constant_placeholder_n = 0;
};

class AnalysisAtom final : public AnalysisElement {
public:
    AnalysisAtom(const AnalysisObjectContext& context,
                 std::size_t atomIndex,
                 PerAtomSubstrateConfig config);

    ~AnalysisAtom() override;

    std::size_t atomIndex() const { return model_index_; }

    void Calculate(std::size_t step) override;
    QJsonObject Truth() const override;

    std::size_t sigmaFolds() const;
    std::size_t relationshipCount() const;
    std::size_t mappedBondCount() const;
    std::size_t mismatchEventCount() const;
    // Compact per-context accumulator counts (replaces the deleted boost grid):
    // total characterized responses, and contexts the atom occupied.
    std::size_t accumulatorResponseCount() const;
    std::size_t accumulatorContextCount() const;
    bool oxygenGatePassed() const;
    std::size_t writeBoundedSigmaRows(QTextStream& out,
                                      const QString& datasetId,
                                      const QString& proteinId) const;
    std::size_t writeClassicalSourceTermRows(QTextStream& out,
                                             const QString& datasetId,
                                             const QString& proteinId) const;
    ClassicalSourceTermRecord classicalSourceTermRecord(const QString& datasetId,
                                                        const QString& proteinId) const;
    std::size_t writeSourceFamilyMatrixRows(QTextStream& out,
                                            const QString& datasetId,
                                            const QString& proteinId) const;
    std::size_t writeSubspaceOverlapRows(QTextStream& out,
                                         const QString& datasetId,
                                         const QString& proteinId) const;
    std::size_t writeEtaByWellRows(QTextStream& out,
                                   const QString& datasetId,
                                   const QString& proteinId) const;

    static void WriteBoundedSigmaHeader(QTextStream& out);
    static void WriteClassicalSourceTermHeader(QTextStream& out);
    static void WriteSourceFamilyMatrixHeader(QTextStream& out);
    static void WriteSubspaceOverlapHeader(QTextStream& out);
    static void WriteEtaByWellHeader(QTextStream& out);
    static bool AssertMolCompOrder(QString* errOut = nullptr);
    static bool AssertPasShapeConvention(QString* errOut = nullptr);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace h5reader::rediscover
