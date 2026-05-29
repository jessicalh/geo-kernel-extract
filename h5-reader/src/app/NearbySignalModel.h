// NearbySignalModel -- focus-neighbourhood candidates for dashboard strips.
//
// This is discovery state, not atom selection. It answers "what nearby
// topology targets could be bound to a strip?" without mutating the
// interactive AtomSelection model.

#pragma once

#include "../model/SignalDictionary.h"

#include <QAbstractTableModel>
#include <QString>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::model {
class Conformation;
class QtProtein;
}

namespace h5reader::app {

class NearbySignalModel final : public QAbstractTableModel {
    Q_OBJECT

public:
    enum class CandidateKind {
        Atom,
        Residue,
        Bond,
        // Per-bond-vector candidate, populated by walking the iRED +
        // ReorientationalDynamics identity tables on the active
        // trajectory. Carries a BondVectorAnchor{residue, kind} —
        // (residue, kind=1) for N-H, kind=2 for Cα-Hα, kind=3 for C=O.
        BondVector,
        Ring,
        AromaticRing,
        SaturatedRing,
        RingMembership,
    };

    struct Candidate {
        CandidateKind kind = CandidateKind::Atom;
        model::SignalAnchor anchor = model::NoneAnchor{};
        std::optional<std::size_t> residueContext;
        std::optional<std::size_t> absoluteRing;
        double distanceAngstrom = 0.0;
        QString label;
    };

    enum Roles {
        KindRole = Qt::UserRole + 1,
        AtomIndexRole,
        ResidueIndexRole,
        BondIndexRole,
        RingIndexRole,
        RingMembershipIndexRole,
        DistanceRole,
        AnchorKindRole,
    };

    explicit NearbySignalModel(QObject* parent = nullptr);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    int columnCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
    QHash<int, QByteArray> roleNames() const override;

    void setContext(const model::QtProtein* protein, model::Conformation* conformation);
    void setRadiusAngstrom(double radius);
    void setAnchor(std::size_t atom, int frame);
    void clear();

    double radiusAngstrom() const { return radiusAngstrom_; }
    std::optional<std::size_t> anchorAtom() const { return anchorAtom_; }
    int frame() const { return frame_; }
    const Candidate* candidateAt(const QModelIndex& index) const;

private:
    QString residueLabel(std::size_t residueIdx) const;
    QString atomLabel(std::size_t atomIdx) const;
    QString bondLabel(std::size_t bondIdx) const;
    QString ringLabel(std::size_t ringIdx) const;
    QString ringMembershipLabel(std::size_t membershipIdx) const;
    void rebuild();

    const model::QtProtein* protein_ = nullptr;
    model::Conformation* conformation_ = nullptr;
    double radiusAngstrom_ = 5.0;
    std::optional<std::size_t> anchorAtom_;
    int frame_ = 0;
    std::vector<Candidate> candidates_;
};

}  // namespace h5reader::app
