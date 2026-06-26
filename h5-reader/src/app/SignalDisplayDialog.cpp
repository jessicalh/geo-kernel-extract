#include "SignalDisplayDialog.h"

#include "DashboardSelectionController.h"
#include "NearbySignalModel.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/DashboardLogging.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/DisplayModeCapability.h"
#include "../model/MetricTaxonomy.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/VisualizationRegistry.h"

#include <QAbstractItemView>
#include <QAbstractTableModel>
#include <QBrush>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFont>
#include <QGroupBox>
#include <QHash>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QItemSelectionModel>
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QPointer>
#include <QSet>
#include <QSignalBlocker>
#include <QSortFilterProxyModel>
#include <QSplitter>
#include <QTableView>
#include <QTreeView>
#include <QUuid>
#include <QVBoxLayout>

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <utility>

namespace h5reader::app {

namespace {

struct ModeControl {
    model::VisualizationType type = model::VisualizationType::TemporalStrip;
    QCheckBox* box = nullptr;
};

struct DescriptorRecord {
    model::SignalDescriptor descriptor;
    QString source;
    QString axis;
    QString requiredAnchor;
    QString valueShape;
    QString units;
    QString residency;
    QStringList displayModes;

    QString searchText() const {
        return QStringList{
                   descriptor.id,
                   descriptor.conceptKey,
                   descriptor.importSet,
                   descriptor.family,
                   descriptor.label,
                   descriptor.description,
                   descriptor.storagePath,
                   source,
                   axis,
                   requiredAnchor,
                   valueShape,
                   units,
                   residency,
                   displayModes.join(QLatin1Char(' ')),
                   descriptor.tags.join(QLatin1Char(' ')),
               }.join(QLatin1Char(' '));
    }
};

QString visualizationTypeLabel(model::VisualizationType type) {
    switch (type) {
    case model::VisualizationType::TemporalStrip:
        return QStringLiteral("Strip");
    case model::VisualizationType::SequenceBar:
        return QStringLiteral("Bar (sequence)");
    case model::VisualizationType::LagCurve:
        return QStringLiteral("Curve (lag)");
    case model::VisualizationType::ChordCoupling:
        return QStringLiteral("Chord (coupling)");
    case model::VisualizationType::FixedFrequency:
        return QStringLiteral("Fixed-freq J(w)");
    case model::VisualizationType::PowerSpectrum:
        return QStringLiteral("Power spectrum");
    case model::VisualizationType::TensorGlyph:
        return QStringLiteral("Glyph / overlay");
    case model::VisualizationType::Newman:
        return QStringLiteral("Newman");
    }
    return {};
}

QString visualizationTypeKey(model::VisualizationType type) {
    switch (type) {
    case model::VisualizationType::TemporalStrip:
        return QStringLiteral("strip");
    case model::VisualizationType::SequenceBar:
        return QStringLiteral("barSequence");
    case model::VisualizationType::LagCurve:
        return QStringLiteral("curveLag");
    case model::VisualizationType::ChordCoupling:
        return QStringLiteral("chordCoupling");
    case model::VisualizationType::FixedFrequency:
        return QStringLiteral("fixedFreq");
    case model::VisualizationType::PowerSpectrum:
        return QStringLiteral("powerSpectrum");
    case model::VisualizationType::TensorGlyph:
        return QStringLiteral("tensorGlyph");
    case model::VisualizationType::Newman:
        return QStringLiteral("newman");
    }
    return {};
}

std::optional<model::VisualizationType> visualizationTypeForKey(const QString& key) {
    if (key == QStringLiteral("strip"))
        return model::VisualizationType::TemporalStrip;
    if (key == QStringLiteral("barSequence"))
        return model::VisualizationType::SequenceBar;
    if (key == QStringLiteral("curveLag"))
        return model::VisualizationType::LagCurve;
    if (key == QStringLiteral("chordCoupling"))
        return model::VisualizationType::ChordCoupling;
    if (key == QStringLiteral("fixedFreq"))
        return model::VisualizationType::FixedFrequency;
    if (key == QStringLiteral("powerSpectrum"))
        return model::VisualizationType::PowerSpectrum;
    if (key == QStringLiteral("tensorGlyph"))
        return model::VisualizationType::TensorGlyph;
    if (key == QStringLiteral("newman"))
        return model::VisualizationType::Newman;
    return std::nullopt;
}

bool modeMatchesType(const QString& modeId, model::VisualizationType type) {
    const QString lower = modeId.toLower();
    switch (type) {
    case model::VisualizationType::TemporalStrip:
        return lower.startsWith(QStringLiteral("strip."))
            && !lower.contains(QStringLiteral("spectrum"))
            && !lower.contains(QStringLiteral("fft"));
    case model::VisualizationType::SequenceBar:
    case model::VisualizationType::LagCurve:
    case model::VisualizationType::ChordCoupling:
    case model::VisualizationType::FixedFrequency:
    case model::VisualizationType::PowerSpectrum:
    case model::VisualizationType::TensorGlyph:
    case model::VisualizationType::Newman:
        if (const model::VisualizationDefinition* definition =
                model::VisualizationRegistry::instance().definitionForMode(modeId)) {
            return definition->type() == type;
        }
        return false;
    }
    return false;
}

QString modeForType(const QStringList& modes, model::VisualizationType type) {
    if (type == model::VisualizationType::TemporalStrip && modes.contains(QStringLiteral("strip.scalar")))
        return QStringLiteral("strip.scalar");
    for (const QString& mode : modes) {
        if (modeMatchesType(mode, type))
            return mode;
    }
    return {};
}

bool modeListContainsType(const QStringList& modes, model::VisualizationType type) {
    return !modeForType(modes, type).isEmpty();
}

const model::VisualizationDefinition*
definitionForType(const QVector<const model::VisualizationDefinition*>& definitions,
                  model::VisualizationType type) {
    for (const model::VisualizationDefinition* definition : definitions) {
        if (definition && definition->type() == type)
            return definition;
    }
    return nullptr;
}

}  // namespace

bool DashboardModeHasVisibleSurface(const QString& modeId) {
    return model::DisplayModeCapabilityFor(modeId).hasVisibleSurface;
}

namespace {

QString modeSummary(const QStringList& displayModes) {
    QStringList labels;
    for (model::VisualizationType type : {model::VisualizationType::TemporalStrip,
                                          model::VisualizationType::SequenceBar,
                                          model::VisualizationType::LagCurve,
                                          model::VisualizationType::ChordCoupling,
                                          model::VisualizationType::FixedFrequency}) {
        if (modeListContainsType(displayModes, type))
            labels.push_back(visualizationTypeLabel(type));
    }
    return labels.isEmpty() ? QStringLiteral("None") : labels.join(QStringLiteral(", "));
}

QString unitsLabel(const model::SignalDescriptor& descriptor) {
    if (!descriptor.defaultDisplayUnits.displaySymbol.isEmpty())
        return descriptor.defaultDisplayUnits.displaySymbol;
    if (!descriptor.sourceUnits.displaySymbol.isEmpty())
        return descriptor.sourceUnits.displaySymbol;
    return descriptor.sourceUnits.sourceSymbol;
}

DescriptorRecord recordFromDescriptor(const model::SignalDescriptor& descriptor) {
    DescriptorRecord record;
    record.descriptor = descriptor;
    record.source = model::ToString(descriptor.sourceKind);
    record.axis = model::ToString(descriptor.nativeAxis);
    record.requiredAnchor = model::ToString(descriptor.requiredAnchor);
    record.valueShape = model::ToString(descriptor.valueShape);
    record.units = unitsLabel(descriptor);
    record.residency = model::ToString(descriptor.residency);
    record.displayModes = model::AllDisplayModes(descriptor);
    return record;
}

model::SignalAxis axisForCandidate(const NearbySignalModel::Candidate& candidate) {
    switch (candidate.kind) {
    case NearbySignalModel::CandidateKind::Atom:
        return model::SignalAxis::Atom;
    case NearbySignalModel::CandidateKind::Residue:
        return model::SignalAxis::Residue;
    case NearbySignalModel::CandidateKind::Bond:
        return model::SignalAxis::Bond;
    case NearbySignalModel::CandidateKind::BondVector:
        return model::SignalAxis::BondVector;
    case NearbySignalModel::CandidateKind::Ring:
        return model::SignalAxis::Ring;
    case NearbySignalModel::CandidateKind::AromaticRing:
        return model::SignalAxis::AromaticRing;
    case NearbySignalModel::CandidateKind::SaturatedRing:
        return model::SignalAxis::SaturatedRing;
    case NearbySignalModel::CandidateKind::RingMembership:
        return model::SignalAxis::RingMembership;
    }
    return model::SignalAxis::None;
}

model::SignalAnchor anchorForCandidate(const NearbySignalModel::Candidate& candidate) {
    return candidate.anchor;
}

QJsonObject signalAnchorToJson(const model::SignalAnchor& anchor) {
    QJsonObject out;
    if (std::holds_alternative<model::NoneAnchor>(anchor)) {
        out["kind"] = "none";
    } else if (const auto* a = std::get_if<model::AtomAnchor>(&anchor)) {
        out["kind"] = "atom"; out["atom"] = static_cast<qint64>(a->atom);
    } else if (const auto* r = std::get_if<model::ResidueAnchor>(&anchor)) {
        out["kind"] = "residue"; out["residue"] = static_cast<qint64>(r->residue);
    } else if (const auto* t = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        QJsonArray atoms;
        for (auto a : t->atoms) atoms.append(static_cast<qint64>(a));
        out["kind"] = "atom_tuple"; out["atoms"] = atoms;
    } else if (const auto* b = std::get_if<model::BondAnchor>(&anchor)) {
        out["kind"] = "bond"; out["bond"] = static_cast<qint64>(b->bond);
    } else if (const auto* v = std::get_if<model::BondVectorAnchor>(&anchor)) {
        out["kind"] = "bond_vector";
        out["residue"] = static_cast<qint64>(v->residue);
        out["kind_id"] = static_cast<qint64>(v->kind);
    } else if (const auto* r = std::get_if<model::RingAnchor>(&anchor)) {
        out["kind"] = "ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* r = std::get_if<model::AromaticRingAnchor>(&anchor)) {
        out["kind"] = "aromatic_ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* r = std::get_if<model::SaturatedRingAnchor>(&anchor)) {
        out["kind"] = "saturated_ring"; out["ring"] = static_cast<qint64>(r->ring);
    } else if (const auto* p = std::get_if<model::RingContributionPairAnchor>(&anchor)) {
        out["kind"] = "ring_contribution_pair"; out["pair"] = static_cast<qint64>(p->pair);
    } else if (const auto* m = std::get_if<model::RingMembershipAnchor>(&anchor)) {
        out["kind"] = "ring_membership"; out["membership"] = static_cast<qint64>(m->membership);
    } else if (const auto* p = std::get_if<model::MutationMatchPairAnchor>(&anchor)) {
        out["kind"] = "mutation_match_pair"; out["pair"] = static_cast<qint64>(p->pair);
    } else if (std::holds_alternative<model::ProteinAnchor>(anchor)) {
        out["kind"] = "protein";
    } else if (std::holds_alternative<model::SystemAnchor>(anchor)) {
        out["kind"] = "system";
    } else if (std::holds_alternative<model::EventAnchor>(anchor)) {
        out["kind"] = "event";
    }
    return out;
}

QString candidateModeDisabledReason(const QString& modeId,
                                    bool structurallySupported,
                                    bool descriptorAvailable) {
    if (modeId.isEmpty())
        return QStringLiteral("This descriptor does not offer that display mode.");
    if (!structurallySupported)
        return QStringLiteral("This display mode has no implemented renderer yet.");
    if (!descriptorAvailable)
        return QStringLiteral("The data for this descriptor is not available in this run.");
    return {};
}

// Delegates to model::AxisCanSatisfy so the dialog's widening rules
// stay in lockstep with the controller-side AnchorMatchesAxis. Earlier
// versions of this dialog kept a local copy that missed the
// Residue → BondVector widening, so iRED / Reorient descriptors were
// silently filtered out of the candidate table when a residue was the
// active anchor (Codex NOW-3, 2026-05-29).
bool anchorAxisCanSatisfy(model::SignalAxis selectedAxis, model::SignalAxis requiredAxis) {
    return model::AxisCanSatisfy(selectedAxis, requiredAxis);
}

std::array<model::VisualizationType, 6> allVisualizationTypes() {
    return {model::VisualizationType::TemporalStrip,
            model::VisualizationType::SequenceBar,
            model::VisualizationType::LagCurve,
            model::VisualizationType::ChordCoupling,
            model::VisualizationType::FixedFrequency,
            model::VisualizationType::Newman};
}

void configureTable(QTableView* view) {
    view->setSelectionBehavior(QAbstractItemView::SelectRows);
    view->setSelectionMode(QAbstractItemView::SingleSelection);
    view->setAlternatingRowColors(true);
    view->setSortingEnabled(true);
    view->setShowGrid(false);
    view->verticalHeader()->hide();
    view->verticalHeader()->setDefaultSectionSize(22);
    view->horizontalHeader()->setStretchLastSection(false);
}

void configureTree(QTreeView* view) {
    view->setSelectionBehavior(QAbstractItemView::SelectRows);
    view->setSelectionMode(QAbstractItemView::SingleSelection);
    view->setUniformRowHeights(true);
    view->setAllColumnsShowFocus(true);
    view->setSortingEnabled(false);
    view->setAnimated(false);
    view->setExpandsOnDoubleClick(true);
    view->header()->setStretchLastSection(false);
    view->header()->setDefaultSectionSize(110);
}

void refillCombo(QComboBox* combo, const QString& allLabel, const QStringList& values) {
    const QVariant previous = combo->currentData();
    QSignalBlocker blocker(combo);
    combo->clear();
    combo->addItem(allLabel, QString());
    for (const QString& value : values)
        combo->addItem(value, value);
    const int previousIndex = combo->findData(previous);
    combo->setCurrentIndex(previousIndex >= 0 ? previousIndex : 0);
}

// ---- Mechanism -> concept -> form tree over the catalog ----------------------
// The candidate list is a MetricTaxonomy tree, not a flat table: the ~188
// descriptors fold to ~35 base concepts in ~11 mechanism groups (the four
// shielding-kernel hypotheses first, then the DFT reference, conditioning
// inputs, dynamics, scaffold). Leaves carry an index into a flat record vector
// so selection / filtering / mode-offering reuse the SAME DescriptorRecord the
// old table model used. Electrostatic is sub-grouped by charge model.

QString groupTitle(model::MetricGroup group) {
    switch (group) {
    case model::MetricGroup::RingCurrent:    return QStringLiteral("Ring current");
    case model::MetricGroup::BondAnisotropy: return QStringLiteral("Bond anisotropy (McConnell)");
    case model::MetricGroup::Electrostatic:  return QStringLiteral("Electrostatic (E-field / EFG)");
    case model::MetricGroup::HBond:          return QStringLiteral("H-bond");
    case model::MetricGroup::DftReference:   return QStringLiteral("DFT reference");
    case model::MetricGroup::Charges:        return QStringLiteral("Charges & electronic structure");
    case model::MetricGroup::Solvation:      return QStringLiteral("Solvation & water");
    case model::MetricGroup::Structure:      return QStringLiteral("Structure & geometry");
    case model::MetricGroup::Dynamics:       return QStringLiteral("Dynamics");
    case model::MetricGroup::Scaffold:       return QStringLiteral("Scaffold & bookkeeping");
    case model::MetricGroup::Mutation:       return QStringLiteral("Mutation (WT / mutant / delta)");
    case model::MetricGroup::Other:          return QStringLiteral("Other");
    }
    return QStringLiteral("Other");
}

QString roleTitle(model::MetricRole role) {
    switch (role) {
    case model::MetricRole::Hypothesis: return QStringLiteral("contribution");
    case model::MetricRole::Reference:  return QStringLiteral("reference");
    case model::MetricRole::Input:      return QStringLiteral("input");
    case model::MetricRole::Dynamics:   return QStringLiteral("dynamics");
    case model::MetricRole::Scaffold:   return QStringLiteral("scaffold");
    }
    return {};
}

QString formTitle(model::MetricForm form) {
    switch (form) {
    case model::MetricForm::Snapshot:   return QStringLiteral("snapshot");
    case model::MetricForm::Series:     return QStringLiteral("series");
    case model::MetricForm::Rollup:     return QStringLiteral("rollup");
    case model::MetricForm::Dynamics:   return QStringLiteral("autocorr");
    case model::MetricForm::Transition: return QStringLiteral("transition");
    case model::MetricForm::Reference:  return QStringLiteral("DFT");
    case model::MetricForm::Spine:      return QStringLiteral("topology");
    case model::MetricForm::Derived:    return QStringLiteral("derived");
    case model::MetricForm::Other:      return QStringLiteral("-");
    }
    return QStringLiteral("-");
}

struct DescriptorTreeNode {
    enum class Kind : std::uint8_t { Group, ChargeModel, Concept, Form };
    Kind kind = Kind::Group;
    QString name;         // column 0 label (group / charge-model / concept / form)
    QString detail;       // group: role; concept or leaf: forms summary
    QString searchExtra;  // ancestor labels, appended to a leaf's search text
    int recordIndex = -1; // >= 0 marks a selectable leaf (index into records_)
    DescriptorTreeNode* parent = nullptr;
    int row = 0;
    QVector<DescriptorTreeNode*> children;

    ~DescriptorTreeNode() { qDeleteAll(children); }
};

DescriptorTreeNode* appendChild(DescriptorTreeNode* parent, DescriptorTreeNode::Kind kind,
                                const QString& name) {
    auto* node = new DescriptorTreeNode;
    node->kind = kind;
    node->name = name;
    node->parent = parent;
    node->row = static_cast<int>(parent->children.size());
    parent->children.push_back(node);
    return node;
}

class DescriptorTreeModel final : public QAbstractItemModel {
public:
    enum Column : std::uint8_t {
        NameColumn,
        FormColumn,
        ShapeColumn,
        DisplaysColumn,
        UnitsColumn,
        ColumnCount,
    };

    enum Role : std::uint16_t {
        DisplayModesRole = Qt::UserRole + 1,
        SearchTextRole,
        SourceRole,
        AxisRole,
        ShapeRole,
        RequiredAnchorAxisRole,
        RecordIndexRole,
        KindRole,
    };

    explicit DescriptorTreeModel(QObject* parent = nullptr)
        : QAbstractItemModel(parent)
        , root_(new DescriptorTreeNode)
    {
        CENSUS_REGISTER(this);
        setObjectName(QStringLiteral("SignalDescriptorTreeModel"));
    }

    ~DescriptorTreeModel() override { delete root_; }

    QModelIndex index(int row, int column, const QModelIndex& parent = QModelIndex()) const override {
        if (!hasIndex(row, column, parent))
            return {};
        DescriptorTreeNode* parentNode = nodeFor(parent);
        if (!parentNode || row < 0 || row >= parentNode->children.size())
            return {};
        return createIndex(row, column, parentNode->children.at(row));
    }

    QModelIndex parent(const QModelIndex& index) const override {
        if (!index.isValid())
            return {};
        auto* node = static_cast<DescriptorTreeNode*>(index.internalPointer());
        DescriptorTreeNode* parentNode = node ? node->parent : nullptr;
        if (!parentNode || parentNode == root_)
            return {};
        return createIndex(parentNode->row, 0, parentNode);
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const override {
        if (parent.column() > 0)
            return 0;
        DescriptorTreeNode* node = nodeFor(parent);
        return node ? static_cast<int>(node->children.size()) : 0;
    }

    int columnCount(const QModelIndex& = QModelIndex()) const override {
        return ColumnCount;
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override {
        DescriptorTreeNode* node = nodeFor(index);
        if (!node || node == root_)
            return {};
        const DescriptorRecord* record =
            (node->recordIndex >= 0 && node->recordIndex < records_.size())
                ? &records_.at(node->recordIndex)
                : nullptr;

        if (role == Qt::DisplayRole) {
            switch (index.column()) {
            case NameColumn:
                return node->name;
            case FormColumn:
                return node->detail;
            case ShapeColumn:
                return record ? record->valueShape : QString();
            case DisplaysColumn:
                return record ? modeSummary(record->displayModes) : QString();
            case UnitsColumn:
                return record ? record->units : QString();
            default:
                return {};
            }
        }
        if (role == Qt::ToolTipRole) {
            if (record) {
                return QStringLiteral("%1\nid: %2\nconcept: %3\nmodes: %4")
                    .arg(record->descriptor.label,
                         record->descriptor.id,
                         record->descriptor.conceptKey.isEmpty() ? QStringLiteral("-")
                                                                 : record->descriptor.conceptKey,
                         record->displayModes.join(QStringLiteral(", ")));
            }
            return node->detail.isEmpty()
                       ? node->name
                       : QStringLiteral("%1  (%2)").arg(node->name, node->detail);
        }
        if (role == Qt::FontRole && node->kind == DescriptorTreeNode::Kind::Group) {
            QFont font;
            font.setBold(true);
            return font;
        }
        if (role == Qt::ForegroundRole && record && record->displayModes.isEmpty())
            return QBrush(Qt::darkGray);
        if (role == RecordIndexRole)
            return node->recordIndex;
        if (role == KindRole)
            return static_cast<int>(node->kind);

        if (record) {
            switch (role) {
            case DisplayModesRole:
                return record->displayModes;
            case SearchTextRole:
                return node->searchExtra.isEmpty()
                           ? record->searchText()
                           : record->searchText() + QLatin1Char(' ') + node->searchExtra;
            case SourceRole:
                return record->source;
            case AxisRole:
                return record->axis;
            case ShapeRole:
                return record->valueShape;
            case RequiredAnchorAxisRole:
                return static_cast<int>(record->descriptor.requiredAnchor);
            default:
                return {};
            }
        }
        if (role == SearchTextRole)
            return node->name + QLatin1Char(' ') + node->detail;
        if (role == RequiredAnchorAxisRole)
            return -1;
        return {};
    }

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const override {
        if (orientation != Qt::Horizontal || role != Qt::DisplayRole)
            return {};
        switch (section) {
        case NameColumn:
            return QStringLiteral("Metric");
        case FormColumn:
            return QStringLiteral("Form");
        case ShapeColumn:
            return QStringLiteral("Shape");
        case DisplaysColumn:
            return QStringLiteral("Displays");
        case UnitsColumn:
            return QStringLiteral("Units");
        default:
            return {};
        }
    }

    void setDescriptors(const QVector<model::SignalDescriptor>& descriptors) {
        beginResetModel();
        delete root_;
        root_ = new DescriptorTreeNode;
        records_.clear();
        records_.reserve(descriptors.size());
        QHash<QString, int> indexById;
        indexById.reserve(descriptors.size());
        for (const model::SignalDescriptor& descriptor : descriptors) {
            indexById.insert(descriptor.id, static_cast<int>(records_.size()));
            records_.push_back(recordFromDescriptor(descriptor));
        }

        const auto conceptLabel = [](const model::MetricConceptNode& concept) {
            return concept.label.isEmpty() ? concept.baseConcept : concept.label;
        };
        const auto formsSummary = [](const model::MetricConceptNode& concept) {
            QStringList names;
            for (const model::MetricFormEntry& form : concept.forms)
                names.push_back(formTitle(form.form));
            names.removeDuplicates();
            return names.join(QStringLiteral(", "));
        };

        QVector<model::MetricGroupNode> groups = model::GroupCatalog(descriptors);
        // Present mechanism groups in the curated hypothesis-first order the
        // MetricGroup enum is authored in (the ring-current / bond / electrostatic
        // / H-bond kernels, then the DFT reference, then inputs / dynamics /
        // scaffold) rather than catalog-insertion order.
        std::sort(groups.begin(), groups.end(),
                  [](const model::MetricGroupNode& a, const model::MetricGroupNode& b) {
                      return static_cast<int>(a.group) < static_cast<int>(b.group);
                  });
        for (const model::MetricGroupNode& group : groups) {
            DescriptorTreeNode* groupNode =
                appendChild(root_, DescriptorTreeNode::Kind::Group, groupTitle(group.group));
            groupNode->detail = roleTitle(group.role);

            const auto addConcept = [&](DescriptorTreeNode* parent,
                                        const model::MetricConceptNode& concept,
                                        const QString& parentLabel) {
                const QString ancestor =
                    (groupNode->name + QLatin1Char(' ') + parentLabel + QLatin1Char(' ')
                     + concept.baseConcept).simplified();
                if (concept.forms.size() == 1) {
                    DescriptorTreeNode* leaf =
                        appendChild(parent, DescriptorTreeNode::Kind::Concept, conceptLabel(concept));
                    leaf->recordIndex = indexById.value(concept.forms.front().descriptorId, -1);
                    leaf->detail = formTitle(concept.forms.front().form);
                    leaf->searchExtra = ancestor;
                    return;
                }
                DescriptorTreeNode* conceptNode =
                    appendChild(parent, DescriptorTreeNode::Kind::Concept, conceptLabel(concept));
                conceptNode->detail = formsSummary(concept);
                for (const model::MetricFormEntry& form : concept.forms) {
                    DescriptorTreeNode* leaf =
                        appendChild(conceptNode, DescriptorTreeNode::Kind::Form, formTitle(form.form));
                    leaf->recordIndex = indexById.value(form.descriptorId, -1);
                    leaf->detail = formTitle(form.form);
                    leaf->searchExtra = ancestor + QLatin1Char(' ') + conceptNode->name;
                }
            };

            if (group.group == model::MetricGroup::Electrostatic) {
                QHash<QString, DescriptorTreeNode*> byChargeModel;
                for (const model::MetricConceptNode& concept : group.concepts) {
                    const QString chargeModel =
                        concept.chargeModel.isEmpty() ? QStringLiteral("Other") : concept.chargeModel;
                    DescriptorTreeNode* bucket = byChargeModel.value(chargeModel, nullptr);
                    if (!bucket) {
                        bucket = appendChild(groupNode, DescriptorTreeNode::Kind::ChargeModel, chargeModel);
                        byChargeModel.insert(chargeModel, bucket);
                    }
                    addConcept(bucket, concept, chargeModel);
                }
            } else {
                for (const model::MetricConceptNode& concept : group.concepts)
                    addConcept(groupNode, concept, QString());
            }
        }
        endResetModel();
    }

    const DescriptorRecord* recordForIndex(const QModelIndex& index) const {
        DescriptorTreeNode* node = nodeFor(index);
        if (!node || node->recordIndex < 0 || node->recordIndex >= records_.size())
            return nullptr;
        return &records_.at(node->recordIndex);
    }

    QStringList uniqueValues(int role) const {
        QSet<QString> seen;
        for (const DescriptorRecord& record : records_) {
            if (role == SourceRole && !record.source.isEmpty())
                seen.insert(record.source);
            else if (role == AxisRole && !record.axis.isEmpty())
                seen.insert(record.axis);
            else if (role == ShapeRole && !record.valueShape.isEmpty())
                seen.insert(record.valueShape);
        }
        QStringList values(seen.cbegin(), seen.cend());
        values.sort(Qt::CaseInsensitive);
        return values;
    }

private:
    DescriptorTreeNode* nodeFor(const QModelIndex& index) const {
        if (!index.isValid())
            return root_;
        return static_cast<DescriptorTreeNode*>(index.internalPointer());
    }

    DescriptorTreeNode* root_ = nullptr;
    QVector<DescriptorRecord> records_;
};

// Depth-first first selectable leaf in a (possibly filtered) candidate model.
// With requireDisplayable, skip leaves that declare no display modes (the
// honestly-non-displayable rows) so auto-selection lands on an ACTIONABLE
// metric rather than a greyed one; callers fall back without the constraint.
QModelIndex firstCandidateLeaf(QAbstractItemModel* model, const QModelIndex& parent,
                               bool requireDisplayable) {
    if (!model)
        return {};
    const int rows = model->rowCount(parent);
    for (int row = 0; row < rows; ++row) {
        const QModelIndex index = model->index(row, 0, parent);
        if (model->data(index, DescriptorTreeModel::RecordIndexRole).toInt() >= 0) {
            if (!requireDisplayable
                || !model->data(index, DescriptorTreeModel::DisplayModesRole).toStringList().isEmpty())
                return index;
            continue;
        }
        const QModelIndex deep = firstCandidateLeaf(model, index, requireDisplayable);
        if (deep.isValid())
            return deep;
    }
    return {};
}

// Default expansion: mechanism groups and charge-model buckets are always open;
// multi-form concept nodes stay collapsed unless a search/filter is active (so
// the matches under them are revealed). Operates on the proxy the view is bound to.
void expandCandidateTree(QTreeView* view, QAbstractItemModel* model,
                         const QModelIndex& parent, bool expandConcepts) {
    if (!view || !model)
        return;
    const int rows = model->rowCount(parent);
    for (int row = 0; row < rows; ++row) {
        const QModelIndex index = model->index(row, 0, parent);
        if (model->data(index, DescriptorTreeModel::RecordIndexRole).toInt() >= 0)
            continue;  // leaf
        const int kind = model->data(index, DescriptorTreeModel::KindRole).toInt();
        const bool structural = kind == static_cast<int>(DescriptorTreeNode::Kind::Group)
                                || kind == static_cast<int>(DescriptorTreeNode::Kind::ChargeModel);
        if (structural || expandConcepts)
            view->expand(index);
        else
            view->collapse(index);
        expandCandidateTree(view, model, index, expandConcepts);
    }
}

class DescriptorFilterProxyModel final : public QSortFilterProxyModel {
public:
    explicit DescriptorFilterProxyModel(QObject* parent = nullptr)
        : QSortFilterProxyModel(parent)
    {
        CENSUS_REGISTER(this);
        setObjectName(QStringLiteral("SignalDescriptorFilterProxyModel"));
        setFilterCaseSensitivity(Qt::CaseInsensitive);
        setSortCaseSensitivity(Qt::CaseInsensitive);
        setRecursiveFilteringEnabled(true);
    }

    void setSearchText(const QString& text) {
        if (searchText_ == text.trimmed())
            return;
        searchText_ = text.trimmed();
        invalidateFilter();
    }

    void setSourceFilter(const QString& value) {
        if (sourceFilter_ == value)
            return;
        sourceFilter_ = value;
        invalidateFilter();
    }

    void setAxisFilter(const QString& value) {
        if (axisFilter_ == value)
            return;
        axisFilter_ = value;
        invalidateFilter();
    }

    void setShapeFilter(const QString& value) {
        if (shapeFilter_ == value)
            return;
        shapeFilter_ = value;
        invalidateFilter();
    }

    void setModeKindFilter(const QString& value) {
        if (modeKindFilter_ == value)
            return;
        modeKindFilter_ = value;
        invalidateFilter();
    }

    void setRequiredAnchorFilter(model::SignalAxis axis, bool enabled) {
        const int value = enabled ? static_cast<int>(axis) : -1;
        if (requiredAnchorFilter_ == value)
            return;
        requiredAnchorFilter_ = value;
        invalidateFilter();
    }

protected:
    bool filterAcceptsRow(int sourceRow, const QModelIndex& sourceParent) const override {
        const QAbstractItemModel* model = sourceModel();
        if (!model)
            return false;
        const QModelIndex index = model->index(sourceRow, 0, sourceParent);
        // Header rows (groups / charge models / multi-form concepts) carry no
        // record: reject them here and let recursive filtering pull them back in
        // whenever a descendant leaf is accepted. Only leaves run the predicates.
        if (model->data(index, DescriptorTreeModel::RecordIndexRole).toInt() < 0)
            return false;
        const auto textForRole = [&](int role) {
            return model->data(index, role).toString();
        };
        if (!sourceFilter_.isEmpty() && textForRole(DescriptorTreeModel::SourceRole) != sourceFilter_)
            return false;
        if (!axisFilter_.isEmpty() && textForRole(DescriptorTreeModel::AxisRole) != axisFilter_)
            return false;
        if (!shapeFilter_.isEmpty() && textForRole(DescriptorTreeModel::ShapeRole) != shapeFilter_)
            return false;
        if (requiredAnchorFilter_ >= 0) {
            const auto selectedAxis = static_cast<model::SignalAxis>(requiredAnchorFilter_);
            const auto requiredAxis =
                static_cast<model::SignalAxis>(
                    model->data(index, DescriptorTreeModel::RequiredAnchorAxisRole).toInt());
            if (!anchorAxisCanSatisfy(selectedAxis, requiredAxis))
                return false;
        }
        if (!modeKindFilter_.isEmpty()) {
            const std::optional<model::VisualizationType> type = visualizationTypeForKey(modeKindFilter_);
            if (!type)
                return false;
            const QStringList modes = model->data(index, DescriptorTreeModel::DisplayModesRole).toStringList();
            if (!modeListContainsType(modes, *type))
                return false;
        }
        if (!searchText_.isEmpty()) {
            const QString haystack = model->data(index, DescriptorTreeModel::SearchTextRole).toString();
            const QStringList tokens = searchText_.split(QLatin1Char(' '), Qt::SkipEmptyParts);
            for (const QString& token : tokens) {
                if (!haystack.contains(token, Qt::CaseInsensitive))
                    return false;
            }
        }
        return true;
    }

private:
    QString searchText_;
    QString sourceFilter_;
    QString axisFilter_;
    QString shapeFilter_;
    QString modeKindFilter_;
    int requiredAnchorFilter_ = -1;
};

}  // namespace

struct SignalDisplayDialog::Impl {
    QPointer<model::TrajectorySignalCatalog> catalog;
    QPointer<model::DashboardSignalModel> activeModel;
    QPointer<model::DashboardPanelModel> panelModel;
    QPointer<DashboardSelectionController> selectionController;
    QPointer<model::AtomSelection> selection;
    const model::QtProtein* protein = nullptr;
    QPointer<model::Conformation> conformation;
    model::VisualizationContext visualizationContext;

    NearbySignalModel* anchorModel = nullptr;
    DescriptorTreeModel* descriptorModel = nullptr;
    DescriptorFilterProxyModel* descriptorProxy = nullptr;
    QSortFilterProxyModel* activeProxy = nullptr;

    QLabel* focusLabel = nullptr;
    QCheckBox* liveBox = nullptr;
    QDoubleSpinBox* radiusSpin = nullptr;
    QTableView* anchorView = nullptr;

    QLineEdit* candidateSearch = nullptr;
    QComboBox* sourceFilter = nullptr;
    QComboBox* axisFilter = nullptr;
    QComboBox* shapeFilter = nullptr;
    QComboBox* modeFilter = nullptr;
    QTreeView* candidateView = nullptr;
    QVector<ModeControl> candidateModes;
    QComboBox* panelCombo = nullptr;
    QLineEdit* newPanelEdit = nullptr;
    QPushButton* addButton = nullptr;

    QLineEdit* activeSearch = nullptr;
    QTableView* activeView = nullptr;
    QVector<ModeControl> activeModes;
    QPushButton* removeButton = nullptr;

    QLabel* statusLabel = nullptr;
    std::optional<std::size_t> focusAtom;
    int frame = 0;
};

SignalDisplayDialog::SignalDisplayDialog(QWidget* parent)
    : QDialog(parent)
    , d_(std::make_unique<Impl>())
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("SignalDisplayDialog"));
    setWindowTitle(QStringLiteral("Metric Picker"));
    resize(1280, 720);
    setMinimumSize(960, 560);

    d_->anchorModel = new NearbySignalModel(this);
    d_->descriptorModel = new DescriptorTreeModel(this);
    d_->descriptorProxy = new DescriptorFilterProxyModel(this);
    d_->descriptorProxy->setSourceModel(d_->descriptorModel);

    d_->activeProxy = new QSortFilterProxyModel(this);
    d_->activeProxy->setFilterCaseSensitivity(Qt::CaseInsensitive);
    d_->activeProxy->setSortCaseSensitivity(Qt::CaseInsensitive);
    d_->activeProxy->setFilterKeyColumn(-1);

    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(6, 6, 6, 6);
    root->setSpacing(6);

    auto* splitter = new QSplitter(Qt::Horizontal, this);
    root->addWidget(splitter, 1);

    auto* candidatesPanel = new QWidget(splitter);
    auto* candidatesLayout = new QVBoxLayout(candidatesPanel);
    candidatesLayout->setContentsMargins(0, 0, 0, 0);
    candidatesLayout->setSpacing(4);

    auto* contextGroup = new QGroupBox(QStringLiteral("Selection context"), candidatesPanel);
    auto* contextLayout = new QVBoxLayout(contextGroup);
    contextLayout->setContentsMargins(6, 4, 6, 6);
    contextLayout->setSpacing(4);

    auto* contextRow = new QHBoxLayout;
    d_->focusLabel = new QLabel(QStringLiteral("Focus: none"), contextGroup);
    d_->focusLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    contextRow->addWidget(d_->focusLabel, 1);
    d_->liveBox = new QCheckBox(QStringLiteral("Live"), contextGroup);
    d_->liveBox->setChecked(true);
    d_->liveBox->setToolTip(QStringLiteral("Keep this picker attached to the current atom focus and playback frame."));
    contextRow->addWidget(d_->liveBox);
    contextRow->addWidget(new QLabel(QStringLiteral("Radius"), contextGroup));
    d_->radiusSpin = new QDoubleSpinBox(contextGroup);
    d_->radiusSpin->setRange(1.0, 30.0);
    d_->radiusSpin->setDecimals(1);
    d_->radiusSpin->setSingleStep(0.5);
    d_->radiusSpin->setSuffix(QStringLiteral(" Angstrom"));
    d_->radiusSpin->setValue(5.0);
    contextRow->addWidget(d_->radiusSpin);
    contextLayout->addLayout(contextRow);

    d_->anchorView = new QTableView(contextGroup);
    configureTable(d_->anchorView);
    d_->anchorView->setModel(d_->anchorModel);
    d_->anchorView->setSortingEnabled(false);
    d_->anchorView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    d_->anchorView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    d_->anchorView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    d_->anchorView->setMaximumHeight(180);
    contextLayout->addWidget(d_->anchorView);
    candidatesLayout->addWidget(contextGroup);

    d_->candidateSearch = new QLineEdit(candidatesPanel);
    d_->candidateSearch->setClearButtonEnabled(true);
    d_->candidateSearch->setPlaceholderText(QStringLiteral("Search metrics for selected atom/residue"));
    candidatesLayout->addWidget(d_->candidateSearch);

    auto* filterRow = new QHBoxLayout;
    filterRow->setSpacing(4);
    d_->sourceFilter = new QComboBox(candidatesPanel);
    d_->axisFilter = new QComboBox(candidatesPanel);
    d_->shapeFilter = new QComboBox(candidatesPanel);
    d_->modeFilter = new QComboBox(candidatesPanel);
    for (QComboBox* combo : {d_->sourceFilter, d_->axisFilter, d_->shapeFilter, d_->modeFilter}) {
        combo->setMinimumContentsLength(10);
        combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
        filterRow->addWidget(combo, 1);
    }
    candidatesLayout->addLayout(filterRow);

    d_->candidateView = new QTreeView(candidatesPanel);
    configureTree(d_->candidateView);
    d_->candidateView->setModel(d_->descriptorProxy);
    d_->candidateView->header()->setSectionResizeMode(DescriptorTreeModel::NameColumn, QHeaderView::Stretch);
    d_->candidateView->header()->setSectionResizeMode(DescriptorTreeModel::FormColumn, QHeaderView::ResizeToContents);
    d_->candidateView->header()->setSectionResizeMode(DescriptorTreeModel::ShapeColumn, QHeaderView::ResizeToContents);
    d_->candidateView->header()->setSectionResizeMode(DescriptorTreeModel::DisplaysColumn, QHeaderView::ResizeToContents);
    d_->candidateView->header()->setSectionResizeMode(DescriptorTreeModel::UnitsColumn, QHeaderView::ResizeToContents);
    candidatesLayout->addWidget(d_->candidateView, 1);

    auto* addGroup = new QGroupBox(QStringLiteral("Add selected descriptor as"), candidatesPanel);
    auto* addLayout = new QHBoxLayout(addGroup);
    addLayout->setContentsMargins(6, 4, 6, 4);
    addLayout->setSpacing(6);
    for (model::VisualizationType type : allVisualizationTypes()) {
        auto* box = new QCheckBox(visualizationTypeLabel(type), addGroup);
        box->setEnabled(false);
        d_->candidateModes.push_back(ModeControl{type, box});
        addLayout->addWidget(box);
    }
    addLayout->addStretch(1);
    addLayout->addWidget(new QLabel(QStringLiteral("Panel"), addGroup));
    d_->panelCombo = new QComboBox(addGroup);
    d_->panelCombo->setMinimumContentsLength(12);
    d_->panelCombo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    addLayout->addWidget(d_->panelCombo);
    d_->newPanelEdit = new QLineEdit(addGroup);
    d_->newPanelEdit->setClearButtonEnabled(true);
    d_->newPanelEdit->setPlaceholderText(QStringLiteral("New panel name"));
    d_->newPanelEdit->setMaximumWidth(180);
    addLayout->addWidget(d_->newPanelEdit);
    d_->addButton = new QPushButton(QStringLiteral("Add Signal"), addGroup);
    d_->addButton->setEnabled(false);
    addLayout->addWidget(d_->addButton);
    candidatesLayout->addWidget(addGroup);

    auto* activePanel = new QWidget(splitter);
    auto* activeLayout = new QVBoxLayout(activePanel);
    activeLayout->setContentsMargins(0, 0, 0, 0);
    activeLayout->setSpacing(4);

    d_->activeSearch = new QLineEdit(activePanel);
    d_->activeSearch->setClearButtonEnabled(true);
    d_->activeSearch->setPlaceholderText(QStringLiteral("Search active signals"));
    activeLayout->addWidget(d_->activeSearch);

    d_->activeView = new QTableView(activePanel);
    configureTable(d_->activeView);
    d_->activeView->setModel(d_->activeProxy);
    d_->activeView->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    activeLayout->addWidget(d_->activeView, 1);

    auto* displayGroup = new QGroupBox(QStringLiteral("Display modes for active signal"), activePanel);
    auto* displayLayout = new QHBoxLayout(displayGroup);
    displayLayout->setContentsMargins(6, 4, 6, 4);
    displayLayout->setSpacing(6);
    for (model::VisualizationType type : allVisualizationTypes()) {
        auto* box = new QCheckBox(visualizationTypeLabel(type), displayGroup);
        box->setEnabled(false);
        d_->activeModes.push_back(ModeControl{type, box});
        displayLayout->addWidget(box);
    }
    displayLayout->addStretch(1);
    d_->removeButton = new QPushButton(QStringLiteral("Remove"), displayGroup);
    d_->removeButton->setEnabled(false);
    displayLayout->addWidget(d_->removeButton);
    activeLayout->addWidget(displayGroup);

    splitter->addWidget(candidatesPanel);
    splitter->addWidget(activePanel);
    splitter->setStretchFactor(0, 3);
    splitter->setStretchFactor(1, 2);

    d_->statusLabel = new QLabel(this);
    d_->statusLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    root->addWidget(d_->statusLabel);

    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
    root->addWidget(buttons);

    refillCombo(d_->sourceFilter, QStringLiteral("All sources"), {});
    refillCombo(d_->axisFilter, QStringLiteral("All axes"), {});
    refillCombo(d_->shapeFilter, QStringLiteral("All shapes"), {});
    d_->modeFilter->addItem(QStringLiteral("Any display"), QString());
    for (model::VisualizationType type : allVisualizationTypes())
        d_->modeFilter->addItem(visualizationTypeLabel(type), visualizationTypeKey(type));

    ACONNECT(d_->liveBox, &QCheckBox::toggled, this, &SignalDisplayDialog::onLiveToggled);
    ACONNECT(d_->radiusSpin, qOverload<double>(&QDoubleSpinBox::valueChanged),
             this, &SignalDisplayDialog::onRadiusChanged);
    ACONNECT(d_->anchorView->selectionModel(),
             &QItemSelectionModel::currentRowChanged,
             this,
             [this](const QModelIndex&, const QModelIndex&) { onAnchorSelectionChanged(); });
    ACONNECT(d_->anchorModel, &QAbstractItemModel::modelReset, this, [this]() {
        if (d_->anchorView && d_->anchorModel->rowCount() > 0 && !d_->anchorView->currentIndex().isValid())
            d_->anchorView->selectRow(0);
        onAnchorSelectionChanged();
    });

    ACONNECT(d_->descriptorModel, &QAbstractItemModel::modelReset, this, [this]() {
        refillCombo(d_->sourceFilter, QStringLiteral("All sources"),
                    d_->descriptorModel->uniqueValues(DescriptorTreeModel::SourceRole));
        refillCombo(d_->axisFilter, QStringLiteral("All axes"),
                    d_->descriptorModel->uniqueValues(DescriptorTreeModel::AxisRole));
        refillCombo(d_->shapeFilter, QStringLiteral("All shapes"),
                    d_->descriptorModel->uniqueValues(DescriptorTreeModel::ShapeRole));
    });
    ACONNECT(d_->candidateSearch, &QLineEdit::textChanged, this, [this](const QString& text) {
        d_->descriptorProxy->setSearchText(text);
        refreshCandidateTree();
    });
    ACONNECT(d_->sourceFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setSourceFilter(d_->sourceFilter->currentData().toString());
        refreshCandidateTree();
    });
    ACONNECT(d_->axisFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setAxisFilter(d_->axisFilter->currentData().toString());
        refreshCandidateTree();
    });
    ACONNECT(d_->shapeFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setShapeFilter(d_->shapeFilter->currentData().toString());
        refreshCandidateTree();
    });
    ACONNECT(d_->modeFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setModeKindFilter(d_->modeFilter->currentData().toString());
        refreshCandidateTree();
    });
    ACONNECT(d_->candidateView->selectionModel(),
             &QItemSelectionModel::currentRowChanged,
             this,
             [this](const QModelIndex&, const QModelIndex&) { onCandidateSelectionChanged(); });
    ACONNECT(d_->activeView->selectionModel(),
             &QItemSelectionModel::currentRowChanged,
             this,
             [this](const QModelIndex&, const QModelIndex&) { onActiveSelectionChanged(); });
    for (const ModeControl& control : d_->candidateModes)
        ACONNECT(control.box, &QCheckBox::toggled, this, &SignalDisplayDialog::onCandidateModeChanged);
    for (const ModeControl& control : d_->activeModes)
        ACONNECT(control.box, &QCheckBox::toggled, this, &SignalDisplayDialog::onActiveModeToggled);

    ACONNECT(d_->activeSearch, &QLineEdit::textChanged, this, [this](const QString& text) {
        d_->activeProxy->setFilterFixedString(text);
        onActiveSelectionChanged();
    });
    ACONNECT(d_->addButton, &QPushButton::clicked, this, &SignalDisplayDialog::onAddSelected);
    ACONNECT(d_->removeButton, &QPushButton::clicked, this, &SignalDisplayDialog::onRemoveActive);
    ACONNECT(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    refreshPanelTargets();
}

SignalDisplayDialog::~SignalDisplayDialog() = default;

void SignalDisplayDialog::setTrajectorySignalCatalog(model::TrajectorySignalCatalog* catalog) {
    ASSERT_THREAD(this);
    d_->catalog = catalog;
    refreshCatalog();
}

void SignalDisplayDialog::setDashboardSignalModel(model::DashboardSignalModel* model) {
    ASSERT_THREAD(this);
    d_->activeModel = model;
    d_->activeProxy->setSourceModel(model);
    d_->activeView->resizeColumnsToContents();
    onActiveSelectionChanged();
}

void SignalDisplayDialog::setDashboardPanelModel(model::DashboardPanelModel* panelModel) {
    ASSERT_THREAD(this);
    if (d_->panelModel)
        disconnect(d_->panelModel, nullptr, this, nullptr);
    d_->panelModel = panelModel;
    if (d_->panelModel) {
        ACONNECT(d_->panelModel.data(), &QAbstractItemModel::rowsInserted,
                 this, &SignalDisplayDialog::refreshPanelTargets);
        ACONNECT(d_->panelModel.data(), &QAbstractItemModel::rowsRemoved,
                 this, &SignalDisplayDialog::refreshPanelTargets);
        ACONNECT(d_->panelModel.data(), &QAbstractItemModel::modelReset,
                 this, &SignalDisplayDialog::refreshPanelTargets);
        ACONNECT(d_->panelModel.data(), &QAbstractItemModel::dataChanged,
                 this, [this](const QModelIndex&, const QModelIndex&, const QList<int>&) {
                     refreshPanelTargets();
                 });
        ACONNECT(d_->panelModel.data(), &model::DashboardPanelModel::activePanelChanged,
                 this, &SignalDisplayDialog::refreshPanelTargets);
    }
    refreshPanelTargets();
}

void SignalDisplayDialog::setDashboardSelectionController(DashboardSelectionController* controller) {
    ASSERT_THREAD(this);
    d_->selectionController = controller;
    onCandidateModeChanged();
}

void SignalDisplayDialog::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    d_->protein = protein;
    d_->conformation = conformation;
    if (d_->anchorModel)
        d_->anchorModel->setContext(protein, conformation);
}

void SignalDisplayDialog::setVisualizationContext(const model::VisualizationContext& ctx) {
    ASSERT_THREAD(this);
    d_->visualizationContext = ctx;
    onCandidateSelectionChanged();
    onActiveSelectionChanged();
}

void SignalDisplayDialog::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    if (d_->selection)
        disconnect(d_->selection, nullptr, this, nullptr);
    d_->selection = selection;
    d_->focusAtom.reset();
    if (d_->selection) {
        ACONNECT(d_->selection.data(), &model::AtomSelection::focusChanged,
                 this, &SignalDisplayDialog::onFocusChanged);
        ACONNECT(d_->selection.data(), &model::AtomSelection::cleared,
                 this, &SignalDisplayDialog::onSelectionCleared);
        if (d_->selection->hasFocus())
            onFocusChanged(d_->selection->focus());
        else
            onSelectionCleared();
    } else {
        onSelectionCleared();
    }
}

model::TrajectorySignalCatalog* SignalDisplayDialog::trajectorySignalCatalog() const {
    return d_->catalog.data();
}

model::DashboardSignalModel* SignalDisplayDialog::dashboardSignalModel() const {
    return d_->activeModel.data();
}

QJsonObject SignalDisplayDialog::pickerState() const {
    ASSERT_THREAD(this);
    QJsonObject out;
    out["open"] = isVisible();

    const QModelIndex proxyIndex = d_->candidateView ? d_->candidateView->currentIndex() : QModelIndex();
    const QModelIndex sourceIndex = proxyIndex.isValid() ? d_->descriptorProxy->mapToSource(proxyIndex)
                                                        : QModelIndex();
    const DescriptorRecord* candidateRecord = d_->descriptorModel->recordForIndex(sourceIndex);
    out["candidate_row"] = candidateRecord ? QJsonValue(proxyIndex.row()) : QJsonValue(QJsonValue::Null);
    out["candidate_descriptor"] = candidateRecord ? QJsonValue(candidateRecord->descriptor.id)
                                                   : QJsonValue(QJsonValue::Null);

    const QModelIndex anchorIndex = d_->anchorView ? d_->anchorView->currentIndex() : QModelIndex();
    const NearbySignalModel::Candidate* candidate = d_->anchorModel->candidateAt(anchorIndex);
    QJsonObject anchor = candidate ? signalAnchorToJson(candidate->anchor)
                                   : QJsonObject{{"kind", "none"}};
    anchor["row"] = candidate ? QJsonValue(anchorIndex.row()) : QJsonValue(QJsonValue::Null);
    anchor["label"] = candidate ? candidate->label : QString();
    anchor["axis"] = candidate ? model::ToString(axisForCandidate(*candidate))
                               : model::ToString(model::SignalAxis::None);
    anchor["distance_angstrom"] = candidate ? QJsonValue(candidate->distanceAngstrom)
                                            : QJsonValue(QJsonValue::Null);
    out["anchor"] = anchor;

    QJsonArray modes;
    for (const ModeControl& control : d_->candidateModes) {
        if (!control.box)
            continue;
        modes.append(QJsonObject{
            {"kind", visualizationTypeKey(control.type)},
            {"mode_id", control.box->property("modeId").toString()},
            {"enabled", control.box->isEnabled()},
            {"checked", control.box->isChecked()},
            {"reason", control.box->toolTip()},
        });
    }
    out["modes"] = modes;
    return out;
}

void SignalDisplayDialog::refreshCandidateTree() {
    ASSERT_THREAD(this);
    if (!d_->candidateView || !d_->descriptorProxy)
        return;
    // A typed search or an explicit source/axis/shape/display filter narrows the
    // set enough that revealing every matching leaf helps; the anchor-axis filter
    // is the normal resting state, so it does NOT force the concepts open.
    const bool expandConcepts =
        (d_->candidateSearch && !d_->candidateSearch->text().trimmed().isEmpty())
        || (d_->sourceFilter && d_->sourceFilter->currentIndex() > 0)
        || (d_->axisFilter && d_->axisFilter->currentIndex() > 0)
        || (d_->shapeFilter && d_->shapeFilter->currentIndex() > 0)
        || (d_->modeFilter && d_->modeFilter->currentIndex() > 0);
    expandCandidateTree(d_->candidateView, d_->descriptorProxy, QModelIndex(), expandConcepts);
    ensureCandidateRowSelected();
    onCandidateSelectionChanged();
}

bool SignalDisplayDialog::ensureCandidateRowSelected() {
    ASSERT_THREAD(this);
    if (!d_->candidateView || !d_->descriptorProxy)
        return false;

    auto* selection = d_->candidateView->selectionModel();
    const QModelIndex current = d_->candidateView->currentIndex();
    if (current.isValid()
        && d_->descriptorProxy->data(current, DescriptorTreeModel::RecordIndexRole).toInt() >= 0)
        return false;  // a real leaf is already current

    QModelIndex leaf = firstCandidateLeaf(d_->descriptorProxy, QModelIndex(), /*requireDisplayable=*/true);
    if (!leaf.isValid())
        leaf = firstCandidateLeaf(d_->descriptorProxy, QModelIndex(), /*requireDisplayable=*/false);
    if (!leaf.isValid()) {
        if (selection) {
            const QSignalBlocker blocker(selection);
            selection->clearSelection();
            selection->clearCurrentIndex();
        }
        return false;
    }

    for (QModelIndex ancestor = leaf.parent(); ancestor.isValid(); ancestor = ancestor.parent())
        d_->candidateView->expand(ancestor);
    if (selection) {
        const QSignalBlocker blocker(selection);
        selection->setCurrentIndex(leaf, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
    } else {
        d_->candidateView->setCurrentIndex(leaf);
    }
    d_->candidateView->scrollTo(leaf, QAbstractItemView::EnsureVisible);
    return true;
}

void SignalDisplayDialog::refreshCatalog() {
    ASSERT_THREAD(this);
    const QVector<model::SignalDescriptor> descriptors = d_->catalog ? d_->catalog->descriptorList()
                                                                     : QVector<model::SignalDescriptor>{};
    d_->descriptorModel->setDescriptors(descriptors);
    refreshCandidateTree();
    if (d_->statusLabel) {
        d_->statusLabel->setText(d_->catalog
                                     ? QStringLiteral("%1 catalog descriptors.").arg(descriptors.size())
                                     : QStringLiteral("No TrajectorySignalCatalog is connected."));
    }
}

void SignalDisplayDialog::setFrame(int frame) {
    ASSERT_THREAD(this);
    d_->frame = std::max(0, frame);
    if (d_->liveBox && d_->liveBox->isChecked() && d_->focusAtom.has_value()) {
        const std::size_t focusAtom = d_->focusAtom.value_or(0);
        d_->anchorModel->setAnchor(focusAtom, d_->frame);
    }
}

void SignalDisplayDialog::onFocusChanged(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    d_->focusAtom = atomIdx;
    QString focusText = QStringLiteral("Focus: atom %1").arg(static_cast<qulonglong>(atomIdx));
    if (d_->protein && atomIdx < d_->protein->atomCount()) {
        const auto& atom = d_->protein->atom(atomIdx);
        const QString atomName = d_->protein->atomNames(atomIdx).iupac;
        if (atom.residueIndex >= 0
            && static_cast<std::size_t>(atom.residueIndex) < d_->protein->residueCount()) {
            const std::size_t residueIndex = static_cast<std::size_t>(atom.residueIndex);
            const auto& residue = d_->protein->residue(residueIndex);
            const QString chain = residue.address.chainId.isEmpty()
                                      ? QString()
                                      : QStringLiteral("%1:").arg(residue.address.chainId);
            const QString residueName = d_->protein->residueLabel(residueIndex,
                                                                  model::NamingConvention::Iupac,
                                                                  model::NamingSource::Verbatim);
            focusText = QStringLiteral("Focus: %1%2%3:%4")
                            .arg(chain, residueName)
                            .arg(residue.address.residueNumber)
                            .arg(atomName);
        } else {
            focusText = QStringLiteral("Focus: atom %1:%2")
                            .arg(static_cast<qulonglong>(atomIdx))
                            .arg(atomName);
        }
    }
    if (d_->focusLabel)
        d_->focusLabel->setText(focusText);
    if (!d_->liveBox || d_->liveBox->isChecked())
        d_->anchorModel->setAnchor(atomIdx, d_->frame);
}

void SignalDisplayDialog::onSelectionCleared() {
    ASSERT_THREAD(this);
    d_->focusAtom.reset();
    if (d_->focusLabel)
        d_->focusLabel->setText(QStringLiteral("Focus: none"));
    if (d_->liveBox && d_->liveBox->isChecked())
        d_->anchorModel->clear();
    onAnchorSelectionChanged();
}

void SignalDisplayDialog::onLiveToggled(bool live) {
    ASSERT_THREAD(this);
    if (live && d_->focusAtom.has_value()) {
        const std::size_t focusAtom = d_->focusAtom.value_or(0);
        d_->anchorModel->setAnchor(focusAtom, d_->frame);
    }
}

void SignalDisplayDialog::onRadiusChanged(double radius) {
    ASSERT_THREAD(this);
    d_->anchorModel->setRadiusAngstrom(radius);
}

void SignalDisplayDialog::onAnchorSelectionChanged() {
    ASSERT_THREAD(this);
    const QModelIndex anchorIndex = d_->anchorView ? d_->anchorView->currentIndex() : QModelIndex();
    const NearbySignalModel::Candidate* candidate = d_->anchorModel->candidateAt(anchorIndex);
    d_->descriptorProxy->setRequiredAnchorFilter(candidate ? axisForCandidate(*candidate) : model::SignalAxis::None,
                                                 candidate != nullptr);
    refreshCandidateTree();
}

void SignalDisplayDialog::onCandidateSelectionChanged() {
    ASSERT_THREAD(this);
    const QModelIndex proxyIndex = d_->candidateView->currentIndex();
    const QModelIndex sourceIndex = proxyIndex.isValid() ? d_->descriptorProxy->mapToSource(proxyIndex) : QModelIndex();
    const DescriptorRecord* record = d_->descriptorModel->recordForIndex(sourceIndex);

    bool checkedOne = false;
    const model::VisualizationRegistry& registry = model::VisualizationRegistry::instance();
    for (const ModeControl& control : d_->candidateModes) {
        QSignalBlocker blocker(control.box);
        const QString modeId = record ? modeForType(record->displayModes, control.type) : QString();
        const QVector<const model::VisualizationDefinition*> structural =
            record ? registry.supporting(record->descriptor)
                   : QVector<const model::VisualizationDefinition*>{};
        const model::VisualizationDefinition* structuralDefinition =
            definitionForType(structural, control.type);
        const bool structurallySupported = !modeId.isEmpty() && structuralDefinition;
        const bool descriptorAvailable =
            !structuralDefinition || !record
            || structuralDefinition->isAvailable(d_->visualizationContext, record->descriptor);
        const QVector<const model::VisualizationDefinition*> offerable =
            record ? registry.visibleOfferable(d_->visualizationContext, record->descriptor)
                   : QVector<const model::VisualizationDefinition*>{};
        const bool supported = !modeId.isEmpty()
            && definitionForType(offerable, control.type) != nullptr;
        control.box->setProperty("modeId", modeId);
        control.box->setEnabled(supported);
        control.box->setChecked(supported && control.type == model::VisualizationType::TemporalStrip);
        control.box->setToolTip(supported
                                    ? QStringLiteral("Add display mode id '%1'.").arg(modeId)
                                    : candidateModeDisabledReason(modeId,
                                                                  structurallySupported,
                                                                  descriptorAvailable));
        checkedOne = checkedOne || (supported && control.box->isChecked());
    }
    if (!checkedOne && record) {
        for (const ModeControl& control : d_->candidateModes) {
            if (control.box->isEnabled()) {
                QSignalBlocker blocker(control.box);
                control.box->setChecked(true);
                break;
            }
        }
    }
    onCandidateModeChanged();

    QStringList enabledKinds;
    QStringList disabledKinds;
    for (const ModeControl& control : d_->candidateModes) {
        const QString kind = visualizationTypeKey(control.type);
        if (control.box->isEnabled()) {
            enabledKinds.push_back(kind);
        } else {
            disabledKinds.push_back(QStringLiteral("%1:%2").arg(kind, control.box->toolTip()));
        }
    }
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=picker_selector_state current_row=%1 enabled=[%2] disabled=[%3]")
               .arg(record ? QString::number(proxyIndex.row()) : QStringLiteral("none"),
                    enabledKinds.join(QStringLiteral(",")),
                    disabledKinds.join(QStringLiteral(",")));
}

void SignalDisplayDialog::onCandidateModeChanged() {
    ASSERT_THREAD(this);
    const QModelIndex anchorIndex = d_->anchorView ? d_->anchorView->currentIndex() : QModelIndex();
    const bool hasAnchor = d_->anchorModel->candidateAt(anchorIndex) != nullptr;
    bool hasMode = false;
    for (const ModeControl& control : d_->candidateModes) {
        if (control.box->isEnabled() && control.box->isChecked()
            && !control.box->property("modeId").toString().isEmpty()) {
            hasMode = true;
            break;
        }
    }
    d_->addButton->setEnabled(d_->activeModel && d_->selectionController && hasAnchor
                              && d_->candidateView->currentIndex().isValid() && hasMode);
}

void SignalDisplayDialog::onAddSelected() {
    ASSERT_THREAD(this);
    const QModelIndex proxyIndex = d_->candidateView->currentIndex();
    const QModelIndex sourceIndex = proxyIndex.isValid() ? d_->descriptorProxy->mapToSource(proxyIndex) : QModelIndex();
    const DescriptorRecord* record = d_->descriptorModel->recordForIndex(sourceIndex);
    if (!record || !d_->activeModel || !d_->selectionController)
        return;
    const QModelIndex anchorIndex = d_->anchorView ? d_->anchorView->currentIndex() : QModelIndex();
    const NearbySignalModel::Candidate* candidate = d_->anchorModel->candidateAt(anchorIndex);
    if (!candidate)
        return;

    const model::SignalAnchor anchor = anchorForCandidate(*candidate);

    QStringList displayModes;
    for (const ModeControl& control : d_->candidateModes) {
        const QString modeId = control.box->property("modeId").toString();
        if (control.box->isEnabled() && control.box->isChecked() && !modeId.isEmpty()) {
            model::DisplaySignalBinding binding;
            binding.sourceKind = record->descriptor.sourceKind;
            binding.descriptorId = record->descriptor.id;
            binding.conceptKey = record->descriptor.conceptKey;
            binding.displayModeId = modeId;
            binding.anchor = anchor;
            binding.followsFocus = false;
            if (d_->catalog && d_->catalog->canBind(binding))
                displayModes.push_back(modeId);
        }
    }
    displayModes.removeDuplicates();
    if (displayModes.isEmpty())
        return;

    const QString newPanelName = d_->newPanelEdit ? d_->newPanelEdit->text().trimmed() : QString();
    DashboardSelectionController::PanelTarget target;
    target.makeActive = true;
    if (d_->panelModel) {
        if (!newPanelName.isEmpty()) {
            target.newPanelName = newPanelName;
        } else if (d_->panelCombo && d_->panelCombo->currentIndex() >= 0) {
            target.panelId = d_->panelCombo->currentData().toUuid();
        }
    }

    const QString label = QStringLiteral("%1 - %2").arg(record->descriptor.label, candidate->label);
    int addedRefs = 0;
    const QUuid id = d_->selectionController->addMetric(record->descriptor,
                                                        anchor,
                                                        displayModes,
                                                        target,
                                                        false,
                                                        label,
                                                        &addedRefs);
    if (id.isNull())
        return;
    if (!newPanelName.isEmpty() && d_->newPanelEdit)
        d_->newPanelEdit->clear();
    if (d_->statusLabel) {
        d_->statusLabel->setText(QStringLiteral("Added '%1' on %2 with %3 (%4 display%5).")
                                     .arg(record->descriptor.label,
                                          candidate->label,
                                          displayModes.join(QStringLiteral(", ")),
                                          QString::number(addedRefs),
                                          addedRefs == 1 ? QString() : QStringLiteral("s")));
    }
    const QModelIndex sourceActive = d_->activeModel->indexForId(id);
    if (sourceActive.isValid()) {
        const QModelIndex proxyActive = d_->activeProxy->mapFromSource(sourceActive);
        if (proxyActive.isValid()) {
            d_->activeView->selectRow(proxyActive.row());
            d_->activeView->scrollTo(proxyActive);
        }
    }
}

void SignalDisplayDialog::refreshPanelTargets() {
    ASSERT_THREAD(this);
    if (!d_->panelCombo)
        return;
    const QUuid previous = d_->panelCombo->currentData().toUuid();
    const QUuid active = d_->panelModel ? d_->panelModel->activePanelId() : QUuid{};

    const QSignalBlocker block(d_->panelCombo);
    d_->panelCombo->clear();
    if (!d_->panelModel) {
        d_->panelCombo->addItem(QStringLiteral("Dashboard"), QUuid{});
        d_->panelCombo->setEnabled(false);
        return;
    }

    d_->panelCombo->setEnabled(true);
    int selectRow = -1;
    for (int row = 0; row < d_->panelModel->rowCount(); ++row) {
        const QModelIndex index = d_->panelModel->index(row, 0);
        const QString name = d_->panelModel->data(index, model::DashboardPanelModel::NameRole).toString();
        const QUuid id = d_->panelModel->data(index, model::DashboardPanelModel::UuidRole).toUuid();
        d_->panelCombo->addItem(name, id);
        if ((!previous.isNull() && id == previous) || (previous.isNull() && id == active))
            selectRow = row;
    }
    if (selectRow < 0)
        selectRow = std::max(0, d_->panelModel->activePanelRow());
    if (selectRow >= 0 && selectRow < d_->panelCombo->count())
        d_->panelCombo->setCurrentIndex(selectRow);
}

QModelIndex currentActiveSourceIndex(QTableView* view, QSortFilterProxyModel* proxy) {
    const QModelIndex proxyIndex = view ? view->currentIndex() : QModelIndex();
    if (!proxyIndex.isValid() || !proxy)
        return {};
    return proxy->mapToSource(proxyIndex);
}

const model::SignalDescriptor* descriptorForActiveSignal(model::TrajectorySignalCatalog* catalog,
                                                         model::DashboardSignalModel* model,
                                                         const QModelIndex& sourceIndex) {
    if (!catalog || !model || !sourceIndex.isValid())
        return nullptr;
    const QString descriptorId = model->data(model->index(sourceIndex.row(), 0),
                                             model::DashboardSignalModel::DescriptorIdRole).toString();
    return catalog->findDescriptor(descriptorId);
}

QUuid signalIdForActiveRow(model::DashboardSignalModel* model, const QModelIndex& sourceIndex) {
    if (!model || !sourceIndex.isValid())
        return {};
    return model->data(model->index(sourceIndex.row(), 0), model::DashboardSignalModel::UuidRole).toUuid();
}

QStringList displayModesForActiveRow(model::DashboardSignalModel* model, const QModelIndex& sourceIndex) {
    if (!model || !sourceIndex.isValid())
        return {};
    return model->data(model->index(sourceIndex.row(), 0),
                       model::DashboardSignalModel::DisplayModesRole).toStringList();
}

void SignalDisplayDialog::onActiveSelectionChanged() {
    ASSERT_THREAD(this);
    const QModelIndex sourceIndex = currentActiveSourceIndex(d_->activeView, d_->activeProxy);
    const bool hasActive = sourceIndex.isValid() && d_->activeModel;
    const model::SignalDescriptor* descriptor =
        descriptorForActiveSignal(d_->catalog, d_->activeModel, sourceIndex);
    const QStringList supportedModes = descriptor ? model::AllDisplayModes(*descriptor) : QStringList{};
    const QStringList enabledModes = displayModesForActiveRow(d_->activeModel, sourceIndex);
    const model::VisualizationRegistry& registry = model::VisualizationRegistry::instance();
    const QVector<const model::VisualizationDefinition*> offerable =
        descriptor ? registry.visibleOfferable(d_->visualizationContext, *descriptor)
                   : QVector<const model::VisualizationDefinition*>{};

    for (const ModeControl& control : d_->activeModes) {
        QSignalBlocker blocker(control.box);
        QString modeId = modeForType(supportedModes, control.type);
        if (modeId.isEmpty())
            modeId = modeForType(enabledModes, control.type);

        const bool supported = hasActive && !modeId.isEmpty()
            && definitionForType(offerable, control.type) != nullptr;
        control.box->setProperty("modeId", modeId);
        control.box->setEnabled(supported);
        control.box->setChecked(supported && modeListContainsType(enabledModes, control.type));
        control.box->setToolTip(supported
                                    ? QStringLiteral("Toggle display mode id '%1'.").arg(modeId)
                                    : QStringLiteral("This display mode does not have an implemented visible renderer."));
    }
    d_->removeButton->setEnabled(hasActive);
}

void SignalDisplayDialog::onActiveModeToggled(bool checked) {
    ASSERT_THREAD(this);
    auto* box = qobject_cast<QCheckBox*>(sender());
    if (!box || !d_->activeModel || !d_->selectionController)
        return;
    const QModelIndex sourceIndex = currentActiveSourceIndex(d_->activeView, d_->activeProxy);
    const QUuid id = signalIdForActiveRow(d_->activeModel, sourceIndex);
    const QString modeId = box->property("modeId").toString();
    if (id.isNull() || modeId.isEmpty())
        return;

    if (!d_->selectionController->setMetricMode(id, modeId, checked)) {
        QSignalBlocker blocker(box);
        box->setChecked(!checked);
        return;
    }
    if (d_->statusLabel)
        d_->statusLabel->setText(QStringLiteral("Updated display mode '%1'.").arg(modeId));
}

void SignalDisplayDialog::onRemoveActive() {
    ASSERT_THREAD(this);
    if (!d_->activeModel || !d_->selectionController)
        return;
    const QModelIndex sourceIndex = currentActiveSourceIndex(d_->activeView, d_->activeProxy);
    const QUuid id = signalIdForActiveRow(d_->activeModel, sourceIndex);
    if (id.isNull())
        return;
    const bool removed = d_->selectionController->removeMetric(id);
    if (d_->statusLabel)
        d_->statusLabel->setText(removed ? QStringLiteral("Removed active signal.")
                                         : QStringLiteral("Could not remove active signal."));
}

}  // namespace h5reader::app
