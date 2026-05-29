#include "SignalDisplayDialog.h"

#include "NearbySignalModel.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/TrajectorySignalCatalog.h"

#include <QAbstractItemView>
#include <QAbstractTableModel>
#include <QBrush>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QItemSelectionModel>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QPointer>
#include <QSet>
#include <QSignalBlocker>
#include <QSortFilterProxyModel>
#include <QSplitter>
#include <QTableView>
#include <QUuid>
#include <QVBoxLayout>

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <utility>

namespace h5reader::app {

namespace {

enum class DisplayModeKind : std::uint8_t {
    Strip,
    Spectrum,
    Table,
    ColorMap,
    TensorGlyph,
    // Phase A-G additions for the new static-panel modes (see
    // h5-reader/notes/SCOPE_NEW_TRS_2026-05-29.md). Each maps to a
    // specific panel subclass via the controller's panel-build
    // dispatch (isPanelMode in DashboardDisplayController.cpp).
    BarSequence,   // SequenceBarPanel — iRED S², Reorient relaxation, Dihedral corr_time
    CurveLag,      // LagDecayPanel    — KernelDynamics ACF, Reorient TCFs, Dihedral ACF
    ChordCoupling, // ChordCouplingPanel — KernelCoherence Pearson matrix
    FixedFreq,     // FixedFreqPanel    — Reorient J(ω) at 5 KTB Larmor combinations (L-3b)
};

struct ModeControl {
    DisplayModeKind kind = DisplayModeKind::Strip;
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

QString modeKindLabel(DisplayModeKind kind) {
    switch (kind) {
    case DisplayModeKind::Strip:
        return QStringLiteral("Strip");
    case DisplayModeKind::Spectrum:
        return QStringLiteral("Spectrum");
    case DisplayModeKind::Table:
        return QStringLiteral("Table");
    case DisplayModeKind::ColorMap:
        return QStringLiteral("Color map");
    case DisplayModeKind::TensorGlyph:
        return QStringLiteral("Glyph / overlay");
    case DisplayModeKind::BarSequence:
        return QStringLiteral("Bar (sequence)");
    case DisplayModeKind::CurveLag:
        return QStringLiteral("Curve (lag)");
    case DisplayModeKind::ChordCoupling:
        return QStringLiteral("Chord (coupling)");
    case DisplayModeKind::FixedFreq:
        return QStringLiteral("Fixed-freq J(ω)");
    }
    return {};
}

QString modeKindKey(DisplayModeKind kind) {
    switch (kind) {
    case DisplayModeKind::Strip:
        return QStringLiteral("strip");
    case DisplayModeKind::Spectrum:
        return QStringLiteral("spectrum");
    case DisplayModeKind::Table:
        return QStringLiteral("table");
    case DisplayModeKind::ColorMap:
        return QStringLiteral("colorMap");
    case DisplayModeKind::TensorGlyph:
        return QStringLiteral("tensorGlyph");
    case DisplayModeKind::BarSequence:
        return QStringLiteral("barSequence");
    case DisplayModeKind::CurveLag:
        return QStringLiteral("curveLag");
    case DisplayModeKind::ChordCoupling:
        return QStringLiteral("chordCoupling");
    case DisplayModeKind::FixedFreq:
        return QStringLiteral("fixedFreq");
    }
    return {};
}

QString canonicalModeId(DisplayModeKind kind) {
    switch (kind) {
    case DisplayModeKind::Strip:
        return QStringLiteral("strip.scalar");
    case DisplayModeKind::Spectrum:
        return QStringLiteral("strip.spectrum");
    case DisplayModeKind::Table:
        return QStringLiteral("static.table");
    case DisplayModeKind::ColorMap:
        return QStringLiteral("static.atomColor");
    case DisplayModeKind::TensorGlyph:
        return QStringLiteral("static.tensor");
    case DisplayModeKind::BarSequence:
        return QStringLiteral("static.bar.sequence");
    case DisplayModeKind::CurveLag:
        return QStringLiteral("static.curve.lag.animated");
    case DisplayModeKind::ChordCoupling:
        return QStringLiteral("static.chord.coupling");
    case DisplayModeKind::FixedFreq:
        return QStringLiteral("static.fixed_freq");
    }
    return {};
}

bool modeMatchesKind(const QString& modeId, DisplayModeKind kind) {
    const QString lower = modeId.toLower();
    switch (kind) {
    case DisplayModeKind::Strip:
        return lower.startsWith(QStringLiteral("strip."))
            && !lower.contains(QStringLiteral("spectrum"))
            && !lower.contains(QStringLiteral("fft"));
    case DisplayModeKind::Spectrum:
        return lower.contains(QStringLiteral("spectrum")) || lower.contains(QStringLiteral("fft"));
    case DisplayModeKind::Table:
        return lower.contains(QStringLiteral("table"));
    case DisplayModeKind::ColorMap:
        return lower.contains(QStringLiteral("color")) || lower.contains(QStringLiteral("map"))
            || lower.contains(QStringLiteral("band"));
    case DisplayModeKind::TensorGlyph:
        return lower.contains(QStringLiteral("glyph")) || lower.contains(QStringLiteral("overlay"))
            || (lower.startsWith(QStringLiteral("static.")) && lower.contains(QStringLiteral("tensor")));
    case DisplayModeKind::BarSequence:
        return lower == QStringLiteral("static.bar.sequence");
    case DisplayModeKind::CurveLag:
        return lower == QStringLiteral("static.curve.lag.animated");
    case DisplayModeKind::ChordCoupling:
        return lower == QStringLiteral("static.chord.coupling");
    case DisplayModeKind::FixedFreq:
        return lower == QStringLiteral("static.fixed_freq");
    }
    return false;
}

QString modeForKind(const QStringList& modes, DisplayModeKind kind) {
    QString preferred = canonicalModeId(kind);
    if (modes.contains(preferred))
        return preferred;
    for (const QString& mode : modes) {
        if (modeMatchesKind(mode, kind))
            return mode;
    }
    return {};
}

bool modeListContainsKind(const QStringList& modes, DisplayModeKind kind) {
    return !modeForKind(modes, kind).isEmpty();
}

QString modeSummary(const QStringList& displayModes) {
    QStringList labels;
    // modeSummary uses the same Table-last order as allModeKinds so
    // the summary labels group concrete renderers before fallbacks.
    for (DisplayModeKind kind : {DisplayModeKind::Strip,
                                 DisplayModeKind::Spectrum,
                                 DisplayModeKind::TensorGlyph,
                                 DisplayModeKind::BarSequence,
                                 DisplayModeKind::CurveLag,
                                 DisplayModeKind::ChordCoupling,
                                 DisplayModeKind::FixedFreq,
                                 DisplayModeKind::ColorMap,
                                 DisplayModeKind::Table}) {
        if (modeListContainsKind(displayModes, kind))
            labels.push_back(modeKindLabel(kind));
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

// Delegates to model::AxisCanSatisfy so the dialog's widening rules
// stay in lockstep with the controller-side AnchorMatchesAxis. Earlier
// versions of this dialog kept a local copy that missed the
// Residue → BondVector widening, so iRED / Reorient descriptors were
// silently filtered out of the candidate table when a residue was the
// active anchor (Codex NOW-3, 2026-05-29).
bool anchorAxisCanSatisfy(model::SignalAxis selectedAxis, model::SignalAxis requiredAxis) {
    return model::AxisCanSatisfy(selectedAxis, requiredAxis);
}

std::array<DisplayModeKind, 9> allModeKinds() {
    // Codex NOW-3 (2026-05-29) reorder: Table moves LAST so the
    // dialog's default-checked-mode walk
    // (onCandidateSelectionChanged) doesn't pick the table fallback
    // ahead of a concrete renderer. Strip + Spectrum stay first
    // because strip-capable descriptors should default to strip
    // mode; the concrete static renderers come next; ColorMap +
    // Table land last as the visualisation-of-last-resort fallbacks.
    return {DisplayModeKind::Strip,
            DisplayModeKind::Spectrum,
            DisplayModeKind::TensorGlyph,
            DisplayModeKind::BarSequence,
            DisplayModeKind::CurveLag,
            DisplayModeKind::ChordCoupling,
            DisplayModeKind::FixedFreq,
            DisplayModeKind::ColorMap,
            DisplayModeKind::Table};
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

class DescriptorTableModel final : public QAbstractTableModel {
public:
    enum Column : std::uint8_t {
        SourceColumn,
        AxisColumn,
        ShapeColumn,
        SignalColumn,
        ModesColumn,
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
    };

    explicit DescriptorTableModel(QObject* parent = nullptr)
        : QAbstractTableModel(parent)
    {
        CENSUS_REGISTER(this);
        setObjectName(QStringLiteral("SignalDescriptorTableModel"));
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const override {
        return parent.isValid() ? 0 : static_cast<int>(records_.size());
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const override {
        return parent.isValid() ? 0 : ColumnCount;
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override {
        if (!index.isValid() || index.row() < 0 || index.row() >= records_.size())
            return {};
        const DescriptorRecord& record = records_.at(index.row());
        if (role == Qt::DisplayRole || role == Qt::EditRole) {
            switch (index.column()) {
            case SourceColumn:
                return record.source;
            case AxisColumn:
                return record.axis;
            case ShapeColumn:
                return record.valueShape;
            case SignalColumn:
                return record.descriptor.label;
            case ModesColumn:
                return modeSummary(record.displayModes);
            case UnitsColumn:
                return record.units;
            default:
                return {};
            }
        }
        if (role == Qt::ToolTipRole) {
            return QStringLiteral("%1\nid: %2\nconcept: %3\nmodes: %4")
                .arg(record.descriptor.label,
                     record.descriptor.id,
                     record.descriptor.conceptKey.isEmpty() ? QStringLiteral("-") : record.descriptor.conceptKey,
                     record.displayModes.join(QStringLiteral(", ")));
        }
        if (role == Qt::ForegroundRole && record.displayModes.isEmpty())
            return QBrush(Qt::darkGray);
        if (role == DisplayModesRole)
            return record.displayModes;
        if (role == SearchTextRole)
            return record.searchText();
        if (role == SourceRole)
            return record.source;
        if (role == AxisRole)
            return record.axis;
        if (role == ShapeRole)
            return record.valueShape;
        if (role == RequiredAnchorAxisRole)
            return static_cast<int>(record.descriptor.requiredAnchor);
        return {};
    }

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override {
        if (orientation != Qt::Horizontal || role != Qt::DisplayRole)
            return {};
        switch (section) {
        case SourceColumn:
            return QStringLiteral("Source");
        case AxisColumn:
            return QStringLiteral("Axis");
        case ShapeColumn:
            return QStringLiteral("Shape");
        case SignalColumn:
            return QStringLiteral("Descriptor");
        case ModesColumn:
            return QStringLiteral("Displays");
        case UnitsColumn:
            return QStringLiteral("Units");
        default:
            return {};
        }
    }

    void setDescriptors(const QVector<model::SignalDescriptor>& descriptors) {
        QVector<DescriptorRecord> records;
        records.reserve(descriptors.size());
        for (const model::SignalDescriptor& descriptor : descriptors)
            records.push_back(recordFromDescriptor(descriptor));

        beginResetModel();
        records_ = std::move(records);
        endResetModel();
    }

    const DescriptorRecord* recordAt(const QModelIndex& index) const {
        if (!index.isValid() || index.row() < 0 || index.row() >= records_.size())
            return nullptr;
        return &records_.at(index.row());
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
        QStringList values;
        values.reserve(seen.size());
        for (const QString& value : seen)
            values.push_back(value);
        values.sort(Qt::CaseInsensitive);
        return values;
    }

private:
    QVector<DescriptorRecord> records_;
};

class DescriptorFilterProxyModel final : public QSortFilterProxyModel {
public:
    explicit DescriptorFilterProxyModel(QObject* parent = nullptr)
        : QSortFilterProxyModel(parent)
    {
        CENSUS_REGISTER(this);
        setObjectName(QStringLiteral("SignalDescriptorFilterProxyModel"));
        setFilterCaseSensitivity(Qt::CaseInsensitive);
        setSortCaseSensitivity(Qt::CaseInsensitive);
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
        const auto textForRole = [&](int role) {
            return model->data(index, role).toString();
        };
        if (!sourceFilter_.isEmpty() && textForRole(DescriptorTableModel::SourceRole) != sourceFilter_)
            return false;
        if (!axisFilter_.isEmpty() && textForRole(DescriptorTableModel::AxisRole) != axisFilter_)
            return false;
        if (!shapeFilter_.isEmpty() && textForRole(DescriptorTableModel::ShapeRole) != shapeFilter_)
            return false;
        if (requiredAnchorFilter_ >= 0) {
            const auto selectedAxis = static_cast<model::SignalAxis>(requiredAnchorFilter_);
            const auto requiredAxis =
                static_cast<model::SignalAxis>(
                    model->data(index, DescriptorTableModel::RequiredAnchorAxisRole).toInt());
            if (!anchorAxisCanSatisfy(selectedAxis, requiredAxis))
                return false;
        }
        if (!modeKindFilter_.isEmpty()) {
            const QStringList modes = model->data(index, DescriptorTableModel::DisplayModesRole).toStringList();
            bool matched = false;
            for (const QString& mode : modes) {
                if ((modeKindFilter_ == QStringLiteral("strip") && modeMatchesKind(mode, DisplayModeKind::Strip))
                    || (modeKindFilter_ == QStringLiteral("spectrum") && modeMatchesKind(mode, DisplayModeKind::Spectrum))
                    || (modeKindFilter_ == QStringLiteral("table") && modeMatchesKind(mode, DisplayModeKind::Table))
                    || (modeKindFilter_ == QStringLiteral("colorMap") && modeMatchesKind(mode, DisplayModeKind::ColorMap))
                    || (modeKindFilter_ == QStringLiteral("tensorGlyph") && modeMatchesKind(mode, DisplayModeKind::TensorGlyph))
                    || (modeKindFilter_ == QStringLiteral("barSequence") && modeMatchesKind(mode, DisplayModeKind::BarSequence))
                    || (modeKindFilter_ == QStringLiteral("curveLag") && modeMatchesKind(mode, DisplayModeKind::CurveLag))
                    || (modeKindFilter_ == QStringLiteral("chordCoupling") && modeMatchesKind(mode, DisplayModeKind::ChordCoupling))
                    // Codex NOW-3 (2026-05-29): fixedFreq filter arm
                    // was missing; the combo entry existed but the
                    // filter dropped instead of matching.
                    || (modeKindFilter_ == QStringLiteral("fixedFreq") && modeMatchesKind(mode, DisplayModeKind::FixedFreq))) {
                    matched = true;
                    break;
                }
            }
            if (!matched)
                return false;
        }
        if (!searchText_.isEmpty()) {
            const QString haystack = model->data(index, DescriptorTableModel::SearchTextRole).toString();
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
    QPointer<model::AtomSelection> selection;
    const model::QtProtein* protein = nullptr;
    QPointer<model::Conformation> conformation;

    NearbySignalModel* anchorModel = nullptr;
    DescriptorTableModel* descriptorModel = nullptr;
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
    QTableView* candidateView = nullptr;
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
    d_->descriptorModel = new DescriptorTableModel(this);
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

    d_->candidateView = new QTableView(candidatesPanel);
    configureTable(d_->candidateView);
    d_->candidateView->setModel(d_->descriptorProxy);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::SourceColumn, QHeaderView::ResizeToContents);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::AxisColumn, QHeaderView::ResizeToContents);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::ShapeColumn, QHeaderView::ResizeToContents);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::SignalColumn, QHeaderView::Stretch);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::ModesColumn, QHeaderView::ResizeToContents);
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::UnitsColumn, QHeaderView::ResizeToContents);
    d_->candidateView->sortByColumn(DescriptorTableModel::SourceColumn, Qt::AscendingOrder);
    candidatesLayout->addWidget(d_->candidateView, 1);

    auto* addGroup = new QGroupBox(QStringLiteral("Add selected descriptor as"), candidatesPanel);
    auto* addLayout = new QHBoxLayout(addGroup);
    addLayout->setContentsMargins(6, 4, 6, 4);
    addLayout->setSpacing(6);
    for (DisplayModeKind kind : allModeKinds()) {
        auto* box = new QCheckBox(modeKindLabel(kind), addGroup);
        box->setEnabled(false);
        d_->candidateModes.push_back(ModeControl{kind, box});
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
    for (DisplayModeKind kind : allModeKinds()) {
        auto* box = new QCheckBox(modeKindLabel(kind), displayGroup);
        box->setEnabled(false);
        d_->activeModes.push_back(ModeControl{kind, box});
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
    for (DisplayModeKind kind : allModeKinds())
        d_->modeFilter->addItem(modeKindLabel(kind), modeKindKey(kind));

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
                    d_->descriptorModel->uniqueValues(DescriptorTableModel::SourceRole));
        refillCombo(d_->axisFilter, QStringLiteral("All axes"),
                    d_->descriptorModel->uniqueValues(DescriptorTableModel::AxisRole));
        refillCombo(d_->shapeFilter, QStringLiteral("All shapes"),
                    d_->descriptorModel->uniqueValues(DescriptorTableModel::ShapeRole));
        onCandidateSelectionChanged();
    });
    ACONNECT(d_->candidateSearch, &QLineEdit::textChanged, this, [this](const QString& text) {
        d_->descriptorProxy->setSearchText(text);
        onCandidateSelectionChanged();
    });
    ACONNECT(d_->sourceFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setSourceFilter(d_->sourceFilter->currentData().toString());
        onCandidateSelectionChanged();
    });
    ACONNECT(d_->axisFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setAxisFilter(d_->axisFilter->currentData().toString());
        onCandidateSelectionChanged();
    });
    ACONNECT(d_->shapeFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setShapeFilter(d_->shapeFilter->currentData().toString());
        onCandidateSelectionChanged();
    });
    ACONNECT(d_->modeFilter, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        d_->descriptorProxy->setModeKindFilter(d_->modeFilter->currentData().toString());
        onCandidateSelectionChanged();
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

void SignalDisplayDialog::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    d_->protein = protein;
    d_->conformation = conformation;
    if (d_->anchorModel)
        d_->anchorModel->setContext(protein, conformation);
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

void SignalDisplayDialog::refreshCatalog() {
    ASSERT_THREAD(this);
    const QVector<model::SignalDescriptor> descriptors = d_->catalog ? d_->catalog->descriptorList()
                                                                     : QVector<model::SignalDescriptor>{};
    d_->descriptorModel->setDescriptors(descriptors);
    d_->candidateView->resizeColumnsToContents();
    d_->candidateView->horizontalHeader()->setSectionResizeMode(DescriptorTableModel::SignalColumn, QHeaderView::Stretch);
    if (d_->descriptorProxy->rowCount() > 0)
        d_->candidateView->selectRow(0);
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
    if (d_->descriptorProxy->rowCount() > 0) {
        const QModelIndex current = d_->candidateView->currentIndex();
        if (!current.isValid() || !d_->descriptorProxy->mapToSource(current).isValid())
            d_->candidateView->selectRow(0);
    }
    onCandidateSelectionChanged();
}

void SignalDisplayDialog::onCandidateSelectionChanged() {
    ASSERT_THREAD(this);
    const QModelIndex proxyIndex = d_->candidateView->currentIndex();
    const QModelIndex sourceIndex = proxyIndex.isValid() ? d_->descriptorProxy->mapToSource(proxyIndex) : QModelIndex();
    const DescriptorRecord* record = d_->descriptorModel->recordAt(sourceIndex);

    bool checkedOne = false;
    for (const ModeControl& control : d_->candidateModes) {
        QSignalBlocker blocker(control.box);
        const QString modeId = record ? modeForKind(record->displayModes, control.kind) : QString();
        const bool supported = !modeId.isEmpty();
        control.box->setProperty("modeId", modeId);
        control.box->setEnabled(supported);
        control.box->setChecked(supported && control.kind == DisplayModeKind::Strip);
        control.box->setToolTip(supported
                                    ? QStringLiteral("Add display mode id '%1'.").arg(modeId)
                                    : QStringLiteral("This descriptor does not advertise that display mode."));
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
    d_->addButton->setEnabled(d_->activeModel && hasAnchor && d_->candidateView->currentIndex().isValid() && hasMode);
}

void SignalDisplayDialog::onAddSelected() {
    ASSERT_THREAD(this);
    const QModelIndex proxyIndex = d_->candidateView->currentIndex();
    const QModelIndex sourceIndex = proxyIndex.isValid() ? d_->descriptorProxy->mapToSource(proxyIndex) : QModelIndex();
    const DescriptorRecord* record = d_->descriptorModel->recordAt(sourceIndex);
    if (!record || !d_->activeModel)
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

    QUuid panelId;
    const QString newPanelName = d_->newPanelEdit ? d_->newPanelEdit->text().trimmed() : QString();
    if (d_->panelModel) {
        if (!newPanelName.isEmpty()) {
            panelId = d_->panelModel->addPanel(newPanelName);
            d_->panelModel->setActivePanel(panelId);
            if (d_->newPanelEdit)
                d_->newPanelEdit->clear();
        } else if (d_->panelCombo && d_->panelCombo->currentIndex() >= 0) {
            panelId = d_->panelCombo->currentData().toUuid();
        }
        if (panelId.isNull())
            panelId = d_->panelModel->activePanelId();
    }

    const QString label = QStringLiteral("%1 - %2").arg(record->descriptor.label, candidate->label);
    const QUuid id = d_->activeModel->addSignal(record->descriptor,
                                                anchor,
                                                QString(),
                                                displayModes,
                                                false,
                                                label);
    int addedRefs = 0;
    if (d_->panelModel && !panelId.isNull()) {
        const QVector<model::DashboardDisplayRef> refs =
            model::DisplayRefsForSignal(id, record->descriptor, displayModes);
        for (const model::DashboardDisplayRef& ref : refs) {
            if (d_->panelModel->addDisplayRef(panelId, ref))
                ++addedRefs;
        }
        d_->panelModel->setActivePanel(panelId);
    }
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

    for (const ModeControl& control : d_->activeModes) {
        QSignalBlocker blocker(control.box);
        QString modeId = modeForKind(supportedModes, control.kind);
        if (modeId.isEmpty())
            modeId = modeForKind(enabledModes, control.kind);
        if (modeId.isEmpty())
            modeId = canonicalModeId(control.kind);

        const bool supported = hasActive && (supportedModes.isEmpty()
                                             || modeListContainsKind(supportedModes, control.kind));
        control.box->setProperty("modeId", modeId);
        control.box->setEnabled(supported);
        control.box->setChecked(hasActive && modeListContainsKind(enabledModes, control.kind));
        control.box->setToolTip(supported
                                    ? QStringLiteral("Toggle display mode id '%1'.").arg(modeId)
                                    : QStringLiteral("The selected signal does not advertise this display mode."));
    }
    d_->removeButton->setEnabled(hasActive);
}

void SignalDisplayDialog::onActiveModeToggled(bool checked) {
    ASSERT_THREAD(this);
    auto* box = qobject_cast<QCheckBox*>(sender());
    if (!box || !d_->activeModel)
        return;
    const QModelIndex sourceIndex = currentActiveSourceIndex(d_->activeView, d_->activeProxy);
    const QUuid id = signalIdForActiveRow(d_->activeModel, sourceIndex);
    const QString modeId = box->property("modeId").toString();
    if (id.isNull() || modeId.isEmpty())
        return;

    if (!d_->activeModel->toggleDisplayMode(id, modeId, checked)) {
        QSignalBlocker blocker(box);
        box->setChecked(!checked);
        return;
    }
    if (d_->panelModel) {
        const model::SignalDescriptor* descriptor =
            descriptorForActiveSignal(d_->catalog, d_->activeModel, sourceIndex);
        const QUuid panelId = d_->panelModel->activePanelId();
        if (descriptor && !panelId.isNull()) {
            const QVector<model::DashboardDisplayRef> refs =
                model::DisplayRefsForSignal(id, *descriptor, {modeId});
            if (checked) {
                d_->panelModel->addDisplayRefs(panelId, refs);
            } else {
                d_->panelModel->removeDisplayRefsForSignalMode(id, modeId);
            }
        }
    }
    if (d_->statusLabel)
        d_->statusLabel->setText(QStringLiteral("Updated display mode '%1'.").arg(modeId));
}

void SignalDisplayDialog::onRemoveActive() {
    ASSERT_THREAD(this);
    if (!d_->activeModel)
        return;
    const QModelIndex sourceIndex = currentActiveSourceIndex(d_->activeView, d_->activeProxy);
    const QUuid id = signalIdForActiveRow(d_->activeModel, sourceIndex);
    if (id.isNull())
        return;
    const bool removed = d_->activeModel->removeSignal(id);
    if (d_->statusLabel)
        d_->statusLabel->setText(removed ? QStringLiteral("Removed active signal.")
                                         : QStringLiteral("Could not remove active signal."));
}

}  // namespace h5reader::app
