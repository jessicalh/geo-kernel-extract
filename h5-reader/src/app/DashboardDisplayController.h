#pragma once

#include "../model/DashboardSignal.h"
#include "../model/SignalDictionary.h"
#include "../model/SignalTimeSeries.h"
#include "AbstractStripPanel.h"

#include <QObject>
#include <QColor>
#include <QPointer>
#include <QString>
#include <QVector>

#include <functional>
#include <memory>
#include <vector>
#include <cstddef>

namespace h5reader::model {
class AtomSelection;
class Conformation;
class DashboardPanelModel;
class DashboardSignalModel;
class DftShieldingStore;
class QtProtein;
class TrajectorySignalCatalog;
}

namespace h5reader::app {

class SceneRevealOverlay;

struct DashboardSmokeSummary {
    struct SeriesSparseness {
        QString signalLabel;
        QString descriptorId;
        QString conceptKey;
        QString sourceKind;
        QString storagePath;
        QString displayModeId;
        QString channelId;
        QString channelLabel;
        long long samples = 0;
        long long validSamples = 0;
        long long gapSamples = 0;
        long long invalidSamples = 0;
        long long pendingGapSamples = 0;
        long long sourceAbsentGapSamples = 0;
        long long frameSourceAbsentGapSamples = 0;
        long long sourceMaskOffGapSamples = 0;
        long long anchorUnavailableGapSamples = 0;
        long long notApplicableGapSamples = 0;
        long long nanSentinelGapSamples = 0;
        long long malformedSourceGapSamples = 0;
        int firstValidFrame = -1;
        int lastValidFrame = -1;
        int longestValidRun = 0;
        int longestGapRun = 0;
    };

    int seriesCount = 0;
    int seriesWithSamples = 0;
    int seriesWithValidSamples = 0;
    int seriesPendingOnly = 0;
    int denseSeries = 0;
    int sparseSeries = 0;
    int allGapSeries = 0;
    int seriesWithFrameSourceAbsentGaps = 0;
    int frameNpySeriesWithFrameSourceAbsentGaps = 0;
    int orcaDftSeriesWithFrameSourceAbsentGaps = 0;
    int seriesWithSourceAbsentGaps = 0;
    int seriesWithAnchorUnavailableGaps = 0;
    int seriesWithMismatchedBuffers = 0;
    int maxLongestGapRun = 0;
    long long samples = 0;
    long long channelValues = 0;
    long long channelValidity = 0;
    long long validSamples = 0;
    long long gapSamples = 0;
    long long pendingGapSamples = 0;
    long long sourceAbsentGapSamples = 0;
    long long frameSourceAbsentGapSamples = 0;
    long long frameNpyFrameSourceAbsentGapSamples = 0;
    long long orcaDftFrameSourceAbsentGapSamples = 0;
    long long anchorUnavailableGapSamples = 0;
    long long invalidSamples = 0;
    QVector<SeriesSparseness> seriesSparseness;
};

class DashboardDisplayController final : public QObject {
    Q_OBJECT

public:
    struct StripTrack {
        const model::ChannelBuffer* buffer = nullptr;
        QColor color;
        bool hasBinding = false;
        model::SignalBinding binding;
    };

    explicit DashboardDisplayController(QObject* parent = nullptr);
    ~DashboardDisplayController() override = default;

    void setContext(const model::QtProtein* protein, model::Conformation* conformation);
    void setSignalModels(model::TrajectorySignalCatalog* catalog,
                         model::DashboardSignalModel* activeModel);
    void setPanelModel(model::DashboardPanelModel* panelModel);
    void setSelection(model::AtomSelection* selection);
    void setDftStore(model::DftShieldingStore* store);
    // L-3a (2026-05-29): the scene overlay receives tensor-glyph
    // reveals when the user activates an h5:reorient_orientation_tensor
    // signal with static.tensor mode. Optional — when null, the
    // controller silently skips the tensor dispatch in rebuild()
    // (the descriptor still works for static.table inspection).
    void setSceneOverlay(SceneRevealOverlay* overlay);

    QVector<StripTrack> stripTracks() const;

    // Move-out the panels built during rebuild() for static-display
    // signals (SequenceBarPanel, future ChordCouplingPanel /
    // PowerSpectrumPanel / LagDecayPanel). The dock forwards them to
    // StripStackWidget::setOwnedPanels(). Returns an empty vector after
    // the move; the controller rebuilds on the next rebuild().
    std::vector<std::unique_ptr<AbstractStripPanel>> takeOwnedPanels();

    DashboardSmokeSummary smokeSummary() const;
    DashboardSmokeSummary smokeSummary(int firstFrame, int lastFrame) const;
    QString statusText() const { return statusText_; }

public slots:
    void setFrame(int frame);
    void setScrubActive(bool active);
    void rebuild();

signals:
    void stripTracksChanged();
    // Emitted ONLY from rebuild() (NOT from setFrame() ticks). Owned
    // panels are static-display artifacts that don't change on time
    // advance; emitting from every tick would drain ownedPanels_ via
    // the dock's setOwnedPanels(controller->takeOwnedPanels()) wiring.
    void ownedPanelsChanged();

private:
    struct ActiveSeries {
        model::DashboardSignal signal;
        model::SignalDescriptor descriptor;
        model::ChannelDescriptor channel;
        QString displayModeId;
        model::SignalBuffer buffer;
        std::function<model::FrameSignalSample(std::size_t frame)> sample;
        QColor color;
        bool hasBinding = false;
        model::SignalBinding binding;
        bool needsFrameSnapshot = false;
        bool needsDftFrame = false;
    };

    void buildGenericTracks(const model::DashboardSignal& signal,
                            const model::SignalDescriptor& descriptor,
                            QVector<ActiveSeries>& series) const;

    // Per-TR panel builders. One per static-display TR landing in
    // Phases C-G; each returns a fully-constructed AbstractStripPanel
    // (or nullptr if the underlying H5 buffer is absent / malformed).
    std::unique_ptr<AbstractStripPanel>
        buildIRedSequenceBarPanel(const model::DashboardSignal& signal,
                                  const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildKernelDynamicsPowerSpectrumPanel(const model::DashboardSignal& signal,
                                              const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildKernelDynamicsLagDecayPanel(const model::DashboardSignal& signal,
                                         const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildReorientSequenceBarPanel(const model::DashboardSignal& signal,
                                      const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildReorientLagDecayPanel(const model::DashboardSignal& signal,
                                   const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildDihedralSequenceBarPanel(const model::DashboardSignal& signal,
                                       const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildDihedralLagDecayPanel(const model::DashboardSignal& signal,
                                    const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildKernelCoherenceChordPanel(const model::DashboardSignal& signal,
                                        const model::SignalDescriptor& descriptor) const;
    std::unique_ptr<AbstractStripPanel>
        buildReorientFixedFreqPanel(const model::DashboardSignal& signal,
                                     const model::SignalDescriptor& descriptor) const;
    // L-4 (2026-05-29): auto-compose. When 2+ Reorient scalar signals
    // (s2/tau_e/r1/r2/noe) are active with static.bar.sequence mode
    // in the same panel, the per-signal loop is skipped for them and
    // this builder folds the group into ONE SequenceBarPanel with
    // the first signal as the primary series + the rest as overlays.
    // Twin-y kicks in automatically when overlay units differ from
    // the primary.
    std::unique_ptr<AbstractStripPanel>
        buildReorientCompositeBarPanel(const QVector<model::DashboardSignal>& group) const;
    QVector<model::ChannelDescriptor> channelsForMode(const model::SignalDescriptor& descriptor,
                                                      const QString& displayModeId) const;
    model::SignalAnchor resolvedAnchorForSignal(const model::DashboardSignal& signal,
                                                const model::SignalDescriptor& descriptor) const;
    QString channelLabel(const model::DashboardSignal& signal,
                         const model::SignalDescriptor& descriptor,
                         const model::ChannelDescriptor& channel,
                         const QString& displayModeId) const;
    QString unitsLabel(const model::SignalDescriptor& descriptor,
                       const model::ChannelDescriptor& channel) const;
    bool seriesIsVisibleInActivePanel(const ActiveSeries& series) const;
    int activePanelSeriesCount() const;
    void refreshPanelVisibility();
    void updateStatusText();
    void extendToFrame(int frame);
    QColor colorForIndex(int index) const;

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    QPointer<model::TrajectorySignalCatalog> catalog_;
    QPointer<model::DashboardSignalModel> activeModel_;
    QPointer<model::DashboardPanelModel> panelModel_;
    QPointer<model::AtomSelection> selection_;
    QPointer<model::DftShieldingStore> dftStore_;
    QPointer<SceneRevealOverlay> sceneOverlay_;

    QVector<ActiveSeries> series_;

    // Owned panels for static-display signals (built during rebuild(),
    // moved out via takeOwnedPanels()).
    std::vector<std::unique_ptr<AbstractStripPanel>> ownedPanels_;

    QString statusText_ = QStringLiteral("No active strip signals.");
    int frame_ = 0;
    bool scrubActive_ = false;
    bool scrubReleasePending_ = false;
    int activeStripSignalCount_ = 0;
    // Number of static-display panels (SequenceBarPanel,
    // PowerSpectrumPanel, etc.) built in the most recent rebuild().
    // Tracked separately from the strip count so updateStatusText() can
    // honestly report iRED / Reorient / etc. as "M panel signal(s)"
    // when no temporal strips are active.
    int activeOwnedPanelCount_ = 0;
};

}  // namespace h5reader::app
