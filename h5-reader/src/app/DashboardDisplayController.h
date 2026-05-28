#pragma once

#include "../model/DashboardSignal.h"
#include "../model/SignalDictionary.h"
#include "../model/SignalTimeSeries.h"

#include <QObject>
#include <QColor>
#include <QPointer>
#include <QString>
#include <QVector>

#include <functional>
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

struct DashboardSmokeSummary {
    int seriesCount = 0;
    int seriesWithSamples = 0;
    int seriesWithValidSamples = 0;
    int seriesPendingOnly = 0;
    int seriesWithMismatchedBuffers = 0;
    long long samples = 0;
    long long channelValues = 0;
    long long channelValidity = 0;
    long long validSamples = 0;
    long long gapSamples = 0;
    long long pendingGapSamples = 0;
    long long sourceAbsentGapSamples = 0;
    long long frameSourceAbsentGapSamples = 0;
    long long anchorUnavailableGapSamples = 0;
    long long invalidSamples = 0;
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

    QVector<StripTrack> stripTracks() const;
    DashboardSmokeSummary smokeSummary() const;
    QString statusText() const { return statusText_; }

public slots:
    void setFrame(int frame);
    void rebuild();

signals:
    void stripTracksChanged();

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
    QVector<model::ChannelDescriptor> channelsForMode(const model::SignalDescriptor& descriptor,
                                                      const QString& displayModeId) const;
    model::SignalBinding revealBindingForSignal(const model::DashboardSignal& signal,
                                                const model::SignalDescriptor& descriptor) const;
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
    void extendToFrame(int frame);
    QColor colorForIndex(int index) const;

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    QPointer<model::TrajectorySignalCatalog> catalog_;
    QPointer<model::DashboardSignalModel> activeModel_;
    QPointer<model::DashboardPanelModel> panelModel_;
    QPointer<model::AtomSelection> selection_;
    QPointer<model::DftShieldingStore> dftStore_;

    QVector<ActiveSeries> series_;
    QString statusText_ = QStringLiteral("No active strip signals.");
    int frame_ = 0;
};

}  // namespace h5reader::app
