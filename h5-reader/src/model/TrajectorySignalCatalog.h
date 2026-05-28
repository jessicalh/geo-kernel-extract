#pragma once

#include "DashboardSignal.h"

#include <QObject>
#include <QVector>

#include <optional>

namespace h5reader::model {

struct SignalDescriptorFilter {
    std::optional<SignalSourceKind> sourceKind;
    std::optional<SignalAxis> axis;
    std::optional<SignalValueShape> valueShape;
    QString displayModeId;
    QString text;
    bool includePending = true;
    bool includeTemporal = true;
    bool includeStatic = true;
};

class TrajectorySignalCatalog final : public QObject {
    Q_OBJECT

public:
    explicit TrajectorySignalCatalog(QObject* parent = nullptr);

    QVector<SignalDescriptor> descriptors() const { return descriptors_; }
    const QVector<SignalDescriptor>& descriptorList() const { return descriptors_; }
    const SignalDescriptor* findDescriptor(const QString& descriptorId) const;
    std::optional<SignalDescriptor> descriptorById(const QString& descriptorId) const;

    QVector<SignalDescriptor> filterDescriptors(const SignalDescriptorFilter& filter) const;
    QVector<SignalDescriptor> descriptorsForSource(SignalSourceKind sourceKind) const;
    QVector<SignalDescriptor> descriptorsForAxis(SignalAxis axis) const;
    QVector<SignalDescriptor> descriptorsForValueShape(SignalValueShape valueShape) const;
    QVector<SignalDescriptor> descriptorsForDisplayMode(const QString& displayModeId) const;
    QVector<SignalDescriptor> search(const QString& text) const;

    bool supportsDisplayMode(const QString& descriptorId, const QString& displayModeId) const;
    bool canBind(const DisplaySignalBinding& binding) const;
    bool canSample(const DisplaySignalBinding& binding) const;

private:
    static QVector<SignalDescriptor> BuildDescriptorCatalog();

    QVector<SignalDescriptor> descriptors_;
};

}  // namespace h5reader::model
