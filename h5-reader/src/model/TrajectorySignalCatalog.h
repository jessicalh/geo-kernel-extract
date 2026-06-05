#pragma once

#include "DashboardSignal.h"
#include "TrajectoryFieldAvailability.h"

#include <QObject>
#include <QVector>

#include <memory>
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
    const QVector<SignalDescriptor>& allDescriptorList() const { return allDescriptors_; }
    void setFieldAvailability(std::shared_ptr<const TrajectoryFieldAvailability> availability);
    std::shared_ptr<const TrajectoryFieldAvailability> fieldAvailability() const { return availability_; }
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
    void rebuildVisibleDescriptors();

    QVector<SignalDescriptor> allDescriptors_;
    QVector<SignalDescriptor> descriptors_;
    std::shared_ptr<const TrajectoryFieldAvailability> availability_;
};

}  // namespace h5reader::model
