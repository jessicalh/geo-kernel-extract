// Conformation — base implementation: topology delegation + current snapshot.

#include "Conformation.h"

#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLoggingCategory>
#include <QPointer>
#include <QThread>
#include <algorithm>
#include <exception>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cSnapshots, "h5reader.snapshots")
}

Conformation::Conformation(const QtProtein* protein, Conformation* snapshotSource)
    : protein_(protein), snapshotSource_(snapshotSource) {
    CENSUS_REGISTER(this);
    if (snapshotSource_) {
        connect(snapshotSource_, &Conformation::snapshotReady, this, &Conformation::snapshotReady);
        connect(snapshotSource_, &Conformation::becameIdle, this, &Conformation::becameIdle);
    }
}

Conformation::~Conformation() {
    if (snapshotWorker_) {
        snapshotWorker_->disconnect(this);
        snapshotWorker_->wait();
        delete snapshotWorker_;
    }
}

std::size_t Conformation::ringCount() const {
    return protein_ ? protein_->ringCount() : 0;
}

std::shared_ptr<const QtConformationSnapshot> Conformation::snapshot(std::size_t frame) const {
    ASSERT_THREAD(this);
    if (snapshotSource_)
        return snapshotSource_->snapshot(frame);
    if (!residentSnapshotFrame_ || *residentSnapshotFrame_ != frame)
        return nullptr;
    return residentSnapshot_;
}

void Conformation::requestSnapshot(std::size_t frame) {
    ASSERT_THREAD(this);
    if (snapshotSource_) {
        snapshotSource_->requestSnapshot(frame);
        return;
    }

    if (!residentSnapshotFrame_ || *residentSnapshotFrame_ != frame) {
        auto reader = snapshotReader(frame);
        residentSnapshot_ = reader ? reader() : nullptr;
        residentSnapshotFrame_ = frame;
    }

    emit snapshotReady(frame);
}

Conformation::SnapshotReader Conformation::snapshotReader(std::size_t) const {
    return {};
}

bool Conformation::isBusy() const {
    ASSERT_THREAD(this);
    return snapshotSource_ ? snapshotSource_->isBusy()
                           : snapshotWorker_ || publishingSnapshot_ || !pendingSnapshotFrames_.empty();
}

void Conformation::requestSnapshotAsync(std::size_t frame) {
    ASSERT_THREAD(this);
    if (snapshotSource_) {
        snapshotSource_->requestSnapshotAsync(frame);
        return;
    }
    if (cancellingSnapshots_)
        return;
    if (residentSnapshotFrame_ == frame) {
        emit snapshotReady(frame);
        return;
    }
    if (activeSnapshotFrame_ == frame ||
        std::find(pendingSnapshotFrames_.begin(), pendingSnapshotFrames_.end(), frame) !=
            pendingSnapshotFrames_.end())
        return;
    pendingSnapshotFrames_.push_back(frame);
    startNextSnapshot();
}

void Conformation::startNextSnapshot() {
    if (snapshotWorker_ || publishingSnapshot_ || pendingSnapshotFrames_.empty() || cancellingSnapshots_)
        return;
    const auto frame = pendingSnapshotFrames_.front();
    pendingSnapshotFrames_.pop_front();
    activeSnapshotFrame_ = frame;
    struct ReadResult {
        std::shared_ptr<const QtConformationSnapshot> snapshot;
        QString error;
    };
    auto result = std::make_shared<ReadResult>();
    snapshotWorker_ = QThread::create([reader = snapshotReader(frame), result] {
        try {
            if (reader)
                result->snapshot = reader();
        } catch (const std::exception& exception) {
            result->error = QString::fromUtf8(exception.what());
        } catch (...) {
            result->error = QStringLiteral("Unknown snapshot read failure");
        }
    });
    snapshotWorker_->setParent(this);
    connect(snapshotWorker_, &QThread::finished, this, [this, frame, result] {
        snapshotWorker_->deleteLater();
        snapshotWorker_ = nullptr;
        activeSnapshotFrame_.reset();
        publishingSnapshot_ = true;
        if (!result->error.isEmpty())
            qCWarning(cSnapshots) << "Frame" << frame << result->error;
        if (!cancellingSnapshots_) {
            residentSnapshotFrame_ = frame;
            residentSnapshot_ = std::move(result->snapshot);
            const QPointer<Conformation> self(this);
            emit snapshotReady(frame);
            if (!self)
                return;
        }
        publishingSnapshot_ = false;
        cancellingSnapshots_ = false;
        startNextSnapshot();
        if (!isBusy())
            emit becameIdle();
    });
    snapshotWorker_->start();
}

void Conformation::cancelPending() {
    ASSERT_THREAD(this);
    if (snapshotSource_) {
        snapshotSource_->cancelPending();
        return;
    }
    pendingSnapshotFrames_.clear();
    cancellingSnapshots_ = snapshotWorker_ != nullptr || publishingSnapshot_;
}

}  // namespace h5reader::model
