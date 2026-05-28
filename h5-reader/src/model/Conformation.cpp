// Conformation — base implementation: topology delegation + current snapshot.

#include "Conformation.h"

#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

namespace h5reader::model {

Conformation::Conformation(const QtProtein* protein)
    : protein_(protein) {
    CENSUS_REGISTER(this);
}

Conformation::~Conformation() = default;

std::size_t Conformation::ringCount() const {
    return protein_ ? protein_->ringCount() : 0;
}

std::shared_ptr<const QtConformationSnapshot> Conformation::snapshot(std::size_t frame) const {
    std::lock_guard<std::mutex> lk(snapshotMutex_);
    if (!residentSnapshotFrame_ || *residentSnapshotFrame_ != frame)
        return nullptr;
    return residentSnapshot_;
}

void Conformation::requestSnapshot(std::size_t frame) {
    ASSERT_THREAD(this);  // v1: the synchronous load runs on the GUI thread

    bool resident = false;
    {
        std::lock_guard<std::mutex> lk(snapshotMutex_);
        if (residentSnapshotFrame_ && *residentSnapshotFrame_ == frame && residentSnapshot_) {
            resident = true;
        }
    }

    if (!resident) {
        // loadSnapshot() does file IO — call it OUTSIDE the lock. v1 is
        // synchronous; the prefetch increment will run this on a worker and
        // hand back a shared_ptr<const> (the pile has no thread affinity).
        std::shared_ptr<const QtConformationSnapshot> snap = loadSnapshot(frame);
        if (!snap)
            return;  // failure already reported at the loader seam; no signal

        std::lock_guard<std::mutex> lk(snapshotMutex_);
        residentSnapshotFrame_ = frame;
        residentSnapshot_ = std::move(snap);
    }

    // Emit OUTSIDE the lock — a connected slot may call back into snapshot().
    emit snapshotReady(frame);
}

}  // namespace h5reader::model
