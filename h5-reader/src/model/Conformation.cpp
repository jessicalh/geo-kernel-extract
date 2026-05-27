// Conformation — base implementation: topology delegation + the snapshot LRU.

#include "Conformation.h"

#include "QtConformationSnapshot.h"
#include "QtProtein.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <algorithm>

namespace h5reader::model {

Conformation::Conformation(const QtProtein* protein, std::size_t cacheCapacity)
    : protein_(protein), cacheCapacity_(cacheCapacity == 0 ? 1 : cacheCapacity) {
    CENSUS_REGISTER(this);
}

Conformation::~Conformation() = default;

std::size_t Conformation::ringCount() const {
    return protein_ ? protein_->ringCount() : 0;
}

std::shared_ptr<const QtConformationSnapshot> Conformation::snapshot(std::size_t frame) const {
    std::lock_guard<std::mutex> lk(cacheMutex_);
    const auto it = cache_.find(frame);
    return it != cache_.end() ? it->second : nullptr;
}

void Conformation::requestSnapshot(std::size_t frame) {
    ASSERT_THREAD(this);  // v1: the synchronous load runs on the GUI thread

    bool resident = false;
    {
        std::lock_guard<std::mutex> lk(cacheMutex_);
        if (cache_.find(frame) != cache_.end()) {
            touch_(frame);
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

        std::lock_guard<std::mutex> lk(cacheMutex_);
        if (cache_.emplace(frame, std::move(snap)).second) {
            lru_.push_back(frame);
            while (cache_.size() > cacheCapacity_) {
                cache_.erase(lru_.front());
                lru_.pop_front();
            }
        }
    }

    // Emit OUTSIDE the lock — a connected slot may call back into snapshot().
    emit snapshotReady(frame);
}

void Conformation::touch_(std::size_t frame) {
    const auto it = std::find(lru_.begin(), lru_.end(), frame);
    if (it != lru_.end()) {
        lru_.erase(it);
        lru_.push_back(frame);
    }
}

}  // namespace h5reader::model
