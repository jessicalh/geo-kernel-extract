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
    ASSERT_THREAD(this);
    if (!residentSnapshotFrame_ || *residentSnapshotFrame_ != frame)
        return nullptr;
    return residentSnapshot_;
}

void Conformation::requestSnapshot(std::size_t frame) {
    ASSERT_THREAD(this);

    if (!residentSnapshotFrame_ || *residentSnapshotFrame_ != frame || !residentSnapshot_) {
        std::shared_ptr<const QtConformationSnapshot> snap = loadSnapshot(frame);
        if (!snap)
            return;

        residentSnapshotFrame_ = frame;
        residentSnapshot_ = std::move(snap);
    }

    emit snapshotReady(frame);
}

}  // namespace h5reader::model
