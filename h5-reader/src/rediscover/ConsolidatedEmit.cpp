#include "ConsolidatedEmit.h"

#include "ScopedProducerCatalog.h"

#include <QDir>
#include <QFileInfo>

namespace h5reader::rediscover {

bool RunConsolidatedEmit(const ConsolidatedEmitOptions& options,
                         ConsolidatedEmitStats* stats,
                         QString* err_out) {
    if (options.root720.isEmpty()) {
        if (err_out) *err_out = QStringLiteral("consolidated_emit requires --root720");
        return false;
    }
    if (options.run1p9j.isEmpty()) {
        if (err_out) *err_out = QStringLiteral("consolidated_emit requires --run for the 1P9J LGS/calcset");
        return false;
    }
    if (options.outDir.isEmpty()) {
        if (err_out) *err_out = QStringLiteral("consolidated_emit requires --out");
        return false;
    }
    if (!QFileInfo(options.root720).isDir()) {
        if (err_out) *err_out = QStringLiteral("720 root does not exist: %1").arg(options.root720);
        return false;
    }
    if (!QFileInfo::exists(options.run1p9j)) {
        if (err_out) *err_out = QStringLiteral("1P9J run path does not exist: %1").arg(options.run1p9j);
        return false;
    }
    QDir().mkpath(options.outDir);

    if (stats) {
        stats->rows = 0;
        stats->scopedFieldCount = ScopedProducerCatalog().size();
    }
    if (err_out) {
        *err_out = QStringLiteral("consolidated_emit implementation is not wired yet");
    }
    return false;
}

}  // namespace h5reader::rediscover
