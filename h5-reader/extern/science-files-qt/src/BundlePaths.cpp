#include "BundlePaths.h"
#include <QStringList>

namespace sciencefiles {
namespace {

bool safeComponent(const QString &name) {
    if (name.isEmpty() || name == "." || name == ".." || name.endsWith('.') || name.endsWith(' '))
        return false;
    for (QChar c : name) {
        if (c.unicode() < 32 || QStringLiteral("/\\:<>\"|?*").contains(c))
            return false;
    }
    const QString stem = name.section('.', 0, 0).toUpper();
    if (stem == "CON" || stem == "PRN" || stem == "AUX" || stem == "NUL")
        return false;
    if (stem.size() == 4 && (stem.startsWith("COM") || stem.startsWith("LPT")) &&
        QStringLiteral("123456789\u00b9\u00b2\u00b3").contains(stem.back()))
        return false;
    return true;
}
} // namespace

bool safeBundleName(const QString &name) {
    return !name.startsWith('.') && safeComponent(name);
}

bool safeRelativePath(const QString &path) {
    const auto parts = path.split('/');
    for (const auto &part : parts) {
        if (!safeComponent(part))
            return false;
    }
    return true;
}

} // namespace sciencefiles
