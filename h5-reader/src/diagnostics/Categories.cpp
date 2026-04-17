// Categories.cpp — out-of-line definitions for the Q_LOGGING_CATEGORY
// symbols declared in ThreadGuard.h and ConnectionAuditor.h.
//
// Qt's logging-category machinery splits declaration and definition:
// Q_DECLARE_LOGGING_CATEGORY in a header gives every translation unit
// the extern function, while Q_LOGGING_CATEGORY in exactly one .cpp
// provides the storage. Forgetting the .cpp half is a linker error
// that points at qloggingcategory.h rather than our code, so we keep
// the definitions together here for visibility.

#include "ConnectionAuditor.h"
#include "ThreadGuard.h"

namespace h5reader::diagnostics {

Q_LOGGING_CATEGORY(cThreadGuard, "h5reader.threadguard")
Q_LOGGING_CATEGORY(cConnections, "h5reader.connections")

}  // namespace h5reader::diagnostics
