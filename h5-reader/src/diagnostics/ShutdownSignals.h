// Bridge terminal signals onto the GUI thread and close through Qt's normal
// window-close path. POSIX uses a self-pipe; Windows handles Ctrl-C/Break.

#pragma once

namespace h5reader::diagnostics {

void InstallShutdownSignalHandlers();

}  // namespace h5reader::diagnostics
