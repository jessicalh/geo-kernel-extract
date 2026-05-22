// Viewer entry point.
//
// Command line parsed via JobSpec (shared with nmr_extract).
// Qt-specific options (--rest-port) handled separately.
//
// Startup sequence:
//   1. QApplication (must be first)
//   2. Library logging (OperationLog — the ONE logging system)
//   3. Parse command line (JobSpec + Qt extras)
//   4. MainWindow
//   5. show()
//   6. Deferred auto-load via QTimer::singleShot (event loop must be running)
//   7. app.exec()

#include <QApplication>
#include <QSurfaceFormat>
#include <QMetaType>
#include <QTimer>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QVTKOpenGLNativeWidget.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "CalculatorConfig.h"
#include "ComputeWorker.h"
#include "JobSpec.h"
#include "MainWindow.h"
#include "OperationLog.h"
#include "RestServer.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "errors.h"

// udp_log shim: legacy viewer code calls this. Forwards to OperationLog.
#include <cstdarg>
extern "C" void udp_log(const char* fmt, ...) {
    char buf[2048];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    nmr::OperationLog::Info(nmr::LogViewer, "Viewer", buf);
}


// ============================================================================
// Extract --rest-port from argv before JobSpec parsing.
// JobSpec ignores unknown flags, but we need to pull this out for Qt.
// ============================================================================

static quint16 ExtractRestPort(int argc, char* argv[]) {
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--rest-port") == 0)
            return static_cast<quint16>(std::atoi(argv[i + 1]));
    }
    return 9147;  // default
}


// ============================================================================
// LoggingConfig: parse [logging] from ~/.nmr_tools.toml.
//
// Lives in main_viewer.cpp (not in OperationLog) because MainWindow's
// QUdpSocket bind needs the host string to call joinMulticastGroup() on
// a multicast destination, and OperationLog has no public getter for
// its configured host (private static host_). Parsing once in main and
// handing (host, port) to MainWindow via constructor keeps the TOML as
// the single source of truth for both sender (OperationLog) and
// receiver (MainWindow's log dock).
//
// Mirrors the OperationLog::LoadChannelConfig parser shape; intentional
// duplication per PATTERNS §17 (duplicate over chain) — the schema is
// trivial and the two sites have orthogonal responsibilities.
// ============================================================================

struct LoggingConfig {
    std::string udp_host = "127.0.0.1";  // unicast fallback if TOML absent
    int udp_port = 9998;
    std::string file_path;
    uint32_t channel_mask = nmr::LogAll;
};

static LoggingConfig ReadLoggingConfig() {
    LoggingConfig cfg;
    const char* home = std::getenv("HOME");
    if (!home)
        return cfg;
    std::string path = std::string(home) + "/.nmr_tools.toml";
    std::ifstream in(path);
    if (!in.is_open())
        return cfg;

    auto trim = [](std::string& s) {
        while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '"' || s.back() == '\r'))
            s.pop_back();
        while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '"'))
            s = s.substr(1);
    };

    std::string line;
    bool in_logging = false;
    while (std::getline(in, line)) {
        if (line.find("[logging]") != std::string::npos) {
            in_logging = true;
            continue;
        }
        if (line.find('[') != std::string::npos && line.find("[logging]") == std::string::npos)
            in_logging = false;
        if (!in_logging)
            continue;

        auto eq = line.find('=');
        if (eq == std::string::npos)
            continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        trim(key);
        trim(val);

        if (key == "udp_host")
            cfg.udp_host = val;
        else if (key == "udp_port") {
            try {
                cfg.udp_port = std::stoi(val);
            } catch (...) {
            }
        } else if (key == "file")
            cfg.file_path = val;
        else if (key == "channels") {
            try {
                cfg.channel_mask = static_cast<uint32_t>(std::stoul(val, nullptr, 0));
            } catch (...) {
            }
        }
    }
    return cfg;
}


int main(int argc, char* argv[]) {
    // 1. QApplication — must be first
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc, argv);
    app.setApplicationName("Protein Tensor Viewer");
    app.setApplicationVersion("2.0");

    // 2. Library environment + logging via Session.
    //
    // nmr::Session is the named entity that holds process-wide
    // resources: RuntimeEnvironment, OperationLog channel config (UDP +
    // file from [logging] in ~/.nmr_tools.toml), the loaded AIMNet2
    // model, the libpq TripeptideDftTable connection, and the
    // LarsenHBondGrid. nmr_extract uses the same Session object; the
    // viewer adopts it here so the two binaries share the same
    // resource-management surface.
    //
    // LoadFromToml does RuntimeEnvironment::Load + OperationLog
    // ::LoadChannelConfig (which configures UDP / file / channel mask
    // from TOML — multicast destinations are auto-detected in the
    // sender, see OperationLog.cpp:36-63). MainWindow still needs the
    // UDP host/port directly so its Log-dock socket can
    // joinMulticastGroup, so we parse the same TOML keys here.
    nmr::Session session;
    if (session.LoadFromToml() != nmr::kOk) {
        fprintf(stderr, "ERROR: Session::LoadFromToml: %s\n", session.LastError().c_str());
        return 1;
    }
    // Session::LoadFromToml already emits LogSessionStart internally
    // (see src/Session.cpp:35); do not call it a second time here.
    const LoggingConfig log_cfg = ReadLoggingConfig();

    // 3. Parse command line via JobSpec
    //    --output is optional for the viewer (no output = just visualise).
    //    --rest-port is viewer-specific, extracted separately.
    auto spec = nmr::ParseJobSpec(argc, argv);
    quint16 restPort = ExtractRestPort(argc, argv);

    // Handle --help: print both JobSpec usage and viewer extras
    if (spec.mode == nmr::JobMode::None && spec.error.empty()) {
        nmr::PrintJobSpecUsage(argv[0]);
        fprintf(stderr,
            "\nViewer options:\n"
            "  --rest-port PORT   REST API port (default 9147, 0 to disable)\n"
            "\n"
            "  The viewer accepts all the same modes as nmr_extract.\n"
            "  --output is optional: omit it to just visualise.\n");
        return 0;
    }

    // Parse error
    if (!spec.Ok()) {
        fprintf(stderr, "ERROR: %s\n\n", spec.error.c_str());
        nmr::PrintJobSpecUsage(argv[0]);
        return 1;
    }

    // Load CalculatorConfig BEFORE ValidateJobSpec so the validation can
    // resolve the AIMNet2 TOML default (aimnet2_model_path) when --aimnet2
    // is not on the CLI. Mirrors the nmr_extract main() pattern; the
    // AIMNet2 gate in ValidateJobSpec is documented as requiring this
    // ordering (see JobSpec.cpp:308-337 and
    // feedback_aimnet2_required_no_weasel).
    if (!spec.config_path.empty()) {
        nmr::CalculatorConfig::Load(spec.config_path);
    } else {
        std::string default_config = std::string(NMR_DATA_DIR) + "/calculator_params.toml";
        if (std::filesystem::exists(default_config))
            nmr::CalculatorConfig::Load(default_config);
    }

    // Validate file existence
    if (!nmr::ValidateJobSpec(spec)) {
        QMessageBox::critical(nullptr, "Validation Error",
            QString::fromStdString(spec.error));
        return 1;
    }

    // Print warnings
    for (const auto& w : spec.warnings)
        fprintf(stderr, "WARNING: %s\n", w.c_str());

    // 4. Load Session-owned process resources. AIMNet2 is mandatory
    //    (gate already passed in ValidateJobSpec). TripeptideDftTable
    //    and LarsenHBondGrid are optional: empty DSN / empty grid dir
    //    means the related calculators are skipped silently. A configured-
    //    but-failed load is fatal — we'd rather see the error than
    //    silently drop expected results.
    if (session.LoadAimnet2Model(spec.aimnet2_model_path) != nmr::kOk) {
        fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }
    if (session.LoadTripeptideDftTable() != nmr::kOk) {
        fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }
    if (session.LoadLarsenHBondGrid() != nmr::kOk) {
        fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }

    // 5. Register metatypes for cross-thread signals
    qRegisterMetaType<nmr::JobSpec>("nmr::JobSpec");
    qRegisterMetaType<ComputeResult>("ComputeResult");

    // 6. Create and show window. Session passes through to ComputeWorker
    //    so OperationRunner::Run sees the AIMNet2 model, Tripeptide DFT
    //    table, and Larsen H-bond grid. Host/port come from the TOML so
    //    the Log dock binds the same UDP destination OperationLog is
    //    sending to (multicast join in MainWindow::setupUI).
    MainWindow window(session, QString::fromStdString(log_cfg.udp_host), static_cast<quint16>(log_cfg.udp_port));
    QObject::connect(&app, &QCoreApplication::aboutToQuit, &window, &MainWindow::shutdown);
    window.show();

    // 5b. REST API server
    RestServer* restServer = nullptr;
    if (restPort > 0) {
        restServer = new RestServer(&window, restPort, &app);
        if (restServer->actualPort() > 0) {
            window.setWindowTitle(QString("Protein Tensor Viewer — REST :%1")
                .arg(restServer->actualPort()));
        }
    }

    // 6. Deferred auto-load from JobSpec — all modes supported
    if (spec.mode != nmr::JobMode::None) {
        QTimer::singleShot(0, &window, [&window, spec]() {
            nmr::OperationLog::Info(nmr::LogViewer, "main",
                "auto-loading from JobSpec");
            window.loadFromJobSpec(spec);
        });
    }

    nmr::OperationLog::Info(nmr::LogViewer, "main", "entering event loop");
    return app.exec();
}
