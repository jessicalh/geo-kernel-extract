#include "ReaderMainWindow.h"

#include "CameraAnchorHelper.h"
#include "CameraComposer.h"
#include "CameraInputFilter.h"
#include "CameraMode.h"
#include "OrientationPolicy.h"
#include "MoleculeScene.h"
#include "TensorGlyphActor.h"
#include "TensorGlyphMath.h"
#include "NearbySignalModel.h"
#include "MeasurementsDock.h"
#include "QtAtomInspectorDock.h"
#include "QtAtomPicker.h"
#include "RestServer.h"
#include "QtBackboneRibbonOverlay.h"
#include "QtAtomTrajectoryOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtPlaybackController.h"
#include "TimeViewportController.h"
#include "MeasurementOverlay.h"
#include "QtRingPolygonOverlay.h"
#include "DashboardStripDock.h"
#include "DashboardDisplayController.h"
#include "DashboardSelectionController.h"
#include "SignalDisplayDialog.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/DashboardLogging.h"
#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/StructuredLogger.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtProteinLoader.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/DftShieldingStore.h"
#include "../model/ExperimentalShieldingMlStore.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"
#include "../model/QtBondVectorBuffers.h"
#include "../io/QtTrajectoryH5.h"
#include "../model/TrajectoryFieldAvailability.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"
#include "../model/VisualizationRegistry.h"
#include "../model/CsaShape.h"
#include "../model/MolecularFrame.h"
#include "../model/MolecularFrameSelect.h"
#include "../model/ConformationGeometry.h"
#include "../physics/SphericalBasis.h"
#include "CsaTensorOverlay.h"

#include <QDockWidget>

#include <QDir>
#include <QFile>
#include <QFileInfo>

#include <QAction>
#include <QActionGroup>
#include <QtGlobal>
#include <QApplication>
#include <QCloseEvent>
#include <QComboBox>
#include <QCoreApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFont>
#include <QFormLayout>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QKeySequence>
#include <QLabel>
#include <QLoggingCategory>
#include <QMenu>
#include <QMenuBar>
#include <QMouseEvent>
#include <QMessageBox>
#include <QRegion>
#include <QSettings>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStackedLayout>
#include <QStringList>
#include <QStatusBar>
#include <QColor>
#include <QIcon>
#include <QPainter>
#include <QPalette>
#include <QPixmap>
#include <QPolygonF>
#include <QPushButton>
#include <QStyle>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QUuid>
#include <QVariant>
#include <QWidget>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkRendererCollection.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cWindow, "h5reader.window")

// QSettings — versioned state blob policy. Bump on dock-object
// additions or any layout-invalidating change so old blobs are
// silently discarded by QMainWindow::restoreState. Schema-evolution
// safe per ROBUSTNESS_BACKLOG_2026-05-30.md item 7.
constexpr int kSettingsVersion = 2;   // bumped: property docks now start hidden
constexpr int kMaxRecentFiles  = 10;

bool fileExistsInDir(const QDir& dir, const QString& fileName) {
    return QFileInfo(dir.filePath(fileName)).isFile();
}

QJsonValue jsonNull() {
    return QJsonValue(QJsonValue::Null);
}

QString diagnosticSeverityName(h5reader::diagnostics::Severity severity) {
    using h5reader::diagnostics::Severity;
    switch (severity) {
        case Severity::Info: return QStringLiteral("Info");
        case Severity::Warning: return QStringLiteral("Warning");
        case Severity::Error: return QStringLiteral("Error");
        case Severity::Fatal: return QStringLiteral("Fatal");
    }
    return QStringLiteral("Unknown");
}

QString diagnosticSeverityStyle(h5reader::diagnostics::Severity severity) {
    using h5reader::diagnostics::Severity;
    switch (severity) {
        case Severity::Info:
            return QStringLiteral("color: #475569;");
        case Severity::Warning:
            return QStringLiteral("color: #7a4b00; font-weight: 600;");
        case Severity::Error:
        case Severity::Fatal:
            return QStringLiteral("color: #9f1239; font-weight: 700;");
    }
    return QString();
}

QJsonArray stringListJson(const QStringList& values) {
    QJsonArray out;
    for (const QString& value : values)
        out.append(value);
    return out;
}

QDir installedExperimentalShieldingMlDir() {
    return QDir(QDir(QCoreApplication::applicationDirPath())
                    .filePath(QStringLiteral("ml/experimental_shielding_ml")));
}

QStringList experimentalShieldingMlRequiredFiles() {
    return {
        QStringLiteral("model.ts"),
        QStringLiteral("manifest.json"),
        QStringLiteral("infer.exe"),
        QStringLiteral("c10.dll"),
        QStringLiteral("torch.dll"),
        QStringLiteral("torch_cpu.dll"),
        QStringLiteral("torch_global_deps.dll"),
        QStringLiteral("libiomp5md.dll"),
        QStringLiteral("libiompstubs5md.dll"),
        QStringLiteral("uv.dll"),
    };
}

QStringList experimentalShieldingMlRocmRequiredFiles() {
    return {
        QStringLiteral("infer.exe"),
        QStringLiteral("c10.dll"),
        QStringLiteral("c10_hip.dll"),
        QStringLiteral("torch_cpu.dll"),
        QStringLiteral("torch_hip.dll"),
        QStringLiteral("caffe2_nvrtc.dll"),
        QStringLiteral("aotriton_v2.dll"),
        QStringLiteral("amdhip64_7.dll"),
        QStringLiteral("amd_comgr0702.dll"),
        QStringLiteral("hiprtc0702.dll"),
        QStringLiteral("MIOpen.dll"),
        QStringLiteral("rocblas.dll"),
        QStringLiteral("libhipblaslt.dll"),
    };
}

QStringList experimentalShieldingMlRocmMissing(const QDir& runtimeDir) {
    QStringList missing;
    for (const QString& fileName : experimentalShieldingMlRocmRequiredFiles()) {
        if (!fileExistsInDir(runtimeDir, fileName))
            missing.append(fileName);
    }
    for (const QString& directory :
         {QStringLiteral("rocblas/library"), QStringLiteral("hipblaslt/library")}) {
        if (!QDir(runtimeDir.filePath(directory)).exists())
            missing.append(directory);
    }
    return missing;
}

bool experimentalShieldingMlRocmRuntimeAvailable(const QString& helperPath) {
    return QFileInfo(helperPath).isFile()
           && experimentalShieldingMlRocmMissing(QFileInfo(helperPath).dir()).isEmpty();
}

QString experimentalShieldingMlDevicePreference() {
    const QString value =
        qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_DEVICE")
            .trimmed()
            .toLower();
    if (value == QStringLiteral("cpu") || value == QStringLiteral("rocm"))
        return value;
    return QStringLiteral("auto");
}

QString developmentExperimentalShieldingMlRocmHelper(const QString& modelPath) {
    const QString explicitPath =
        qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_ROCM_HELPER");
    if (!explicitPath.isEmpty())
        return explicitPath;
    return QFileInfo(modelPath).dir().filePath(QStringLiteral("rocm/infer.exe"));
}

bool installedExperimentalShieldingMlRuntimeAvailable() {
    const QDir mlDir = installedExperimentalShieldingMlDir();
    const QStringList requiredFiles = experimentalShieldingMlRequiredFiles();
    for (const QString& fileName : requiredFiles) {
        if (!fileExistsInDir(mlDir, fileName))
            return false;
    }
    return model::ExperimentalShieldingMlStore::ManifestHasInferenceSchema(
        mlDir.filePath(QStringLiteral("manifest.json")));
}

QStringList experimentalShieldingMlDevMissingFiles(const QString& modelPath,
                                                   const QString& manifestPath,
                                                   const QString& helperPath) {
    QStringList missing;
    if (!QFileInfo(modelPath).isFile())
        missing.append(QStringLiteral("model.ts"));
    if (!QFileInfo(manifestPath).isFile())
        missing.append(QStringLiteral("manifest.json"));
    else if (!model::ExperimentalShieldingMlStore::ManifestHasInferenceSchema(manifestPath))
        missing.append(QStringLiteral("manifest.inference_schema"));
    if (!QFileInfo(helperPath).isFile())
        missing.append(QStringLiteral("infer.exe"));

    const QDir helperDir = QFileInfo(helperPath).dir();
    for (const QString& fileName :
         {QStringLiteral("c10.dll"),
          QStringLiteral("torch.dll"),
          QStringLiteral("torch_cpu.dll"),
          QStringLiteral("torch_global_deps.dll"),
          QStringLiteral("libiomp5md.dll"),
          QStringLiteral("libiompstubs5md.dll"),
          QStringLiteral("uv.dll")}) {
        if (!fileExistsInDir(helperDir, fileName))
            missing.append(fileName);
    }
    return missing;
}

QJsonObject summarizeExperimentalShieldingMlTraining(const QJsonObject& training) {
    QJsonObject out;
    const QJsonValue fold = training.value(QStringLiteral("fold"));
    if (fold.isString())
        out.insert(QStringLiteral("fold"), fold.toString());
    const QJsonValue claim = training.value(QStringLiteral("training_claim"));
    if (claim.isString())
        out.insert(QStringLiteral("claim"), claim.toString());
    const QJsonValue bestEpoch = training.value(QStringLiteral("best_epoch"));
    if (bestEpoch.isDouble())
        out.insert(QStringLiteral("bestEpoch"), bestEpoch.toInt());
    const QJsonValue bestVal = training.value(QStringLiteral("best_val"));
    if (bestVal.isDouble())
        out.insert(QStringLiteral("bestVal"), bestVal.toDouble());
    const QJsonValue vocabPolicy = training.value(QStringLiteral("label_vocab_policy"));
    if (vocabPolicy.isString())
        out.insert(QStringLiteral("labelVocabPolicy"), vocabPolicy.toString());
    const QJsonValue run = training.value(QStringLiteral("run"));
    if (run.isString())
        out.insert(QStringLiteral("run"), run.toString());
    const QJsonValue producerCommit = training.value(QStringLiteral("producer_commit"));
    if (producerCommit.isString())
        out.insert(QStringLiteral("producerCommit"), producerCommit.toString());
    const QJsonObject args = training.value(QStringLiteral("args")).toObject();
    const QJsonValue inputPreset = args.value(QStringLiteral("input_preset"));
    if (inputPreset.isString())
        out.insert(QStringLiteral("inputPreset"), inputPreset.toString());
    return out;
}

QJsonObject summarizeExperimentalShieldingMlModel(const QJsonObject& model) {
    QJsonObject out;
    const auto copyString = [&model, &out](const QString& key, const QString& outKey) {
        const QJsonValue value = model.value(key);
        if (value.isString())
            out.insert(outKey, value.toString());
    };
    copyString(QStringLiteral("id"), QStringLiteral("id"));
    copyString(QStringLiteral("label"), QStringLiteral("label"));
    copyString(QStringLiteral("role"), QStringLiteral("role"));
    copyString(QStringLiteral("model_file"), QStringLiteral("modelFile"));
    copyString(QStringLiteral("input_preset"), QStringLiteral("inputPreset"));
    copyString(QStringLiteral("selection"), QStringLiteral("selection"));
    const QJsonObject training = model.value(QStringLiteral("training")).toObject();
    if (!training.isEmpty())
        out.insert(QStringLiteral("training"), summarizeExperimentalShieldingMlTraining(training));
    return out;
}

QJsonObject readExperimentalShieldingMlManifestSummary(const QString& path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        return {};
    QJsonParseError error{};
    const QJsonDocument doc = QJsonDocument::fromJson(file.readAll(), &error);
    if (error.error != QJsonParseError::NoError || !doc.isObject())
        return {};
    const QJsonObject manifest = doc.object();
    QJsonObject out;
    const auto copyString = [&manifest, &out](const char* key, const QString& outKey) {
        const QJsonValue value = manifest.value(QLatin1String(key));
        if (value.isString())
            out.insert(outKey, value.toString());
    };
    copyString("name", QStringLiteral("name"));
    copyString("bundle_version", QStringLiteral("bundleVersion"));
    copyString("bundle_date", QStringLiteral("bundleDate"));
    copyString("source_checkpoint", QStringLiteral("sourceCheckpoint"));
    copyString("target", QStringLiteral("target"));
    const QJsonObject training = manifest.value(QStringLiteral("training")).toObject();
    if (!training.isEmpty()) {
        out.insert(QStringLiteral("training"), summarizeExperimentalShieldingMlTraining(training));
    }
    const QJsonArray models = manifest.value(QStringLiteral("models")).toArray();
    if (!models.isEmpty()) {
        QJsonArray outModels;
        for (const QJsonValue& value : models) {
            if (value.isObject())
                outModels.append(summarizeExperimentalShieldingMlModel(value.toObject()));
        }
        out.insert(QStringLiteral("models"), outModels);
    }
    const QJsonObject inferenceSchema =
        manifest.value(QStringLiteral("inference_schema")).toObject();
    if (!inferenceSchema.isEmpty())
        out.insert(QStringLiteral("inferenceSchemaVersion"),
                   inferenceSchema.value(QStringLiteral("version")).toInt());
    return out;
}

QJsonObject experimentalShieldingMlInputProfileJson(
    const model::TrajectoryFieldAvailability*,
    bool loaded) {
    QJsonObject profile;
    profile.insert(QStringLiteral("loaded"), loaded);
    profile.insert(QStringLiteral("contract"),
                   QStringLiteral("july_full720_f003_no_mopac_common_sense_v1"));
    return profile;
}

QJsonObject selectedExperimentalShieldingMlModelJson(const QJsonObject& inputProfile,
                                                     bool runtimeAvailable) {
    QJsonObject selected;
    if (!runtimeAvailable) {
        selected.insert(QStringLiteral("id"), jsonNull());
        selected.insert(QStringLiteral("reason"), QStringLiteral("runtime_missing"));
        return selected;
    }
    if (!inputProfile.value(QStringLiteral("loaded")).toBool()) {
        selected.insert(QStringLiteral("id"), jsonNull());
        selected.insert(QStringLiteral("reason"), QStringLiteral("no_loaded_run"));
        return selected;
    }
    selected.insert(QStringLiteral("id"), QStringLiteral("f003_r004"));
    selected.insert(QStringLiteral("modelFile"), QStringLiteral("model.ts"));
    selected.insert(QStringLiteral("inputPreset"),
                    QStringLiteral("f003_no_mopac_common_sense"));
    selected.insert(QStringLiteral("reason"), QStringLiteral("july_contract_loaded"));
    return selected;
}

struct ExperimentalShieldingMlRuntimePaths {
    QString model;
    QString manifest;
    QString helper;
    QString device;
    QString fallbackHelper;
};

std::optional<ExperimentalShieldingMlRuntimePaths>
resolveExperimentalShieldingMlRuntime() {
    QString model = qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL");
    QString manifest = qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST");
    QString cpuHelper = qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER");
    QString rocmHelper = developmentExperimentalShieldingMlRocmHelper(model);
    if (!experimentalShieldingMlDevMissingFiles(model, manifest, cpuHelper).isEmpty()) {
        const QDir installed = installedExperimentalShieldingMlDir();
        if (!installedExperimentalShieldingMlRuntimeAvailable())
            return std::nullopt;
        model = installed.filePath(QStringLiteral("model.ts"));
        manifest = installed.filePath(QStringLiteral("manifest.json"));
        cpuHelper = installed.filePath(QStringLiteral("infer.exe"));
        rocmHelper = installed.filePath(QStringLiteral("rocm/infer.exe"));
    }

    ExperimentalShieldingMlRuntimePaths paths;
    paths.manifest = manifest;
    const QString preference = experimentalShieldingMlDevicePreference();
    const bool rocmAvailable = experimentalShieldingMlRocmRuntimeAvailable(rocmHelper);
    if (preference == QStringLiteral("rocm") && !rocmAvailable)
        return std::nullopt;
    if (preference != QStringLiteral("cpu") && rocmAvailable) {
        paths.helper = rocmHelper;
        paths.device = QStringLiteral("rocm");
        if (preference == QStringLiteral("auto"))
            paths.fallbackHelper = cpuHelper;
    } else {
        paths.helper = cpuHelper;
        paths.device = QStringLiteral("cpu");
    }
    paths.model = model;
    return paths;
}

QJsonObject experimentalShieldingMlRuntimeJson(
    const model::TrajectoryFieldAvailability* availability,
    bool loaded) {
    QJsonObject out;
    out.insert(QStringLiteral("available"), false);
    out.insert(QStringLiteral("runtime"), QStringLiteral("missing"));
    const QJsonObject inputProfile =
        experimentalShieldingMlInputProfileJson(availability, loaded);
    out.insert(QStringLiteral("inputProfile"), inputProfile);

    const QString modelPath =
        qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL");
    const QString manifestPath =
        qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST");
    const QString helperPath =
        qEnvironmentVariable("H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER");
    const bool devRuntimeRequested =
        !modelPath.isEmpty() || !manifestPath.isEmpty() || !helperPath.isEmpty();
    QStringList devMissing;
    if (devRuntimeRequested) {
        devMissing = experimentalShieldingMlDevMissingFiles(modelPath, manifestPath, helperPath);
        if (devMissing.isEmpty()) {
            out.insert(QStringLiteral("available"), true);
            out.insert(QStringLiteral("runtime"), QStringLiteral("development"));
            if (QFileInfo(manifestPath).isFile()) {
                out.insert(QStringLiteral("manifest"),
                           readExperimentalShieldingMlManifestSummary(manifestPath));
            }
            out.insert(QStringLiteral("selectedModel"),
                       selectedExperimentalShieldingMlModelJson(inputProfile, true));
            return out;
        }
    }

    const QDir mlDir = installedExperimentalShieldingMlDir();
    QStringList missing;
    for (const QString& fileName : experimentalShieldingMlRequiredFiles()) {
        if (!fileExistsInDir(mlDir, fileName))
            missing.append(fileName);
    }
    const QString installedManifest = mlDir.filePath(QStringLiteral("manifest.json"));
    if (QFileInfo(installedManifest).isFile()
        && !model::ExperimentalShieldingMlStore::ManifestHasInferenceSchema(installedManifest)) {
        missing.append(QStringLiteral("manifest.inference_schema"));
    }

    if (!missing.isEmpty() && devRuntimeRequested) {
        out.insert(QStringLiteral("runtime"), QStringLiteral("development"));
        out.insert(QStringLiteral("missing"), stringListJson(devMissing));
        out.insert(QStringLiteral("installedMissing"), stringListJson(missing));
        out.insert(QStringLiteral("selectedModel"),
                   selectedExperimentalShieldingMlModelJson(inputProfile, false));
        return out;
    }

    out.insert(QStringLiteral("available"), missing.isEmpty());
    out.insert(QStringLiteral("runtime"), QStringLiteral("installed"));
    if (!missing.isEmpty())
        out.insert(QStringLiteral("missing"), stringListJson(missing));
    if (devRuntimeRequested)
        out.insert(QStringLiteral("developmentMissing"), stringListJson(devMissing));
    if (QFileInfo(installedManifest).isFile())
        out.insert(QStringLiteral("manifest"), readExperimentalShieldingMlManifestSummary(installedManifest));
    out.insert(QStringLiteral("selectedModel"),
               selectedExperimentalShieldingMlModelJson(inputProfile, missing.isEmpty()));
    return out;
}

QString fitModeToolTip() {
    return QStringLiteral(
        "Stabilisation mode — click to switch.\n"
        "Locked backbone: Kabsch fit of the backbone (industry standard) — removes global tumbling; the backbone holds still while sidechains/loops move.\n"
        "Kabsch with give: all-atom fit — removes tumbling but lets real internal motion show.");
}

// Note: locateDftJobsDir was deleted as part of the 2026-05-31 SIMPLIFY
// pass; the DFT campaign now comes from the `.LGS` `dft.frames[]` array
// (see CalcsetManifest + DftShieldingStore).

}  // namespace

ReaderMainWindow::ReaderMainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("ReaderMainWindow"));

    qCInfo(cWindow).noquote() << "ctor entered";

    buildUi();
    buildToolbar();
    buildStatusBar();
    ACONNECT(h5reader::diagnostics::ErrorBus::Instance(),
             &h5reader::diagnostics::ErrorBus::errorReported,
             this,
             &ReaderMainWindow::handleErrorBusReport);
    buildDocks();

    // Default size — wide enough for the playback + camera + transform +
    // metrics + overlays + panel controls to fit in one
    // toolbar without Qt's overflow chevron. QSettings restore overrides
    // this on later launches.
    resize(1600, 900);
    setWindowTitle(QStringLiteral("h5-reader"));

    // QSettings restore — geometry, dock state, log mask, recent menu.
    // Tolerant: missing or version-mismatched blobs leave the ctor's
    // explicit defaults intact. Runs AFTER all docks/toolbars exist so
    // restoreState has named docks to bind to.
    restoreAllSettings();
    setEmptyState();

    qCInfo(cWindow).noquote() << "ctor done";
}

bool ReaderMainWindow::loadRunPath(const QString& path) {
    ASSERT_THREAD(this);
    lastLoadError_.clear();
    lastDiagnosticSeverity_.clear();
    lastDiagnosticSource_.clear();
    lastDiagnosticMessage_.clear();
    lastDiagnosticValues_.clear();
    if (diagnosticLabel_) {
        diagnosticLabel_->clear();
        diagnosticLabel_->setToolTip(QString());
        diagnosticLabel_->setVisible(false);
    }
    if (path.isEmpty()) {
        lastLoadError_ = QStringLiteral("No calcset path was provided.");
        return false;
    }

    qCInfo(cWindow).noquote() << "loading" << path;
    auto loaded = h5reader::io::QtProteinLoader::LoadRunPath(path);
    if (!loaded.ok) {
        lastLoadError_ = loaded.error.isEmpty()
            ? QStringLiteral("Load failed for %1").arg(path)
            : loaded.error;
        qCCritical(cWindow).noquote() << "Load failed:" << lastLoadError_;
        return false;
    }
    if (loaded.decodeWarnings > 0) {
        qCWarning(cWindow).noquote()
            << "Decode completed with" << loaded.decodeWarnings << "warnings";
    }

    installLoadedRun(std::move(loaded));
    lastLoadError_.clear();
    return true;
}

void ReaderMainWindow::installLoadedRun(h5reader::io::QtLoadResult&& loaded) {
    ASSERT_THREAD(this);
    clearLoadedRun();
    loaded_ = std::make_unique<h5reader::io::QtLoadResult>(std::move(loaded));

    // Upstream data-transform layer (feedback_viewer_two_layers_transform_and_camera).
    // Wraps the loader's Conformation so consumers (scene, picker, overlays,
    // REST /positions) read positions through a runtime-switchable rigid-body
    // display transform. Startup mode is backbone fit with the iterative mean
    // seeded/anchored at frame 0 so the reader opens stationary.
    transformed_ = new h5reader::model::TransformedConformation(loaded_->conformation.get(), this);
    const auto backboneSubset =
        h5reader::model::TransformedConformation::BackboneSubset(*loaded_->protein);
    using TMode = h5reader::model::TransformedConformation::Mode;
    if (backboneSubset.size() >= 3) {
        transformed_->setMode(TMode::FitSubset, 0, backboneSubset);
    } else {
        qCWarning(cWindow).noquote()
            << "backbone fit unavailable at startup; falling back to all-atom fit";
        transformed_->setMode(TMode::FitReference, 0);
    }
    ACONNECT(transformed_, &h5reader::model::TransformedConformation::transformChanged,
             this, [this]() {
                 updateFitModeLabel();
                 if (scene_) scene_->refreshCurrentFrame();
             });
    updateFitModeLabel();

    // Scene binds to the VTK widget's render window. The scene reads
    // positions through the wrapped conformation so transform mode
    // changes are visible immediately.
    scene_ = new MoleculeScene(vtkWidget_, renderWindow_, this);
    scene_->Build(*loaded_->protein, *transformed_);
    applyOverlayActionState();
    scene_->refreshCurrentFrame();
    scene_->ResetCamera();
    ACONNECT(scene_, &MoleculeScene::cameraPlaneLockChanged,
             this, [this](bool) { updateCameraModeActions(); });
    if (scene_->cameraComposer()) {
        ACONNECT(scene_->cameraComposer(), &CameraComposer::modeChanged,
                 this, [this]() { updateCameraModeActions(); });
    }

    // Playback controller — frameChanged drives the scene, which drives
    // the render. Toolbar controls drive the playback.
    const int T = static_cast<int>(loaded_->conformation->frameCount());
    playback_ = new QtPlaybackController(T, this);
    timeViewport_ = new TimeViewportController(T, this);

    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             scene_,    &MoleculeScene::setFrame);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             this,      &ReaderMainWindow::onFrameChanged);
    ACONNECT(playback_,     &QtPlaybackController::frameChanged,
             timeViewport_, &TimeViewportController::setCurrentFrame);
    ACONNECT(timeViewport_, &TimeViewportController::playbackFrameRequested,
             playback_,     &QtPlaybackController::setFrame);
    ACONNECT(playback_, &QtPlaybackController::playingChanged,
             this,      [this](bool) { refreshControlStates(); });

    if (frameSlider_) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setRange(0, std::max(0, T - 1));
        frameSlider_->setValue(0);
    }
    if (fpsSpinner_) {
        const QSignalBlocker block(fpsSpinner_);
        fpsSpinner_->setRange(1, 60);
        fpsSpinner_->setValue(playback_->fps());
    }

    // Atom picker — event filter on the VTK widget. Emits atomPicked(idx,
    // modifiers) on double-click. It stays dumb; AtomSelection interprets
    // the gesture and fans typed changes out.
    auto* firstRenderer = renderWindow_->GetRenderers()->GetFirstRenderer();
    picker_ = new QtAtomPicker(vtkWidget_, firstRenderer,
                                loaded_->protein.get(),
                                transformed_,
                                playback_, this);

    // Camera input filter — installed AFTER the picker so Qt's filter
    // chain runs THIS first. Double-click events fall through to the picker.
    cameraInputFilter_ = new CameraInputFilter(vtkWidget_, scene_,
                                                 scene_->cameraComposer(), this);

    // Click on empty space (no atom hit, no drag) stops/restarts the animation.
    ACONNECT(cameraInputFilter_, &CameraInputFilter::viewportClicked,
             this, [this](QPointF pos) {
                 if (!picker_ || !playback_) return;
                 const auto hit = picker_->atomAt(
                     static_cast<int>(pos.x()), static_cast<int>(pos.y()));
                 if (!hit) playback_->togglePlayPause();
             });

    inspectorDock_->setContext(loaded_->protein.get(), transformed_);
    ACONNECT(playback_,  &QtPlaybackController::frameChanged,
             inspectorDock_, &QtAtomInspectorDock::setFrame);

    if (measurementsDock_) {
        measurementsDock_->setContext(loaded_->protein.get(), loaded_->conformation.get());
        ACONNECT(playback_, &QtPlaybackController::frameChanged,
                 measurementsDock_, &MeasurementsDock::setFrame);
    }

    // ---- Selection model — the single source of selection truth ----------
    selection_ = new model::AtomSelection(loaded_->protein.get(), this);

    // Nearby-residue source for the Filter checklist. Prefer the transformed
    // (displayed) conformation; fall back to the base conformation if no
    // backbone fit is active so the checklist still populates. Lazily created
    // once and re-pointed on each load (parented to this; never per-load leak).
    if (!filterNearby_)
        filterNearby_ = new NearbySignalModel(this);
    filterNearby_->setContext(
        loaded_->protein.get(),
        transformed_ ? static_cast<model::Conformation*>(transformed_)
                     : loaded_->conformation.get());
    filterNearby_->setRadiusAngstrom(5.0);
    filterResidues_.clear();

    signalCatalog_ = new model::TrajectorySignalCatalog(this);
    if (const auto runtime = resolveExperimentalShieldingMlRuntime()) {
        experimentalMlStore_ = new model::ExperimentalShieldingMlStore(
            loaded_->protein.get(),
            loaded_->conformation.get(),
            runtime->model,
            runtime->manifest,
            loaded_->extractionManifestPath,
            runtime->helper,
            runtime->device,
            runtime->fallbackHelper,
            this);
    }

    // The availability gate needs the startup-loaded topology spine sizes and
    // the DFT job count to classify Topology + live ORCA descriptors honestly
    // (neither travels through the per-frame NPY path the gate probes). The DFT
    // count comes from the manifest (== DftShieldingStore::jobCount()). The ML
    // source is live only when both its runtime and this run's required inputs
    // passed the store's contract check.
    const model::TrajectoryFieldAvailability::TopologyExtent topologyExtent{
        static_cast<qsizetype>(loaded_->protein->atomCount()),
        static_cast<qsizetype>(loaded_->protein->bondCount()),
        static_cast<qsizetype>(loaded_->protein->residueCount()),
        static_cast<qsizetype>(loaded_->protein->ringCount()),
        static_cast<qsizetype>(loaded_->protein->ringMembershipCount())};
    const std::size_t dftJobCount =
        loaded_->manifest.dft.has_value() ? loaded_->manifest.dft->frames.size() : 0;
    fieldAvailability_ = std::make_shared<model::TrajectoryFieldAvailability>(
        model::TrajectoryFieldAvailability::Build(loaded_->conformation.get(),
                                                  topologyExtent, dftJobCount,
                                                  experimentalMlStore_
                                                      && experimentalMlStore_->isReady(),
                                                  signalCatalog_->allDescriptorList()));
    signalCatalog_->setFieldAvailability(fieldAvailability_);
    visualizationContext_ = {};
    visualizationContext_.availability = fieldAvailability_.get();
    visualizationContext_.hasTrajectory = loaded_->conformation
        && loaded_->conformation->asTrajectory() != nullptr;
    visualizationContext_.tensorGlyphGestureEnabled =
        scene_ && scene_->csaOverlay()
        && experimentalMlStore_ && experimentalMlStore_->isReady();
    const QStringList unresolvedModes =
        model::VisualizationRegistry::instance().unresolvedStaticModes(*signalCatalog_);
    for (const QString& mode : unresolvedModes) {
        qCWarning(diagnostics::cDash).noquote()
            << QStringLiteral("event=viz_unresolved_static_mode mode=%1").arg(mode);
    }
    Q_ASSERT(unresolvedModes.isEmpty());
    inspectorDock_->setFieldAvailability(fieldAvailability_);
    dashboardSignals_ = new model::DashboardSignalModel(this);
    dashboardSignals_->setFieldAvailability(fieldAvailability_);
    dashboardPanels_ = new model::DashboardPanelModel(this);
    dashboardSelectionController_ =
        new DashboardSelectionController(signalCatalog_, dashboardSignals_, dashboardPanels_, this);
    lastDashboardSelectedCount_ = dashboardSelectionController_->selectedCount();

    signalDisplayDialog_ = new SignalDisplayDialog(this);
    signalDisplayDialog_->setTrajectorySignalCatalog(signalCatalog_);
    signalDisplayDialog_->setDashboardSignalModel(dashboardSignals_);
    signalDisplayDialog_->setDashboardPanelModel(dashboardPanels_);
    signalDisplayDialog_->setDashboardSelectionController(dashboardSelectionController_.data());
    signalDisplayDialog_->setContext(loaded_->protein.get(), transformed_);
    signalDisplayDialog_->setVisualizationContext(visualizationContext_);
    signalDisplayDialog_->setSelection(selection_);
    ACONNECT(playback_, &QtPlaybackController::frameChanged,
             signalDisplayDialog_, &SignalDisplayDialog::setFrame);

    ACONNECT(picker_,    &QtAtomPicker::atomPicked,
             selection_, &model::AtomSelection::applyPick);
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             scene_,  &MoleculeScene::clearReveal);
    // Tag the render scheduler so the EndEvent observer logs source=picker
    // for the render that follows. selection_->applyPick triggers
    // refreshCurrentFrame which itself calls requestRender(Timer);
    // tagging Picker afterward overrides the source.
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             this,   [this](std::size_t, Qt::KeyboardModifiers) {
                 if (scene_) scene_->requestRender(
                     MoleculeScene::RenderSource::Picker);
             });

    // Reveal-on-pick: picking an atom brings up its Inspector (it starts hidden;
    // this is the dock's reveal path now that the Panels menu is gone).
    ACONNECT(picker_, &QtAtomPicker::atomPicked,
             this,   [this](std::size_t, Qt::KeyboardModifiers) {
                 if (inspectorDock_ && !inspectorDock_->isVisible())
                     revealDockQueued(inspectorDock_);
             });

    ACONNECT(selection_, &model::AtomSelection::focusChanged,
             inspectorDock_, &QtAtomInspectorDock::setPickedAtom);
    ACONNECT(selection_, &model::AtomSelection::cleared,
             inspectorDock_, &QtAtomInspectorDock::clearSelection);

    if (measurementsDock_) {
        // The whole ORDERED tuple drives a measurement (not just focus), so this
        // tracks AtomSelection::changed; it reveals only once a 2+ atom geometry
        // exists (a single pick is the Inspector's job, not a measurement).
        ACONNECT(selection_, &model::AtomSelection::changed, this, [this]() {
            if (!measurementsDock_)
                return;
            measurementsDock_->setAtoms(selection_->atoms());
            if (selection_->atoms().size() >= 2 && !measurementsDock_->isVisible())
                revealDockQueued(measurementsDock_);
        });
        ACONNECT(selection_, &model::AtomSelection::cleared,
                 measurementsDock_, &MeasurementsDock::clear);
    }

    // Atom Info is the DEFAULT front tab on a single-atom focus: the Measurements
    // / strip docks reveal alongside as tabs but must not steal it. Covers REST
    // picks too (the picker-signal reveal above is GUI-only). Deferred via
    // revealDockQueued so this raise wins; gated on a single atom so a 2-4 atom
    // geometry instead raises the Measurements tab (handled above).
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
             [this](std::size_t) {
                 if (inspectorDock_ && selection_ && selection_->atoms().size() < 2)
                     revealDockQueued(inspectorDock_);
             });
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
             [this](std::size_t) { refreshControlStates(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this,
             [this]() { refreshControlStates(); });

    // CSA tensor glyph (mode-2): focus + frame driven; honest gap on a missing
    // DFT frame; raw->display alignment via the molecular frame.
    ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
             [this](std::size_t) { updateCsaGlyph(true); updateOrientationTensorGlyph(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, [this]() {
        updateCsaGlyph(true);
        if (scene_ && scene_->orientationGlyph())
            scene_->orientationGlyph()->clear();
        if (scene_)
            scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    });
    ACONNECT(playback_, &QtPlaybackController::frameChanged, this,
             [this](int) { updateCsaGlyph(false); updateOrientationTensorGlyph(); });
    ACONNECT(playback_, &QtPlaybackController::playingChanged, this,
             [this](bool playing) {
        if (!playing)
            updateCsaGlyph(true);
    });

    if (auto* meas = scene_->measurementOverlay()) {
        meas->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::changed,
                 meas,       &MeasurementOverlay::onSelectionChanged);
    }

    // Trajectory overlay: focus-driven (NOT changed() -- focus-only, so a
    // plain pick doesn't rebuild twice) + transform-driven (a fit-mode change
    // moves every aligned position, so the path is stale). The overlay
    // self-guards single-pose / rigid atoms. A render is requested after
    // focus/clear rebuilds (coalesced); transformChanged is already followed
    // by refreshCurrentFrame above.
    if (auto* traj = scene_->atomTrajectoryOverlay()) {
        traj->setSelection(selection_);
        ACONNECT(selection_, &model::AtomSelection::focusChanged,
                 traj,       &QtAtomTrajectoryOverlay::onFocusChanged);
        ACONNECT(selection_, &model::AtomSelection::cleared,
                 traj,       &QtAtomTrajectoryOverlay::onSelectionCleared);
        ACONNECT(transformed_, &model::TransformedConformation::transformChanged,
                 traj,         &QtAtomTrajectoryOverlay::onTransformChanged);
        ACONNECT(traj, &QtAtomTrajectoryOverlay::rebuildStarted,
                 this, [this](int frames) {
                     QApplication::setOverrideCursor(Qt::WaitCursor);
                     statusBar()->showMessage(
                         QStringLiteral("Loading trajectory envelope (%1 frames)...").arg(frames));
                 });
        ACONNECT(traj, &QtAtomTrajectoryOverlay::rebuildFinished,
                 this, [this](int frames, int /*dftSamples*/, int loadMs) {
                     QApplication::restoreOverrideCursor();
                     statusBar()->showMessage(
                         QStringLiteral("Trajectory envelope: %1 frames (%2 ms)")
                             .arg(frames).arg(loadMs));
                 });
        ACONNECT(selection_, &model::AtomSelection::focusChanged, this,
                 [this](std::size_t) {
                     if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
                 });
        ACONNECT(selection_, &model::AtomSelection::cleared, this,
                 [this]() {
                     if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
                 });
    }
    ACONNECT(selection_, &model::AtomSelection::changed,
             this, [this]() {
                 if (scene_ && scene_->cameraComposer()
                     && scene_->cameraComposer()->mode().kind
                            == CameraMode::Kind::Plane) {
                     const std::size_t t = playback_
                         ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
                     scene_->cameraComposer()->setMode(FreeMode(), FreePolicy(), t);
                 }
                 refreshControlStates();
                 if (scene_) scene_->refreshCurrentFrame();
             });

    // Selection summary in the status bar (count + measurement kind).
    ACONNECT(selection_, &model::AtomSelection::changed, this, [this]() { updateSelectionStatus(); });
    ACONNECT(selection_, &model::AtomSelection::cleared, this, [this]() { updateSelectionStatus(); });
    updateSelectionStatus();

    dashboardStripDock_->setContext(loaded_->protein.get(), transformed_);
    dashboardStripDock_->setSignalModels(signalCatalog_, dashboardSignals_);
    dashboardStripDock_->setPanelModel(dashboardPanels_);
    dashboardStripDock_->setSelectionController(dashboardSelectionController_.data());
    dashboardStripDock_->setSelection(selection_);
    dashboardStripDock_->setTimeViewport(timeViewport_);
    dashboardController_ = dashboardStripDock_->displayController();
    if (dashboardController_) {
        dashboardController_->setVisualizationContext(visualizationContext_);
        ACONNECT(dashboardController_.data(),
                 &DashboardDisplayController::sceneTensorBindingChanged,
                 this,
                 [this](const QString& descriptorId, qint64 atom) {
                     activeExperimentalMlTensorDescriptor_.clear();
                     activeExperimentalMlTensorAtom_.reset();
                     if (descriptorId
                             == QStringLiteral("ml:experimental_shielding_t2")
                         && atom >= 0
                         && loaded_ && loaded_->protein
                         && static_cast<std::size_t>(atom)
                                < loaded_->protein->atomCount()) {
                         activeExperimentalMlTensorDescriptor_ = descriptorId;
                         activeExperimentalMlTensorAtom_ =
                             static_cast<std::size_t>(atom);
                     }
                     updateCsaGlyph(true);
                 });
    }

    ACONNECT(dashboardSelectionController_.data(),
             &DashboardSelectionController::selectedCountChanged,
             this,
             [this](int count) {
                 ASSERT_THREAD(this);
                 const bool added = count > lastDashboardSelectedCount_;
                 lastDashboardSelectedCount_ = count;
                 if (!added || !dashboardStripDock_ || dashboardStripDock_->isVisible())
                     return;
                 qCInfo(diagnostics::cDash).noquote()
                     << QStringLiteral("event=dock_reveal_on_add count=%1").arg(count);
                 revealDockQueued(dashboardStripDock_);
             });
    ACONNECT(playback_,           &QtPlaybackController::frameChanged,
             dashboardStripDock_, &DashboardStripDock::setFrame);

    // Expose the shared scene overlay to dashboard visualizations that explicitly
    // support 3-D geometry, including the common tensor glyph.
    if (scene_->revealOverlay()) {
        dashboardStripDock_->setSceneOverlay(scene_->revealOverlay());
        visualizationContext_.hasSceneOverlay = true;
        if (dashboardController_)
            dashboardController_->setVisualizationContext(visualizationContext_);
        if (signalDisplayDialog_)
            signalDisplayDialog_->setVisualizationContext(visualizationContext_);
    }

    if (experimentalMlStore_ && experimentalMlStore_->isReady()) {
        dashboardStripDock_->setExperimentalShieldingMlStore(experimentalMlStore_);
        ACONNECT(experimentalMlStore_,
                 &model::ExperimentalShieldingMlStore::frameReady,
                 this,
                 [this](std::size_t frame) {
                     const int current = playback_ ? playback_->currentFrame() : 0;
                     if (frame
                         == static_cast<std::size_t>(std::max(0, current))) {
                         updateCsaGlyph(false);
                     }
                 });
        qCInfo(cWindow).noquote()
            << QStringLiteral("Experimental Shielding ML store wired | model=%1 device=%2")
                   .arg(experimentalMlStore_->modelId())
                   .arg(experimentalMlStore_->device());
    }

    // DFT shielding campaign (optional): make the frame-local source
    // available to descriptor-family samplers. The `.LGS` carries the
    // typed `dft.frames[]` map — frame_index → meta.json — so the store
    // builds straight from it.
    if (loaded_->manifest.dft.has_value()) {
        const auto& dft = *loaded_->manifest.dft;
        dftStore_ = new model::DftShieldingStore(loaded_->protein.get(), dft.frames, this);
        if (scene_ && scene_->atomTrajectoryOverlay())
            scene_->atomTrajectoryOverlay()->setDftStore(dftStore_);
        ACONNECT(dftStore_, &model::DftShieldingStore::frameReady,
                 this, [this](std::size_t originalIndex) {
            if (!loaded_ || !loaded_->conformation || !playback_)
                return;
            const int frameI = playback_->currentFrame();
            const std::size_t frame = static_cast<std::size_t>(frameI < 0 ? 0 : frameI);
            if (loaded_->conformation->originalFrameIndex(frame) == originalIndex)
                updateCsaGlyph(false);
        });
        dashboardStripDock_->setDftStore(dftStore_);
        visualizationContext_.hasDftStore = true;
        if (dashboardController_)
            dashboardController_->setVisualizationContext(visualizationContext_);
        if (signalDisplayDialog_)
            signalDisplayDialog_->setVisualizationContext(visualizationContext_);
        qCInfo(cWindow).noquote() << "DFT shielding store wired from .LGS |"
                                  << "method=" << dft.method
                                  << "| frames=" << dftStore_->jobCount()
                                  << "| campaign_target=" << dft.campaign_target_frames;
    }

    QString alternatePath;
    if (loaded_->manifest.kind == h5reader::io::CalcsetManifest::Kind::MutantPair
        && loaded_->manifest.mutant_pair.has_value()) {
        alternatePath = loaded_->manifest.mutant_pair->ala_lgs_abspath;
    }
    updateMutantAlternateAction(alternatePath);

    refreshControlStates();
    onFrameChanged(0);
    resetDashboardStateForRunLoad();
    if (centralContainer_) {
        if (auto* stack = qobject_cast<QStackedLayout*>(centralContainer_->layout()))
            stack->setCurrentWidget(vtkWidget_);
    }
    if (emptyPlaceholder_)
        emptyPlaceholder_->setVisible(false);
    setWindowTitle(QStringLiteral("h5-reader — %1").arg(loaded_->proteinId));
    if (proteinLabel_)
        proteinLabel_->setText(loaded_->proteinId);
    if (!loaded_->runPath.isEmpty())
        addToRecentFiles(QDir(loaded_->runPath).absolutePath());
    syncRestServerContext();

    qCInfo(cWindow).noquote()
        << "run loaded | protein=" << loaded_->proteinId
        << "| path=" << loaded_->runPath;
}

void ReaderMainWindow::updateCsaGlyph(bool requestMissingDft) {
    ASSERT_THREAD(this);
    CsaTensorOverlay* overlay = scene_ ? scene_->csaOverlay() : nullptr;
    if (!overlay)
        return;
    auto redraw = [this] {
        if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    };
    auto hide = [&] {
        overlay->clear();
        if (inspectorDock_) inspectorDock_->clearCsaTensor();
        experimentalMlTensorDisplayed_ = false;
        experimentalMlTensorDisplayedFrame_.reset();
        redraw();
    };
    if (trajectoryOverlayActive()) {
        hide();
        return;
    }

    if (!loaded_ || !transformed_) {
        hide();
        return;
    }

    const int frameI = playback_ ? playback_->currentFrame() : 0;
    const std::size_t frame = static_cast<std::size_t>(frameI < 0 ? 0 : frameI);

    // A dashboard-selected F003 tensor owns the shared shielding glyph while
    // active. The network emits its equivariant tensor in the raw coordinate
    // frame; apply the exact display Kabsch rotation used by atomPosition().
    if (activeExperimentalMlTensorDescriptor_
            == QStringLiteral("ml:experimental_shielding_t2")
        && activeExperimentalMlTensorAtom_.has_value()) {
        const std::size_t atom = *activeExperimentalMlTensorAtom_;
        if (!experimentalMlStore_ || !experimentalMlStore_->isReady()
            || !loaded_->protein || atom >= loaded_->protein->atomCount()
            || frame >= transformed_->frameCount()) {
            hide();
            return;
        }
        const auto values = experimentalMlStore_->tensor(frame, atom);
        if (!values) {
            if (requestMissingDft && (!playback_ || !playback_->isPlaying()))
                experimentalMlStore_->requestFrame(frame);
            hide();
            return;
        }

        const std::array<double, 5> t2{
            (*values)[1], (*values)[2], (*values)[3], (*values)[4], (*values)[5]};
        const model::Mat3 rawTensor =
            physics::ReconstructLibraryT2Matrix((*values)[0], t2);
        const model::Mat3 rotation = transformed_->displayRotation(frame);
        const model::Mat3 displayTensor =
            rotation * rawTensor * rotation.transpose();
        const model::CsaShape shape = model::ComputeCsaShape(displayTensor);
        if (!shape.valid) {
            hide();
            return;
        }

        const model::Vec3 atomPos = transformed_->atomPosition(frame, atom);
        overlay->show(atomPos, shape, std::nullopt);
        experimentalMlTensorDisplayed_ = true;
        experimentalMlTensorDisplayedFrame_ = frame;
        if (inspectorDock_) {
            CsaTensorInfo info;
            info.sourceLabel = QStringLiteral("Experimental Shielding ML");
            info.sourceDetail = experimentalMlStore_->modelId();
            info.frameKind =
                QStringLiteral("equivariant output, display-aligned");
            info.sigmaIso = shape.sigma_iso;
            info.span = shape.span;
            info.skew = shape.skew;
            info.eta = shape.eta;
            info.sigma11 = shape.principal_values[0];
            info.sigma22 = shape.principal_values[1];
            info.sigma33 = shape.principal_values[2];
            inspectorDock_->setCsaTensor(atom, info);
        }
        qCDebug(cWindow).noquote()
            << "Experimental shielding ML glyph | atom=" << atom
            << "| frame=" << frame
            << "| model=" << experimentalMlStore_->modelId()
            << "| iso=" << shape.sigma_iso
            << "| eta=" << shape.eta
            << "| span=" << shape.span;
        redraw();
        return;
    }

    experimentalMlTensorDisplayed_ = false;
    experimentalMlTensorDisplayedFrame_.reset();
    if (!dftStore_ || !selection_ || !selection_->hasFocus()) {
        hide();
        return;
    }
    const std::size_t atom = selection_->focus();
    const model::QtProtein* protein = loaded_->protein.get();
    model::Conformation* rawConf = loaded_->conformation.get();
    if (!protein || !rawConf) {
        hide();
        return;
    }
    const std::size_t original = rawConf->originalFrameIndex(frame);
    if (!dftStore_->hasJob(original)) {
        hide();
        return;
    }
    if (!dftStore_->frame(original)) {
        if (requestMissingDft && (!playback_ || !playback_->isPlaying()))
            dftStore_->requestFrameAsync(original);
        hide();
        return;
    }

    const model::AtomCsaResult r =
        model::ComputeAtomCsa(*protein, *rawConf, *transformed_, *dftStore_,
                              atom, frame, false);
    if (!r.valid) {
        hide();
        return;
    }
    qCDebug(cWindow).noquote()
        << "CSA glyph | atom=" << atom << "| framed=" << r.framed
        << "| kind=" << model::MolecularFrameKindName(r.frameKind)
        << "| iso=" << r.shape.sigma_iso << "| eta=" << r.shape.eta
        << "| span=" << r.shape.span;
    overlay->show(r.atomPos, r.shape, r.molecularAxes);
    if (inspectorDock_) {
        CsaTensorInfo info;
        info.framed = r.framed;
        info.sourceLabel = QStringLiteral("ORCA DFT");
        info.sourceDetail =
            loaded_->manifest.dft.has_value()
                ? loaded_->manifest.dft->method
                : QStringLiteral("ORCA");
        info.frameKind = r.framed
                             ? QString::fromLatin1(model::MolecularFrameKindName(r.frameKind))
                             : QStringLiteral("unframed (raw PAS)");
        info.sigmaIso = r.shape.sigma_iso;
        info.span = r.shape.span;
        info.skew = r.shape.skew;
        info.eta = r.shape.eta;
        info.sigma11 = r.shape.principal_values[0];
        info.sigma22 = r.shape.principal_values[1];
        info.sigma33 = r.shape.principal_values[2];
        inspectorDock_->setCsaTensor(atom, info);
    }
    redraw();
}

void ReaderMainWindow::updateOrientationTensorGlyph() {
    ASSERT_THREAD(this);
    TensorGlyphActor* glyph = scene_ ? scene_->orientationGlyph() : nullptr;
    if (!glyph)
        return;
    auto hide = [&] {
        glyph->clear();
        if (inspectorDock_) inspectorDock_->clearOrientationTensor();
        if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    };
    if (trajectoryOverlayActive()) {
        hide();
        return;
    }
    if (!loaded_ || !transformed_ || !selection_ || !selection_->hasFocus()) {
        hide();
        return;
    }
    const model::QtProtein* protein = loaded_->protein.get();
    model::Conformation* conf = loaded_->conformation.get();
    if (!protein || !conf) {
        hide();
        return;
    }
    const auto* traj = conf->asTrajectory();
    const auto* h5 = traj ? traj->h5() : nullptr;
    const auto* rd = h5 ? h5->reorientationalDynamics() : nullptr;
    if (!rd || rd->identity.n_vectors == 0) {
        hide();
        return;
    }

    // The first reorient bond vector whose tail or head IS the focused atom.
    const std::size_t atom = selection_->focus();
    std::optional<std::size_t> row;
    for (std::size_t i = 0; i < rd->identity.n_vectors; ++i) {
        const std::int32_t t = rd->identity.tail_atom[i];
        const std::int32_t h = rd->identity.head_atom[i];
        if ((t >= 0 && static_cast<std::size_t>(t) == atom)
            || (h >= 0 && static_cast<std::size_t>(h) == atom)) {
            row = i;
            break;
        }
    }
    if (!row || (*row + 1) * 9 > rd->orientation_tensor.data.size()) {
        hide();
        return;
    }
    const std::size_t tailAtom = static_cast<std::size_t>(rd->identity.tail_atom[*row]);
    const std::size_t headAtom = static_cast<std::size_t>(rd->identity.head_atom[*row]);
    if (tailAtom >= protein->atomCount() || headAtom >= protein->atomCount()) {
        hide();
        return;
    }

    std::array<double, 9> flat{};
    for (std::size_t k = 0; k < 9; ++k)
        flat[k] = rd->orientation_tensor.data[*row * 9 + k];
    const math::TensorEllipsoid eig = math::decomposeSymmetric3x3(flat);

    const int frameI = playback_ ? playback_->currentFrame() : 0;
    const std::size_t frame = static_cast<std::size_t>(frameI < 0 ? 0 : frameI);
    const model::Vec3 tailPos = transformed_->atomPosition(frame, tailAtom);
    const model::Vec3 headPos = transformed_->atomPosition(frame, headAtom);
    const model::Vec3 mid = (tailPos + headPos) * 0.5;
    const model::Vec3 bondVec = headPos - tailPos;
    if (bondVec.norm() < 1e-9) {
        hide();
        return;
    }
    const model::Vec3 bondDir = bondVec.normalized();

    // Body->lab: the bond_orientation_tensor accumulates in the producer's
    // frame-0-aligned BODY frame. Align its dominant eigenvector to the CURRENT
    // displayed bond direction (the same approximation the old reorient ellipsoid
    // used) so the glyph tracks the molecule as it tumbles.
    const Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(eig.axes[0], bondDir);
    model::Mat3 labAxes;
    std::array<double, 3> pv{};
    for (int i = 0; i < 3; ++i) {
        const std::size_t ai = static_cast<std::size_t>(i);
        const model::Vec3 lab = q * eig.axes[ai];
        labAxes.col(i) = lab;
        const double r = eig.radii[ai];  // decomposeSymmetric3x3 returns sqrt(lambda)
        pv[ai] = r * r;
    }
    const double iso = (pv[0] + pv[1] + pv[2]) / 3.0;

    glyph->show(mid, pv, labAxes, iso);

    // Mirror the numbers into the Atom Info panel (text -> window), exactly as the
    // CSA glyph does, so picture and numbers agree and the scene stays geometry.
    if (inspectorDock_) {
        OrientationTensorInfo info;
        QString bondName;
        switch (rd->identity.kind[*row]) {
        case 1: bondName = QStringLiteral("N-H"); break;
        case 2: bondName = QStringLiteral("CA-HA"); break;
        case 3: bondName = QStringLiteral("C-O"); break;
        default: bondName = QStringLiteral("bond"); break;
        }
        info.bond = QStringLiteral("%1 (residue %2)")
                        .arg(bondName)
                        .arg(rd->identity.residue_index[*row]);
        // Lipari-Szabo order parameter straight from the order-tensor
        // eigenvalues: S^2 = (3*sum(lambda^2) - 1)/2. Universal across NH/CaHa/CO,
        // unlike rd->s2 (NH-only -> NaN on a C=O or CA-HA bond).
        info.s2 = (3.0 * (pv[0] * pv[0] + pv[1] * pv[1] + pv[2] * pv[2]) - 1.0) / 2.0;
        info.lambda1 = pv[0];
        info.lambda2 = pv[1];
        info.lambda3 = pv[2];
        inspectorDock_->setOrientationTensor(atom, info);
    }

    if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
}

// Compute one atom's CSA result for the current frame -- the SAME orchestration
// the glyph uses, exposed so REST /csa vets exactly what is drawn.
model::AtomCsaResult ReaderMainWindow::probeAtomCsa(std::size_t atom, int requestedFrame) {
    if (!loaded_ || !dftStore_ || !transformed_)
        return {};
    const model::QtProtein* protein = loaded_->protein.get();
    model::Conformation* rawConf = loaded_->conformation.get();
    if (!protein || !rawConf)
        return {};
    // requestedFrame < 0 -> live frame; >= 0 -> that frame, read directly from
    // the DFT store (ComputeAtomCsa already takes an explicit frame), so probing
    // a past frame never moves playback.
    const int frameI = requestedFrame >= 0
                           ? requestedFrame
                           : (playback_ ? playback_->currentFrame() : 0);
    const std::size_t frame = static_cast<std::size_t>(frameI < 0 ? 0 : frameI);
    return model::ComputeAtomCsa(*protein, *rawConf, *transformed_, *dftStore_, atom, frame);
}

void ReaderMainWindow::clearLoadedRun() {
    ASSERT_THREAD(this);

    if (playback_)
        playback_->pause();

    if (signalDisplayDialog_) {
        signalDisplayDialog_->hide();
        delete signalDisplayDialog_;
        signalDisplayDialog_ = nullptr;
    }

    if (dashboardStripDock_) {
        dashboardStripDock_->setSceneOverlay(nullptr);
        dashboardStripDock_->setSelection(nullptr);
        dashboardStripDock_->setSelectionController(nullptr);
        dashboardStripDock_->setTimeViewport(nullptr);
        dashboardStripDock_->setDftStore(nullptr);
        dashboardStripDock_->setExperimentalShieldingMlStore(nullptr);
        dashboardStripDock_->setPanelModel(nullptr);
        dashboardStripDock_->setSignalModels(nullptr, nullptr);
        dashboardStripDock_->setContext(nullptr, nullptr);
    }
    if (dashboardController_)
        dashboardController_->setVisualizationContext({});
    if (selectionLabel_)
        selectionLabel_->setText(QStringLiteral("no selection"));
    if (inspectorDock_) {
        inspectorDock_->setContext(nullptr, nullptr);
        inspectorDock_->setFieldAvailability({});
        inspectorDock_->clearSelection();
    }
    if (measurementsDock_)
        measurementsDock_->setContext(nullptr, nullptr);

    delete cameraInputFilter_;
    cameraInputFilter_ = nullptr;
    delete picker_;
    picker_ = nullptr;
    delete scene_;
    scene_ = nullptr;

    delete playback_;
    playback_ = nullptr;
    delete timeViewport_;
    timeViewport_ = nullptr;

    delete dftStore_;
    dftStore_ = nullptr;
    delete experimentalMlStore_;
    experimentalMlStore_ = nullptr;
    activeExperimentalMlTensorDescriptor_.clear();
    activeExperimentalMlTensorAtom_.reset();
    experimentalMlTensorDisplayed_ = false;
    experimentalMlTensorDisplayedFrame_.reset();
    delete dashboardSelectionController_;
    dashboardSelectionController_ = nullptr;
    delete dashboardPanels_;
    dashboardPanels_ = nullptr;
    delete dashboardSignals_;
    dashboardSignals_ = nullptr;
    delete signalCatalog_;
    signalCatalog_ = nullptr;
    delete selection_;
    selection_ = nullptr;

    delete transformed_;
    transformed_ = nullptr;

    visualizationContext_ = {};
    fieldAvailability_.reset();
    lastDashboardSelectedCount_ = 0;
    loaded_.reset();
    updateMutantAlternateAction({});

    setEmptyState();
    syncRestServerContext();
}

void ReaderMainWindow::setEmptyState() {
    ASSERT_THREAD(this);
    refreshControlStates();
    updateFitModeLabel();
    if (emptyPlaceholder_)
        emptyPlaceholder_->setVisible(true);
    if (centralContainer_) {
        if (auto* stack = qobject_cast<QStackedLayout*>(centralContainer_->layout()))
            stack->setCurrentWidget(emptyPlaceholder_);
    }
    if (proteinLabel_)
        proteinLabel_->setText(QStringLiteral("No calcset loaded"));
    if (frameLabel_)
        frameLabel_->setText(QStringLiteral("frame —"));
    if (timeLabel_)
        timeLabel_->setText(QStringLiteral("t=— ps"));
    if (frameSlider_) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setRange(0, 0);
        frameSlider_->setValue(0);
    }
    if (fpsSpinner_) {
        const QSignalBlocker block(fpsSpinner_);
        fpsSpinner_->setRange(1, 60);
        fpsSpinner_->setValue(5);
    }
    setWindowTitle(QStringLiteral("h5-reader"));
}

void ReaderMainWindow::refreshControlStates() {
    ASSERT_THREAD(this);
    // Derive the whole operating state once — single source of truth.
    const bool loaded   = scene_ != nullptr;
    const int  frames   = (loaded && loaded_ && loaded_->conformation)
                            ? static_cast<int>(loaded_->conformation->frameCount()) : 0;
    const bool playable = loaded && frames > 1;            // single pose: nothing to play
    const bool playing  = playback_ && playback_->isPlaying();
    const bool traj     = loaded && visualizationContext_.hasTrajectory;
    const bool hasFocus = selection_ && selection_->hasFocus();
    const bool hasRings = loaded && loaded_ && loaded_->protein
                            && loaded_->protein->ringCount() > 0;
    // Filter (isolation) mode hides the whole-molecule overlays, so their
    // toggles can't do anything while it's on — keep them on the toolbar
    // (consistent layout) but disabled.
    const bool filtered = loaded && scene_->atomFilterActive();

    const auto en = [](QAction* a, bool e) { if (a) a->setEnabled(e); };

    // Transport — only meaningful for a multi-frame trajectory.
    en(playBackAction_,    playable);
    en(stepBackAction_,    playable);
    en(stepForwardAction_, playable);
    en(playForwardAction_, playable);
    en(stopAction_,        playable && playing);   // stop is inert when stopped
    if (frameSlider_) frameSlider_->setEnabled(playable);
    if (fpsSpinner_)  fpsSpinner_->setEnabled(playable);

    // Analysis controls.
    en(transformFitAction_,   loaded && transformed_ != nullptr);
    en(goToAtomAction_,       loaded && loaded_ && loaded_->protein);
    en(signalDisplaysAction_, loaded && hasFocus);

    // Overlays — gated on the data that makes each one mean something.
    en(showRibbonAction_,    loaded   && !filtered);
    en(showRingsAction_,     hasRings && !filtered);
    en(showButterflyAction_, hasRings && !filtered);
    en(showNullConeAction_,  hasRings && !filtered);
    en(showBFieldAction_,    hasRings && !filtered);
    en(showTrajectoryAction_, traj);   // selected-atom path; needs trajectory

    // Camera modes (enable gating + checked-state sync).
    updateCameraModeActions();

    // Filter button: usable whenever a run is loaded (so an active filter can
    // always be cleared) — its dropdown content adapts to the current focus.
    if (filterButton_) filterButton_->setEnabled(loaded);
    updateFilterButton();
}

QJsonObject ReaderMainWindow::uiStateJson() const {
    const bool loaded  = scene_ != nullptr;
    const int  frames  = (loaded && loaded_ && loaded_->conformation)
                           ? static_cast<int>(loaded_->conformation->frameCount()) : 0;
    const bool playing = playback_ && playback_->isPlaying();

    const auto ctl = [](const QAction* a) {
        QJsonObject o;
        o[QStringLiteral("present")] = (a != nullptr);
        o[QStringLiteral("enabled")] = (a != nullptr && a->isEnabled());
        if (a && a->isCheckable())
            o[QStringLiteral("checked")] = a->isChecked();
        return o;
    };

    QJsonObject controls;
    controls[QStringLiteral("playBack")]     = ctl(playBackAction_);
    controls[QStringLiteral("stepBack")]     = ctl(stepBackAction_);
    controls[QStringLiteral("stop")]         = ctl(stopAction_);
    controls[QStringLiteral("stepForward")]  = ctl(stepForwardAction_);
    controls[QStringLiteral("playForward")]  = ctl(playForwardAction_);
    controls[QStringLiteral("focus")]        = ctl(focusAction_);
    controls[QStringLiteral("transformFit")] = ctl(transformFitAction_);
    controls[QStringLiteral("goToAtom")]     = ctl(goToAtomAction_);
    controls[QStringLiteral("metrics")]      = ctl(signalDisplaysAction_);
    controls[QStringLiteral("ribbon")]       = ctl(showRibbonAction_);
    controls[QStringLiteral("rings")]        = ctl(showRingsAction_);
    controls[QStringLiteral("butterfly")]    = ctl(showButterflyAction_);
    controls[QStringLiteral("nullcone")]     = ctl(showNullConeAction_);
    controls[QStringLiteral("bfield")]       = ctl(showBFieldAction_);
    controls[QStringLiteral("trajectory")]   = ctl(showTrajectoryAction_);

    // Filter is a QToolButton (not a QAction), so report it explicitly:
    // present/enabled like the others, plus whether isolation is active and
    // how many residues are pinned. Keeps /ui/state an honest mirror of the UI.
    QJsonObject filter;
    filter[QStringLiteral("present")] = (filterButton_ != nullptr);
    filter[QStringLiteral("enabled")] = (filterButton_ && filterButton_->isEnabled());
    filter[QStringLiteral("active")]  = (scene_ && scene_->atomFilterActive());
    filter[QStringLiteral("count")]   = static_cast<int>(filterResidues_.size());
    controls[QStringLiteral("filter")] = filter;

    QJsonObject sel;
    sel[QStringLiteral("count")] = static_cast<int>(selection_ ? selection_->count() : 0);
    sel[QStringLiteral("focus")] = (selection_ && selection_->hasFocus());

    QJsonObject diagnostic;
    diagnostic[QStringLiteral("present")] = !lastDiagnosticMessage_.isEmpty();
    diagnostic[QStringLiteral("severity")] = lastDiagnosticSeverity_;
    diagnostic[QStringLiteral("source")] = lastDiagnosticSource_;
    diagnostic[QStringLiteral("message")] = lastDiagnosticMessage_;
    diagnostic[QStringLiteral("values")] = lastDiagnosticValues_;

    QJsonObject out;
    out[QStringLiteral("loaded")]        = loaded;
    out[QStringLiteral("protein")]       = (loaded && loaded_) ? loaded_->proteinId : QString();
    out[QStringLiteral("frames")]        = frames;
    out[QStringLiteral("currentFrame")]  = playback_ ? playback_->currentFrame() : 0;
    out[QStringLiteral("playing")]       = playing;
    out[QStringLiteral("playDirection")] = playback_ ? playback_->direction() : 1;
    out[QStringLiteral("selection")]     = sel;
    out[QStringLiteral("cameraMode")]    =
        (loaded && scene_ && scene_->cameraComposer())
            ? QString::fromLatin1(NameFor(scene_->cameraComposer()->mode().kind))
            : QStringLiteral("none");
    out[QStringLiteral("controls")]      = controls;
    out[QStringLiteral("diagnostic")]    = diagnostic;
    QJsonObject experimentalMl =
        experimentalShieldingMlRuntimeJson(fieldAvailability_.get(), loaded);
    const auto configuredExperimentalMl =
        resolveExperimentalShieldingMlRuntime();
    experimentalMl.insert(QStringLiteral("devicePreference"),
                          experimentalShieldingMlDevicePreference());
    experimentalMl.insert(
        QStringLiteral("configuredDevice"),
        configuredExperimentalMl ? QJsonValue(configuredExperimentalMl->device) : jsonNull());
    experimentalMl.insert(
        QStringLiteral("cpuFallbackConfigured"),
        configuredExperimentalMl && !configuredExperimentalMl->fallbackHelper.isEmpty());
    experimentalMl.insert(QStringLiteral("inferenceReady"),
                          experimentalMlStore_ && experimentalMlStore_->isReady());
    experimentalMl.insert(QStringLiteral("inferenceRunning"),
                          experimentalMlStore_ && experimentalMlStore_->isRunning());
    if (experimentalMlStore_) {
        experimentalMl.insert(QStringLiteral("activeModelId"), experimentalMlStore_->modelId());
        experimentalMl.insert(QStringLiteral("activeDevice"), experimentalMlStore_->device());
        experimentalMl.insert(QStringLiteral("usingCpuFallback"),
                              experimentalMlStore_->usingFallback());
        if (!experimentalMlStore_->errorReason().isEmpty())
            experimentalMl.insert(QStringLiteral("inferenceError"), experimentalMlStore_->errorReason());
    }
    QJsonObject tensorDisplay;
    const bool tensorActive =
        activeExperimentalMlTensorDescriptor_
            == QStringLiteral("ml:experimental_shielding_t2")
        && activeExperimentalMlTensorAtom_.has_value();
    const int tensorFrameI = playback_ ? playback_->currentFrame() : 0;
    const std::size_t tensorFrame =
        static_cast<std::size_t>(std::max(0, tensorFrameI));
    tensorDisplay.insert(QStringLiteral("active"), tensorActive);
    tensorDisplay.insert(QStringLiteral("descriptorId"),
                         tensorActive
                             ? QJsonValue(activeExperimentalMlTensorDescriptor_)
                             : jsonNull());
    tensorDisplay.insert(
        QStringLiteral("atom"),
        tensorActive
            ? QJsonValue(static_cast<qint64>(*activeExperimentalMlTensorAtom_))
            : jsonNull());
    tensorDisplay.insert(QStringLiteral("frame"),
                         tensorActive ? QJsonValue(static_cast<qint64>(tensorFrame))
                                      : jsonNull());
    tensorDisplay.insert(QStringLiteral("source"),
                         tensorActive
                             ? QJsonValue(QStringLiteral("Experimental Shielding ML"))
                             : jsonNull());
    tensorDisplay.insert(QStringLiteral("modelId"),
                         tensorActive && experimentalMlStore_
                             ? QJsonValue(experimentalMlStore_->modelId())
                             : jsonNull());
    const auto tensorValues =
        tensorActive && experimentalMlStore_
            ? experimentalMlStore_->tensor(
                  tensorFrame, *activeExperimentalMlTensorAtom_)
            : std::optional<std::array<double, 6>>{};
    tensorDisplay.insert(QStringLiteral("resident"), tensorValues.has_value());
    tensorDisplay.insert(
        QStringLiteral("displayed"),
        experimentalMlTensorDisplayed_
            && experimentalMlTensorDisplayedFrame_.has_value()
            && *experimentalMlTensorDisplayedFrame_ == tensorFrame);
    if (tensorValues) {
        tensorDisplay.insert(QStringLiteral("sigmaIsoPpm"), (*tensorValues)[0]);
        QJsonArray t2;
        for (std::size_t i = 1; i < tensorValues->size(); ++i)
            t2.append((*tensorValues)[i]);
        tensorDisplay.insert(QStringLiteral("t2"), t2);
    }
    experimentalMl.insert(QStringLiteral("tensorDisplay"), tensorDisplay);
    out[QStringLiteral("experimentalShieldingMl")] = experimentalMl;
    return out;
}

void ReaderMainWindow::applyOverlayActionState() {
    ASSERT_THREAD(this);
    if (!scene_)
        return;
    if (showRibbonAction_ && scene_->ribbonOverlay())
        scene_->ribbonOverlay()->setVisible(showRibbonAction_->isChecked());
    if (showRingsAction_ && scene_->ringPolygonOverlay())
        scene_->ringPolygonOverlay()->setVisible(showRingsAction_->isChecked());
    if (showButterflyAction_ && scene_->fieldGridOverlay())
        scene_->fieldGridOverlay()->setVisible(showButterflyAction_->isChecked());
    if (showNullConeAction_ && scene_->fieldGridOverlay())
        scene_->fieldGridOverlay()->setNullConeVisible(showNullConeAction_->isChecked());
    if (showBFieldAction_ && scene_->bfieldStreamOverlay())
        scene_->bfieldStreamOverlay()->setVisible(showBFieldAction_->isChecked());
    if (showTrajectoryAction_ && scene_->atomTrajectoryOverlay())
        scene_->atomTrajectoryOverlay()->setVisible(showTrajectoryAction_->isChecked());
}

void ReaderMainWindow::updateMutantAlternateAction(const QString& alternatePath) {
    ASSERT_THREAD(this);
    if (mutantAlternateAction_) {
        if (fileMenu_)
            fileMenu_->removeAction(mutantAlternateAction_);
        delete mutantAlternateAction_.data();
        mutantAlternateAction_ = nullptr;
    }
    if (alternatePath.isEmpty() || !fileMenu_)
        return;

    mutantAlternateAction_ = fileMenu_->addAction(
        QStringLiteral("Open mutant alternate (ALA)…"));
    mutantAlternateAction_->setToolTip(QStringLiteral(
        "This run is a mutant pair; WT is opened in this window. "
        "Click to load the ALA pose in this window: %1").arg(alternatePath));
    ACONNECT(mutantAlternateAction_.data(), &QAction::triggered, this, [this, alternatePath]() {
        if (!loadRunPath(alternatePath)) {
            QMessageBox::critical(this,
                                  QStringLiteral("Open calcset failed"),
                                  lastLoadError());
        }
    });
}

void ReaderMainWindow::syncRestServerContext() {
    ASSERT_THREAD(this);
    if (!restServer_)
        return;
    restServer_->setContext(scene_,
                            selection_,
                            dashboardSignals_,
                            dashboardPanels_,
                            dashboardSelectionController_.data(),
                            dashboardController_.data(),
                            signalCatalog_,
                            playback_,
                            loaded_.get(),
                            this,
                            this,
                            transformed_);
}

void ReaderMainWindow::showEvent(QShowEvent* event) {
    QMainWindow::showEvent(event);
    if (glInfoLogged_) return;
    glInfoLogged_ = true;

    // The VTK widget hasn't painted yet at showEvent time, so the GL context
    // isn't current and ReportCapabilities() is empty. Hook the render window's
    // first EndEvent — the actual moment the context is current and a frame has
    // rendered — and read caps there. One-shot: the callback removes its own
    // observer. Event-driven, not a guessed "next tick". Mirrors
    // MoleculeScene's EndEvent observer.
    if (!renderWindow_)
        return;
    auto capsCb = vtkSmartPointer<vtkCallbackCommand>::New();
    capsCb->SetClientData(this);
    capsCb->SetCallback(
        [](vtkObject* /*caller*/, unsigned long, void* clientData, void*) {
            auto* self = static_cast<ReaderMainWindow*>(clientData);
            if (!self || !self->renderWindow_)
                return;
            if (self->glCapsObserverTag_ != 0) {
                self->renderWindow_->RemoveObserver(self->glCapsObserverTag_);
                self->glCapsObserverTag_ = 0;
            }
            const QString caps =
                QString::fromUtf8(self->renderWindow_->ReportCapabilities());
            const QStringList wanted = {
                QStringLiteral("OpenGL vendor"),
                QStringLiteral("OpenGL renderer"),
                QStringLiteral("OpenGL version"),
                QStringLiteral("OpenGL vendor-specific"),
            };
            for (const QString& line : caps.split(QChar('\n'))) {
                for (const QString& key : wanted) {
                    if (line.contains(key, Qt::CaseInsensitive)) {
                        qCInfo(cWindow).noquote() << "GL:" << line.trimmed();
                        break;
                    }
                }
            }
        });
    glCapsObserverTag_ = renderWindow_->AddObserver(vtkCommand::EndEvent, capsCb);
}

ReaderMainWindow::~ReaderMainWindow() {
    // Most cleanup runs in shutdown(). The destructor only handles the
    // pathological case where shutdown() was never called (e.g. window
    // deleted outside the normal quit flow).
    if (!shutdownDone_) {
        qCWarning(cWindow).noquote()
            << "destructor called without prior shutdown(); running now";
        shutdown();
    }
}

QJsonArray ReaderMainWindow::inspectorTreeJson() const {
    return inspectorDock_ ? inspectorDock_->dumpTree() : QJsonArray();
}


quint16 ReaderMainWindow::startRestServer(quint16 port) {
    ASSERT_THREAD(this);
    if (!loaded_ || !loaded_->protein || !loaded_->conformation) {
        qCCritical(cWindow).noquote() << "REST start refused: loader result not wired";
        return 0;
    }
    if (restServer_) {
        qCWarning(cWindow).noquote() << "REST server already running; ignoring re-start";
        return 0;
    }
    restServer_ = new RestServer(this);
    restServer_->setContext(scene_,
                            selection_,
                            dashboardSignals_,
                            dashboardPanels_,
                            dashboardSelectionController_.data(),
                            dashboardController_.data(),
                            signalCatalog_,
                            playback_,
                            loaded_.get(),
                            this,
                            this,
                            transformed_);
    const quint16 bound = restServer_->listen(port);
    if (bound == 0) {
        qCCritical(cWindow).noquote() << "REST server failed to bind port" << port;
        restServer_->deleteLater();
        restServer_ = nullptr;
    }
    return bound;
}

bool ReaderMainWindow::setOverlayVisible(const QString& name, bool on) {
    ASSERT_THREAD(this);
    // Map a stable automation key → the toolbar QAction, then setChecked()
    // so the already-connected QAction::toggled handler runs the real overlay
    // logic (setVisible + refreshCurrentFrame for the per-frame kernel
    // overlays). This keeps REST control and the toolbar UI on one code path.
    const QString key = name.toLower();
    QPointer<QAction> a;
    if (key == QStringLiteral("ribbon"))
        a = showRibbonAction_;
    else if (key == QStringLiteral("rings"))
        a = showRingsAction_;
    else if (key == QStringLiteral("butterfly") || key == QStringLiteral("fieldgrid")
             || key == QStringLiteral("field") || key == QStringLiteral("isosurface"))
        a = showButterflyAction_;
    else if (key == QStringLiteral("nullcone") || key == QStringLiteral("null_cone")
             || key == QStringLiteral("ring_null") || key == QStringLiteral("ringnull")
             || key == QStringLiteral("magiccone") || key == QStringLiteral("magic_cone"))
        a = showNullConeAction_;
    else if (key == QStringLiteral("bfield") || key == QStringLiteral("streamlines")
             || key == QStringLiteral("stream"))
        a = showBFieldAction_;
    else if (key == QStringLiteral("trajectory") || key == QStringLiteral("path"))
        a = showTrajectoryAction_;
    if (!a)
        return false;
    if (a->isChecked() != on)
        a->setChecked(on);   // fires QAction::toggled → overlay setVisible + refresh
    return true;
}

bool ReaderMainWindow::setFieldThreshold(double ppm) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    scene_->fieldGridOverlay()->setThresholdPpm(ppm);
    // Only the isovalue changed (scalar field unchanged), so a plain render
    // re-runs the contour filters — no need to re-evaluate the kernel grid.
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

bool ReaderMainWindow::setFieldExtent(double sigmaA) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    // setGaussianExtent re-evaluates the scalar grid (the taper reshapes the
    // field); a render then re-contours it.
    scene_->fieldGridOverlay()->setGaussianExtent(sigmaA);
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

bool ReaderMainWindow::setFieldPeak(double amplitude) {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->fieldGridOverlay())
        return false;
    scene_->fieldGridOverlay()->setGaussianPeak(amplitude);
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    return true;
}

void ReaderMainWindow::setResidueFilter(const std::vector<std::size_t>& residues) {
    ASSERT_THREAD(this);
    if (!scene_ || !loaded_ || !loaded_->protein)
        return;

    // The whole-molecule overlays (ribbon, rings, field) would keep
    // drawing the entire structure and defeat the isolation, so hide them while
    // filtered and restore them (per the toolbar toggles) when the filter clears.
    const auto restoreOverlays = [this]() { applyOverlayActionState(); };
    const auto hideOverlays = [this]() {
        if (scene_->ribbonOverlay())          scene_->ribbonOverlay()->setVisible(false);
        if (scene_->ringPolygonOverlay())     scene_->ringPolygonOverlay()->setVisible(false);
        if (scene_->fieldGridOverlay()) {
            scene_->fieldGridOverlay()->setVisible(false);
            scene_->fieldGridOverlay()->setNullConeVisible(false);
        }
        if (scene_->bfieldStreamOverlay())    scene_->bfieldStreamOverlay()->setVisible(false);
    };

    const auto* protein = loaded_->protein.get();
    std::vector<std::size_t> atoms;
    if (!residues.empty()) {
        std::vector<char> keep(protein->residueCount(), 0);
        for (std::size_t r : residues)
            if (r < keep.size()) keep[r] = 1;
        for (std::size_t i = 0; i < protein->atomCount(); ++i) {
            const auto& atom = protein->atom(i);
            if (atom.residueIndex >= 0
                && static_cast<std::size_t>(atom.residueIndex) < keep.size()
                && keep[static_cast<std::size_t>(atom.residueIndex)])
                atoms.push_back(i);
        }
    }

    if (atoms.empty()) {
        filterResidues_.clear();
        restoreOverlays();
        scene_->clearAtomFilter();
    } else {
        filterResidues_ = residues;   // single source of truth: keeps the button
        hideOverlays();               // label + checklist consistent whether the
        scene_->setAtomFilter(atoms); // filter was driven by the menu or REST
    }
    refreshControlStates();   // grey overlay toggles while filtered; restore on clear
}

void ReaderMainWindow::setDocksVisible(bool visible) {
    ASSERT_THREAD(this);

    // Hide path: stash each dock's pre-hide visibility so a later restore
    // can return individually-user-hidden docks to their hidden state.
    // No-op if already hidden (don't double-stash).
    if (!visible) {
        if (docksHidden_)
            return;
        stashedDockVisibility_.clear();
        const std::vector<QDockWidget*> docks = {
            inspectorDock_, dashboardStripDock_
        };
        for (QDockWidget* d : docks) {
            if (!d) continue;
            stashedDockVisibility_.push_back({QPointer<QDockWidget>(d), d->isVisible()});
            d->setVisible(false);
        }
        docksHidden_ = true;
        qCInfo(cWindow).noquote()
            << "docks hidden | count=" << stashedDockVisibility_.size();
        return;
    }

    // Restore path: walk the stash; QPointer-safe iteration in case any
    // dock was destroyed since the hide. Each dock returns to its stashed
    // visibility, so a dock that was already hidden before the harness
    // requested hide stays hidden.
    if (!docksHidden_)
        return;
    for (const DockVis& dv : stashedDockVisibility_) {
        if (dv.dock)
            dv.dock->setVisible(dv.wasVisible);
    }
    qCInfo(cWindow).noquote()
        << "docks restored | count=" << stashedDockVisibility_.size();
    stashedDockVisibility_.clear();
    docksHidden_ = false;
}

void ReaderMainWindow::setDashboardDockVisible(bool visible) {
    ASSERT_THREAD(this);
    if (!dashboardStripDock_)
        return;
    if (visible) {
        revealDockQueued(dashboardStripDock_);
    } else {
        dashboardStripDock_->setVisible(false);
    }
}

bool ReaderMainWindow::dashboardDockVisible() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ && dashboardStripDock_->isVisible();
}

int ReaderMainWindow::dashboardDockWidth() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->width() : 0;
}

bool ReaderMainWindow::dashboardDockRaised() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ && dashboardStripDock_->isVisible()
        && !dashboardStripDock_->visibleRegion().isEmpty();
}

int ReaderMainWindow::dashboardOwnedPanelCount() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->ownedPanelCount() : 0;
}

QJsonArray ReaderMainWindow::dashboardPanelManifest() const {
    ASSERT_THREAD(this);
    QJsonArray out;
    if (!dashboardStripDock_) return out;
    for (const PanelDisplayData& d : dashboardStripDock_->ownedPanelDisplayData()) {
        QJsonObject o{
            {QStringLiteral("kind"), d.kind},
            {QStringLiteral("title"), d.title},
            {QStringLiteral("series_count"), d.seriesCount},
            {QStringLiteral("point_count"), d.pointCount},
            {QStringLiteral("finite_count"), d.finiteCount},
            {QStringLiteral("nan_count"), d.nanCount},
            {QStringLiteral("empty"), d.empty()},
        };
        if (d.finiteCount > 0) {
            o.insert(QStringLiteral("min"), d.minValue);
            o.insert(QStringLiteral("max"), d.maxValue);
        }
        if (!d.descriptorId.isEmpty())
            o.insert(QStringLiteral("descriptor_id"), d.descriptorId);
        out.append(o);
    }
    return out;
}

int ReaderMainWindow::dashboardStripTrackCount() const {
    ASSERT_THREAD(this);
    return dashboardStripDock_ ? dashboardStripDock_->stripTrackCount() : 0;
}

bool ReaderMainWindow::openSignalDisplayPicker(QString* blockedReason) {
    ASSERT_THREAD(this);
    if (blockedReason)
        blockedReason->clear();
    if (!signalDisplayDialog_) {
        if (blockedReason)
            *blockedReason = QStringLiteral("signal display dialog not wired");
        return false;
    }
    if (!selection_ || !selection_->hasFocus()) {
        if (blockedReason)
            *blockedReason = QStringLiteral("Metrics action is disabled because AtomSelection has no focus.");
        return false;
    }
    if (signalDisplaysAction_ && !signalDisplaysAction_->isEnabled()) {
        if (blockedReason)
            *blockedReason = QStringLiteral("Metrics action is disabled.");
        return false;
    }
    if (playback_)
        signalDisplayDialog_->setFrame(playback_->currentFrame());
    signalDisplayDialog_->refreshCatalog();
    signalDisplayDialog_->show();
    signalDisplayDialog_->raise();
    signalDisplayDialog_->activateWindow();
    return true;
}

QJsonObject ReaderMainWindow::signalDisplayPickerState() const {
    ASSERT_THREAD(this);
    if (!signalDisplayDialog_)
        return QJsonObject{{"open", false}};
    return signalDisplayDialog_->pickerState();
}

QJsonObject ReaderMainWindow::addSelectedSignalFromPicker() {
    ASSERT_THREAD(this);
    if (signalDisplayDialog_)
        signalDisplayDialog_->onAddSelected();
    return signalDisplayPickerState();
}

void ReaderMainWindow::revealDockQueued(QDockWidget* dock) {
    ASSERT_THREAD(this);
    if (!dock)
        return;
    dock->setVisible(true);
    // Queue the resize behind the show/relayout events already in this window's
    // queue (deterministic FIFO), instead of guessing a tick on the timer wheel.
    // Same queued-invoke pattern as MoleculeScene::requestRender.
    QMetaObject::invokeMethod(this, [this, dock = QPointer<QDockWidget>(dock)]() {
        if (!dock)
            return;
        resizeDocks({dock.data()}, {360}, Qt::Horizontal);
        dock->raise();
    }, Qt::QueuedConnection);
}

void ReaderMainWindow::updateSelectionStatus() {
    ASSERT_THREAD(this);
    if (!selectionLabel_)
        return;
    if (!selection_ || selection_->empty()) {
        selectionLabel_->setText(QStringLiteral("no selection"));
        return;
    }
    const std::size_t        n = selection_->count();
    const model::GeometryKind k = selection_->geometryKind();
    if (k == model::GeometryKind::None)
        selectionLabel_->setText(QStringLiteral("%1 atom selected").arg(n));
    else
        selectionLabel_->setText(QStringLiteral("%1 atoms · %2")
            .arg(n).arg(QString::fromLatin1(model::NameForGeometryKind(k))));
}

void ReaderMainWindow::buildFilterMenu() {
    ASSERT_THREAD(this);
    if (!filterMenu_)
        return;
    filterMenu_->clear();

    // Leaving filter mode is always offered, enabled only while a filter is on.
    QAction* showAll = filterMenu_->addAction(QStringLiteral("Show whole structure"));
    showAll->setEnabled(scene_ && scene_->atomFilterActive());
    ACONNECT(showAll, &QAction::triggered, this, [this]() {
        setResidueFilter({});      // empty restores the full structure + overlays
    });
    filterMenu_->addSeparator();

    // The checklist is the residues near the focused atom. No focus → nothing
    // to list; say so rather than present an empty menu.
    if (!selection_ || !selection_->hasFocus() || !filterNearby_) {
        QAction* hint = filterMenu_->addAction(
            QStringLiteral("Select an atom to list nearby residues"));
        hint->setEnabled(false);
        return;
    }

    // Rebuild the neighbourhood for the current focus + frame, then offer one
    // checkable row per nearby residue, nearest first.
    filterNearby_->setAnchor({selection_->focus(),
                              playback_ ? playback_->currentFrame() : 0});

    struct Row { std::size_t residue; double dist; QString label; };
    std::vector<Row> rows;
    const int rowN = filterNearby_->rowCount();
    for (int r = 0; r < rowN; ++r) {
        const NearbySignalModel::Candidate* c =
            filterNearby_->candidateAt(filterNearby_->index(r, 0));
        if (!c || c->kind != NearbySignalModel::CandidateKind::Residue
              || !c->residueContext)
            continue;
        rows.push_back({*c->residueContext, c->distanceAngstrom, c->label});
    }
    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.dist < b.dist; });

    if (rows.empty()) {
        QAction* none = filterMenu_->addAction(
            QStringLiteral("No residues within %1 Å")
                .arg(filterNearby_->radiusAngstrom(), 0, 'f', 1));
        none->setEnabled(false);
        return;
    }

    for (const Row& row : rows) {
        QAction* a = filterMenu_->addAction(
            QStringLiteral("%1 · %2 Å").arg(row.label).arg(row.dist, 0, 'f', 1));
        a->setCheckable(true);
        a->setChecked(std::find(filterResidues_.begin(), filterResidues_.end(),
                                row.residue) != filterResidues_.end());
        const std::size_t residue = row.residue;
        ACONNECT(a, &QAction::toggled, this, [this, residue](bool on) {
            onFilterResidueToggled(residue, on);
        });
    }
}

void ReaderMainWindow::onFilterResidueToggled(std::size_t residue, bool on) {
    ASSERT_THREAD(this);
    // Build the desired set locally; setResidueFilter is the sole writer of
    // filterResidues_, so it (not this handler) keeps the button label synced.
    std::vector<std::size_t> next = filterResidues_;
    const auto it = std::find(next.begin(), next.end(), residue);
    if (on) {
        if (it == next.end())
            next.push_back(residue);
    } else if (it != next.end()) {
        next.erase(it);
    }
    setResidueFilter(next);   // empty set restores the whole structure
}

void ReaderMainWindow::updateFilterButton() {
    ASSERT_THREAD(this);
    if (!filterButton_)
        return;
    const bool active = scene_ && scene_->atomFilterActive();
    filterButton_->setText(active
        ? QStringLiteral("Filter (%1)").arg(filterResidues_.size())
        : QStringLiteral("Filter"));
}

void ReaderMainWindow::resetDashboardStateForRunLoad() {
    ASSERT_THREAD(this);
    if (dashboardSelectionController_) {
        dashboardSelectionController_->clearAllMetrics();
        lastDashboardSelectedCount_ = dashboardSelectionController_->selectedCount();
    }
    if (dashboardStripDock_)
        dashboardStripDock_->setVisible(false);
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_reset_on_load count=%1 dock_visible=%2")
               .arg(dashboardSelectionController_ ? dashboardSelectionController_->selectedCount() : 0)
               .arg(dashboardStripDock_ && dashboardStripDock_->isVisible() ? 1 : 0);
}

void ReaderMainWindow::handleErrorBusReport(h5reader::diagnostics::Severity severity,
                                            const QString& source,
                                            const QString& message,
                                            const QString& values) {
    ASSERT_THREAD(this);
    lastDiagnosticSeverity_ = diagnosticSeverityName(severity);
    lastDiagnosticSource_ = source;
    lastDiagnosticMessage_ = message;
    lastDiagnosticValues_ = values;

    QString body = message;
    if (!source.isEmpty())
        body = QStringLiteral("%1: %2").arg(source, message);
    const QString display = QStringLiteral("%1: %2")
        .arg(lastDiagnosticSeverity_, body);

    if (diagnosticLabel_) {
        diagnosticLabel_->setText(display);
        diagnosticLabel_->setToolTip(values.isEmpty()
            ? display
            : QStringLiteral("%1\n%2").arg(display, values));
        diagnosticLabel_->setStyleSheet(diagnosticSeverityStyle(severity));
        diagnosticLabel_->setVisible(true);
    }
    statusBar()->showMessage(display, 10000);
}

void ReaderMainWindow::shutdown() {
    ASSERT_THREAD(this);
    if (shutdownDone_) return;
    shutdownDone_ = true;

    qCInfo(cWindow).noquote() << "shutdown entered";

    // Per spec/viewport_pipeline_2026-05-30.md §4.4:
    //
    // 1. Stop the REST server SYNCHRONOUSLY. The /shutdown endpoint
    //    fires from a request handler; the server needs to drain
    //    before timers stop so a follow-up request can't trigger a
    //    race with timer teardown.
    if (restServer_) {
        // RestServer doesn't expose stopListening(); the QHttpServer
        // owned by it tears down when the RestServer is deleted, but
        // deleteLater on shutdown is enough for this path because
        // aboutToQuit drains the event loop afterwards. We do hold a
        // direct pointer; do a synchronous delete here.
        delete restServer_;
        restServer_ = nullptr;
    }

    // 2. Stop every timer owned by us or our children. The generic
    //    findChildren sweep catches QtPlaybackController's timer too.
    const auto timers = findChildren<QTimer*>();
    for (auto* timer : timers) {
        if (timer->isActive()) timer->stop();
    }

    // 3. Detach the render window from the widget BEFORE dropping our
    //    smart pointer. setRenderWindow(nullptr) makes the context
    //    current and calls Finalize on the old render window via the
    //    adapter's destructor (QVTKRenderWindowAdapter.cxx:150-166).
    //    The explicit renderWindow_->Finalize() that used to live here
    //    is gone — doing it AFTER detaching the widget left the adapter
    //    holding a destroyed window for the brief moment between the
    //    two calls.
    if (vtkWidget_) {
        vtkWidget_->setRenderWindow(static_cast<vtkGenericOpenGLRenderWindow*>(nullptr));
    }

    qCInfo(cWindow).noquote() << "shutdown done";
}

void ReaderMainWindow::buildUi() {
    centralContainer_ = new QWidget(this);
    auto* stack = new QStackedLayout(centralContainer_);
    stack->setContentsMargins(0, 0, 0, 0);
    stack->setStackingMode(QStackedLayout::StackAll);

    vtkWidget_    = new QVTKOpenGLNativeWidget(centralContainer_);
    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWidget_->setRenderWindow(renderWindow_);
    stack->addWidget(vtkWidget_);

    emptyPlaceholder_ = new QLabel(QStringLiteral("Open a calcset (.LGS) to begin."), centralContainer_);
    emptyPlaceholder_->setAlignment(Qt::AlignCenter);
    emptyPlaceholder_->setWordWrap(true);
    emptyPlaceholder_->setAttribute(Qt::WA_TransparentForMouseEvents);
    emptyPlaceholder_->setStyleSheet(QStringLiteral(
        "QLabel { color: #5f6872; background: #fafafa; font-size: 18px; }"));
    stack->addWidget(emptyPlaceholder_);
    stack->setCurrentWidget(emptyPlaceholder_);
    setCentralWidget(centralContainer_);

    // File ▸ Open… loads a calcset into this window.
    fileMenu_ = menuBar()->addMenu(QStringLiteral("&File"));
    auto* openFileAct = fileMenu_->addAction(QStringLiteral("Open…"));
    openFileAct->setShortcut(QKeySequence::Open);  // Ctrl+O — pick a .LGS file with the mouse
    ACONNECT(openFileAct, &QAction::triggered, this, &ReaderMainWindow::onOpenFile);

    auto* openDirAct = fileMenu_->addAction(QStringLiteral("Open Directory…"));
    ACONNECT(openDirAct, &QAction::triggered, this, &ReaderMainWindow::onOpenDirectory);

    // File ▸ Recent — populated from QSettings during restoreAllSettings.
    // Empty until then; each entry loads into this window on click.
    recentMenu_ = fileMenu_->addMenu(QStringLiteral("&Recent"));
    recentMenu_->setObjectName(QStringLiteral("RecentMenu"));
}

namespace {
// Transport-control glyphs drawn in a given colour so they stay legible on any
// palette — Qt's SP_Media* standard icons render near-black on the dark Fusion
// toolbar. Centred in a 32 px square.
enum class TransportGlyph { PlayForward, PlayBackward, StepForward, StepBackward, Stop };

QIcon makeTransportIcon(TransportGlyph kind, const QColor& color) {
    constexpr int s = 32;
    QPixmap pm(s, s);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setPen(Qt::NoPen);
    p.setBrush(color);
    const double m   = s * 0.30;             // margin
    const double x0  = m, y0 = m, x1 = s - m, y1 = s - m, yc = s * 0.5;
    const double bar = (x1 - x0) * 0.32;     // bar thickness for the step glyphs
    switch (kind) {
        case TransportGlyph::PlayForward: {
            QPolygonF t; t << QPointF(x0, y0) << QPointF(x1, yc) << QPointF(x0, y1);
            p.drawPolygon(t);
            break;
        }
        case TransportGlyph::PlayBackward: {
            QPolygonF t; t << QPointF(x1, y0) << QPointF(x0, yc) << QPointF(x1, y1);
            p.drawPolygon(t);
            break;
        }
        case TransportGlyph::StepForward: {
            QPolygonF t; t << QPointF(x0, y0) << QPointF(x1 - bar, yc) << QPointF(x0, y1);
            p.drawPolygon(t);
            p.drawRect(QRectF(x1 - bar, y0, bar, y1 - y0));
            break;
        }
        case TransportGlyph::StepBackward: {
            QPolygonF t; t << QPointF(x1, y0) << QPointF(x0 + bar, yc) << QPointF(x1, y1);
            p.drawPolygon(t);
            p.drawRect(QRectF(x0, y0, bar, y1 - y0));
            break;
        }
        case TransportGlyph::Stop: {
            p.drawRect(QRectF(x0, y0, x1 - x0, y1 - y0));
            break;
        }
    }
    p.end();
    return QIcon(pm);
}

// A QMenu that does NOT close when a checkable item is clicked, so the Filter
// checklist lets you tick several residues in one go.
class FilterMenu final : public QMenu {
public:
    using QMenu::QMenu;
protected:
    void mouseReleaseEvent(QMouseEvent* e) override {
        QAction* a = activeAction();
        if (a && a->isEnabled() && a->isCheckable()) {
            a->trigger();   // toggle + fire; keep the menu open
            return;
        }
        QMenu::mouseReleaseEvent(e);
    }
};
}  // namespace

void ReaderMainWindow::buildToolbar() {
    auto* tb = addToolBar(QStringLiteral("Playback"));
    tb->setObjectName(QStringLiteral("PlaybackToolbar"));
    tb->setMovable(false);
    playbackToolbar_ = tb;
    QFont toolbarFont = tb->font();
    if (toolbarFont.pointSize() > 8)
        toolbarFont.setPointSize(toolbarFont.pointSize() - 1);
    else if (toolbarFont.pixelSize() > 10)
        toolbarFont.setPixelSize(toolbarFont.pixelSize() - 1);
    tb->setFont(toolbarFont);

    // Transport: ⏪ play-back · ◀ step-back · ■ stop · ▶ step-fwd · ⏩ play-fwd.
    // Custom glyphs in the palette text colour (legible on the dark toolbar).
    const QColor glyph = palette().color(QPalette::ButtonText);

    playBackAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::PlayBackward, glyph),
        QStringLiteral("Play backward"));
    playBackAction_->setToolTip(QStringLiteral("Play continuously, backward in time."));
    ACONNECT(playBackAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->playBackward(); });

    stepBackAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::StepBackward, glyph),
        QStringLiteral("Step back"));
    stepBackAction_->setToolTip(QStringLiteral("Step one frame back."));
    ACONNECT(stepBackAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->stepBackward(); });

    stopAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::Stop, glyph),
        QStringLiteral("Stop"));
    stopAction_->setToolTip(QStringLiteral("Stop the animation."));
    ACONNECT(stopAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->pause(); });

    stepForwardAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::StepForward, glyph),
        QStringLiteral("Step forward"));
    stepForwardAction_->setToolTip(QStringLiteral("Step one frame forward."));
    ACONNECT(stepForwardAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->stepForward(); });

    playForwardAction_ = tb->addAction(
        makeTransportIcon(TransportGlyph::PlayForward, glyph),
        QStringLiteral("Play forward"));
    playForwardAction_->setToolTip(QStringLiteral("Play continuously, forward in time."));
    ACONNECT(playForwardAction_.data(), &QAction::triggered,
             this, [this]() { if (playback_) playback_->playForward(); });

    tb->addSeparator();

    frameSlider_ = new QSlider(Qt::Horizontal, tb);
    frameSlider_->setMinimumWidth(400);
    tb->addWidget(frameSlider_);
    ACONNECT(frameSlider_.data(), &QSlider::valueChanged,
             this, [this](int frame) {
                 if (playback_) playback_->setFrame(frame);
             });

    tb->addSeparator();
    tb->addWidget(new QLabel(QStringLiteral("fps"), tb));
    fpsSpinner_ = new QSpinBox(tb);
    fpsSpinner_->setSuffix(QStringLiteral(" /s"));
    tb->addWidget(fpsSpinner_);
    ACONNECT(fpsSpinner_.data(), qOverload<int>(&QSpinBox::valueChanged),
             this, [this](int fps) {
                 if (playback_) playback_->setFps(fps);
             });

    addToolBarBreak();
    tb = addToolBar(QStringLiteral("Tools"));
    tb->setObjectName(QStringLiteral("ToolsToolbar"));
    tb->setMovable(false);
    toolsToolbar_ = tb;
    tb->setFont(toolbarFont);

    // Camera-mode cluster removed: Focus is now a single de-emphasized toggle
    // at the toolbar tail (built below). The Newman / Plane-lock / Free buttons
    // are gone; the composer keeps every mode for the dashboard-reveal + REST
    // paths, which are their real consumers.
    transformFitAction_ = tb->addAction(QStringLiteral("Mode: Locked backbone  ⇄"));
    transformFitAction_->setToolTip(fitModeToolTip());
    ACONNECT(transformFitAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onTransformFitClicked);

    tb->addSeparator();

    goToAtomAction_ = tb->addAction(QStringLiteral("Go..."));
    goToAtomAction_->setEnabled(false);
    goToAtomAction_->setToolTip(QStringLiteral(
        "Jump to a residue number, atom, and frame."));
    ACONNECT(goToAtomAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onGoToAtomTriggered);

    signalDisplaysAction_ = tb->addAction(QStringLiteral("Metrics..."));
    signalDisplaysAction_->setEnabled(false);
    signalDisplaysAction_->setToolTip(QStringLiteral("Select a nearby atom or residue and add a metric display."));
    ACONNECT(signalDisplaysAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onOpenSignalDisplays);

    // Display isolation ("Filter"): a dropdown checklist of residues near the
    // focused atom. Ticking residues enters filter mode (only those render);
    // "Show whole structure" leaves it. Disabled until a run is loaded; the
    // dropdown is (re)built lazily on aboutToShow so it tracks the current
    // focus + frame. First QToolButton-with-menu on this toolbar.
    filterMenu_ = new FilterMenu(this);
    filterButton_ = new QToolButton(this);
    filterButton_->setText(QStringLiteral("Filter"));
    filterButton_->setToolButtonStyle(Qt::ToolButtonTextOnly);
    filterButton_->setPopupMode(QToolButton::InstantPopup);
    filterButton_->setMenu(filterMenu_);
    filterButton_->setEnabled(false);
    filterButton_->setToolTip(QStringLiteral(
        "Show only chosen residues near the selected atom and step through "
        "frames isolated. Select an atom first."));
    tb->addWidget(filterButton_);
    ACONNECT(filterMenu_.data(), &QMenu::aboutToShow,
             this, &ReaderMainWindow::buildFilterMenu);

    tb->addSeparator();

    // Overlay toggles — inert until a scene has been loaded.
    showRibbonAction_ = tb->addAction(QStringLiteral("Ribbon"));
    showRibbonAction_->setCheckable(true);
    showRibbonAction_->setChecked(true);
    showRibbonAction_->setToolTip(QStringLiteral(
        "Backbone ribbon; secondary structure driven by per-frame DSSP."));

    showRingsAction_ = tb->addAction(QStringLiteral("Rings"));
    showRingsAction_->setCheckable(true);
    showRingsAction_->setChecked(true);
    showRingsAction_->setToolTip(QStringLiteral(
        "Aromatic ring polygons + normal arrows (per-frame ring_geometry)."));

    showButterflyAction_ = tb->addAction(QStringLiteral("Butterfly"));
    showButterflyAction_->setCheckable(true);
    showButterflyAction_->setChecked(false);   // off by default — expensive
    showButterflyAction_->setToolTip(QStringLiteral(
        "BS / HM volumetric isosurfaces around each aromatic ring. "
        "Re-evaluates closed-form kernel per frame on a 20³ grid."));

    showNullConeAction_ = tb->addAction(QStringLiteral("Ring null"));
    showNullConeAction_->setCheckable(true);
    showNullConeAction_->setChecked(false);
    showNullConeAction_->setToolTip(QStringLiteral(
        "Transparent ring null surface: the magic-angle boundary used for "
        "operational ring-null-crossing tests."));

    showBFieldAction_ = tb->addAction(QStringLiteral("B-field"));
    showBFieldAction_->setCheckable(true);
    showBFieldAction_->setChecked(false);   // off by default — expensive
    showBFieldAction_->setToolTip(QStringLiteral(
        "Biot-Savart B-field streamlines around each aromatic ring, "
        "seeded on a circle at 1.5× ring radius, coloured by |B|."));

    showTrajectoryAction_ = tb->addAction(QStringLiteral("Trajectory"));
    showTrajectoryAction_->setCheckable(true);
    showTrajectoryAction_->setChecked(false);   // off by default
    showTrajectoryAction_->setToolTip(QStringLiteral(
        "Focused-atom trajectory envelope across the loaded trajectory. Frame "
        "changes move the atom through the same shell."));

    ACONNECT(showRibbonAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->ribbonOverlay()) return;
                 scene_->ribbonOverlay()->setVisible(on);
                 scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showRingsAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->ringPolygonOverlay()) return;
                 scene_->ringPolygonOverlay()->setVisible(on);
                 scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showButterflyAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->fieldGridOverlay()) return;
                 scene_->fieldGridOverlay()->setVisible(on);
                 if (on) scene_->refreshCurrentFrame();
                 else    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showNullConeAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->fieldGridOverlay()) return;
                 scene_->fieldGridOverlay()->setNullConeVisible(on);
                 if (on) scene_->refreshCurrentFrame();
                 else    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showBFieldAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                 if (!scene_ || !scene_->bfieldStreamOverlay()) return;
                 scene_->bfieldStreamOverlay()->setVisible(on);
                 if (on) scene_->refreshCurrentFrame();
                 else    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
             });
    ACONNECT(showTrajectoryAction_.data(), &QAction::toggled,
             this, [this](bool on) {
                  if (!scene_ || !scene_->atomTrajectoryOverlay()) return;
                  scene_->atomTrajectoryOverlay()->setVisible(on);
                  if (on) {
                      updateCsaGlyph(false);
                      updateOrientationTensorGlyph();
                  } else {
                      updateCsaGlyph(true);
                      updateOrientationTensorGlyph();
                  }
                  scene_->requestRender(MoleculeScene::RenderSource::Overlay);
              });

    // Focus — a self-contained toggle at the toolbar tail (deliberately
    // de-emphasized; the Filter feature now covers the common "see one atom +
    // neighbours" need). Checked = track the focused atom, keeping it centred
    // as frames play; unchecked = manual mouse camera. No separate Free/lock
    // mode is surfaced — the button's own state is the only indicator, and it
    // also releases a lock a dashboard reveal engaged.
    tb->addSeparator();
    focusAction_ = tb->addAction(QStringLiteral("Focus"));
    focusAction_->setCheckable(true);
    focusAction_->setEnabled(false);
    focusAction_->setToolTip(QStringLiteral(
        "Track the focused atom — keep it centred as frames play. "
        "Toggle off for free mouse control. Requires a focused atom."));
    ACONNECT(focusAction_.data(), &QAction::triggered,
             this, &ReaderMainWindow::onFocusCameraTriggered);
}

void ReaderMainWindow::buildStatusBar() {
    selectionLabel_ = new QLabel(QStringLiteral("no selection"), this);
    proteinLabel_   = new QLabel(QStringLiteral("No calcset loaded"), this);
    frameLabel_     = new QLabel(QStringLiteral("frame —"), this);
    timeLabel_      = new QLabel(QStringLiteral("t=— ps"), this);

    diagnosticLabel_ = new QLabel(this);
    diagnosticLabel_->setVisible(false);
    diagnosticLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);

    // Selection summary on the LEFT; identity/frame/time pinned on the right.
    statusBar()->addWidget(selectionLabel_);
    statusBar()->addWidget(diagnosticLabel_, 1);
    statusBar()->addPermanentWidget(proteinLabel_);
    statusBar()->addPermanentWidget(frameLabel_);
    statusBar()->addPermanentWidget(timeLabel_);
}

void ReaderMainWindow::buildDocks() {
    // Atom Info dock — tabified on the LEFT alongside Selection + Strip.
    // It is constructed before load and starts with its own placeholder.
    inspectorDock_ = new QtAtomInspectorDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, inspectorDock_);

    // Dashboard strips — the dock and controller are stable chrome; run-backed
    // models are swapped in by installLoadedRun().
    dashboardStripDock_ = new DashboardStripDock(this);
    dashboardStripDock_->setContext(nullptr, nullptr);
    dashboardStripDock_->setSignalModels(nullptr, nullptr);
    dashboardStripDock_->setPanelModel(nullptr);
    dashboardStripDock_->setSelectionController(nullptr);
    dashboardStripDock_->setSelection(nullptr);
    dashboardStripDock_->setTimeViewport(nullptr);
    dashboardController_ = dashboardStripDock_->displayController();
    if (dashboardController_)
        dashboardController_->setVisualizationContext({});
    addDockWidget(Qt::LeftDockWidgetArea, dashboardStripDock_);
    tabifyDockWidget(inspectorDock_, dashboardStripDock_);

    // Measurements dock -- tabified alongside; reveals when a 2-4 atom geometry
    // is selected, and re-reads the distance/angle/dihedral live as frames play.
    // The value lives here, not as floating text on the molecule.
    measurementsDock_ = new MeasurementsDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, measurementsDock_);
    tabifyDockWidget(inspectorDock_, measurementsDock_);
    measurementsDock_->setVisible(false);

    inspectorDock_->raise();
    resizeDocks({inspectorDock_}, {360}, Qt::Horizontal);

    // Start clean — docks hidden on launch. The Inspector reveals on atom pick,
    // the Strip dock when a metric is added. QSettings restore can override this
    // for users who intentionally left a dock visible.
    inspectorDock_->setVisible(false);
    dashboardStripDock_->setVisible(false);

    ACONNECT(dashboardStripDock_, &QDockWidget::visibilityChanged,
             this, [this](bool visible) {
                 qCInfo(diagnostics::cDash).noquote()
                     << QStringLiteral("event=dock_visibility_changed visible=%1 width=%2")
                            .arg(visible ? 1 : 0)
                            .arg(dashboardDockWidth());
             });
    ACONNECT(dashboardStripDock_, &DashboardStripDock::revealRequested,
             this, [this](const model::SignalBinding& binding) {
                 if (scene_)
                     scene_->revealBinding(binding);
             });
    ACONNECT(dashboardStripDock_, &DashboardStripDock::metricPickerRequested,
             this, &ReaderMainWindow::onOpenSignalDisplays);

    // The "Panels" menu/toolbar button was removed: it exposed dock toggles that
    // greyed out with no working route. Docks reveal themselves where it makes
    // sense — the Strip dock when a metric is added, the Inspector on atom pick;
    // the Selection dock was retired (redundant with the in-scene measurements).

    if (frameSlider_) {
        ACONNECT(frameSlider_.data(), &QSlider::sliderPressed,
                 this, [this]() {
                     if (dashboardController_)
                         dashboardController_->setScrubActive(true);
                 });
        ACONNECT(frameSlider_.data(), &QSlider::sliderReleased,
                 this, [this]() {
                     if (dashboardController_)
                         dashboardController_->setScrubActive(false);
                     updateCsaGlyph(true);
                 });
    }
}

void ReaderMainWindow::onFrameChanged(int t) {
    ASSERT_THREAD(this);
    if (!loaded_ || !loaded_->conformation) {
        if (frameLabel_)
            frameLabel_->setText(QStringLiteral("frame —"));
        if (timeLabel_)
            timeLabel_->setText(QStringLiteral("t=— ps"));
        if (frameSlider_ && frameSlider_->value() != 0) {
            const QSignalBlocker block(frameSlider_);
            frameSlider_->setValue(0);
        }
        return;
    }

    const int T = static_cast<int>(loaded_->conformation->frameCount());
    const double t_ps = loaded_->conformation->timePicoseconds(
        static_cast<size_t>(std::clamp(t, 0, T - 1)));

    if (frameLabel_) {
        frameLabel_->setText(QStringLiteral("frame %1 / %2").arg(t + 1).arg(T));
    }
    if (timeLabel_) {
        timeLabel_->setText(QStringLiteral("t=%1 ps").arg(t_ps, 0, 'f', 1));
    }
    if (frameSlider_ && frameSlider_->value() != t) {
        const QSignalBlocker block(frameSlider_);
        frameSlider_->setValue(t);
    }
}

void ReaderMainWindow::updateCameraModeActions() {
    // Focus is the only camera control now. Enabled when there's a focused atom
    // to track OR the composer is already tracking (so a lock can always be
    // released); checked iff the composer is tracking anything — any non-Free
    // mode, including one a dashboard reveal engaged. Signal-blocked because we
    // connect via QAction::triggered (user-only), not toggled.
    if (!focusAction_)
        return;
    const bool hasFocus = selection_ && selection_->hasFocus();
    const bool tracking = scene_ && scene_->cameraComposer()
        && scene_->cameraComposer()->mode().kind != CameraMode::Kind::Free;
    focusAction_->setEnabled(scene_ && (hasFocus || tracking));
    const QSignalBlocker block(focusAction_);
    focusAction_->setChecked(tracking);
}

void ReaderMainWindow::updateFitModeLabel() {
    if (!transformFitAction_)
        return;
    if (!transformed_) {
        transformFitAction_->setText(QStringLiteral("Mode: no calcset"));
        transformFitAction_->setToolTip(fitModeToolTip());
        return;
    }

    using TMode = h5reader::model::TransformedConformation::Mode;
    transformFitAction_->setText(transformed_->mode() == TMode::FitSubset
        ? QStringLiteral("Mode: Locked backbone  ⇄")
        : QStringLiteral("Mode: Kabsch with give  ⇄"));
    transformFitAction_->setToolTip(fitModeToolTip());
}

bool ReaderMainWindow::trajectoryOverlayActive() const {
    return showTrajectoryAction_ && showTrajectoryAction_->isChecked();
}

void ReaderMainWindow::onFocusCameraTriggered() {
    ASSERT_THREAD(this);
    if (!scene_ || !scene_->cameraComposer()) {
        updateCameraModeActions();
        return;
    }
    const std::size_t t = playback_ ? static_cast<std::size_t>(playback_->currentFrame()) : 0u;
    // Toggle: unchecked → release to manual (Free); modeChanged re-syncs the
    // button. This also releases a lock a dashboard reveal engaged, so the one
    // toggle is the universal "stop tracking" control.
    if (focusAction_ && !focusAction_->isChecked()) {
        scene_->cameraComposer()->setMode(FreeMode(), FreePolicy(), t);
        return;
    }
    if (!selection_ || !selection_->hasFocus() || !loaded_ || !loaded_->protein) {
        updateCameraModeActions();
        return;
    }
    auto result = h5reader::app::DeriveFocusAnchor(*loaded_->protein,
                                                    selection_->focus(),
                                                    FocusAnchorKind::Plane);
    if (result.outcome != FocusAnchorOutcome::Ok) {
        qCWarning(cWindow).noquote()
            << "Focus camera: derive failed | atom=" << selection_->focus()
            << "| outcome=" << static_cast<int>(result.outcome);
        updateCameraModeActions();
        return;
    }
    scene_->cameraComposer()->setMode(result.mode, result.policy, t);
    (void)scene_->cameraComposer()->write(t);
    scene_->syncCameraClippingRange();
    scene_->requestRender(MoleculeScene::RenderSource::External);
    updateCameraModeActions();
}

void ReaderMainWindow::onTransformFitClicked() {
    ASSERT_THREAD(this);
    using TMode = h5reader::model::TransformedConformation::Mode;
    if (!transformed_ || !loaded_ || !loaded_->protein)
        return;

    if (transformed_->mode() == TMode::FitReference) {
        auto subset = h5reader::model::TransformedConformation::BackboneSubset(*loaded_->protein);
        if (subset.size() < 3) {
            qCWarning(cWindow).noquote()
                << "backbone fit requested but subset has <3 atoms; keeping all-atom fit";
            transformed_->setMode(TMode::FitReference, 0);
            return;
        }
        transformed_->setMode(TMode::FitSubset, 0, std::move(subset));
    } else {
        transformed_->setMode(TMode::FitReference, 0);
    }
}

void ReaderMainWindow::closeEvent(QCloseEvent* event) {
    ASSERT_THREAD(this);
    saveAllSettings();
    // Accept unconditionally — a failed save is logged but not allowed
    // to trap the user inside the application. aboutToQuit fires the
    // existing shutdown() chain after this returns.
    event->accept();
}

void ReaderMainWindow::saveAllSettings() {
    ASSERT_THREAD(this);
    QSettings s;   // org/app names set in main_reader.cpp
    s.setValue(QStringLiteral("viewer/window/geometry"), saveGeometry());
    s.setValue(QStringLiteral("viewer/window/state"),
               saveState(kSettingsVersion));
    s.setValue(QStringLiteral("viewer/log/mask"),
               static_cast<uint>(h5reader::diagnostics::StructuredLogger::CategoryMask()));
    // Recent files list is write-through (addToRecentFiles writes
    // immediately) so no batch write here.
    qCInfo(cWindow).noquote() << "settings saved | mask="
                              << h5reader::diagnostics::StructuredLogger::CategoryMask();
}

void ReaderMainWindow::restoreAllSettings() {
    ASSERT_THREAD(this);
    QSettings s;
    const QByteArray geom = s.value(QStringLiteral("viewer/window/geometry")).toByteArray();
    if (!geom.isEmpty())
        restoreGeometry(geom);
    const QByteArray state = s.value(QStringLiteral("viewer/window/state")).toByteArray();
    if (!state.isEmpty())
        restoreState(state, kSettingsVersion);
    const QVariant maskVar = s.value(QStringLiteral("viewer/log/mask"));
    if (maskVar.isValid()) {
        bool ok = false;
        const uint mask = maskVar.toUInt(&ok);
        if (ok)
            h5reader::diagnostics::StructuredLogger::SetCategoryMask(mask);
    }
    const QStringList recent = s.value(QStringLiteral("viewer/recent/files")).toStringList();
    rebuildRecentFilesMenu(recent);
}

void ReaderMainWindow::addToRecentFiles(const QString& path) {
    ASSERT_THREAD(this);
    if (path.isEmpty()) return;
    QSettings s;
    QStringList recent = s.value(QStringLiteral("viewer/recent/files")).toStringList();
    recent.removeAll(path);
    recent.prepend(path);
    while (recent.size() > kMaxRecentFiles)
        recent.removeLast();
    s.setValue(QStringLiteral("viewer/recent/files"), recent);
    rebuildRecentFilesMenu(recent);
}

void ReaderMainWindow::rebuildRecentFilesMenu(const QStringList& paths) {
    ASSERT_THREAD(this);
    if (!recentMenu_) return;
    recentMenu_->clear();
    if (paths.isEmpty()) {
        QAction* empty = recentMenu_->addAction(QStringLiteral("(none)"));
        empty->setEnabled(false);
        return;
    }
    for (const QString& path : paths) {
        QAction* a = recentMenu_->addAction(path);
        ACONNECT(a, &QAction::triggered, this, [this, path]() {
            openRecentPath(path);
        });
    }
}

void ReaderMainWindow::openRecentPath(const QString& path) {
    ASSERT_THREAD(this);
    if (!loadRunPath(path)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onPlayPauseClicked() {
    ASSERT_THREAD(this);
    if (playback_) playback_->togglePlayPause();
}

void ReaderMainWindow::onOpenFile() {
    ASSERT_THREAD(this);
    // Pick a .LGS calcset file directly with the mouse. CalcsetManifest::Load
    // (via ResolveLgsPath) accepts a .LGS file path, so load it into this
    // window.
    const QString file = QFileDialog::getOpenFileName(
        this, QStringLiteral("Open a calcset (.LGS) file"),
        qEnvironmentVariable("H5READER_OPEN_DIR"),
        QStringLiteral("Calcset manifest (*.lgs *.LGS);;All files (*)"));
    if (file.isEmpty())
        return;
    if (!loadRunPath(file)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onOpenDirectory() {
    ASSERT_THREAD(this);
    const QString dir = QFileDialog::getExistingDirectory(
        this, QStringLiteral("Open a run directory (trajectory or single pose)"));
    if (dir.isEmpty())
        return;
    if (!loadRunPath(dir)) {
        QMessageBox::critical(this,
                              QStringLiteral("Open calcset failed"),
                              lastLoadError());
    }
}

void ReaderMainWindow::onGoToAtomTriggered() {
    ASSERT_THREAD(this);
    if (!loaded_ || !loaded_->protein || !loaded_->conformation || !selection_) {
        QMessageBox::information(this,
                                 QStringLiteral("Go to atom"),
                                 QStringLiteral("Load a calcset before jumping to an atom."));
        return;
    }

    const model::QtProtein& protein = *loaded_->protein;
    if (protein.residueCount() == 0 || protein.atomCount() == 0) {
        QMessageBox::information(this,
                                 QStringLiteral("Go to atom"),
                                 QStringLiteral("This calcset has no residues or atoms to select."));
        return;
    }

    int minResidueNumber = protein.residue(0).address.residueNumber;
    int maxResidueNumber = minResidueNumber;
    for (std::size_t i = 1; i < protein.residueCount(); ++i) {
        const int n = protein.residue(i).address.residueNumber;
        minResidueNumber = std::min(minResidueNumber, n);
        maxResidueNumber = std::max(maxResidueNumber, n);
    }

    int initialResidueNumber = minResidueNumber;
    std::optional<std::size_t> preferredAtom;
    if (selection_->hasFocus()) {
        const std::size_t atom = selection_->focus();
        if (atom < protein.atomCount()) {
            const model::QtAtom& picked = protein.atom(atom);
            if (picked.residueIndex >= 0
                && static_cast<std::size_t>(picked.residueIndex) < protein.residueCount()) {
                initialResidueNumber =
                    protein.residue(static_cast<std::size_t>(picked.residueIndex))
                        .address.residueNumber;
                preferredAtom = atom;
            }
        }
    }

    QDialog dialog(this);
    dialog.setWindowTitle(QStringLiteral("Go to atom"));
    dialog.setModal(true);

    auto* form = new QFormLayout(&dialog);

    auto* residueSpin = new QSpinBox(&dialog);
    residueSpin->setRange(minResidueNumber, maxResidueNumber);
    residueSpin->setValue(initialResidueNumber);
    residueSpin->setToolTip(QStringLiteral("Visible residue sequence number."));
    form->addRow(QStringLiteral("Residue"), residueSpin);

    auto* atomCombo = new QComboBox(&dialog);
    atomCombo->setMinimumContentsLength(18);
    atomCombo->setToolTip(QStringLiteral("Atoms in matching residues."));
    form->addRow(QStringLiteral("Atom"), atomCombo);

    auto* frameSpin = new QSpinBox(&dialog);
    const int frameCount = static_cast<int>(loaded_->conformation->frameCount());
    frameSpin->setRange(0, std::max(0, frameCount - 1));
    frameSpin->setValue(playback_ ? playback_->currentFrame() : 0);
    frameSpin->setToolTip(QStringLiteral("Zero-based frame index."));
    form->addRow(QStringLiteral("Frame"), frameSpin);

    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
                                         &dialog);
    form->addRow(buttons);

    auto residueLabel = [&protein](std::size_t residueIndex) {
        const model::QtResidue& residue = protein.residue(residueIndex);
        QString label = QStringLiteral("%1%2")
            .arg(protein.residueLabel(residueIndex,
                                      model::NamingConvention::Amber,
                                      model::NamingSource::Verbatim))
            .arg(residue.address.residueNumber);
        if (!residue.address.chainId.isEmpty())
            label += QStringLiteral(" chain %1").arg(residue.address.chainId);
        if (!residue.address.insertionCode.isEmpty())
            label += QStringLiteral(" ins %1").arg(residue.address.insertionCode);
        return label;
    };

    auto rebuildAtomChoices = [&]() {
        const QSignalBlocker block(atomCombo);
        atomCombo->clear();
        const int residueNumber = residueSpin->value();
        int preferredIndex = -1;
        for (std::size_t r = 0; r < protein.residueCount(); ++r) {
            const model::QtResidue& residue = protein.residue(r);
            if (residue.address.residueNumber != residueNumber)
                continue;
            const QString residueText = residueLabel(r);
            for (int32_t atomIndex : residue.atomIndices) {
                if (atomIndex < 0)
                    continue;
                const std::size_t atom = static_cast<std::size_t>(atomIndex);
                if (atom >= protein.atomCount())
                    continue;
                const QString atomName = protein.atomNames(atom).amber;
                atomCombo->addItem(QStringLiteral("%1:%2  #%3")
                                       .arg(residueText)
                                       .arg(atomName)
                                       .arg(atom),
                                   QVariant::fromValue<qulonglong>(
                                       static_cast<qulonglong>(atom)));
                if (preferredAtom && *preferredAtom == atom)
                    preferredIndex = atomCombo->count() - 1;
            }
        }
        if (preferredIndex >= 0)
            atomCombo->setCurrentIndex(preferredIndex);
        if (QPushButton* ok = buttons->button(QDialogButtonBox::Ok))
            ok->setEnabled(atomCombo->count() > 0);
    };

    ACONNECT(residueSpin, qOverload<int>(&QSpinBox::valueChanged),
             &dialog,    [rebuildAtomChoices](int) mutable { rebuildAtomChoices(); });
    ACONNECT(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    ACONNECT(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);

    rebuildAtomChoices();
    if (dialog.exec() != QDialog::Accepted)
        return;

    bool ok = false;
    const qulonglong atomValue = atomCombo->currentData().toULongLong(&ok);
    if (!ok || atomValue >= static_cast<qulonglong>(protein.atomCount()))
        return;

    if (playback_)
        playback_->setFrame(frameSpin->value());
    selection_->bulkSet({static_cast<std::size_t>(atomValue)});
    if (scene_)
        scene_->clearReveal();
    statusBar()->showMessage(QStringLiteral("Jumped to %1 at frame %2")
                                 .arg(atomCombo->currentText())
                                 .arg(frameSpin->value()));
}

void ReaderMainWindow::onOpenSignalDisplays() {
    ASSERT_THREAD(this);
    openSignalDisplayPicker();
}

}  // namespace h5reader::app
