#include "MainWindow.h"
#include "BackboneRibbonOverlay.h"
#include "RingCurrentOverlay.h"
#include "PeptideBondOverlay.h"
#include "TensorGlyph.h"
#include "EllipsoidGlyph.h"
#include "FieldOverlay.h"
#include "IsosurfaceOverlay.h"
#include "ButterflyOverlay.h"
#include "FieldGridOverlay.h"
#include "ComputeWorker.h"

// Library headers — the viewer reads the object model directly
#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Atom.h"
#include "Residue.h"
#include "Ring.h"
#include "Bond.h"
#include "DsspResult.h"
#include "MopacResult.h"
#include "MopacCoulombResult.h"
#include "MopacMcConnellResult.h"
#include "GeometryChoice.h"
#include "ConformationResult.h"

#include <filesystem>

#include <QMenuBar>
#include <QToolBar>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QSlider>
#include <QComboBox>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QScrollArea>
#include <QStatusBar>
#include <QThread>
#include <QElapsedTimer>
#include <QTimer>
#include <QProgressDialog>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNew.h>
#include <vtkCamera.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QMouseEvent>
#include <QApplication>
#include <cmath>
#include <QDir>
#include <QFileInfo>
#include <QPlainTextEdit>
#include <QUdpSocket>
#include <QNetworkDatagram>
#include <QFont>
#include <vtkSphereSource.h>
#include <vtkCellPicker.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLookupTable.h>
#include <vtkDoubleArray.h>

using namespace nmr;

extern "C" void udp_log(const char* fmt, ...);

// GAP: Types.h has SymbolForElement, ThreeLetterCodeForAminoAcid, etc. but no
// string conversion for AtomRole or BondCategory. These should be added to
// the library alongside the existing enum→string functions in Types.h.
static const char* NameForAtomRole(AtomRole r) {
    switch (r) {
        case AtomRole::BackboneN:   return "BackboneN";
        case AtomRole::BackboneCA:  return "BackboneCA";
        case AtomRole::BackboneC:   return "BackboneC";
        case AtomRole::BackboneO:   return "BackboneO";
        case AtomRole::SidechainC:  return "SidechainC";
        case AtomRole::SidechainN:  return "SidechainN";
        case AtomRole::SidechainO:  return "SidechainO";
        case AtomRole::SidechainS:  return "SidechainS";
        case AtomRole::AromaticC:   return "AromaticC";
        case AtomRole::AromaticN:   return "AromaticN";
        case AtomRole::AmideH:      return "AmideH";
        case AtomRole::AlphaH:      return "AlphaH";
        case AtomRole::MethylH:     return "MethylH";
        case AtomRole::AromaticH:   return "AromaticH";
        case AtomRole::HydroxylH:   return "HydroxylH";
        case AtomRole::OtherH:      return "OtherH";
        case AtomRole::Unknown:     return "Unknown";
    }
    return "?";
}

static const char* NameForBondCategory(BondCategory c) {
    switch (c) {
        case BondCategory::PeptideCO:      return "PeptideCO";
        case BondCategory::PeptideCN:      return "PeptideCN";
        case BondCategory::BackboneOther:  return "BackboneOther";
        case BondCategory::SidechainCO:    return "SidechainCO";
        case BondCategory::Aromatic:       return "Aromatic";
        case BondCategory::Disulfide:      return "Disulfide";
        case BondCategory::SidechainOther: return "SidechainOther";
        case BondCategory::Unknown:        return "Unknown";
    }
    return "?";
}

MainWindow::MainWindow(const QString& initialDir, QWidget* parent)
    : QMainWindow(parent)
    , initialDir_(initialDir)
    , ringOverlay_(nullptr)
    , peptideBondOverlay_(nullptr)
    , tensorGlyph_(nullptr)
    , ellipsoidGlyph_(nullptr)
    , fieldOverlay_(nullptr)
    , isosurfaceOverlay_(nullptr)
    , isosurfaceOverlayPass_(nullptr)
    , butterflyOverlay_(nullptr)
    , fieldGridOverlay_(nullptr)
{
    udp_log("[lifecycle] MainWindow constructor entered\n");
    setWindowTitle("NMR Shielding Tensor Viewer");
    resize(1400, 900);

    setupUI();
    udp_log("[lifecycle] setupUI done\n");
    setupMenuBar();

    // Double-click on the 3D view picks an atom
    vtkWidget_->installEventFilter(this);

    udp_log("[lifecycle] MainWindow constructor done\n");
}

MainWindow::~MainWindow() {
    // Destructor may fire after QApplication is gone — minimal work only.
    // Heavy cleanup happens in shutdown().
}

void MainWindow::shutdown() {
    udp_log("[lifecycle] shutdown() entered\n");

    // 1. Stop all timers
    const auto timers = findChildren<QTimer*>();
    for (auto* timer : timers)
        timer->stop();

    // 2. Stop async workers
    cancelCompute();

    // 3. Finalize VTK before Qt tears down the GL context
    if (renderWindow_)
        renderWindow_->Finalize();

    udp_log("[lifecycle] shutdown() done\n");
}

void MainWindow::setupMenuBar() {
    auto* fileMenu = menuBar()->addMenu("&File");

    auto* screenshotAct = fileMenu->addAction("&Save Screenshot...");
    connect(screenshotAct, &QAction::triggered, this, &MainWindow::saveScreenshot);

    exportFeaturesAct_ = fileMenu->addAction("&Export Features...");
    exportFeaturesAct_->setEnabled(false);
    connect(exportFeaturesAct_, &QAction::triggered, this, &MainWindow::exportFeatures);

    fileMenu->addSeparator();
    auto* quitAct = fileMenu->addAction("&Quit");
    connect(quitAct, &QAction::triggered, this, &QWidget::close);
}

void MainWindow::exportFeatures() {
    if (!protein_) return;

    QString dir = QFileDialog::getExistingDirectory(this, "Export Features",
        QDir::currentPath());
    if (dir.isEmpty()) return;

    std::string outDir = dir.toStdString();
    int totalArrays = 0;

    if (pendingSpec_.mode == nmr::JobMode::Fleet) {
        for (size_t i = 0; i < protein_->ConformationCount(); ++i) {
            auto& conf = protein_->ConformationAt(i);
            std::string frameDir = outDir + "/frame_" +
                std::to_string(i + 1);
            std::filesystem::create_directories(frameDir);
            totalArrays += nmr::ConformationResult::WriteAllFeatures(
                conf, frameDir);
        }
    } else {
        std::filesystem::create_directories(outDir);
        auto& conf = protein_->Conformation();
        totalArrays = nmr::ConformationResult::WriteAllFeatures(
            conf, outDir);
    }

    statusLabel_->setText(QString("Exported %1 arrays to %2")
        .arg(totalArrays).arg(dir));
}

void MainWindow::setupUI() {
    // Central VTK widget
    vtkWidget_ = new QVTKOpenGLNativeWidget(this);
    setCentralWidget(vtkWidget_);

    // Match pdbviewer-v1 init order: setRenderWindow -> AddRenderer -> peeling -> alpha
    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWidget_->setRenderWindow(renderWindow_);

    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderWindow_->AddRenderer(renderer_);

    renderer_->SetBackground(1.0, 1.0, 1.0);
    renderer_->SetUseFXAA(true);
    renderer_->SetUseDepthPeeling(0);
    renderWindow_->SetAlphaBitPlanes(1);
    renderWindow_->SetMultiSamples(0);  // MSAA off — conflicts with translucency; FXAA handles AA

    vtkNew<vtkInteractorStyleTrackballCamera> style;
    renderWindow_->GetInteractor()->SetInteractorStyle(style);

    // Sidebar dock
    auto* dock = new QDockWidget(this);
    dock->setTitleBarWidget(new QWidget());
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    dock->setFeatures(QDockWidget::DockWidgetMovable);

    auto* scrollArea = new QScrollArea(dock);
    scrollArea->setWidgetResizable(true);
    scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scrollArea->setFrameShape(QFrame::NoFrame);
    auto* sidebarWidget = new QWidget();
    auto* sidebar = new QVBoxLayout(sidebarWidget);
    sidebar->setContentsMargins(2, 2, 2, 2);
    sidebar->setSpacing(0);

    auto makeSection = [&sidebar](const QString& title, bool startOpen) {
        auto* header = new QPushButton(startOpen ? QString("- %1").arg(title)
                                                 : QString("+ %1").arg(title));
        header->setFlat(true);
        header->setStyleSheet(
            "QPushButton { text-align: left; font-weight: bold; padding: 4px 6px;"
            "  background: palette(mid); border: none; }"
            "QPushButton:hover { background: palette(dark); }");
        auto* content = new QWidget();
        content->setVisible(startOpen);
        auto* contentLayout = new QVBoxLayout(content);
        contentLayout->setContentsMargins(6, 4, 6, 4);
        contentLayout->setSpacing(3);

        QObject::connect(header, &QPushButton::clicked, [header, content, title]() {
            bool show = !content->isVisible();
            content->setVisible(show);
            header->setText(show ? QString("- %1").arg(title)
                                 : QString("+ %1").arg(title));
        });

        sidebar->addWidget(header);
        sidebar->addWidget(content);
        return contentLayout;
    };

    // --- Rendering ---
    {
        auto* lay = makeSection("Rendering", true);

        renderModeCombo_ = new QComboBox();
        renderModeCombo_->addItems({"Ball & Stick", "Liquorice"});
        connect(renderModeCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &MainWindow::onRenderModeChanged);
        lay->addWidget(renderModeCombo_);

        showRibbonCheck_ = new QCheckBox("Backbone ribbon");
        showRibbonCheck_->setChecked(false);
        connect(showRibbonCheck_, &QCheckBox::toggled, this, [this](bool checked) {
            if (ribbonOverlay_) { ribbonOverlay_->setVisible(checked); renderWindow_->Render(); }
        });
        lay->addWidget(showRibbonCheck_);
    }

    // --- Ring Currents ---
    {
        auto* lay = makeSection("Ring Currents", true);
        showRingsCheck_ = new QCheckBox("Ring outlines");
        showRingsCheck_->setChecked(true);
        connect(showRingsCheck_, &QCheckBox::toggled, this, &MainWindow::onShowRingsToggled);
        lay->addWidget(showRingsCheck_);

        showFieldGridCheck_ = new QCheckBox("Shielded lobes (blue)");
        showFieldGridCheck_->setChecked(true);
        connect(showFieldGridCheck_, &QCheckBox::toggled, this, [this](bool checked) {
            if (fieldGridOverlay_) { fieldGridOverlay_->setShieldedVisible(checked); renderWindow_->Render(); }
        });
        lay->addWidget(showFieldGridCheck_);

        showDeshieldedCheck_ = new QCheckBox("Deshielded lobes (coral)");
        showDeshieldedCheck_->setChecked(true);
        connect(showDeshieldedCheck_, &QCheckBox::toggled, this, [this](bool checked) {
            if (fieldGridOverlay_) { fieldGridOverlay_->setDeshieldedVisible(checked); renderWindow_->Render(); }
        });
        lay->addWidget(showDeshieldedCheck_);

        showButterflyCheck_ = new QCheckBox("B-field streamlines");
        showButterflyCheck_->setChecked(false);
        connect(showButterflyCheck_, &QCheckBox::toggled, this, &MainWindow::onShowButterflyToggled);
        lay->addWidget(showButterflyCheck_);

        isoThresholdLabel_ = new QLabel("Iso threshold: 0.10 ppm");
        lay->addWidget(isoThresholdLabel_);
        isoThresholdSlider_ = new QSlider(Qt::Horizontal);
        isoThresholdSlider_->setRange(1, 100);  // 0.01 to 1.00 ppm
        isoThresholdSlider_->setValue(10);       // default 0.10 ppm
        connect(isoThresholdSlider_, &QSlider::valueChanged, this, [this](int value) {
            isoThresholdLabel_->setText(QString("Iso threshold: %1 ppm").arg(value / 100.0, 0, 'f', 2));
        });
        connect(isoThresholdSlider_, &QSlider::sliderReleased, this, &MainWindow::onIsoThresholdChanged);
        lay->addWidget(isoThresholdSlider_);

        lay->addWidget(new QLabel("Iso opacity:"));
        currentScaleSlider_ = new QSlider(Qt::Horizontal);
        currentScaleSlider_->setRange(5, 80);  // 0.05 to 0.80
        currentScaleSlider_->setValue(35);      // default 0.35
        connect(currentScaleSlider_, &QSlider::sliderReleased, this, [this]() {
            if (fieldGridOverlay_ && !fieldGrids_.empty()) {
                double opacity = currentScaleSlider_->value() / 100.0;
                double threshold = isoThresholdSlider_->value() / 100.0;
                fieldGridOverlay_->setData(fieldGrids_, threshold, opacity, 0);
                renderWindow_->Render();
            }
        });
        lay->addWidget(currentScaleSlider_);
    }

    // --- Bond Anisotropy ---
    {
        auto* lay = makeSection("Bond Anisotropy", true);
        showPeptideBondsCheck_ = new QCheckBox("Peptide bond tubes");
        showPeptideBondsCheck_->setChecked(false);
        connect(showPeptideBondsCheck_, &QCheckBox::toggled, this, &MainWindow::onShowPeptideBondsToggled);
        lay->addWidget(showPeptideBondsCheck_);
    }

    // --- MOPAC Electronic ---
    {
        auto* lay = makeSection("MOPAC Electronic", true);
        showBondOrderCheck_ = new QCheckBox("Bond order coloring");
        showBondOrderCheck_->setChecked(false);
        connect(showBondOrderCheck_, &QCheckBox::toggled, this, &MainWindow::onShowBondOrderToggled);
        lay->addWidget(showBondOrderCheck_);
    }

    // --- Tensor Display (global settings for future per-calculator glyphs) ---
    {
        auto* lay = makeSection("Tensor Display", false);
        lay->addWidget(new QLabel("Scale:"));
        glyphScaleSlider_ = new QSlider(Qt::Horizontal);
        glyphScaleSlider_->setRange(1, 200);
        glyphScaleSlider_->setValue(50);
        lay->addWidget(glyphScaleSlider_);
        lay->addWidget(new QLabel("Opacity:"));
        opacitySlider_ = new QSlider(Qt::Horizontal);
        opacitySlider_->setRange(0, 100);
        opacitySlider_->setValue(70);
        lay->addWidget(opacitySlider_);
    }

    sidebar->addStretch();
    sidebarWidget->setLayout(sidebar);
    scrollArea->setWidget(sidebarWidget);
    dock->setWidget(scrollArea);
    addDockWidget(Qt::RightDockWidgetArea, dock);

    // Status bar
    statusLabel_ = new QLabel("Ready");
    statusBar()->addWidget(statusLabel_);

    // Atom info panel — double-click an atom to inspect its full object model
    atomInfoDock_ = new QDockWidget("Atom Inspector", this);
    atomInfoDock_->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);
    atomInfoDock_->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    atomInfoTree_ = new QTreeWidget();
    atomInfoTree_->setHeaderLabels({"Property", "Value"});
    atomInfoTree_->setColumnWidth(0, 220);
    atomInfoTree_->setAlternatingRowColors(true);
    atomInfoTree_->setIndentation(16);
    atomInfoDock_->setWidget(atomInfoTree_);
    addDockWidget(Qt::BottomDockWidgetArea, atomInfoDock_);

    // Bond browser panel — tree view of all bonds grouped by category
    bondInfoDock_ = new QDockWidget("Bonds", this);
    bondInfoDock_->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);
    bondInfoDock_->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    bondInfoTree_ = new QTreeWidget();
    bondInfoTree_->setHeaderLabels({"Property", "Value"});
    bondInfoTree_->setColumnWidth(0, 220);
    bondInfoTree_->setAlternatingRowColors(true);
    bondInfoTree_->setIndentation(16);
    bondInfoDock_->setWidget(bondInfoTree_);
    tabifyDockWidget(atomInfoDock_, bondInfoDock_);

    // GeometryChoice panel — calculator decisions for picked atom
    gcDock_ = new QDockWidget("Geometry Choices", this);
    gcDock_->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);
    gcDock_->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    gcTree_ = new QTreeWidget();
    gcTree_->setHeaderLabels({"Decision", "Detail"});
    gcTree_->setColumnWidth(0, 250);
    gcTree_->setAlternatingRowColors(true);
    gcTree_->setIndentation(16);
    gcDock_->setWidget(gcTree_);
    tabifyDockWidget(atomInfoDock_, gcDock_);

    // Operations log panel — listens on UDP for library log stream
    logDock_ = new QDockWidget("Operations Log", this);
    logDock_->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);
    logDock_->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    logText_ = new QPlainTextEdit();
    logText_->setReadOnly(true);
    logText_->setMaximumBlockCount(5000);
    logText_->setFont(QFont("Monospace", 9));
    logDock_->setWidget(logText_);
    tabifyDockWidget(atomInfoDock_, logDock_);

    // Bind UDP listener on port 9998 (same port OperationLog sends to)
    logSocket_ = new QUdpSocket(this);
    if (logSocket_->bind(QHostAddress::LocalHost, 9998, QUdpSocket::ShareAddress)) {
        connect(logSocket_, &QUdpSocket::readyRead, this, &MainWindow::onLogDatagramReady);
    } else {
        // Try next ports if 9998 is taken
        for (quint16 p = 9999; p < 10008; ++p) {
            if (logSocket_->bind(QHostAddress::LocalHost, p, QUdpSocket::ShareAddress)) {
                connect(logSocket_, &QUdpSocket::readyRead, this, &MainWindow::onLogDatagramReady);
                break;
            }
        }
    }

    // Start with the log tab visible — it shows computation progress
    logDock_->raise();
}

void MainWindow::loadFromJobSpec(const nmr::JobSpec& spec) {
    pendingSpec_ = spec;

    // Derive a display name from the spec
    switch (spec.mode) {
    case nmr::JobMode::Pdb:
    case nmr::JobMode::ProtonatedPdb:
        currentProteinId_ = QFileInfo(QString::fromStdString(spec.pdb_path)).baseName().toStdString();
        break;
    case nmr::JobMode::Orca:
        currentProteinId_ = QFileInfo(QString::fromStdString(spec.orca_files.xyz_path)).baseName().toStdString();
        break;
    case nmr::JobMode::Mutant:
        currentProteinId_ = QFileInfo(QString::fromStdString(spec.wt_files.xyz_path)).baseName().toStdString();
        break;
    case nmr::JobMode::Fleet:
        currentProteinId_ = QFileInfo(QString::fromStdString(spec.fleet_paths.sampled_poses_dir)).fileName().toStdString();
        break;
    case nmr::JobMode::None:
        return;
    }

    loadMolecule();
}

void MainWindow::loadMolecule() {
    udp_log("[lifecycle] loadMolecule entered: %s\n", currentProteinId_.c_str());
    QElapsedTimer timer;
    timer.start();

    cancelCompute();
    udp_log("[lifecycle] cancelCompute done\n");

    renderer_->RemoveAllViewProps();

    // Clean up old overlays (will be rebuilt in onComputeFinished)
    if (ribbonOverlay_) { delete ribbonOverlay_; ribbonOverlay_ = nullptr; }
    if (ringOverlay_) { delete ringOverlay_; ringOverlay_ = nullptr; }
    if (tensorGlyph_) { delete tensorGlyph_; tensorGlyph_ = nullptr; }
    if (ellipsoidGlyph_) { delete ellipsoidGlyph_; ellipsoidGlyph_ = nullptr; }
    if (peptideBondOverlay_) { delete peptideBondOverlay_; peptideBondOverlay_ = nullptr; }
    if (fieldOverlay_) { delete fieldOverlay_; fieldOverlay_ = nullptr; }
    if (isosurfaceOverlay_) { delete isosurfaceOverlay_; isosurfaceOverlay_ = nullptr; }
    if (isosurfaceOverlayPass_) { delete isosurfaceOverlayPass_; isosurfaceOverlayPass_ = nullptr; }
    if (butterflyOverlay_) { delete butterflyOverlay_; butterflyOverlay_ = nullptr; }

    tensorGlyph_ = new TensorGlyph(renderer_);
    ellipsoidGlyph_ = new EllipsoidGlyph(renderer_);
    fieldOverlay_ = new FieldOverlay(renderer_);
    isosurfaceOverlay_ = new IsosurfaceOverlay(renderer_);
    isosurfaceOverlayPass_ = new IsosurfaceOverlay(renderer_);

    udp_log("[diag] VTK setup: %lld ms\n", (long long)timer.elapsed());
    timer.restart();

    renderer_->ResetCamera();
    udp_log("[diag] ResetCamera done, calling Render()\n");
    renderWindow_->Render();
    udp_log("[diag] First render: %lld ms\n", (long long)timer.elapsed());

    statusLabel_->setText(QString("Loaded %1 — computing features...")
        .arg(QString::fromStdString(currentProteinId_)));

    // Clear stale data
    protein_.reset();
    fieldGrids_.clear();
    butterflyFields_.clear();

    QTimer::singleShot(0, this, &MainWindow::startCompute);
}

void MainWindow::startCompute() {
    udp_log("[lifecycle] startCompute entered\n");
    workerThread_ = new QThread;
    worker_ = new ComputeWorker;
    worker_->moveToThread(workerThread_);

    connect(this, &MainWindow::computeRequested, worker_, &ComputeWorker::computeAll);
    connect(worker_, &ComputeWorker::progress, this, &MainWindow::onComputeProgress);
    connect(worker_, &ComputeWorker::finished, this, &MainWindow::onComputeFinished);

    progressDialog_ = new QProgressDialog("Computing features...", "Cancel", 0, 100, this);
    progressDialog_->setMinimumDuration(0);
    progressDialog_->setValue(0);
    connect(progressDialog_, &QProgressDialog::canceled, this, &MainWindow::cancelCompute);

    workerThread_->start();
    emit computeRequested(pendingSpec_);
}

void MainWindow::cancelCompute() {
    if (worker_)
        worker_->cancel();
    if (workerThread_ && workerThread_->isRunning()) {
        workerThread_->quit();
        workerThread_->wait();
    }
    delete worker_;
    worker_ = nullptr;
    delete workerThread_;
    workerThread_ = nullptr;
    if (progressDialog_) {
        progressDialog_->close();
        progressDialog_->deleteLater();
        progressDialog_ = nullptr;
    }
}

void MainWindow::onComputeProgress(int current, int total, QString phase) {
    if (!progressDialog_) return;
    progressDialog_->setMaximum(total);
    progressDialog_->setValue(std::min(current, total - 1));
    progressDialog_->setLabelText(phase);
}

void MainWindow::onComputeFinished(ComputeResult result) {
    udp_log("[diag] onComputeFinished: protein=%s\n",
            result.protein ? "valid" : "null");

    // Store the protein — the library object model, fully const after Run
    protein_ = std::move(result.protein);
    fieldGrids_ = std::move(result.fieldGrids);
    butterflyFields_ = std::move(result.butterflyFields);

    // Stop worker thread
    if (workerThread_ && workerThread_->isRunning()) {
        workerThread_->quit();
        workerThread_->wait();
    }
    delete worker_;
    worker_ = nullptr;
    delete workerThread_;
    workerThread_ = nullptr;
    udp_log("[diag] worker cleaned up\n");

    if (progressDialog_) {
        udp_log("[diag] closing progress dialog\n");
        progressDialog_->disconnect();
        progressDialog_->close();
        progressDialog_->deleteLater();
        progressDialog_ = nullptr;
    }

    if (!protein_) {
        statusLabel_->setText("Load failed");
        QMessageBox::critical(this, "Load Failed",
            QString("Failed to load protein.\n\n"
                    "Check the command-line arguments and file paths."));
        QApplication::exit(EXIT_FAILURE);
        return;
    }

    exportFeaturesAct_->setEnabled(true);

    const auto& protein = *protein_;
    const auto& conf = protein.Conformation();

    udp_log("[diag] protein: %zu atoms, %zu bonds, %zu rings\n",
            protein.AtomCount(), protein.BondCount(), protein.RingCount());

    // Build vtkMolecule FIRST — RemoveAllViewProps clears the renderer,
    // so this must happen before any overlay actors are added.
    {
        molecule_ = vtkSmartPointer<vtkMolecule>::New();

        for (size_t i = 0; i < protein.AtomCount(); ++i) {
            unsigned short anum = static_cast<unsigned short>(
                AtomicNumberForElement(protein.AtomAt(i).element));
            const Vec3& pos = conf.AtomAt(i).Position();
            molecule_->AppendAtom(anum, pos.x(), pos.y(), pos.z());
        }

        for (size_t i = 0; i < protein.BondCount(); ++i) {
            const auto& bond = protein.BondAt(i);
            unsigned short vtkOrder = 1;
            switch (bond.order) {
                case BondOrder::Double:   vtkOrder = 2; break;
                case BondOrder::Triple:   vtkOrder = 3; break;
                case BondOrder::Aromatic: vtkOrder = 2; break;  // show as double
                case BondOrder::Peptide:  vtkOrder = 1; break;  // partial double, show single
                default: break;
            }
            molecule_->AppendBond(bond.atom_index_a, bond.atom_index_b, vtkOrder);
        }

        molMapper_ = vtkSmartPointer<vtkOpenGLMoleculeMapper>::New();
        molMapper_->SetInputData(molecule_);
        molMapper_->UseBallAndStickSettings();

        renderer_->RemoveAllViewProps();
        molActor_ = vtkSmartPointer<vtkActor>::New();
        molActor_->SetMapper(molMapper_);
        renderer_->AddActor(molActor_);

        udp_log("[diag] vtkMolecule built: %d atoms, %d bonds\n",
                (int)molecule_->GetNumberOfAtoms(), (int)molecule_->GetNumberOfBonds());
    }

    // Overlays — added AFTER the molecule actor so they render on top.
    if (ribbonOverlay_) { delete ribbonOverlay_; ribbonOverlay_ = nullptr; }
    ribbonOverlay_ = new BackboneRibbonOverlay(renderer_, protein, conf);
    ribbonOverlay_->setVisible(showRibbonCheck_->isChecked());

    if (ringOverlay_) { delete ringOverlay_; ringOverlay_ = nullptr; }
    ringOverlay_ = new RingCurrentOverlay(renderer_, protein, conf);

    if (peptideBondOverlay_) { delete peptideBondOverlay_; peptideBondOverlay_ = nullptr; }
    peptideBondOverlay_ = new PeptideBondOverlay(renderer_, protein, conf);
    peptideBondOverlay_->setVisible(showPeptideBondsCheck_->isChecked());

    // Butterfly overlay (B-field streamlines)
    udp_log("[diag] butterfly fields: %zu grids\n", butterflyFields_.size());
    if (!butterflyFields_.empty()) {
        butterflyOverlay_ = new ButterflyOverlay(renderer_);
        butterflyOverlay_->setData(butterflyFields_);
        butterflyOverlay_->setVisible(showButterflyCheck_->isChecked());
    }

    // Field grid overlay — T0 shielding isosurface around each ring
    if (!fieldGrids_.empty()) {
        udp_log("[diag] %zu field grids, T0 range shown in status\n", fieldGrids_.size());
        if (fieldGridOverlay_) { delete fieldGridOverlay_; fieldGridOverlay_ = nullptr; }
        fieldGridOverlay_ = new FieldGridOverlay(renderer_);
        double threshold = isoThresholdSlider_->value() / 100.0;
        fieldGridOverlay_->setData(fieldGrids_, threshold, 0.35, 0);
        fieldGridOverlay_->setVisible(showFieldGridCheck_->isChecked());
    }

    // Bond order color overlay — tubes colored by MOPAC Wiberg order
    if (bondOrderActor_) {
        renderer_->RemoveActor(bondOrderActor_);
        bondOrderActor_ = nullptr;
    }
    if (conf.HasResult<MopacResult>()) {
        const auto& mopac = conf.Result<MopacResult>();
        const auto& topoOrders = mopac.TopologyBondOrders();

        vtkNew<vtkPoints> pts;
        vtkNew<vtkCellArray> lines;
        vtkNew<vtkDoubleArray> scalars;
        scalars->SetName("WibergOrder");
        scalars->SetNumberOfComponents(1);

        for (size_t i = 0; i < protein.BondCount(); ++i) {
            double bo = (i < topoOrders.size()) ? topoOrders[i] : 0.0;
            if (bo < 0.01) continue;  // skip electronically insignificant bonds

            const Bond& bond = protein.BondAt(i);
            Vec3 posA = conf.AtomAt(bond.atom_index_a).Position();
            Vec3 posB = conf.AtomAt(bond.atom_index_b).Position();
            vtkIdType id0 = pts->InsertNextPoint(posA.data());
            vtkIdType id1 = pts->InsertNextPoint(posB.data());
            vtkNew<vtkLine> ln;
            ln->GetPointIds()->SetId(0, id0);
            ln->GetPointIds()->SetId(1, id1);
            lines->InsertNextCell(ln);
            scalars->InsertNextValue(bo);
        }

        if (lines->GetNumberOfCells() > 0) {
            vtkNew<vtkPolyData> pd;
            pd->SetPoints(pts);
            pd->SetLines(lines);
            pd->GetCellData()->SetScalars(scalars);

            vtkNew<vtkTubeFilter> tube;
            tube->SetInputData(pd);
            tube->SetRadius(0.08);
            tube->SetNumberOfSides(8);

            // Cool-warm colormap: 0.8=blue(single) → 1.5=white(aromatic) → 2.0=red(double)
            vtkNew<vtkLookupTable> lut;
            lut->SetNumberOfTableValues(256);
            lut->SetRange(0.8, 2.1);
            lut->SetHueRange(0.65, 0.0);    // blue → red
            lut->SetSaturationRange(0.7, 0.9);
            lut->SetValueRange(0.8, 0.95);
            lut->Build();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(tube->GetOutputPort());
            mapper->SetScalarModeToUseCellData();
            mapper->SetLookupTable(lut);
            mapper->SetScalarRange(0.8, 2.1);

            bondOrderActor_ = vtkSmartPointer<vtkActor>::New();
            bondOrderActor_->SetMapper(mapper);
            bondOrderActor_->GetProperty()->SetOpacity(0.7);
            bondOrderActor_->SetVisibility(0);  // off by default, toggled by checkbox
            renderer_->AddActor(bondOrderActor_);

            udp_log("[diag] Bond order overlay: %lld bonds with BO >= 0.01\n",
                    (long long)lines->GetNumberOfCells());
        }
    }

    renderer_->ResetCamera();
    renderWindow_->Render();

    // Status: read tier classification from library (set by PredictionResult)
    int nReport = 0, nPass = 0, nSilent = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        switch (conf.AtomAt(i).tier) {
            case HeuristicTier::REPORT: nReport++; break;
            case HeuristicTier::PASS:   nPass++;   break;
            case HeuristicTier::SILENT: nSilent++; break;
        }
    }

    udp_log("[diag] setting status text\n");
    QString status = QString("Loaded %1: %2 atoms (%3 REPORT, %4 PASS, %5 SILENT), %6 rings")
        .arg(QString::fromStdString(currentProteinId_))
        .arg(protein.AtomCount())
        .arg(nReport)
        .arg(nPass)
        .arg(nSilent)
        .arg(protein.RingCount());
    statusLabel_->setText(status);

    udp_log("[diag] calling updateOverlay\n");
    updateOverlay();
    udp_log("[diag] onComputeFinished done\n");
}

// Placeholder — per-calculator visualizations will replace the old overlay system

void MainWindow::saveScreenshot() {
    QString path = QFileDialog::getSaveFileName(this, "Save Screenshot",
        "screenshot.png", "PNG Files (*.png)");
    if (path.isEmpty()) return;

    vtkNew<vtkWindowToImageFilter> filter;
    filter->SetInput(renderWindow_);
    filter->SetScale(1);
    filter->SetInputBufferTypeToRGBA();
    filter->ReadFrontBufferOn();
    filter->ShouldRerenderOff();
    filter->Update();

    vtkNew<vtkPNGWriter> writer;
    writer->SetFileName(path.toStdString().c_str());
    writer->SetInputConnection(filter->GetOutputPort());
    writer->Write();

    statusLabel_->setText(QString("Screenshot saved: %1").arg(path));
}

void MainWindow::updateOverlay() {
    // Per-calculator visualizations will replace this. For now, just ensure
    // the glyph infrastructure stays clear when no calculator viz is active.
    // Field grid and butterfly overlays are managed by their own toggles —
    // do NOT clear them here.
    if (!tensorGlyph_ || !ellipsoidGlyph_ || !protein_) return;
    tensorGlyph_->clear();
    ellipsoidGlyph_->clear();
    if (isosurfaceOverlay_) isosurfaceOverlay_->clear();
    if (isosurfaceOverlayPass_) isosurfaceOverlayPass_->clear();
    renderWindow_->Render();
}

void MainWindow::onRenderModeChanged(int index) {
    if (!molMapper_) return;
    switch (index) {
        case 0: molMapper_->UseBallAndStickSettings(); break;
        case 1: molMapper_->UseLiquoriceStickSettings(); break;
    }
    renderWindow_->Render();
}

void MainWindow::onShowRingsToggled(bool checked) {
    if (ringOverlay_) {
        ringOverlay_->setVisible(checked);
        renderWindow_->Render();
    }
}

void MainWindow::onShowPeptideBondsToggled(bool checked) {
    if (peptideBondOverlay_) {
        peptideBondOverlay_->setVisible(checked);
        renderWindow_->Render();
    }
}

void MainWindow::onShowBondOrderToggled(bool checked) {
    if (bondOrderActor_)
        bondOrderActor_->SetVisibility(checked ? 1 : 0);
    renderWindow_->Render();
}

void MainWindow::onShowButterflyToggled(bool checked) {
    if (butterflyOverlay_)
        butterflyOverlay_->setVisible(checked);
    renderWindow_->Render();
}

void MainWindow::onIsoThresholdChanged() {
    if (fieldGridOverlay_ && !fieldGrids_.empty()) {
        double threshold = isoThresholdSlider_->value() / 100.0;
        double opacity = currentScaleSlider_->value() / 100.0;
        fieldGridOverlay_->setData(fieldGrids_, threshold, opacity, 0);
        renderWindow_->Render();
    }
}

// ================================================================
// Operations log
// ================================================================

void MainWindow::onLogDatagramReady() {
    while (logSocket_->hasPendingDatagrams()) {
        QNetworkDatagram dg = logSocket_->receiveDatagram();
        QString msg = QString::fromUtf8(dg.data()).trimmed();
        if (msg.isEmpty()) continue;
        logText_->appendPlainText(msg);
    }
}

// ================================================================
// Atom picking and inspection
// ================================================================

bool MainWindow::eventFilter(QObject* obj, QEvent* event) {
    if (obj == vtkWidget_ && event->type() == QEvent::MouseButtonDblClick) {
        auto* me = static_cast<QMouseEvent*>(event);
        pickAtom(me->pos().x(), me->pos().y());
        return true;
    }
    return QMainWindow::eventFilter(obj, event);
}

void MainWindow::pickAtom(int displayX, int displayY) {
    if (!protein_) return;

    // Original session-4 ray casting, with DPR correction and logging.
    double dpr = vtkWidget_->devicePixelRatioF();
    int vtkX = static_cast<int>(displayX * dpr);
    int vtkY = static_cast<int>((vtkWidget_->height() - displayY) * dpr);

    udp_log("[pick] pickAtom qt=(%d,%d) vtk=(%d,%d) dpr=%.1f winSz=(%d,%d)\n",
            displayX, displayY, vtkX, vtkY, dpr,
            renderWindow_->GetSize()[0], renderWindow_->GetSize()[1]);

    auto* camera = renderer_->GetActiveCamera();
    double camPos[3];
    camera->GetPosition(camPos);
    Vec3 rayOrigin(camPos[0], camPos[1], camPos[2]);

    renderer_->SetDisplayPoint(vtkX, vtkY, 0.0);
    renderer_->DisplayToWorld();
    double worldPt[4];
    renderer_->GetWorldPoint(worldPt);
    Vec3 clickWorld(worldPt[0] / worldPt[3],
                    worldPt[1] / worldPt[3],
                    worldPt[2] / worldPt[3]);

    Vec3 rayDir = (clickWorld - rayOrigin).normalized();

    udp_log("[pick] ray origin=(%.2f,%.2f,%.2f) dir=(%.3f,%.3f,%.3f)\n",
            rayOrigin.x(), rayOrigin.y(), rayOrigin.z(),
            rayDir.x(), rayDir.y(), rayDir.z());

    const auto& conf = protein_->Conformation();
    double bestDist = 1e30;
    int bestAtom = -1;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        Vec3 pos = conf.AtomAt(i).Position();
        Vec3 toAtom = pos - rayOrigin;
        double projLen = toAtom.dot(rayDir);
        if (projLen < 0) continue;
        Vec3 closest = rayOrigin + projLen * rayDir;
        double dist = (pos - closest).norm();
        if (dist < bestDist) {
            bestDist = dist;
            bestAtom = static_cast<int>(i);
        }
    }

    udp_log("[pick] best atom=%d ray-dist=%.3f A\n", bestAtom, bestDist);

    if (bestAtom >= 0 && bestDist < 2.0) {
        size_t atomIndex = static_cast<size_t>(bestAtom);
        const auto& id = protein_->AtomAt(atomIndex);
        const auto& res = protein_->ResidueAt(id.residue_index);
        Vec3 pos = conf.AtomAt(atomIndex).Position();

        udp_log("[pick] → atom %zu: %s %s-%d\n", atomIndex,
                id.pdb_atom_name.c_str(),
                ThreeLetterCodeForAminoAcid(res.type).c_str(),
                res.sequence_number);

        populateAtomInfo(atomIndex);
        populateAtomBonds(atomIndex);
        populateGeometryChoices(atomIndex);
        atomInfoDock_->raise();

        if (selectionActor_) renderer_->RemoveActor(selectionActor_);
        vtkNew<vtkSphereSource> sphere;
        sphere->SetCenter(pos.x(), pos.y(), pos.z());
        sphere->SetRadius(1.0);
        sphere->SetPhiResolution(16);
        sphere->SetThetaResolution(16);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(sphere->GetOutputPort());
        selectionActor_ = vtkSmartPointer<vtkActor>::New();
        selectionActor_->SetMapper(mapper);
        selectionActor_->GetProperty()->SetColor(1.0, 1.0, 0.0);
        selectionActor_->GetProperty()->SetOpacity(0.3);
        renderer_->AddActor(selectionActor_);
        renderWindow_->Render();

        statusLabel_->setText(QString("Atom %1: %2 %3-%4")
            .arg(atomIndex)
            .arg(QString::fromStdString(id.pdb_atom_name))
            .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
            .arg(res.sequence_number));
    }
}

// Helper: format a SphericalTensor as a tree item
static QTreeWidgetItem* stItem(const QString& name, const SphericalTensor& st) {
    auto* item = new QTreeWidgetItem({name,
        QString("T0=%1").arg(st.T0, 0, 'f', 4)});
    item->addChild(new QTreeWidgetItem({"T0", QString::number(st.T0, 'f', 6)}));
    item->addChild(new QTreeWidgetItem({"T1",
        QString("(%1, %2, %3)")
            .arg(st.T1[0], 0, 'f', 5).arg(st.T1[1], 0, 'f', 5).arg(st.T1[2], 0, 'f', 5)}));
    item->addChild(new QTreeWidgetItem({"T2",
        QString("(%1, %2, %3, %4, %5)")
            .arg(st.T2[0], 0, 'f', 5).arg(st.T2[1], 0, 'f', 5).arg(st.T2[2], 0, 'f', 5)
            .arg(st.T2[3], 0, 'f', 5).arg(st.T2[4], 0, 'f', 5)}));
    return item;
}

// Helper: format a Vec3
static QString vec3Str(const Vec3& v) {
    return QString("(%1, %2, %3)").arg(v.x(), 0, 'f', 4).arg(v.y(), 0, 'f', 4).arg(v.z(), 0, 'f', 4);
}

void MainWindow::populateAtomInfo(size_t idx) {
    atomInfoTree_->clear();
    if (!protein_ || idx >= protein_->AtomCount()) return;

    const auto& protein = *protein_;
    const auto& conf = protein.Conformation();
    const auto& id = protein.AtomAt(idx);        // Atom (identity)
    const auto& ca = conf.AtomAt(idx);            // ConformationAtom (computed)
    const auto& res = protein.ResidueAt(id.residue_index);

    // ---- Header ----
    QString header = QString("Atom %1: %2 %3 (%4-%5-%6)")
        .arg(idx)
        .arg(QString::fromStdString(SymbolForElement(id.element)))
        .arg(QString::fromStdString(id.pdb_atom_name))
        .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
        .arg(res.sequence_number)
        .arg(QString::fromStdString(res.chain_id));

    // ---- Identity ----
    auto* identity = new QTreeWidgetItem({header, ""});
    identity->addChild(new QTreeWidgetItem({"Element", QString::fromStdString(SymbolForElement(id.element))
        + QString(" (Z=%1)").arg(AtomicNumberForElement(id.element))}));
    identity->addChild(new QTreeWidgetItem({"PDB name", QString::fromStdString(id.pdb_atom_name)}));
    identity->addChild(new QTreeWidgetItem({"Residue", QString("%1 %2 %3")
        .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
        .arg(res.sequence_number)
        .arg(QString::fromStdString(res.chain_id))}));
    identity->addChild(new QTreeWidgetItem({"Role", QString::fromStdString(NameForAtomRole(ca.role))}));
    identity->addChild(new QTreeWidgetItem({"Backbone", ca.is_backbone ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Amide H", ca.is_amide_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Alpha H", ca.is_alpha_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Aromatic H", ca.is_aromatic_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Methyl", ca.is_methyl ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Position", vec3Str(ca.Position())}));
    if (id.parent_atom_index != SIZE_MAX)
        identity->addChild(new QTreeWidgetItem({"Parent atom", QString::number(id.parent_atom_index)}));
    atomInfoTree_->addTopLevelItem(identity);
    identity->setExpanded(true);

    // ---- Charges & Electronic ----
    auto* charges = new QTreeWidgetItem({"Charges & Electronic", ""});
    charges->addChild(new QTreeWidgetItem({"Partial (ff14SB)", QString::number(ca.partial_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"MOPAC (PM7)", QString::number(ca.mopac_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"Delta (PM7-ff14SB)",
        QString::number(ca.mopac_charge - ca.partial_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"s-orbital pop", QString::number(ca.mopac_s_pop, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"p-orbital pop", QString::number(ca.mopac_p_pop, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"Valency (sum BO)", QString::number(ca.mopac_valency, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"VdW radius", QString::number(ca.vdw_radius, 'f', 3) + " A"}));
    atomInfoTree_->addTopLevelItem(charges);

    // ---- Shielding contributions (8 classical + 2 MOPAC-derived) ----
    auto* shielding = new QTreeWidgetItem({"Shielding contributions", ""});
    shielding->addChild(stItem("Biot-Savart", ca.bs_shielding_contribution));
    shielding->addChild(stItem("Haigh-Mallion", ca.hm_shielding_contribution));
    shielding->addChild(stItem("McConnell", ca.mc_shielding_contribution));
    shielding->addChild(stItem("Dispersion", ca.disp_shielding_contribution));
    shielding->addChild(stItem("Coulomb EFG", ca.coulomb_shielding_contribution));
    shielding->addChild(stItem("Pi-Quadrupole", ca.piquad_shielding_contribution));
    shielding->addChild(stItem("Ring Suscept.", ca.ringchi_shielding_contribution));
    shielding->addChild(stItem("H-Bond", ca.hbond_shielding_contribution));
    shielding->addChild(stItem("MOPAC Coulomb", ca.mopac_coulomb_shielding_contribution));
    shielding->addChild(stItem("MOPAC McConnell", ca.mopac_mc_shielding_contribution));

    // Predicted T0 from library (set by PredictionResult)
    shielding->addChild(new QTreeWidgetItem({"Predicted T0", QString::number(ca.predicted_T0, 'f', 4) + " ppm"}));
    shielding->addChild(new QTreeWidgetItem({"Tier", ca.tier == HeuristicTier::REPORT ? "REPORT" :
                                                     ca.tier == HeuristicTier::PASS ? "PASS" : "SILENT"}));
    atomInfoTree_->addTopLevelItem(shielding);
    shielding->setExpanded(true);

    // ---- Vector fields ----
    auto* fields = new QTreeWidgetItem({"Vector fields", ""});
    fields->addChild(new QTreeWidgetItem({"B field (BS)", vec3Str(ca.total_B_field)}));
    fields->addChild(new QTreeWidgetItem({"E field (ff14SB)", vec3Str(ca.coulomb_E_total)}));
    fields->addChild(new QTreeWidgetItem({"|E| ff14SB", QString::number(ca.coulomb_E_magnitude, 'f', 3) + " V/A"}));
    fields->addChild(new QTreeWidgetItem({"E bb frac (ff14SB)", QString::number(ca.coulomb_E_backbone_frac, 'f', 3)}));
    fields->addChild(new QTreeWidgetItem({"E field (MOPAC)", vec3Str(ca.mopac_coulomb_E_total)}));
    fields->addChild(new QTreeWidgetItem({"|E| MOPAC", QString::number(ca.mopac_coulomb_E_magnitude, 'f', 3) + " V/A"}));
    fields->addChild(new QTreeWidgetItem({"E bb frac (MOPAC)", QString::number(ca.mopac_coulomb_E_backbone_frac, 'f', 3)}));
    if (ca.apbs_efield.norm() > 1e-10)
        fields->addChild(new QTreeWidgetItem({"E field (APBS)", vec3Str(ca.apbs_efield)}));
    atomInfoTree_->addTopLevelItem(fields);

    // ---- Ring neighbours ----
    if (!ca.ring_neighbours.empty()) {
        auto* rings = new QTreeWidgetItem({"Ring neighbours",
            QString::number(ca.ring_neighbours.size())});
        for (const auto& rn : ca.ring_neighbours) {
            const Ring& ring = protein.RingAt(rn.ring_index);
            QString label = QString("%1 ring %2 (d=%3 A)")
                .arg(QString::fromStdString(ring.TypeName()))
                .arg(rn.ring_index)
                .arg(rn.distance_to_center, 0, 'f', 2);
            auto* rnItem = new QTreeWidgetItem({label, ""});
            rnItem->addChild(new QTreeWidgetItem({"rho", QString::number(rn.rho, 'f', 3) + " A"}));
            rnItem->addChild(new QTreeWidgetItem({"z", QString::number(rn.z, 'f', 3) + " A"}));
            rnItem->addChild(new QTreeWidgetItem({"theta", QString::number(rn.theta * 180.0 / M_PI, 'f', 1) + " deg"}));
            rnItem->addChild(stItem("G (BS kernel)", rn.G_spherical));
            rnItem->addChild(stItem("HM raw integral (H)", rn.hm_H_spherical));
            rnItem->addChild(stItem("HM shielding (G)", rn.hm_G_spherical));
            rnItem->addChild(new QTreeWidgetItem({"B field", vec3Str(rn.B_field)}));
            rnItem->addChild(new QTreeWidgetItem({"B cylindrical", vec3Str(rn.B_cylindrical)}));
            if (rn.chi_scalar != 0.0)
                rnItem->addChild(stItem("Chi (suscept.)", rn.chi_spherical));
            if (rn.quad_scalar != 0.0)
                rnItem->addChild(stItem("Quad (pi-quad)", rn.quad_spherical));
            if (rn.disp_contacts > 0) {
                rnItem->addChild(stItem("Dispersion", rn.disp_spherical));
                rnItem->addChild(new QTreeWidgetItem({"Disp contacts", QString::number(rn.disp_contacts)}));
            }
            rings->addChild(rnItem);
        }
        atomInfoTree_->addTopLevelItem(rings);
        rings->setExpanded(true);
    }

    // ---- MOPAC bond orders (direct covalent bonds from QM) ----
    if (!ca.mopac_bond_neighbours.empty()) {
        auto* mopacBonds = new QTreeWidgetItem({"MOPAC bonds (Wiberg order)",
            QString::number(ca.mopac_bond_neighbours.size())});
        for (const auto& mb : ca.mopac_bond_neighbours) {
            const auto& otherId = protein.AtomAt(mb.other_atom);
            const auto& otherRes = protein.ResidueAt(otherId.residue_index);
            QString label = QString("%1 %2-%3 (BO=%4)")
                .arg(QString::fromStdString(SymbolForElement(otherId.element)))
                .arg(QString::fromStdString(otherId.pdb_atom_name))
                .arg(otherRes.sequence_number)
                .arg(mb.wiberg_order, 0, 'f', 3);
            auto* mbItem = new QTreeWidgetItem({label, ""});
            mbItem->addChild(new QTreeWidgetItem({"Atom index", QString::number(mb.other_atom)}));
            if (mb.topology_bond_index != SIZE_MAX) {
                const Bond& bond = protein.BondAt(mb.topology_bond_index);
                mbItem->addChild(new QTreeWidgetItem({"Bond category",
                    QString::fromStdString(NameForBondCategory(bond.category))}));
            }
            mopacBonds->addChild(mbItem);
        }
        atomInfoTree_->addTopLevelItem(mopacBonds);
        mopacBonds->setExpanded(true);
    }

    // ---- McConnell bond neighbours (dipolar, within 10A) ----
    if (!ca.bond_neighbours.empty()) {
        auto* bonds = new QTreeWidgetItem({"McConnell bond neighbours",
            QString::number(ca.bond_neighbours.size())});
        for (const auto& bn : ca.bond_neighbours) {
            const Bond& bond = protein.BondAt(bn.bond_index);
            // Look up MOPAC bond order for this topology bond
            double mopacBO = 0.0;
            if (conf.HasResult<MopacResult>())
                mopacBO = conf.Result<MopacResult>().TopologyBondOrder(bn.bond_index);
            QString label = QString("Bond %1 (%2, d=%3 A)")
                .arg(bn.bond_index)
                .arg(QString::fromStdString(NameForBondCategory(bond.category)))
                .arg(bn.distance_to_midpoint, 0, 'f', 2);
            auto* bnItem = new QTreeWidgetItem({label, ""});
            bnItem->addChild(new QTreeWidgetItem({"McConnell scalar", QString::number(bn.mcconnell_scalar, 'f', 5)}));
            if (mopacBO > 1e-6)
                bnItem->addChild(new QTreeWidgetItem({"Wiberg order", QString::number(mopacBO, 'f', 3)}));
            bnItem->addChild(stItem("Dipolar tensor", bn.dipolar_spherical));
            bonds->addChild(bnItem);
        }
        atomInfoTree_->addTopLevelItem(bonds);
    }

    // ---- H-bond ----
    if (ca.hbond_nearest_dist > 0.01) {
        auto* hbond = new QTreeWidgetItem({"H-bond", ""});
        hbond->addChild(new QTreeWidgetItem({"Nearest dist", QString::number(ca.hbond_nearest_dist, 'f', 3) + " A"}));
        hbond->addChild(new QTreeWidgetItem({"Direction", vec3Str(ca.hbond_nearest_dir)}));
        hbond->addChild(new QTreeWidgetItem({"Count <3.5A", QString::number(ca.hbond_count_within_3_5A)}));
        hbond->addChild(new QTreeWidgetItem({"Is donor", ca.hbond_is_donor ? "yes" : "no"}));
        hbond->addChild(new QTreeWidgetItem({"Is acceptor", ca.hbond_is_acceptor ? "yes" : "no"}));
        hbond->addChild(new QTreeWidgetItem({"Backbone", ca.hbond_is_backbone ? "yes" : "no"}));
        atomInfoTree_->addTopLevelItem(hbond);
    }

    // ---- DSSP ----
    if (conf.HasResult<DsspResult>()) {
        const auto& dssp = conf.Result<DsspResult>();
        auto* dsspItem = new QTreeWidgetItem({"DSSP", ""});
        dsspItem->addChild(new QTreeWidgetItem({"Secondary", QString(QChar(dssp.SecondaryStructure(id.residue_index)))}));
        dsspItem->addChild(new QTreeWidgetItem({"Phi", QString::number(dssp.Phi(id.residue_index), 'f', 1) + " deg"}));
        dsspItem->addChild(new QTreeWidgetItem({"Psi", QString::number(dssp.Psi(id.residue_index), 'f', 1) + " deg"}));
        dsspItem->addChild(new QTreeWidgetItem({"SASA", QString::number(dssp.SASA(id.residue_index), 'f', 1) + " A^2"}));
        atomInfoTree_->addTopLevelItem(dsspItem);
    }

    // ---- ORCA DFT ----
    if (ca.has_orca_shielding) {
        auto* orca = new QTreeWidgetItem({"ORCA DFT shielding", ""});
        orca->addChild(stItem("Total", ca.orca_shielding_total_spherical));
        orca->addChild(stItem("Diamagnetic", ca.orca_shielding_diamagnetic_spherical));
        orca->addChild(stItem("Paramagnetic", ca.orca_shielding_paramagnetic_spherical));
        atomInfoTree_->addTopLevelItem(orca);
        orca->setExpanded(true);
    }

    // ---- McConnell breakdown (unweighted) ----
    if (std::abs(ca.mc_shielding_contribution.T0) > 1e-6) {
        auto* mc = new QTreeWidgetItem({"McConnell (unweighted)", ""});
        mc->addChild(new QTreeWidgetItem({"CO sum", QString::number(ca.mcconnell_co_sum, 'f', 5)}));
        mc->addChild(new QTreeWidgetItem({"CN sum", QString::number(ca.mcconnell_cn_sum, 'f', 5)}));
        mc->addChild(new QTreeWidgetItem({"Sidechain sum", QString::number(ca.mcconnell_sidechain_sum, 'f', 5)}));
        mc->addChild(new QTreeWidgetItem({"Aromatic sum", QString::number(ca.mcconnell_aromatic_sum, 'f', 5)}));
        mc->addChild(new QTreeWidgetItem({"Nearest CO dist", QString::number(ca.nearest_CO_dist, 'f', 3) + " A"}));
        atomInfoTree_->addTopLevelItem(mc);
    }

    // ---- MOPAC McConnell breakdown (bond-order-weighted) ----
    if (std::abs(ca.mopac_mc_shielding_contribution.T0) > 1e-6) {
        auto* mmc = new QTreeWidgetItem({"McConnell (BO-weighted)", ""});
        mmc->addChild(new QTreeWidgetItem({"CO sum", QString::number(ca.mopac_mc_co_sum, 'f', 5)}));
        mmc->addChild(new QTreeWidgetItem({"CN sum", QString::number(ca.mopac_mc_cn_sum, 'f', 5)}));
        mmc->addChild(new QTreeWidgetItem({"Sidechain sum", QString::number(ca.mopac_mc_sidechain_sum, 'f', 5)}));
        mmc->addChild(new QTreeWidgetItem({"Aromatic sum", QString::number(ca.mopac_mc_aromatic_sum, 'f', 5)}));
        mmc->addChild(new QTreeWidgetItem({"Nearest CO (BO-wt)", QString::number(ca.mopac_mc_co_nearest, 'f', 5)}));
        if (ca.mopac_mc_nearest_CO_dist > 0.01)
            mmc->addChild(new QTreeWidgetItem({"Nearest CO dist", QString::number(ca.mopac_mc_nearest_CO_dist, 'f', 3) + " A"}));
        if (ca.mopac_mc_nearest_CN_dist > 0.01)
            mmc->addChild(new QTreeWidgetItem({"Nearest CN dist", QString::number(ca.mopac_mc_nearest_CN_dist, 'f', 3) + " A"}));
        mmc->addChild(stItem("T2 backbone", ca.mopac_mc_T2_backbone_total));
        mmc->addChild(stItem("T2 sidechain", ca.mopac_mc_T2_sidechain_total));
        mmc->addChild(stItem("T2 aromatic", ca.mopac_mc_T2_aromatic_total));
        atomInfoTree_->addTopLevelItem(mmc);
    }
}

// ================================================================
// Bond tab — shows bonds for the currently picked atom
// ================================================================

void MainWindow::populateAtomBonds(size_t idx) {
    bondInfoTree_->clear();
    if (!protein_ || idx >= protein_->AtomCount()) return;

    const auto& protein = *protein_;
    const auto& conf = protein.Conformation();
    const auto& id = protein.AtomAt(idx);
    const auto& ca = conf.AtomAt(idx);

    // ---- Covalent bonds (from topology) ----
    if (!id.bond_indices.empty()) {
        auto* covalent = new QTreeWidgetItem({"Covalent bonds",
            QString::number(id.bond_indices.size())});
        for (size_t bi : id.bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t other = (bond.atom_index_a == idx) ? bond.atom_index_b : bond.atom_index_a;
            const auto& otherId = protein.AtomAt(other);
            const auto& otherRes = protein.ResidueAt(otherId.residue_index);

            QString label = QString("%1 %2-%3 (%4)")
                .arg(QString::fromStdString(otherId.pdb_atom_name))
                .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(otherRes.type)))
                .arg(otherRes.sequence_number)
                .arg(QString::fromStdString(NameForBondCategory(bond.category)));

            auto* bondItem = new QTreeWidgetItem({label,
                QString::number(conf.bond_lengths[bi], 'f', 3) + " A"});
            bondItem->addChild(new QTreeWidgetItem({"Direction", vec3Str(conf.bond_directions[bi])}));

            if (conf.HasResult<MopacResult>()) {
                double bo = conf.Result<MopacResult>().TopologyBondOrder(bi);
                bondItem->addChild(new QTreeWidgetItem({"Wiberg order", QString::number(bo, 'f', 4)}));
            }

            covalent->addChild(bondItem);
        }
        bondInfoTree_->addTopLevelItem(covalent);
        covalent->setExpanded(true);
    }

    // ---- MOPAC bond neighbours (may include non-topology bonds) ----
    if (!ca.mopac_bond_neighbours.empty()) {
        auto* mopacBonds = new QTreeWidgetItem({"MOPAC bonds",
            QString::number(ca.mopac_bond_neighbours.size())});
        for (const auto& mb : ca.mopac_bond_neighbours) {
            const auto& otherId = protein.AtomAt(mb.other_atom);
            const auto& otherRes = protein.ResidueAt(otherId.residue_index);
            QString label = QString("%1 %2-%3 (BO=%4)")
                .arg(QString::fromStdString(otherId.pdb_atom_name))
                .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(otherRes.type)))
                .arg(otherRes.sequence_number)
                .arg(mb.wiberg_order, 0, 'f', 3);
            auto* mbItem = new QTreeWidgetItem({label, ""});
            if (mb.topology_bond_index != SIZE_MAX) {
                const Bond& bond = protein.BondAt(mb.topology_bond_index);
                mbItem->addChild(new QTreeWidgetItem({"Category",
                    QString::fromStdString(NameForBondCategory(bond.category))}));
            }
            mopacBonds->addChild(mbItem);
        }
        bondInfoTree_->addTopLevelItem(mopacBonds);
        mopacBonds->setExpanded(true);
    }

    // ---- McConnell bond neighbours (dipolar, within 10A) ----
    if (!ca.bond_neighbours.empty()) {
        auto* mcBonds = new QTreeWidgetItem({"McConnell neighbours",
            QString::number(ca.bond_neighbours.size())});
        for (const auto& bn : ca.bond_neighbours) {
            const Bond& bond = protein.BondAt(bn.bond_index);
            QString label = QString("Bond %1 (%2, d=%3 A)")
                .arg(bn.bond_index)
                .arg(QString::fromStdString(NameForBondCategory(bond.category)))
                .arg(bn.distance_to_midpoint, 0, 'f', 2);
            auto* bnItem = new QTreeWidgetItem({label, ""});
            bnItem->addChild(new QTreeWidgetItem({"McConnell scalar",
                QString::number(bn.mcconnell_scalar, 'f', 5)}));
            bnItem->addChild(stItem("Dipolar tensor", bn.dipolar_spherical));
            mcBonds->addChild(bnItem);
        }
        bondInfoTree_->addTopLevelItem(mcBonds);
    }
}

// ================================================================
// GeometryChoice tab — calculator decisions affecting picked atom
// ================================================================

// GAP: no NameForCalculatorId in library Types.h yet
static const char* NameForCalculatorId(nmr::CalculatorId id) {
    using nmr::CalculatorId;
    switch (id) {
        case CalculatorId::BiotSavart:        return "Biot-Savart";
        case CalculatorId::HaighMallion:      return "Haigh-Mallion";
        case CalculatorId::McConnell:         return "McConnell";
        case CalculatorId::Coulomb:           return "Coulomb";
        case CalculatorId::PiQuadrupole:      return "Pi-Quadrupole";
        case CalculatorId::RingSusceptibility: return "Ring Susceptibility";
        case CalculatorId::Dispersion:        return "Dispersion";
        case CalculatorId::HBond:             return "H-Bond";
        case CalculatorId::MopacCoulomb:      return "MOPAC Coulomb";
        case CalculatorId::MopacMcConnell:    return "MOPAC McConnell";
        default:                              return "Other";
    }
}

static const char* NameForOutcome(nmr::EntityOutcome o) {
    switch (o) {
        case nmr::EntityOutcome::Included:     return "included";
        case nmr::EntityOutcome::Excluded:     return "EXCLUDED";
        case nmr::EntityOutcome::Triggered:    return "TRIGGERED";
        case nmr::EntityOutcome::NotTriggered: return "not triggered";
    }
    return "?";
}

void MainWindow::populateGeometryChoices(size_t atomIndex) {
    gcTree_->clear();
    if (!protein_ || atomIndex >= protein_->AtomCount()) return;

    const auto& conf = protein_->Conformation();
    const auto& choices = conf.geometry_choices;

    // Group choices by calculator where this atom appears as Target
    using nmr::CalculatorId;
    using nmr::EntityRole;
    using nmr::EntityOutcome;

    // Map: calculator → list of (choice index)
    std::map<CalculatorId, std::vector<size_t>> byCalc;

    for (size_t ci = 0; ci < choices.size(); ++ci) {
        const auto& gc = choices[ci];
        for (const auto& ent : gc.Entities()) {
            if (ent.atom_index == atomIndex && ent.role == EntityRole::Target) {
                byCalc[gc.Calculator()].push_back(ci);
                break;
            }
        }
    }

    if (byCalc.empty()) {
        gcTree_->addTopLevelItem(new QTreeWidgetItem(
            {"No geometry choices reference this atom", ""}));
        return;
    }

    for (const auto& [calcId, indices] : byCalc) {
        // Count included vs excluded
        int nIncluded = 0, nExcluded = 0, nTriggered = 0;
        for (size_t ci : indices) {
            const auto& gc = choices[ci];
            for (const auto& ent : gc.Entities()) {
                if (ent.atom_index != atomIndex) continue;
                if (ent.outcome == EntityOutcome::Included) ++nIncluded;
                else if (ent.outcome == EntityOutcome::Excluded) ++nExcluded;
                else if (ent.outcome == EntityOutcome::Triggered) ++nTriggered;
            }
        }

        QString summary = QString("%1 inc, %2 exc").arg(nIncluded).arg(nExcluded);
        if (nTriggered > 0) summary += QString(", %1 trig").arg(nTriggered);

        auto* calcItem = new QTreeWidgetItem(
            {NameForCalculatorId(calcId), summary});

        for (size_t ci : indices) {
            const auto& gc = choices[ci];

            // Find this atom's outcome
            EntityOutcome atomOutcome = EntityOutcome::Included;
            QString filterName;
            for (const auto& ent : gc.Entities()) {
                if (ent.atom_index == atomIndex) {
                    atomOutcome = ent.outcome;
                    if (!ent.filter_name.empty())
                        filterName = QString::fromStdString(ent.filter_name);
                    break;
                }
            }

            QString label = QString::fromStdString(gc.Label());
            QString detail = NameForOutcome(atomOutcome);
            if (!filterName.isEmpty())
                detail += " (" + filterName + ")";

            auto* choiceItem = new QTreeWidgetItem({label, detail});

            // Source entities (rings/bonds that generated this decision)
            for (const auto& ent : gc.Entities()) {
                if (ent.role != EntityRole::Source) continue;
                if (ent.ring) {
                    choiceItem->addChild(new QTreeWidgetItem(
                        {"source ring", QString::fromStdString(ent.ring->TypeName())
                         + " res " + QString::number(ent.ring->parent_residue_number)}));
                } else if (ent.bond) {
                    choiceItem->addChild(new QTreeWidgetItem(
                        {"source bond", QString::fromStdString(
                            NameForBondCategory(ent.bond->category))}));
                }
            }

            // Named numbers
            for (const auto& nn : gc.Numbers()) {
                choiceItem->addChild(new QTreeWidgetItem(
                    {QString::fromStdString(nn.name),
                     QString::number(nn.value, 'f', 3)
                     + (nn.unit.empty() ? "" : " " + QString::fromStdString(nn.unit))}));
            }

            // Sampler indicator
            if (gc.HasSampler()) {
                choiceItem->addChild(new QTreeWidgetItem(
                    {"sampler", "available (field probe)"}));
            }

            calcItem->addChild(choiceItem);
        }

        calcItem->setExpanded(true);
        gcTree_->addTopLevelItem(calcItem);
    }
}
