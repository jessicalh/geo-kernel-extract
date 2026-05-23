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
#include "AminoAcidType.h"
#include "Atom.h"
#include "Bond.h"
#include "ConformationAtom.h"
#include "ConformationResult.h"
#include "DsspResult.h"
#include "GeometryChoice.h"
#include "LegacyAmberTopology.h"
#include "MopacCoulombResult.h"
#include "MopacMcConnellResult.h"
#include "MopacResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "SemanticEnums.h"
#include "Session.h"

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

// AMBER substrate enum-to-string helpers. These are display-only —
// the runtime uses typed enums for every decision. NameFor* live here
// (not in SemanticEnums.h) because the library is mid-lint and adding
// to a public header is out of ui/ scope (per ui/CLAUDE.md). Same
// local pattern as NameForAtomRole / NameForBondCategory above.

static const char* NameForPlanarGroupKind(PlanarGroupKind k) {
    switch (k) {
    case PlanarGroupKind::None:
        return "—";
    case PlanarGroupKind::PeptideAmide:
        return "PeptideAmide";
    case PlanarGroupKind::SidechainAmide:
        return "SidechainAmide";
    case PlanarGroupKind::Guanidinium:
        return "Guanidinium";
    case PlanarGroupKind::Imidazole:
        return "Imidazole";
    case PlanarGroupKind::Aromatic6Ring:
        return "Aromatic6Ring";
    case PlanarGroupKind::Aromatic5Ring:
        return "Aromatic5Ring";
    case PlanarGroupKind::Carboxylate:
        return "Carboxylate";
    case PlanarGroupKind::AromaticHydroxyl:
        return "AromaticHydroxyl";
    case PlanarGroupKind::AromaticOxide:
        return "AromaticOxide";
    }
    return "?";
}

static const char* NameForPolarHKind(PolarHKind k) {
    switch (k) {
    case PolarHKind::NotPolar:
        return "NotPolar";
    case PolarHKind::BackboneAmide:
        return "BackboneAmide";
    case PolarHKind::SidechainPrimaryAmide:
        return "SidechainPrimaryAmide";
    case PolarHKind::IndoleNH:
        return "IndoleNH";
    case PolarHKind::AmmoniumNH:
        return "AmmoniumNH";
    case PolarHKind::GuanidiniumNH:
        return "GuanidiniumNH";
    case PolarHKind::ImidazoleNH:
        return "ImidazoleNH";
    case PolarHKind::CarboxylOH:
        return "CarboxylOH";
    case PolarHKind::HydroxylOH_Aliphatic:
        return "HydroxylOH (aliphatic)";
    case PolarHKind::HydroxylOH_Aromatic:
        return "HydroxylOH (aromatic)";
    case PolarHKind::ThiolSH:
        return "ThiolSH";
    case PolarHKind::AmineNH:
        return "AmineNH";
    case PolarHKind::OtherPolarH:
        return "OtherPolarH";
    }
    return "?";
}

static const char* NameForProchiralStereo(ProchiralStereo s) {
    switch (s) {
    case ProchiralStereo::NotProchiral:
        return "NotProchiral";
    case ProchiralStereo::ProR:
        return "ProR";
    case ProchiralStereo::ProS:
        return "ProS";
    case ProchiralStereo::Unassigned:
        return "Unassigned";
    }
    return "?";
}

static const char* NameForPseudoatomKind(PseudoatomKind k) {
    switch (k) {
    case PseudoatomKind::None:
        return "None";
    case PseudoatomKind::M:
        return "M (methyl)";
    case PseudoatomKind::Q:
        return "Q (equivalent-H group)";
    case PseudoatomKind::R:
        return "R (ring aggregator)";
    }
    return "?";
}

static const char* NameForRingSystemKind(RingSystemKind k) {
    switch (k) {
    case RingSystemKind::NotInRing:
        return "—";
    case RingSystemKind::Benzene_Phe:
        return "Phe benzene";
    case RingSystemKind::Benzene_Tyr:
        return "Tyr benzene";
    case RingSystemKind::Imidazole_His:
        return "His imidazole";
    case RingSystemKind::Indole_Trp_5:
        return "Trp pyrrole";
    case RingSystemKind::Indole_Trp_6:
        return "Trp benzene";
    case RingSystemKind::Pyrrolidine_Pro:
        return "Pro pyrrolidine";
    case RingSystemKind::Indole_Trp_9:
        return "Trp indole perimeter";
    }
    return "?";
}

static const char* NameForRingPositionLabel(RingPositionLabel p) {
    switch (p) {
    case RingPositionLabel::NotInRing:
        return "—";
    case RingPositionLabel::Ipso:
        return "ipso";
    case RingPositionLabel::Ortho1:
        return "ortho1";
    case RingPositionLabel::Ortho2:
        return "ortho2";
    case RingPositionLabel::Meta1:
        return "meta1";
    case RingPositionLabel::Meta2:
        return "meta2";
    case RingPositionLabel::Para:
        return "para";
    case RingPositionLabel::PyrroleAlpha:
        return "pyrrole α";
    case RingPositionLabel::PyrroleBeta:
        return "pyrrole β";
    case RingPositionLabel::BridgeFusion:
        return "bridge (fused)";
    case RingPositionLabel::Heteroatom_NH:
        return "heteroatom NH";
    case RingPositionLabel::Heteroatom_NoH:
        return "heteroatom (no H)";
    case RingPositionLabel::Heteroatom_OH:
        return "heteroatom OH";
    case RingPositionLabel::Saturated:
        return "saturated";
    case RingPositionLabel::ProRingNitrogen:
        return "Pro ring N";
    case RingPositionLabel::ProRingAlphaCarbon:
        return "Pro ring Cα";
    case RingPositionLabel::ProRingBeta:
        return "Pro ring Cβ";
    case RingPositionLabel::ProRingPuckerPivot:
        return "Pro pucker pivot (Cγ)";
    case RingPositionLabel::ProRingDelta:
        return "Pro ring Cδ";
    case RingPositionLabel::PerimeterMember:
        return "perimeter member";
    }
    return "?";
}

MainWindow::MainWindow(nmr::Session& session,
                       const QString& udpHost,
                       quint16 udpPort,
                       const QString& initialDir,
                       QWidget* parent)
    : QMainWindow(parent)
    , session_(session)
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
    , udpHost_(udpHost)
    , udpPort_(udpPort) {
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
    for (auto* timer : timers) {
        timer->stop();
    }

    // 2. Stop async workers
    cancelCompute();

    // 3. Finalize VTK before Qt tears down the GL context
    if (renderWindow_) {
        renderWindow_->Finalize();
    }

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

    QString const dir = QFileDialog::getExistingDirectory(this, "Export Features", QDir::currentPath());
    if (dir.isEmpty()) return;

    std::string const outDir = dir.toStdString();
    int totalArrays = 0;

    std::filesystem::create_directories(outDir);
    auto& conf = protein_->Conformation();
    totalArrays = nmr::ConformationResult::WriteAllFeatures(conf, outDir);

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

    vtkNew<vtkInteractorStyleTrackballCamera> const style;
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

        // Context object (3rd arg) = header: when the QPushButton is
        // destroyed, Qt auto-disconnects this lambda so it can't fire
        // on dangling header/content captures. clazy -Wclazy-connect-
        // 3arg-lambda flagged the missing context.
        QObject::connect(header, &QPushButton::clicked, header, [header, content, title]() {
            bool const show = !content->isVisible();
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
                double const opacity = currentScaleSlider_->value() / 100.0;
                double const threshold = isoThresholdSlider_->value() / 100.0;
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

    // Bind UDP listener on the configured destination (udpHost_:udpPort_
    // sourced from [logging] in ~/.nmr_tools.toml by main_viewer.cpp).
    //
    // Multicast (239.0.0.0/8): bind AnyIPv4 so the kernel hands us
    // datagrams destined for the group, then join the group via
    // joinMulticastGroup. This lets udp_listen.py and other subscribers
    // co-listen on the same port — Linux unicast UDP would otherwise
    // route each datagram to a single bound socket. ShareAddress +
    // ReuseAddressHint let stale sockets be reclaimed cleanly.
    //
    // Unicast: keep the original host-specific bind (preserves the
    // "127.0.0.1 to local listener" path for offline use).
    //
    // There is intentionally NO fallback to ports 9999+. A fallback
    // would silently bind the wrong port (library always sends to
    // udpPort_), making the tab appear to work while receiving nothing.
    // Fail visibly instead so the cause is obvious.
    logSocket_ = new QUdpSocket(this);

    auto isMulticast = [](const QString& host) {
        const QStringList parts = host.split('.');
        if (parts.size() != 4) {
            return false;
        }
        bool ok = false;
        const int first = parts[0].toInt(&ok);
        return ok && first >= 224 && first <= 239;
    };

    const bool multicast = isMulticast(udpHost_);
    bool bound = false;
    if (multicast) {
        bound = logSocket_->bind(QHostAddress::AnyIPv4, udpPort_, QUdpSocket::ShareAddress | QUdpSocket::ReuseAddressHint);
        if (bound) {
            const QHostAddress group(udpHost_);
            if (!logSocket_->joinMulticastGroup(group)) {
                logText_->appendPlainText(QStringLiteral("** Log tab: bound %1:%2 but failed to join "
                                                         "multicast group; tab will be quiet. **")
                                              .arg(udpHost_)
                                              .arg(udpPort_));
            }
        }
    } else {
        bound = logSocket_->bind(QHostAddress(udpHost_), udpPort_, QUdpSocket::ShareAddress);
    }

    if (bound) {
        connect(logSocket_, &QUdpSocket::readyRead, this, &MainWindow::onLogDatagramReady);
    } else {
        logText_->appendPlainText(QStringLiteral("** Log tab inactive: failed to bind UDP %1:%2 (%3). **\n"
                                                 "Multicast: check that the configured group is in 239.0.0.0/8 "
                                                 "and IGMP isn't blocked locally.\n"
                                                 "Unicast: is another process already holding the port?\n"
                                                 "See ui/launch_viewer.sh + [logging] in ~/.nmr_tools.toml.")
                                      .arg(udpHost_)
                                      .arg(udpPort_)
                                      .arg(multicast ? "multicast" : "unicast"));
    }

    // Start with the log tab visible — it shows computation progress
    logDock_->raise();

    // Time Series (H5) panel — read-only per-atom frame-0 slice of the
    // companion analysis H5 when --analysis-h5 was supplied on the
    // command line. The viewer never writes H5 files.
    timeSeriesDock_ = new QDockWidget("Time Series (H5)", this);
    timeSeriesDock_->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);
    timeSeriesDock_->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    timeSeriesTree_ = new QTreeWidget();
    timeSeriesTree_->setHeaderLabels({"Property", "Value"});
    timeSeriesTree_->setColumnWidth(0, 260);
    timeSeriesTree_->setAlternatingRowColors(true);
    timeSeriesTree_->setIndentation(16);
    timeSeriesDock_->setWidget(timeSeriesTree_);
    tabifyDockWidget(atomInfoDock_, timeSeriesDock_);
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
    case nmr::JobMode::Trajectory:
        // Trajectory mode is CLI-only; not loaded via viewer.
        return;
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

    udp_log("[diag] VTK setup: %lld ms\n", static_cast<long long>(timer.elapsed()));
    timer.restart();

    renderer_->ResetCamera();
    udp_log("[diag] ResetCamera done, calling Render()\n");
    renderWindow_->Render();
    udp_log("[diag] First render: %lld ms\n", static_cast<long long>(timer.elapsed()));

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
    worker_ = new ComputeWorker(session_);
    worker_->moveToThread(workerThread_);

    connect(this, &MainWindow::computeRequested, worker_, &ComputeWorker::computeAll);
    connect(worker_, &ComputeWorker::progress, this, &MainWindow::onComputeProgress);
    connect(worker_, &ComputeWorker::finished, this, &MainWindow::onComputeFinished);

    // Worker / thread lifetime via Qt's canonical pattern. When the
    // worker emits finished:
    //   1. workerThread quits its event loop (processes pending events
    //      including the deleteLater below before exiting),
    //   2. the worker is deleteLater'd on its own thread (event loop
    //      drains and runs the deferred delete cleanly),
    //   3. the thread emits QThread::finished after exec() returns,
    //      triggering its own deleteLater on the main thread (where it
    //      was created).
    // cancelCompute and onComputeFinished still call wait() so shutdown
    // is deterministic; with this wiring they no longer raw-delete the
    // worker or thread objects. The agent flagged the previous explicit
    // delete as safe today but fragile against a future ComputeWorker
    // that grows QObject children (timers, network access manager).
    connect(worker_, &ComputeWorker::finished, workerThread_, &QThread::quit);
    connect(worker_, &ComputeWorker::finished, worker_, &QObject::deleteLater);
    connect(workerThread_, &QThread::finished, workerThread_, &QObject::deleteLater);

    progressDialog_ = new QProgressDialog("Computing features...", "Cancel", 0, 100, this);
    progressDialog_->setMinimumDuration(0);
    progressDialog_->setValue(0);
    connect(progressDialog_, &QProgressDialog::canceled, this, &MainWindow::cancelCompute);

    workerThread_->start();
    emit computeRequested(pendingSpec_);
}

void MainWindow::cancelCompute() {
    if (worker_) {
        worker_->cancel();
    }
    if (workerThread_ && workerThread_->isRunning()) {
        workerThread_->quit();
        workerThread_->wait();
    }
    // worker_ / workerThread_ auto-delete via the deleteLater wiring in
    // startCompute. Just drop our references so we don't reuse stale
    // pointers; Qt processes the deferred deletions on the appropriate
    // event loops.
    worker_ = nullptr;
    workerThread_ = nullptr;
    if (progressDialog_) {
        progressDialog_->close();
        progressDialog_->deleteLater();
        progressDialog_ = nullptr;
    }
}

void MainWindow::onComputeProgress(int current, int total, const QString& phase) {
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
    analysisBinding_ = std::move(result.analysisBinding);  // Valid() false if no H5

    // Stop worker thread. worker_ / workerThread_ auto-delete via the
    // deleteLater wiring in startCompute (see comment there).
    if (workerThread_ && workerThread_->isRunning()) {
        workerThread_->quit();
        workerThread_->wait();
    }
    worker_ = nullptr;
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
            unsigned short const anum = static_cast<unsigned short>(AtomicNumberForElement(protein.AtomAt(i).element));
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
                static_cast<int>(molecule_->GetNumberOfAtoms()),
                static_cast<int>(molecule_->GetNumberOfBonds()));
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
        double const threshold = isoThresholdSlider_->value() / 100.0;
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
            double const bo = (i < topoOrders.size()) ? topoOrders[i] : 0.0;
            if (bo < 0.01) continue;  // skip electronically insignificant bonds

            const Bond& bond = protein.BondAt(i);
            Vec3 posA = conf.AtomAt(bond.atom_index_a).Position();
            Vec3 posB = conf.AtomAt(bond.atom_index_b).Position();
            vtkIdType const id0 = pts->InsertNextPoint(posA.data());
            vtkIdType const id1 = pts->InsertNextPoint(posB.data());
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
                    static_cast<long long>(lines->GetNumberOfCells()));
        }
    }

    renderer_->ResetCamera();
    renderWindow_->Render();

    // Status line: protein name + atom + bond + ring counts. The old
    // HeuristicTier counting (REPORT/PASS/SILENT) was removed — that
    // field is from the pre-kernel-catalogue prediction era (UI_ROADMAP
    // Known Issues #1) and no current ConformationResult writes it.
    udp_log("[diag] setting status text\n");
    QString const status = QString("Loaded %1: %2 atoms, %3 bonds, %4 rings")
                               .arg(QString::fromStdString(currentProteinId_))
                               .arg(protein.AtomCount())
                               .arg(protein.BondCount())
                               .arg(protein.RingCount());
    statusLabel_->setText(status);

    udp_log("[diag] calling updateOverlay\n");
    updateOverlay();
    udp_log("[diag] onComputeFinished done\n");
}

// Placeholder — per-calculator visualizations will replace the old overlay system

void MainWindow::saveScreenshot() {
    QString const path = QFileDialog::getSaveFileName(this, "Save Screenshot", "screenshot.png", "PNG Files (*.png)");
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
    if (bondOrderActor_) {
        bondOrderActor_->SetVisibility(checked ? 1 : 0);
    }
    renderWindow_->Render();
}

void MainWindow::onShowButterflyToggled(bool checked) {
    if (butterflyOverlay_) {
        butterflyOverlay_->setVisible(checked);
    }
    renderWindow_->Render();
}

void MainWindow::onIsoThresholdChanged() {
    if (fieldGridOverlay_ && !fieldGrids_.empty()) {
        double const threshold = isoThresholdSlider_->value() / 100.0;
        double const opacity = currentScaleSlider_->value() / 100.0;
        fieldGridOverlay_->setData(fieldGrids_, threshold, opacity, 0);
        renderWindow_->Render();
    }
}

// ================================================================
// Operations log
// ================================================================

void MainWindow::onLogDatagramReady() {
    while (logSocket_->hasPendingDatagrams()) {
        QNetworkDatagram const dg = logSocket_->receiveDatagram();
        QString const msg = QString::fromUtf8(dg.data()).trimmed();
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
    double const dpr = vtkWidget_->devicePixelRatioF();
    int const vtkX = static_cast<int>(displayX * dpr);
    int const vtkY = static_cast<int>((vtkWidget_->height() - displayY) * dpr);

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
    Vec3 const clickWorld(worldPt[0] / worldPt[3], worldPt[1] / worldPt[3], worldPt[2] / worldPt[3]);

    Vec3 rayDir = (clickWorld - rayOrigin).normalized();

    udp_log("[pick] ray origin=(%.2f,%.2f,%.2f) dir=(%.3f,%.3f,%.3f)\n",
            rayOrigin.x(), rayOrigin.y(), rayOrigin.z(),
            rayDir.x(), rayDir.y(), rayDir.z());

    const auto& conf = protein_->Conformation();
    double bestDist = 1e30;
    int bestAtom = -1;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        Vec3 const pos = conf.AtomAt(i).Position();
        Vec3 const toAtom = pos - rayOrigin;
        double const projLen = toAtom.dot(rayDir);
        if (projLen < 0) continue;
        Vec3 const closest = rayOrigin + projLen * rayDir;
        double const dist = (pos - closest).norm();
        if (dist < bestDist) {
            bestDist = dist;
            bestAtom = static_cast<int>(i);
        }
    }

    udp_log("[pick] best atom=%d ray-dist=%.3f A\n", bestAtom, bestDist);

    if (bestAtom >= 0 && bestDist < 2.0) {
        size_t const atomIndex = static_cast<size_t>(bestAtom);
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
        populateTimeSeries(atomIndex);
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
                                  .arg(QString::number(static_cast<qulonglong>(atomIndex)),
                                       QString::fromStdString(id.pdb_atom_name),
                                       QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)),
                                       QString::number(res.sequence_number)));
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
    QString const header = QString("Atom %1: %2 %3 (%4-%5-%6)")
                               .arg(QString::number(static_cast<qulonglong>(idx)),
                                    QString::fromStdString(SymbolForElement(id.element)),
                                    QString::fromStdString(id.pdb_atom_name),
                                    QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)),
                                    QString::number(res.sequence_number),
                                    QString::fromStdString(res.chain_id));

    // ---- Identity ----
    auto* identity = new QTreeWidgetItem({header, ""});
    identity->addChild(new QTreeWidgetItem({"Element", QString::fromStdString(SymbolForElement(id.element))
        + QString(" (Z=%1)").arg(AtomicNumberForElement(id.element))}));
    identity->addChild(new QTreeWidgetItem({"PDB name", QString::fromStdString(id.pdb_atom_name)}));
    identity->addChild(new QTreeWidgetItem({"Residue", QString("%1 %2 %3")
        .arg(QString::fromStdString(ThreeLetterCodeForAminoAcid(res.type)))
        .arg(res.sequence_number)
        .arg(QString::fromStdString(res.chain_id))}));
    // Terminal state + protonation variant are typed Residue fields
    // populated at load time by LegacyAmberTopology / ProtonationDetection.
    // Variant name comes from AminoAcidType::variants[i].name —
    // distinguishes HID/HIE/HIP, ASP/ASH, GLU/GLH, LYS/LYN, CYS/CYX/CYM.
    identity->addChild(new QTreeWidgetItem({"Terminal state", ResidueTerminalStateName(res.terminal_state)}));
    if (res.protonation_variant_index >= 0) {
        const AminoAcidType& aaType = GetAminoAcidType(res.type);
        if (res.protonation_variant_index < static_cast<int>(aaType.variants.size())) {
            const auto& variant = aaType.variants[res.protonation_variant_index];
            identity->addChild(new QTreeWidgetItem(
                {"Protonation variant",
                 QString("%1 — %2").arg(QString::fromUtf8(variant.name), QString::fromUtf8(variant.description))}));
        }
    }
    identity->addChild(new QTreeWidgetItem({"Role", QString::fromStdString(NameForAtomRole(ca.role))}));
    identity->addChild(new QTreeWidgetItem({"Backbone", ca.is_backbone ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Amide H", ca.is_amide_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Alpha H", ca.is_alpha_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Aromatic H", ca.is_aromatic_H ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Methyl", ca.is_methyl ? "yes" : "no"}));
    identity->addChild(new QTreeWidgetItem({"Position", vec3Str(ca.Position())}));
    if (id.parent_atom_index != SIZE_MAX) {
        identity->addChild(new QTreeWidgetItem({"Parent atom", QString::number(id.parent_atom_index)}));
    }
    atomInfoTree_->addTopLevelItem(identity);
    identity->setExpanded(true);

    // ---- AMBER substrate (typed chemistry) ----
    // LegacyAmberTopology::SemanticAt(i) is the typed atom-chemistry
    // substrate populated at load time by ComposeAtomSemantic
    // (CCD + RDKit perception, baked into src/generated/
    // LegacyAmberSemanticTables.cpp). Every field here is the typed
    // answer to a chemistry question that used to require a string
    // match against pdb_atom_name — backbone_role, planar_group,
    // ring membership, polar-H taxonomy, prochiral CIP assignment,
    // pseudoatom grouping, formal charge.
    //
    // The flags already in Identity (role, is_amide_H, is_methyl, …)
    // are EnrichmentResult-derived and overlap with the substrate.
    // This section shows what the substrate adds that EnrichmentResult
    // does NOT carry: the finer-grained chemistry taxonomy.
    //
    // Gate: HasAtomSemantic() — empty for any load path that didn't
    // compose the substrate. PDB load + AMBER-as-standard path always
    // populates it.
    {
        const auto& topo = protein.LegacyAmber();
        if (topo.HasAtomSemantic()) {
            const auto& sem = topo.SemanticAt(idx);
            auto* asub = new QTreeWidgetItem({"AMBER substrate (typed)", ""});

            asub->addChild(new QTreeWidgetItem({"Planar group", NameForPlanarGroupKind(sem.planar_group)}));

            // Ring position: skip if not in any ring. When present, show
            // the primary ring system + position label. Secondary /
            // tertiary slots populated only for fused Trp atoms.
            if (sem.ring_position.IsInAnyRing()) {
                auto* ring = new QTreeWidgetItem({"Ring position", ""});
                const auto& p = sem.ring_position.primary;
                ring->addChild(new QTreeWidgetItem({"Primary ring", NameForRingSystemKind(p.ring)}));
                ring->addChild(new QTreeWidgetItem({"Primary position", NameForRingPositionLabel(p.position)}));
                if (sem.ring_position.HasSecondaryRing()) {
                    const auto& s = sem.ring_position.secondary;
                    ring->addChild(
                        new QTreeWidgetItem({"Secondary ring",
                                             QString("%1 (%2)").arg(QString::fromUtf8(NameForRingSystemKind(s.ring)),
                                                                    QString::fromUtf8(NameForRingPositionLabel(s.position)))}));
                }
                if (sem.ring_position.HasTertiaryRing()) {
                    const auto& t = sem.ring_position.tertiary;
                    ring->addChild(
                        new QTreeWidgetItem({"Tertiary ring",
                                             QString("%1 (%2)").arg(QString::fromUtf8(NameForRingSystemKind(t.ring)),
                                                                    QString::fromUtf8(NameForRingPositionLabel(t.position)))}));
                }
                asub->addChild(ring);
            }

            // Polar H taxonomy. Only meaningful for H atoms; for non-H
            // it's always NotPolar (skip to reduce visual noise).
            if (sem.polar_h != PolarHKind::NotPolar) {
                asub->addChild(new QTreeWidgetItem({"Polar H kind", NameForPolarHKind(sem.polar_h)}));
            }

            // Prochiral CIP designation (Markley 1998 §2.1.4). Only
            // meaningful for prochiral atoms; skip NotProchiral.
            if (sem.prochiral != ProchiralStereo::NotProchiral) {
                asub->addChild(new QTreeWidgetItem({"Prochiral", NameForProchiralStereo(sem.prochiral)}));
            }

            // Pseudoatom membership (Markley 1998 Table 1). Show only
            // when the atom is in some pseudoatom group.
            if (sem.pseudoatom.IsMember()) {
                QString detail = NameForPseudoatomKind(sem.pseudoatom.kind);
                if (sem.pseudoatom.in_super_group) {
                    detail += " (+ super-group)";
                }
                asub->addChild(new QTreeWidgetItem({"Pseudoatom", detail}));
            }

            // Formal charge — show only when non-zero. Most atoms are 0.
            if (sem.formal_charge != 0) {
                asub->addChild(new QTreeWidgetItem({"Formal charge", QString::number(static_cast<int>(sem.formal_charge))}));
            }

            // Exchangeable flag (H/D exchange chemistry). Show only when true.
            if (sem.is_exchangeable) {
                asub->addChild(new QTreeWidgetItem({"Exchangeable", "yes"}));
            }

            atomInfoTree_->addTopLevelItem(asub);
            asub->setExpanded(true);
        }
    }

    // ---- Charges & Electronic ----
    auto* charges = new QTreeWidgetItem({"Charges & Electronic", ""});
    charges->addChild(new QTreeWidgetItem({"Partial (ff14SB)", QString::number(ca.partial_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"AIMNet2 (Hirshfeld wB97M)", QString::number(ca.aimnet2_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"EEQ (Caldeweyher 2019)", QString::number(ca.eeq_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"EEQ coord. number", QString::number(ca.eeq_cn, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"MOPAC (PM7)", QString::number(ca.mopac_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"Delta (PM7-ff14SB)",
        QString::number(ca.mopac_charge - ca.partial_charge, 'f', 4)}));
    charges->addChild(new QTreeWidgetItem({"s-orbital pop", QString::number(ca.mopac_s_pop, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"p-orbital pop", QString::number(ca.mopac_p_pop, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"Valency (sum BO)", QString::number(ca.mopac_valency, 'f', 3)}));
    charges->addChild(new QTreeWidgetItem({"PB radius", QString::number(ca.pb_radius, 'f', 3) + " A"}));
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
    shielding->addChild(stItem("H-Bond (kernel)", ca.hbond_shielding_contribution));
    shielding->addChild(stItem("Larsen H-Bond (grid)", ca.larsen_hbond_shielding_spherical));
    shielding->addChild(stItem("Tripeptide BB (DFT)", ca.tripeptide_bb_shielding_spherical));
    shielding->addChild(stItem("Tripeptide neighbor", ca.tripeptide_neighbor_shielding_spherical));
    shielding->addChild(stItem("AIMNet2 EFG", ca.aimnet2_shielding_contribution));
    shielding->addChild(stItem("MOPAC Coulomb", ca.mopac_coulomb_shielding_contribution));
    shielding->addChild(stItem("MOPAC McConnell", ca.mopac_mc_shielding_contribution));

    // Predicted T0 / Tier were removed: PredictionResult is pre-kernel-
    // catalogue framing per UI_ROADMAP Known Issues #1; no current
    // ConformationResult writes those fields. Don't surface what no
    // calculator populates.
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
    if (ca.apbs_efield.norm() > 1e-10) {
        fields->addChild(new QTreeWidgetItem({"E field (APBS)", vec3Str(ca.apbs_efield)}));
    }
    // APBS EFG tensor (solvated). Same kernel as Coulomb EFG but solvated.
    if (ca.apbs_efg.norm() > 1e-10) {
        fields->addChild(stItem("EFG (APBS)", ca.apbs_efg_spherical));
    }
    // AIMNet2 EFG family — same kernel as CoulombResult but with the
    // wB97M neural Hirshfeld charges. Total + decomposed by source.
    fields->addChild(stItem("EFG (AIMNet2 total)", ca.aimnet2_EFG_total_spherical));
    fields->addChild(stItem("EFG (AIMNet2 backbone)", ca.aimnet2_EFG_backbone_spherical));
    fields->addChild(stItem("EFG (AIMNet2 aromatic)", ca.aimnet2_EFG_aromatic_spherical));
    atomInfoTree_->addTopLevelItem(fields);

    // ---- Local geometry (planar pyramidalization, SASA) ----
    // pyramidalization is signed out-of-plane displacement (Å) at atoms
    // whose AtomSemanticTable::planar_group != None (PlanarGeometryResult).
    // atom_sasa is Shrake-Rupley SASA (Å²); sasa_normal is the outward
    // surface normal from non-occluded test points (SasaResult). All three
    // are written for every atom every run, so they always appear (even
    // when zero, which is itself informative — sp3 / buried interior).
    {
        auto* geo = new QTreeWidgetItem({"Local geometry", ""});
        geo->addChild(new QTreeWidgetItem({"Pyramidalization", QString::number(ca.pyramidalization, 'f', 4) + " A (signed)"}));
        geo->addChild(new QTreeWidgetItem({"SASA", QString::number(ca.atom_sasa, 'f', 2) + " A^2"}));
        geo->addChild(new QTreeWidgetItem({"SASA normal", vec3Str(ca.sasa_normal)}));
        atomInfoTree_->addTopLevelItem(geo);
    }

    // ---- Graph topology (MolecularGraphResult) ----
    // Through-bond BFS features: distance (in bonds) to the nearest
    // ring atom and to the nearest N / O; electronegativity sums over
    // 1-bond and 2-bond neighbours; pi-bond count within 3 bonds; a
    // bool flag for chains of alternating single/double bonds; an
    // exponential decay e^(-d/BFS_DECAY_LENGTH) over the ring distance
    // (BFS_DECAY_LENGTH = 4.0 from MolecularGraphResult.h).
    //
    // Distances default to -1 when the atom is not reachable from the
    // target class (e.g. an isolated methyl C with no path to any ring
    // atom). Show "—" instead of "-1" so the inspector reads cleanly.
    {
        auto distOrDash = [](int d) { return d < 0 ? QStringLiteral("—") : QString::number(d); };
        auto* graph = new QTreeWidgetItem({"Graph topology", ""});
        graph->addChild(new QTreeWidgetItem({"Bonds to nearest ring atom", distOrDash(ca.graph_dist_ring)}));
        graph->addChild(new QTreeWidgetItem({"Bonds to nearest N", distOrDash(ca.graph_dist_N)}));
        graph->addChild(new QTreeWidgetItem({"Bonds to nearest O", distOrDash(ca.graph_dist_O)}));
        graph->addChild(new QTreeWidgetItem({"e^(-d_ring/4) decay", QString::number(ca.bfs_decay, 'f', 4)}));
        graph->addChild(new QTreeWidgetItem({"Σ electronegativity (1-bond)", QString::number(ca.eneg_sum_1, 'f', 3)}));
        graph->addChild(new QTreeWidgetItem({"Σ electronegativity (2-bond)", QString::number(ca.eneg_sum_2, 'f', 3)}));
        graph->addChild(new QTreeWidgetItem({"π bonds within 3 bonds", QString::number(ca.n_pi_bonds_3)}));
        graph->addChild(new QTreeWidgetItem({"Conjugated", ca.is_conjugated ? "yes" : "no"}));
        atomInfoTree_->addTopLevelItem(graph);
    }

    // ---- AIMNet2 polarisability (charge-response gradient) ----
    // dL/d(r_i) where L = sum_j q_j^2, computed by autograd through the
    // AIMNet2 forward pass (AIMNet2ChargeResponseGradientResult). The
    // scalar is the L2 norm of the vector — charge-weighted per-atom
    // polarisability. Vector direction encodes WHICH way moving the
    // atom changes the network's charge prediction most.
    {
        auto* pol = new QTreeWidgetItem({"AIMNet2 polarisability", ""});
        pol->addChild(new QTreeWidgetItem({"dL/dr", vec3Str(ca.aimnet2_charge_response_gradient_vector)}));
        pol->addChild(new QTreeWidgetItem({"|dL/dr|", QString::number(ca.aimnet2_charge_response_gradient_scalar, 'f', 4)}));
        atomInfoTree_->addTopLevelItem(pol);
    }

    // ---- AIMNet2 aim embedding (256-dim, geometry-dependent) ----
    // The learned electronic-structure embedding from the AIMNet2 wB97M
    // model. 256 dims is too many to show; summarise via L2 norm + the
    // first four components (mirrors the H5 Time Series tab format).
    {
        double n2 = 0.0;
        for (size_t k = 0; k < nmr::AIMNET2_AIM_DIMS; ++k) {
            double const v = static_cast<double>(ca.aimnet2_aim[k]);
            n2 += v * v;
        }
        auto* aim = new QTreeWidgetItem(
            {"AIMNet2 embedding", QString("%1 dims (float32)").arg(static_cast<qulonglong>(nmr::AIMNET2_AIM_DIMS))});
        aim->addChild(new QTreeWidgetItem({"L2 norm^2", QString::number(n2, 'f', 5)}));
        aim->addChild(new QTreeWidgetItem({"aim[0..3]",
                                           QString("%1  %2  %3  %4")
                                               .arg(ca.aimnet2_aim[0], 0, 'f', 4)
                                               .arg(ca.aimnet2_aim[1], 0, 'f', 4)
                                               .arg(ca.aimnet2_aim[2], 0, 'f', 4)
                                               .arg(ca.aimnet2_aim[3], 0, 'f', 4)}));
        atomInfoTree_->addTopLevelItem(aim);
    }

    // ---- Ring neighbours ----
    if (!ca.ring_neighbours.empty()) {
        auto* rings = new QTreeWidgetItem({"Ring neighbours",
            QString::number(ca.ring_neighbours.size())});
        for (const auto& rn : ca.ring_neighbours) {
            const Ring& ring = protein.RingAt(rn.ring_index);
            QString const label = QString("%1 ring %2 (d=%3 A)")
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
            if (rn.chi_scalar != 0.0) {
                rnItem->addChild(stItem("Chi (suscept.)", rn.chi_spherical));
            }
            if (rn.quad_scalar != 0.0) {
                rnItem->addChild(stItem("Quad (pi-quad)", rn.quad_spherical));
            }
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
            QString const label = QString("%1 %2-%3 (BO=%4)")
                                      .arg(QString::fromStdString(SymbolForElement(otherId.element)),
                                           QString::fromStdString(otherId.pdb_atom_name),
                                           QString::number(otherRes.sequence_number),
                                           QString::number(mb.wiberg_order, 'f', 3));
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
            if (conf.HasResult<MopacResult>()) {
                mopacBO = conf.Result<MopacResult>().TopologyBondOrder(bn.bond_index);
            }
            QString const label = QString("Bond %1 (%2, d=%3 A)")
                                      .arg(bn.bond_index)
                                      .arg(QString::fromStdString(NameForBondCategory(bond.category)))
                                      .arg(bn.distance_to_midpoint, 0, 'f', 2);
            auto* bnItem = new QTreeWidgetItem({label, ""});
            bnItem->addChild(new QTreeWidgetItem({"McConnell scalar", QString::number(bn.mcconnell_scalar, 'f', 5)}));
            if (mopacBO > 1e-6) {
                bnItem->addChild(new QTreeWidgetItem({"Wiberg order", QString::number(mopacBO, 'f', 3)}));
            }
            bnItem->addChild(stItem("Dipolar tensor", bn.dipolar_spherical));
            bonds->addChild(bnItem);
        }
        atomInfoTree_->addTopLevelItem(bonds);
    }

    // ---- H-bond (kernel form, HBondResult) ----
    if (ca.hbond_nearest_dist > 0.01) {
        auto* hbond = new QTreeWidgetItem({"H-bond (kernel)", ""});
        hbond->addChild(new QTreeWidgetItem({"Nearest dist", QString::number(ca.hbond_nearest_dist, 'f', 3) + " A"}));
        hbond->addChild(new QTreeWidgetItem({"Direction", vec3Str(ca.hbond_nearest_dir)}));
        hbond->addChild(new QTreeWidgetItem({"Count <3.5A", QString::number(ca.hbond_count_within_3_5A)}));
        hbond->addChild(new QTreeWidgetItem({"Is donor", ca.hbond_is_donor ? "yes" : "no"}));
        hbond->addChild(new QTreeWidgetItem({"Is acceptor", ca.hbond_is_acceptor ? "yes" : "no"}));
        hbond->addChild(new QTreeWidgetItem({"Backbone", ca.hbond_is_backbone ? "yes" : "no"}));
        atomInfoTree_->addTopLevelItem(hbond);
    }

    // ---- Larsen H-bond (grid lookup; LarsenHBondShieldingResult) ----
    // Four contribution classes from the 6-archive ProCS15 DFT scans
    // (Δσ_1°HB + Δσ_2°HB + Δσ_1°HαB + Δσ_2°HαB) plus the Δσ_w water term
    // (2.07 ppm isotropic on amide H atoms with no H-bond pairs).
    // Methods-accumulate sibling of the kernel-form HBondResult — both
    // run side-by-side and produce comparable T0 the user can eyeball.
    // CB diagnostic should be ~0 per Larsen Table 2 (Cβ gets no HB
    // contribution); nonzero CB indicates a parser/loader/rotation bug.
    if (ca.larsen_hbond_n_pairs > 0 || std::abs(ca.larsen_hbond_water_term) > 1e-9
        || std::abs(ca.larsen_hbond_shielding_spherical.T0) > 1e-9) {
        auto* lhb = new QTreeWidgetItem({"Larsen H-bond (grid)", ""});
        lhb->addChild(new QTreeWidgetItem({"Pairs contributing", QString::number(ca.larsen_hbond_n_pairs)}));
        if (std::abs(ca.larsen_hbond_water_term) > 1e-9) {
            lhb->addChild(new QTreeWidgetItem(
                {"Water term Δσ_w", QString::number(ca.larsen_hbond_water_term, 'f', 4) + " ppm (isotropic)"}));
        }
        lhb->addChild(stItem("Sum (all classes)", ca.larsen_hbond_shielding_spherical));
        // Per-class breakdowns only when non-trivial; reduces visual noise on
        // atom types that get only one or two of the four contributions.
        auto addClassIfPresent = [&](const QString& name, const nmr::SphericalTensor& st) {
            if (std::abs(st.T0) > 1e-9) {
                lhb->addChild(stItem(name, st));
            }
        };
        addClassIfPresent("Δσ_1°HB  (donor primary)", ca.larsen_hbond_1pHB_spherical);
        addClassIfPresent("Δσ_2°HB  (donor secondary)", ca.larsen_hbond_2pHB_spherical);
        addClassIfPresent("Δσ_1°HαB (Hα primary)", ca.larsen_hbond_1pHaB_spherical);
        addClassIfPresent("Δσ_2°HαB (Hα secondary)", ca.larsen_hbond_2pHaB_spherical);
        // CB diagnostic — should be ~0 per Larsen Table 2; nonzero =
        // pipeline regression (parser → loader → rotation).
        if (std::abs(ca.larsen_hbond_diagnostic_CB_spherical.T0) > 1e-6) {
            lhb->addChild(stItem("Cβ diagnostic (expected ~0)", ca.larsen_hbond_diagnostic_CB_spherical));
        }
        if (ca.larsen_hbond_any_corner_imputed) {
            lhb->addChild(new QTreeWidgetItem({"Corner imputed", "yes (nearest-neighbour fill in at least one grid lookup)"}));
        }
        atomInfoTree_->addTopLevelItem(lhb);
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

    // ---- Tripeptide DFT shielding (Larsen 2015 σ_BB^i + Δσ_BB^{i±1}) ----
    // Pulled from the local ProCS15 tensorcs15 replica via libpq (BB,
    // central atom) and AXA-scan reuse (neighbor i±1 sidechain effect).
    // Tensor lab-frame: Kabsch backbone alignment + sidechain rotation
    // around CA-CB. method_tag distinguishes the source DFT method
    // (1 = Gaussian OPBE/6-31G(d,p) for 19 residues; 2 = ORCA PBE/6-31G(d,p)
    // for SER) — see project_serine_pbe_discontinuity.
    if (ca.tripeptide_bb_has_match || ca.tripeptide_neighbor_has_match) {
        auto* tp = new QTreeWidgetItem({"Tripeptide DFT (ProCS15)", ""});

        if (ca.tripeptide_bb_has_match) {
            auto* bb = new QTreeWidgetItem({"Backbone (i)", ""});
            bb->addChild(stItem("σ_BB", ca.tripeptide_bb_shielding_spherical));
            bb->addChild(
                new QTreeWidgetItem({"Match distance", QString::number(ca.tripeptide_bb_match_distance, 'f', 3) + " A"}));
            bb->addChild(new QTreeWidgetItem({"Residual (aligned − protein)", vec3Str(ca.tripeptide_bb_residual_vec)}));
            const char* tag = ca.tripeptide_bb_method_tag == 1   ? "Gaussian OPBE/6-31G(d,p)"
                              : ca.tripeptide_bb_method_tag == 2 ? "ORCA PBE/6-31G(d,p) (SER regen)"
                                                                 : "no match";
            bb->addChild(new QTreeWidgetItem({"DFT method", tag}));
            tp->addChild(bb);
        }

        if (ca.tripeptide_neighbor_has_match) {
            auto* nb = new QTreeWidgetItem({"Neighbor Δσ_BB^{i±1}", ""});
            nb->addChild(stItem("σ (sum of i-1 + i+1)", ca.tripeptide_neighbor_shielding_spherical));
            nb->addChild(new QTreeWidgetItem({"Residual prev (i-1)", vec3Str(ca.tripeptide_neighbor_residual_vec_prev)}));
            nb->addChild(new QTreeWidgetItem({"Residual next (i+1)", vec3Str(ca.tripeptide_neighbor_residual_vec_next)}));
            tp->addChild(nb);
        }

        atomInfoTree_->addTopLevelItem(tp);
        tp->setExpanded(true);
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
        if (ca.mopac_mc_nearest_CO_dist > 0.01) {
            mmc->addChild(new QTreeWidgetItem({"Nearest CO dist", QString::number(ca.mopac_mc_nearest_CO_dist, 'f', 3) + " A"}));
        }
        if (ca.mopac_mc_nearest_CN_dist > 0.01) {
            mmc->addChild(new QTreeWidgetItem({"Nearest CN dist", QString::number(ca.mopac_mc_nearest_CN_dist, 'f', 3) + " A"}));
        }
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
        for (size_t const bi : id.bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t const other = (bond.atom_index_a == idx) ? bond.atom_index_b : bond.atom_index_a;
            const auto& otherId = protein.AtomAt(other);
            const auto& otherRes = protein.ResidueAt(otherId.residue_index);

            QString const label = QString("%1 %2-%3 (%4)")
                                      .arg(QString::fromStdString(otherId.pdb_atom_name),
                                           QString::fromStdString(ThreeLetterCodeForAminoAcid(otherRes.type)),
                                           QString::number(otherRes.sequence_number),
                                           QString::fromStdString(NameForBondCategory(bond.category)));

            auto* bondItem = new QTreeWidgetItem({label,
                QString::number(conf.bond_lengths[bi], 'f', 3) + " A"});
            bondItem->addChild(new QTreeWidgetItem({"Direction", vec3Str(conf.bond_directions[bi])}));

            if (conf.HasResult<MopacResult>()) {
                double const bo = conf.Result<MopacResult>().TopologyBondOrder(bi);
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
            QString const label = QString("%1 %2-%3 (BO=%4)")
                                      .arg(QString::fromStdString(otherId.pdb_atom_name),
                                           QString::fromStdString(ThreeLetterCodeForAminoAcid(otherRes.type)),
                                           QString::number(otherRes.sequence_number),
                                           QString::number(mb.wiberg_order, 'f', 3));
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
            QString const label = QString("Bond %1 (%2, d=%3 A)")
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
        int nIncluded = 0;
        int nExcluded = 0;
        int nTriggered = 0;
        for (size_t const ci : indices) {
            const auto& gc = choices[ci];
            for (const auto& ent : gc.Entities()) {
                if (ent.atom_index != atomIndex) continue;
                if (ent.outcome == EntityOutcome::Included) {
                    ++nIncluded;
                } else if (ent.outcome == EntityOutcome::Excluded) {
                    ++nExcluded;
                } else if (ent.outcome == EntityOutcome::Triggered) {
                    ++nTriggered;
                }
            }
        }

        QString summary = QString("%1 inc, %2 exc").arg(nIncluded).arg(nExcluded);
        if (nTriggered > 0) summary += QString(", %1 trig").arg(nTriggered);

        auto* calcItem = new QTreeWidgetItem(
            {NameForCalculatorId(calcId), summary});

        for (size_t const ci : indices) {
            const auto& gc = choices[ci];

            // Find this atom's outcome
            EntityOutcome atomOutcome = EntityOutcome::Included;
            QString filterName;
            for (const auto& ent : gc.Entities()) {
                if (ent.atom_index == atomIndex) {
                    atomOutcome = ent.outcome;
                    if (!ent.filter_name.empty()) {
                        filterName = QString::fromStdString(ent.filter_name);
                    }
                    break;
                }
            }

            QString const label = QString::fromStdString(gc.Label());
            QString detail = NameForOutcome(atomOutcome);
            if (!filterName.isEmpty()) {
                detail += " (" + filterName + ")";
            }

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

// ────────────────────────────────────────────────────────────────────
//  Time Series (H5): read-only per-atom view of the trajectory.h5
//  companion. Reads the typed TrajectoryH5 cached on AnalysisBinding.
//
//  Each section is sparse-tolerant: it appears in the tree only if
//  the underlying group is present in the file. Welford rollups and
//  frame-0 slabs are surfaced side-by-side per physics family — the
//  "live calc at this pose vs. ensemble mean/std" diagnostic.
//
//  Adding a new TR family = one TrajectoryH5 accessor pair + one
//  if-block here. Per-PATTERNS no string dispatch, no generic
//  Get(name) surface: the file is typed at the H5 ingress and the
//  inspector reads typed accessors.
// ────────────────────────────────────────────────────────────────────

void MainWindow::populateTimeSeries(size_t idx) {
    timeSeriesTree_->clear();
    if (!analysisBinding_.Valid()) {
        timeSeriesTree_->addTopLevelItem(new QTreeWidgetItem(
            {QStringLiteral("(no trajectory H5 loaded)"), QStringLiteral("pass --analysis-h5 PATH on the command line")}));
        return;
    }

    const size_t h5idx = analysisBinding_.H5IndexFor(idx);
    const TrajectoryH5& h5 = *analysisBinding_.h5;

    // ── H5 meta + per-atom identity ──────────────────────────────
    {
        auto* g = new QTreeWidgetItem({"H5", QString::fromStdString(h5.ProteinId())});
        g->addChild(new QTreeWidgetItem({"frame", QString("0 of %1").arg(static_cast<qulonglong>(h5.FrameCount()))}));
        g->addChild(new QTreeWidgetItem({"frame time (ps)", QString::number(h5.FrameTimePs(0), 'f', 3)}));
        g->addChild(new QTreeWidgetItem({"n_atoms", QString::number(static_cast<qulonglong>(h5.AtomCount()))}));
        g->addChild(new QTreeWidgetItem({"library atom index",
            QString::number(static_cast<qulonglong>(idx))}));
        g->addChild(new QTreeWidgetItem({"H5 atom index (via H5IndexFor)",
            QString::number(static_cast<qulonglong>(h5idx))}));
        g->addChild(new QTreeWidgetItem({"atom element (Z, from H5)", QString::number(h5.ElementAt(h5idx))}));
        g->addChild(new QTreeWidgetItem({"atom name (H5)", QString::fromStdString(h5.AtomNameAt(h5idx))}));
        QStringList groups;
        for (const auto& gn : h5.GroupsPresent()) {
            groups << QString::fromStdString(gn);
        }
        g->addChild(new QTreeWidgetItem({"groups present", groups.isEmpty() ? QStringLiteral("(none)") : groups.join(", ")}));
        g->addChild(new QTreeWidgetItem({"atom-name mismatches (total)",
                                         QString::number(static_cast<qulonglong>(analysisBinding_.nameMismatches.size()))}));
        g->setExpanded(true);
        timeSeriesTree_->addTopLevelItem(g);
    }

    // Shared renderer for tensor-source shielding sections — same shape
    // across BS / HM / McConnell / PiQuad / RingChi / Disp / HBond.
    // Welford rollup may be present, absent, or accompanied by frame 0;
    // every combination is sparse-tolerant.
    auto addShieldingSection = [&](const QString& title,
                                   const QString& units,
                                   std::optional<TrajectoryH5::ShieldingWelfordRow> w,
                                   std::optional<TrajectoryH5::ShieldingFrame0Row> f) {
        if (!w && !f)
            return;
        auto* g = new QTreeWidgetItem({title, QString()});
        if (f) {
            g->addChild(new QTreeWidgetItem({"T0 frame 0", QString::number(f->T0, 'f', 5) + " " + units}));
            g->addChild(new QTreeWidgetItem({"|T2| frame 0", QString::number(f->T2_magnitude, 'f', 5) + " " + units}));
        }
        if (w) {
            g->addChild(new QTreeWidgetItem(
                {"T0 mean ± std", QString("%1 ± %2 %3").arg(w->t0.mean, 0, 'f', 5).arg(w->t0.std, 0, 'f', 5).arg(units)}));
            g->addChild(new QTreeWidgetItem(
                {"|T2| mean ± std",
                 QString("%1 ± %2 %3").arg(w->t2magnitude.mean, 0, 'f', 5).arg(w->t2magnitude.std, 0, 'f', 5).arg(units)}));
        }
        timeSeriesTree_->addTopLevelItem(g);
    };

    addShieldingSection("Ring Current (BS)", "ppm", h5.BsWelford(h5idx), h5.BsShieldingFrame0(h5idx));
    addShieldingSection("Ring Current (HM)", "Å⁻¹", h5.HmWelford(h5idx), h5.HmShieldingFrame0(h5idx));
    addShieldingSection("Bond Anisotropy (McConnell)", "Å⁻³", h5.McWelford(h5idx), h5.McShieldingFrame0(h5idx));
    addShieldingSection("Pi-quadrupole", "Å⁻⁵", std::nullopt, h5.PiQuadShieldingFrame0(h5idx));
    addShieldingSection("Ring susceptibility", "Å⁻³", std::nullopt, h5.RingChiShieldingFrame0(h5idx));
    addShieldingSection("Dispersion", "Å⁻⁶", std::nullopt, h5.DispShieldingFrame0(h5idx));
    addShieldingSection("H-bond (kernel-form)", "Å⁻³", std::nullopt, h5.HBondShieldingFrame0(h5idx));

    // ── SASA — welford rollup + frame 0 ─────────────────────────
    {
        auto sw = h5.SasaWelford(h5idx);
        auto sf = h5.SasaFrame0(h5idx);
        if (sw || sf) {
            auto* g = new QTreeWidgetItem({"SASA", QString()});
            if (sf) {
                g->addChild(new QTreeWidgetItem({"frame 0", QString::number(*sf, 'f', 3) + " Å²"}));
            }
            if (sw) {
                g->addChild(new QTreeWidgetItem(
                    {"mean ± std", QString("%1 ± %2 Å²").arg(sw->sasa.mean, 0, 'f', 3).arg(sw->sasa.std, 0, 'f', 3)}));
            }
            timeSeriesTree_->addTopLevelItem(g);
        }
    }

    // ── EEQ charge — welford rollup ─────────────────────────────
    if (auto ew = h5.EeqWelford(h5idx)) {
        auto* g = new QTreeWidgetItem({"EEQ charge", QString()});
        g->addChild(new QTreeWidgetItem(
            {"mean ± std", QString("%1 ± %2 e").arg(ew->charge.mean, 0, 'f', 5).arg(ew->charge.std, 0, 'f', 5)}));
        timeSeriesTree_->addTopLevelItem(g);
    }

    // ── AIMNet2 charge — frame 0 ────────────────────────────────
    if (auto ac = h5.Aimnet2ChargeFrame0(h5idx)) {
        auto* g = new QTreeWidgetItem({"AIMNet2 charge", QString()});
        g->addChild(new QTreeWidgetItem({"frame 0", QString::number(*ac, 'f', 5) + " e"}));
        timeSeriesTree_->addTopLevelItem(g);
    }

    // ── H-bond count — welford rollup ───────────────────────────
    if (auto hc = h5.HBondCountWelford(h5idx)) {
        auto* g = new QTreeWidgetItem({"H-bond count (3.5 Å)", QString()});
        g->addChild(new QTreeWidgetItem(
            {"mean ± std", QString("%1 ± %2 pairs").arg(hc->count.mean, 0, 'f', 3).arg(hc->count.std, 0, 'f', 3)}));
        timeSeriesTree_->addTopLevelItem(g);
    }

    // ── APBS efield — Cartesian Vec3 frame 0 ────────────────────
    if (auto ef = h5.ApbsEfieldFrame0(h5idx)) {
        auto* g = new QTreeWidgetItem({"APBS efield", QString()});
        g->addChild(new QTreeWidgetItem(
            {"E frame 0 (V/Å)", QString("(%1, %2, %3)").arg(ef->x, 0, 'f', 5).arg(ef->y, 0, 'f', 5).arg(ef->z, 0, 'f', 5)}));
        const double mag = std::sqrt(ef->x * ef->x + ef->y * ef->y + ef->z * ef->z);
        g->addChild(new QTreeWidgetItem({"|E| frame 0", QString::number(mag, 'f', 5) + " V/Å"}));
        timeSeriesTree_->addTopLevelItem(g);
    }
}
