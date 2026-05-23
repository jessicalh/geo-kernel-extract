#pragma once
#include <QMainWindow>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkOpenGLMoleculeMapper.h>
#include <vtkMolecule.h>
#include "ComputeWorker.h"
// AnalysisBinding (incl. its TrajectoryH5 inner pointer) comes via ComputeWorker.h.

class BackboneRibbonOverlay;
class RingCurrentOverlay;
class PeptideBondOverlay;
class TensorGlyph;
class EllipsoidGlyph;
class FieldOverlay;
class IsosurfaceOverlay;
class ButterflyOverlay;
class FieldGridOverlay;
class QSlider;
class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QGroupBox;
class QThread;
class QProgressDialog;
class QTimer;
class QTreeWidget;
class QTreeWidgetItem;
class QDockWidget;
class QPlainTextEdit;
class QUdpSocket;

namespace nmr {
class Protein;
class ProteinConformation;
class Session;
}

class RestServer;

class MainWindow : public QMainWindow {
    Q_OBJECT
    friend class RestServer;  // REST server needs access to UI controls and VTK renderer
public:
    // session: process-wide resource owner (AIMNet2 model, Tripeptide
    // DFT table libpq connection, Larsen H-bond grids). Constructed in
    // main_viewer.cpp; MainWindow holds a reference for the duration
    // and threads it through to ComputeWorker so OperationRunner::Run
    // can attach the dependent calculators. The reference outlives
    // MainWindow (Session is on main's stack; this widget exits before
    // main returns).
    //
    // udpHost / udpPort: the same UDP destination OperationLog is
    // sending to (sourced from [logging] in ~/.nmr_tools.toml by
    // main_viewer.cpp). When udpHost is an IPv4 multicast address
    // (239.0.0.0/8), the Log-dock socket binds AnyIPv4 and joins the
    // group via joinMulticastGroup, letting udp_listen.py and other
    // subscribers co-listen on the same port. Unicast hosts get a
    // direct host-specific bind (preserving the original behaviour
    // when the TOML is set to 127.0.0.1 for offline use).
    explicit MainWindow(nmr::Session& session,
                        const QString& udpHost,
                        quint16 udpPort,
                        const QString& initialDir = QString(),
                        QWidget* parent = nullptr);
    ~MainWindow() override;

    // Orderly shutdown while QApplication is still alive.
    // Called from aboutToQuit handler — stops timers, workers, VTK.
    void shutdown();

    // Load from a validated JobSpec (all modes: pdb, orca, mutant, fleet)
    void loadFromJobSpec(const nmr::JobSpec& spec);

signals:
    void computeRequested(nmr::JobSpec spec);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private slots:
    void saveScreenshot();
    void onRenderModeChanged(int index);
    void onShowRingsToggled(bool checked);
    void onShowPeptideBondsToggled(bool checked);
    void onShowBondOrderToggled(bool checked);
    void onShowButterflyToggled(bool checked);
    void onIsoThresholdChanged();

    // Async compute slots
    void onComputeProgress(int current, int total, const QString& phase);
    void onComputeFinished(ComputeResult result);

private:
    // Process-wide resource holder. Lifetime owned by main_viewer.cpp;
    // MainWindow holds a reference for the entire widget lifetime and
    // hands it to every ComputeWorker constructed via startCompute().
    nmr::Session& session_;
    QString initialDir_;
    void setupUI();
    void setupMenuBar();
    void exportFeatures();
    void loadMolecule();
    void startCompute();
    void cancelCompute();
    void updateOverlay();

    // Atom picking and inspection
    void pickAtom(int displayX, int displayY);
    void populateAtomInfo(size_t atomIndex);

    // Bond tab — shows bonds for the currently picked atom
    void populateAtomBonds(size_t atomIndex);

    // GeometryChoice tab — shows calculator decisions for picked atom
    void populateGeometryChoices(size_t atomIndex);

    // Time Series tab — per-atom view of the trajectory.h5 companion.
    // Shows Welford rollups (mean ± std) + frame-0 slabs per TR group.
    // Each section is sparse-tolerant: only appears if the group is in
    // the file. No-op (with "no H5 loaded" placeholder) when binding
    // is not Valid().
    void populateTimeSeries(size_t atomIndex);

    // VTK rendering
    QVTKOpenGLNativeWidget* vtkWidget_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkMolecule> molecule_;
    vtkSmartPointer<vtkOpenGLMoleculeMapper> molMapper_;
    vtkSmartPointer<vtkActor> molActor_;

    // Overlays
    BackboneRibbonOverlay* ribbonOverlay_ = nullptr;
    RingCurrentOverlay* ringOverlay_;
    PeptideBondOverlay* peptideBondOverlay_;
    TensorGlyph* tensorGlyph_;
    EllipsoidGlyph* ellipsoidGlyph_;
    FieldOverlay* fieldOverlay_;
    IsosurfaceOverlay* isosurfaceOverlay_;
    IsosurfaceOverlay* isosurfaceOverlayPass_;  // PASS tier (reduced opacity)
    ButterflyOverlay* butterflyOverlay_;
    FieldGridOverlay* fieldGridOverlay_ = nullptr;

    // Selection highlight
    vtkSmartPointer<vtkActor> selectionActor_;

    // Data — the library Protein, fully const after OperationRunner::Run.
    std::shared_ptr<nmr::Protein> protein_;
    std::vector<ViewerFieldGrid> fieldGrids_;
    std::vector<ViewerButterflyData> butterflyFields_;
    std::string currentProteinId_;

    // Companion time-series binding.  The viewer never writes H5 files
    // and never triggers a new extraction run.  Valid() iff --analysis-h5
    // was supplied AND the identity check passed.  All time-series reads
    // route through analysisBinding_.H5IndexFor(libAtomIdx) — one call
    // site to grow if a future producer emits non-identity ordering.
    AnalysisBinding analysisBinding_;

    // Async computation
    QThread* workerThread_ = nullptr;
    ComputeWorker* worker_ = nullptr;
    QProgressDialog* progressDialog_ = nullptr;


    // UI Controls — sidebar
    QComboBox* renderModeCombo_;
    QCheckBox* showRibbonCheck_;
    QSlider* glyphScaleSlider_;
    QSlider* opacitySlider_;
    QSlider* currentScaleSlider_;
    QCheckBox* showRingsCheck_;
    QCheckBox* showPeptideBondsCheck_;
    QCheckBox* showBondOrderCheck_;
    QCheckBox* showButterflyCheck_;
    QCheckBox* showFieldGridCheck_;
    QCheckBox* showDeshieldedCheck_;
    QSlider* isoThresholdSlider_;
    QLabel* isoThresholdLabel_;
    QLabel* statusLabel_;

    // Atom info panel — shows full object model for picked atom
    QDockWidget* atomInfoDock_;
    QTreeWidget* atomInfoTree_;

    // Bond tab — bonds for picked atom
    QDockWidget* bondInfoDock_;
    QTreeWidget* bondInfoTree_;

    // GeometryChoice tab — calculator decisions for picked atom
    QDockWidget* gcDock_;
    QTreeWidget* gcTree_;

    // Time Series tab — frame-0 values from the analysis H5 for picked atom
    QDockWidget* timeSeriesDock_ = nullptr;
    QTreeWidget* timeSeriesTree_ = nullptr;

    // Menu actions
    QAction* exportFeaturesAct_ = nullptr;

    // Bond order color overlay (tubes colored by Wiberg order)
    vtkSmartPointer<vtkActor> bondOrderActor_;

    // Operations log panel — shows library log stream via UDP.
    // udpHost_ / udpPort_ are set from main_viewer.cpp at construction
    // (sourced from [logging] in ~/.nmr_tools.toml). Multicast hosts
    // (239.0.0.0/8) trigger an AnyIPv4 bind + joinMulticastGroup so
    // udp_listen.py and the Log dock co-listen; unicast keeps the
    // host-specific bind.
    QString udpHost_;
    quint16 udpPort_ = 0;
    QDockWidget* logDock_;
    QPlainTextEdit* logText_;
    QUdpSocket* logSocket_;
    void onLogDatagramReady();

    // Pending load state (set by loadFromJobSpec, consumed by startCompute)
    nmr::JobSpec pendingSpec_;
};
