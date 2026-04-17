#pragma once
#include <QMainWindow>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkOpenGLMoleculeMapper.h>
#include <vtkMolecule.h>
#include "ComputeWorker.h"
#include "analysis_file.h"  // read-only time-series companion data

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
}

class RestServer;

class MainWindow : public QMainWindow {
    Q_OBJECT
    friend class RestServer;  // REST server needs access to UI controls and VTK renderer
public:
    explicit MainWindow(const QString& initialDir = QString(), QWidget* parent = nullptr);
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
    void onComputeProgress(int current, int total, QString phase);
    void onComputeFinished(ComputeResult result);

private:
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

    // Time Series tab — per-atom, frame-0 slice of the companion analysis H5.
    // Populates only when analysisFile_ is set AND atomIndex is in range.
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

    // Operations log panel — shows library log stream via UDP
    QDockWidget* logDock_;
    QPlainTextEdit* logText_;
    QUdpSocket* logSocket_;
    void onLogDatagramReady();

    // Pending load state (set by loadFromJobSpec, consumed by startCompute)
    nmr::JobSpec pendingSpec_;
};
