#pragma once
#include <QMainWindow>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkOpenGLMoleculeMapper.h>
#include <vtkMolecule.h>
#include "ComputeWorker.h"

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

    // Load a PDB directly (for command-line usage)
    void loadPdb(const std::string& pdbPath);
    void loadProteinDir(const std::string& dirPath);

signals:
    void computeRequested(ViewerLoadRequest request);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

private slots:
    void saveScreenshot();
    void onRenderModeChanged(int index);
    void onOverlayModeChanged(int index);
    void onGlyphStyleChanged(int index);
    void onGlyphScaleChanged(int value);
    void onOpacityChanged(int value);
    void onCurrentScaleChanged(int value);
    void onShowRingsToggled(bool checked);
    void onShowPeptideBondsToggled(bool checked);
    void onShowBondOrderToggled(bool checked);
    void onShowButterflyToggled(bool checked);
    void onPhysicsCheckChanged();
    void onVizModeChanged(int index);
    void onIsoThresholdChanged(double value);
    void onGaussianRadiusChanged(double value);

    // Async compute slots
    void onComputeProgress(int current, int total, QString phase);
    void onComputeFinished(ComputeResult result);

private:
    QString initialDir_;
    void setupUI();
    void setupMenuBar();
    void loadMolecule(const std::string& pdbPath);
    void startCompute();
    void cancelCompute();
    void updateOverlay();

    // Atom picking and inspection
    void pickAtom(int displayX, int displayY);
    void populateAtomInfo(size_t atomIndex);

    // Bond picking and inspection
    void pickBond(int displayX, int displayY);
    void populateBondInfo(size_t bondIndex);

    // Helpers: sum checked calculator contributions for one atom (by index)
    double checkedCalcT0(size_t atomIndex) const;
    nmr::SphericalTensor checkedCalcST(size_t atomIndex) const;

    // VTK rendering
    QVTKOpenGLNativeWidget* vtkWidget_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkMolecule> molecule_;
    vtkSmartPointer<vtkOpenGLMoleculeMapper> molMapper_;
    vtkSmartPointer<vtkActor> molActor_;

    // Overlays
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

    // Async computation
    QThread* workerThread_ = nullptr;
    ComputeWorker* worker_ = nullptr;
    QProgressDialog* progressDialog_ = nullptr;

    // Slider debounce
    QTimer* sliderDebounce_ = nullptr;

    // UI Controls
    QComboBox* renderModeCombo_;
    QComboBox* overlayModeCombo_;
    QComboBox* glyphStyleCombo_;
    QSlider* glyphScaleSlider_;
    QSlider* opacitySlider_;
    QSlider* currentScaleSlider_;
    QCheckBox* showRingsCheck_;
    QCheckBox* showPeptideBondsCheck_;
    QCheckBox* showBondOrderCheck_;
    QCheckBox* showButterflyCheck_;

    // Physics contribution checkboxes (8 calculators)
    QCheckBox* physicsChecks_[8];  // BS, HM, MC, LD, CE, PQ, RSA, HB
    QComboBox* vizModeCombo_;
    QDoubleSpinBox* isoThreshold_;
    QDoubleSpinBox* gaussianRadius_;
    QLabel* statusLabel_;

    // Atom info panel — shows full object model for picked atom
    QDockWidget* atomInfoDock_;
    QTreeWidget* atomInfoTree_;

    // Bond info panel — shows full bond data including MOPAC order
    QDockWidget* bondInfoDock_;
    QTreeWidget* bondInfoTree_;

    // Bond order color overlay (tubes colored by Wiberg order)
    vtkSmartPointer<vtkActor> bondOrderActor_;

    // Operations log panel — shows library log stream via UDP
    QDockWidget* logDock_;
    QPlainTextEdit* logText_;
    QUdpSocket* logSocket_;
    void onLogDatagramReady();

    // Pending load state (set by loadPdb/loadProteinDir, consumed by startCompute)
    ViewerLoadRequest pendingRequest_;
};
