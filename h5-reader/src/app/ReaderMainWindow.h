// ReaderMainWindow — top-level QMainWindow for h5-reader.
//
// Holds the QVTKOpenGLNativeWidget as the central widget, owns the
// MoleculeScene and QtPlaybackController, wires frame-scrub controls
// in a toolbar and status bar. Designed to accept additional dock
// widgets (atom inspector, time-series tab) in later commits without
// restructuring the central layout.
//
// Shutdown protocol — see feedback_qt_discipline and the library
// viewer's MainWindow::shutdown(). When QApplication is about to
// quit, stop all timers, drop VTK references in order, and
// renderWindow_->Finalize() BEFORE Qt destroys the GL context.

#pragma once

#include <QMainWindow>
#include <QPointer>

#include <vtkSmartPointer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <memory>

class QLabel;
class QSlider;
class QSpinBox;
class QToolBar;
class QVTKOpenGLNativeWidget;

namespace h5reader::io {
struct QtLoadResult;
}

namespace h5reader::app {

class MoleculeScene;
class QtPlaybackController;

class ReaderMainWindow final : public QMainWindow {
    Q_OBJECT

public:
    // Takes the loader's result by rvalue — the window assumes ownership
    // of the protein, conformation, and AnalysisFile. The result must
    // have ok=true; otherwise the caller should have aborted already.
    explicit ReaderMainWindow(h5reader::io::QtLoadResult&& loaded,
                              QWidget* parent = nullptr);
    ~ReaderMainWindow() override;

public slots:
    // Called from aboutToQuit. Stops timers, cancels any workers, and
    // finalises the VTK render window before Qt tears down the GL context.
    void shutdown();

private slots:
    void onFrameChanged(int t);
    void onPlayPauseClicked();
    void onFpsChanged(int fps);

private:
    void buildUi();
    void buildToolbar();
    void buildStatusBar();

    // The loaded model. Owned by the window for its lifetime.
    std::unique_ptr<h5reader::io::QtLoadResult> loaded_;

    // VTK viewport widget.
    QVTKOpenGLNativeWidget*                       vtkWidget_   = nullptr;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;

    // Scene + playback.
    MoleculeScene*        scene_    = nullptr;
    QtPlaybackController* playback_ = nullptr;

    // Atom picker + inspector dock. Picker is an event filter on the
    // VTK widget; inspector is a tabified QDockWidget on the right.
    class QtAtomPicker*          picker_          = nullptr;
    class QtAtomInspectorDock*   inspectorDock_   = nullptr;
    class QtAtomTimeSeriesDock*  timeSeriesDock_  = nullptr;

    // Toolbar controls.
    QPointer<QSlider>  frameSlider_;
    QPointer<QSpinBox> fpsSpinner_;
    QPointer<QAction>  playAction_;
    QPointer<QAction>  showRibbonAction_;
    QPointer<QAction>  showRingsAction_;
    QPointer<QAction>  showButterflyAction_;
    QPointer<QAction>  showBFieldAction_;

    // Status bar labels.
    QPointer<QLabel>   proteinLabel_;
    QPointer<QLabel>   frameLabel_;
    QPointer<QLabel>   timeLabel_;

    bool shutdownDone_ = false;
};

}  // namespace h5reader::app
