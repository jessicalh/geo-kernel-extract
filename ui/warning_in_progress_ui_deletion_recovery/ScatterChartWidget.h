#pragma once
#include <QWidget>
#include <vtkSmartPointer.h>
#include <vector>
#include <utility>

class QVTKOpenGLNativeWidget;
class vtkChartXY;
class vtkContextView;
class vtkGenericOpenGLRenderWindow;

// Embedded VTK 2D scatter chart: DFT vs predicted (BS and/or ML)
// Shows y=x reference line, colored series, axis labels

class ScatterChartWidget : public QWidget {
    Q_OBJECT
public:
    explicit ScatterChartWidget(QWidget* parent = nullptr);

    void setData(const std::vector<std::pair<double,double>>& bsPairs,
                 const std::vector<std::pair<double,double>>& mlPairs = {});
    void clear();

private:
    QVTKOpenGLNativeWidget* vtkWidget_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkContextView> view_;
    vtkSmartPointer<vtkChartXY> chart_;
};
