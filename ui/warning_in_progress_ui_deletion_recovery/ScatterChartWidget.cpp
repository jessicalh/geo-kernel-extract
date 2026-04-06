#include "ScatterChartWidget.h"
#include <QVBoxLayout>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkContextView.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkPlotPoints.h>
#include <vtkPlotLine.h>
#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkAxis.h>
#include <vtkPen.h>
#include <vtkContextScene.h>
#include <vtkNew.h>
#include <cmath>
#include <algorithm>

ScatterChartWidget::ScatterChartWidget(QWidget* parent)
    : QWidget(parent)
{
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    vtkWidget_ = new QVTKOpenGLNativeWidget(this);
    vtkWidget_->setMinimumHeight(250);
    vtkWidget_->setMaximumHeight(300);
    layout->addWidget(vtkWidget_);

    renderWindow_ = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWidget_->setRenderWindow(renderWindow_);

    view_ = vtkSmartPointer<vtkContextView>::New();
    view_->SetRenderWindow(renderWindow_);

    chart_ = vtkSmartPointer<vtkChartXY>::New();
    view_->GetScene()->AddItem(chart_);

    chart_->SetTitle("DFT vs Predicted");
    chart_->GetAxis(vtkAxis::BOTTOM)->SetTitle("DFT delta_iso (ppm)");
    chart_->GetAxis(vtkAxis::LEFT)->SetTitle("Predicted sigma_iso (ppm)");
    chart_->SetShowLegend(true);
}

void ScatterChartWidget::clear() {
    chart_->ClearPlots();
    renderWindow_->Render();
}

void ScatterChartWidget::setData(
    const std::vector<std::pair<double,double>>& bsPairs,
    const std::vector<std::pair<double,double>>& mlPairs)
{
    chart_->ClearPlots();

    double minVal = 1e30, maxVal = -1e30;

    // Biot-Savart scatter (blue circles)
    if (!bsPairs.empty()) {
        vtkNew<vtkTable> table;
        vtkNew<vtkFloatArray> colDFT;
        colDFT->SetName("DFT (BS)");
        vtkNew<vtkFloatArray> colBS;
        colBS->SetName("Biot-Savart");

        for (const auto& [dft, bs] : bsPairs) {
            colDFT->InsertNextValue(static_cast<float>(dft));
            colBS->InsertNextValue(static_cast<float>(bs));
            minVal = std::min({minVal, dft, bs});
            maxVal = std::max({maxVal, dft, bs});
        }

        table->AddColumn(colDFT);
        table->AddColumn(colBS);

        auto* plot = chart_->AddPlot(vtkChart::POINTS);
        plot->SetInputData(table, 0, 1);
        plot->SetColor(51, 102, 204, 180);
        plot->SetWidth(4.0f);
        plot->SetLabel("BS");
    }

    // ML scatter (red circles)
    if (!mlPairs.empty()) {
        vtkNew<vtkTable> table;
        vtkNew<vtkFloatArray> colDFT;
        colDFT->SetName("DFT (ML)");
        vtkNew<vtkFloatArray> colML;
        colML->SetName("ML");

        for (const auto& [dft, ml] : mlPairs) {
            colDFT->InsertNextValue(static_cast<float>(dft));
            colML->InsertNextValue(static_cast<float>(ml));
            minVal = std::min({minVal, dft, ml});
            maxVal = std::max({maxVal, dft, ml});
        }

        table->AddColumn(colDFT);
        table->AddColumn(colML);

        auto* plot = chart_->AddPlot(vtkChart::POINTS);
        plot->SetInputData(table, 0, 1);
        plot->SetColor(204, 51, 51, 180);
        plot->SetWidth(4.0f);
        plot->SetLabel("ML");
    }

    // y=x reference line (gray dashed)
    if (minVal < maxVal) {
        double margin = (maxVal - minVal) * 0.05;
        double lo = minVal - margin;
        double hi = maxVal + margin;

        vtkNew<vtkTable> refTable;
        vtkNew<vtkFloatArray> colX;
        colX->SetName("RefX");
        vtkNew<vtkFloatArray> colY;
        colY->SetName("RefY");
        colX->InsertNextValue(static_cast<float>(lo));
        colY->InsertNextValue(static_cast<float>(lo));
        colX->InsertNextValue(static_cast<float>(hi));
        colY->InsertNextValue(static_cast<float>(hi));
        refTable->AddColumn(colX);
        refTable->AddColumn(colY);

        auto* line = chart_->AddPlot(vtkChart::LINE);
        line->SetInputData(refTable, 0, 1);
        line->SetColor(128, 128, 128, 128);
        line->SetWidth(1.5f);
        line->GetPen()->SetLineType(vtkPen::DASH_LINE);
        line->SetLabel("y=x");
    }

    renderWindow_->Render();
}
