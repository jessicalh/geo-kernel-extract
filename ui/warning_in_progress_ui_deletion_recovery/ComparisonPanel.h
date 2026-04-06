#pragma once
#include "ComputeWorker.h"
#include <QWidget>
#include <QTableWidget>
#include <QLabel>
#include <vector>

class ScatterChartWidget;

// Side panel showing DFT vs predicted comparison statistics
// Includes a scatter chart, table, and summary R^2/RMSE values

class ComparisonPanel : public QWidget {
    Q_OBJECT
public:
    explicit ComparisonPanel(QWidget* parent = nullptr);

    // Update from ViewerAtomResults (atoms that have hasDft=true are compared)
    void update(const std::vector<ViewerAtomResult>& atoms);

    void clear();

private:
    ScatterChartWidget* chart_;
    QTableWidget* table_;
    QLabel* statsLabel_;
};
