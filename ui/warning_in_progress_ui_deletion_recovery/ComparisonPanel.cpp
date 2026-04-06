#include "ComparisonPanel.h"
#include "ScatterChartWidget.h"
#include <QVBoxLayout>
#include <QHeaderView>
#include <cmath>
#include <algorithm>

ComparisonPanel::ComparisonPanel(QWidget* parent)
    : QWidget(parent)
{
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    statsLabel_ = new QLabel("No data loaded");
    statsLabel_->setWordWrap(true);
    layout->addWidget(statsLabel_);

    chart_ = new ScatterChartWidget(this);
    layout->addWidget(chart_);

    table_ = new QTableWidget(this);
    table_->setColumnCount(6);
    table_->setHorizontalHeaderLabels({"Atom", "Elem", "DFT iso", "BS iso", "ML iso", "Resid(BS)"});
    table_->horizontalHeader()->setStretchLastSection(true);
    table_->setMaximumHeight(200);
    layout->addWidget(table_);
}

void ComparisonPanel::clear() {
    table_->setRowCount(0);
    chart_->clear();
    statsLabel_->setText("No data loaded");
}

void ComparisonPanel::update(const std::vector<ViewerAtomResult>& atoms)
{
    // Collect atoms that have DFT comparison data
    struct MatchedAtom {
        std::string element;
        std::string atomName;
        double dftT0;
        double bsT0;
        double totalT0;
    };
    std::vector<MatchedAtom> matches;

    for (const auto& a : atoms) {
        if (!a.hasDft) continue;
        matches.push_back({a.element, a.atomName, a.dftT0, a.bsT0, a.totalT0});
    }

    if (matches.empty()) {
        statsLabel_->setText("No DFT comparison atoms");
        return;
    }

    // Compute BS statistics (predicted BS T0 vs DFT T0)
    double ss_res_bs = 0.0, ss_tot = 0.0, mean_dft = 0.0;
    int n_H = 0;
    double ss_res_bs_H = 0.0, ss_tot_H = 0.0, mean_dft_H = 0.0;

    for (const auto& m : matches) {
        mean_dft += m.dftT0;
        if (m.element == "H") { mean_dft_H += m.dftT0; n_H++; }
    }
    mean_dft /= matches.size();
    if (n_H > 0) mean_dft_H /= n_H;

    for (const auto& m : matches) {
        double res_bs = m.bsT0 - m.dftT0;
        ss_res_bs += res_bs * res_bs;
        ss_tot += (m.dftT0 - mean_dft) * (m.dftT0 - mean_dft);
        if (m.element == "H") {
            ss_res_bs_H += res_bs * res_bs;
            ss_tot_H += (m.dftT0 - mean_dft_H) * (m.dftT0 - mean_dft_H);
        }
    }

    double rmse_bs = std::sqrt(ss_res_bs / matches.size());
    double r2_bs = (ss_tot > 1e-30) ? (1.0 - ss_res_bs / ss_tot) : 0.0;
    double rmse_bs_H = (n_H > 0) ? std::sqrt(ss_res_bs_H / n_H) : 0.0;
    double r2_bs_H = (ss_tot_H > 1e-30) ? (1.0 - ss_res_bs_H / ss_tot_H) : 0.0;

    QString statsText = QString(
        "BS: RMSE=%1 ppm, R2=%2\n"
        "BS(H, %3): RMSE=%4 ppm, R2=%5")
        .arg(rmse_bs, 0, 'f', 3).arg(r2_bs, 0, 'f', 4)
        .arg(n_H)
        .arg(rmse_bs_H, 0, 'f', 3).arg(r2_bs_H, 0, 'f', 4);

    // Total predicted statistics
    {
        double ss_res_total = 0.0, ss_res_total_H = 0.0;
        for (const auto& m : matches) {
            double res = m.totalT0 - m.dftT0;
            ss_res_total += res * res;
            if (m.element == "H")
                ss_res_total_H += res * res;
        }
        double rmse_total = std::sqrt(ss_res_total / matches.size());
        double r2_total = (ss_tot > 1e-30) ? (1.0 - ss_res_total / ss_tot) : 0.0;
        double rmse_total_H = (n_H > 0) ? std::sqrt(ss_res_total_H / n_H) : 0.0;
        double r2_total_H = (ss_tot_H > 1e-30) ? (1.0 - ss_res_total_H / ss_tot_H) : 0.0;
        statsText += QString(
            "\nTotal: RMSE=%1 ppm, R2=%2\n"
            "Total(H, %3): RMSE=%4 ppm, R2=%5")
            .arg(rmse_total, 0, 'f', 3).arg(r2_total, 0, 'f', 4)
            .arg(n_H)
            .arg(rmse_total_H, 0, 'f', 3).arg(r2_total_H, 0, 'f', 4);
    }

    statsLabel_->setText(statsText);

    // Build scatter chart data: DFT vs BS and DFT vs Total
    std::vector<std::pair<double,double>> bsPairs, totalPairs;
    for (const auto& m : matches) {
        bsPairs.push_back({m.dftT0, m.bsT0});
        totalPairs.push_back({m.dftT0, m.totalT0});
    }
    chart_->setData(bsPairs, totalPairs);

    // Sort by |BS residual| descending
    std::sort(matches.begin(), matches.end(), [](const auto& a, const auto& b) {
        return std::abs(a.bsT0 - a.dftT0) > std::abs(b.bsT0 - b.dftT0);
    });

    // Populate table (top 100)
    int nShow = std::min(static_cast<int>(matches.size()), 100);
    table_->setRowCount(nShow);
    for (int i = 0; i < nShow; ++i) {
        const auto& m = matches[i];
        table_->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(m.atomName)));
        table_->setItem(i, 1, new QTableWidgetItem(QString::fromStdString(m.element)));
        table_->setItem(i, 2, new QTableWidgetItem(QString::number(m.dftT0, 'f', 3)));
        table_->setItem(i, 3, new QTableWidgetItem(QString::number(m.bsT0, 'f', 3)));
        table_->setItem(i, 4, new QTableWidgetItem(QString::number(m.totalT0, 'f', 3)));
        table_->setItem(i, 5, new QTableWidgetItem(
            QString::number(m.bsT0 - m.dftT0, 'f', 3)));
    }
}
