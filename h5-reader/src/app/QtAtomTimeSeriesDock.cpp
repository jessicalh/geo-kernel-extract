#include "QtAtomTimeSeriesDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../io/QtNamingRegistry.h"
#include "../model/QtFrame.h"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLoggingCategory>
#include <QPen>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::app {

using model::QtFrame;
using model::SphericalTensor;
using model::Vec3;

namespace {
Q_LOGGING_CATEGORY(cTs, "h5reader.timeseries")

// ---------------------------------------------------------------------------
// Scalar catalogue. Each entry: category group, display label, unit, and a
// function pointer that reads the scalar at (frame, atomIdx). Non-capturing
// lambdas decay to plain function pointers, which keeps the table POD.
// ---------------------------------------------------------------------------

using ExtractFn = double(*)(const QtFrame&, std::size_t);

struct Desc {
    const char* group;
    const char* label;
    const char* unit;
    ExtractFn   extract;
};

double T2Mag(const SphericalTensor& st) {
    double s = 0;
    for (double c : st.T2) s += c * c;
    return std::sqrt(s);
}

const Desc& descAt(int i);
int         descCount();

const std::vector<Desc>& descs() {
    static const std::vector<Desc> kDescs = {
        // Position
        {"Position", "x",       "Å",   [](const QtFrame& f, std::size_t a){ return f.position(a).x(); }},
        {"Position", "y",       "Å",   [](const QtFrame& f, std::size_t a){ return f.position(a).y(); }},
        {"Position", "z",       "Å",   [](const QtFrame& f, std::size_t a){ return f.position(a).z(); }},

        // Ring current
        {"Ring current", "bs_shielding T0",      "ppm",  [](const QtFrame& f, std::size_t a){ return f.bsShielding(a).T0; }},
        {"Ring current", "bs_shielding |T2|",    "ppm",  [](const QtFrame& f, std::size_t a){ return T2Mag(f.bsShielding(a)); }},
        {"Ring current", "hm_shielding T0",      "ppm",  [](const QtFrame& f, std::size_t a){ return f.hmShielding(a).T0; }},
        {"Ring current", "hm_shielding |T2|",    "ppm",  [](const QtFrame& f, std::size_t a){ return T2Mag(f.hmShielding(a)); }},
        {"Ring current", "rs_shielding T0",      "ppm",  [](const QtFrame& f, std::size_t a){ return f.rsShielding(a).T0; }},
        {"Ring current", "|B| total",            "T",    [](const QtFrame& f, std::size_t a){ return f.totalBField(a).norm(); }},
        {"Ring current", "n rings ≤ 3 Å",        "",     [](const QtFrame& f, std::size_t a){ return double(f.nRings3A(a)); }},
        {"Ring current", "n rings ≤ 5 Å",        "",     [](const QtFrame& f, std::size_t a){ return double(f.nRings5A(a)); }},
        {"Ring current", "n rings ≤ 8 Å",        "",     [](const QtFrame& f, std::size_t a){ return double(f.nRings8A(a)); }},
        {"Ring current", "mean ring dist",       "Å",    [](const QtFrame& f, std::size_t a){ return f.meanRingDist(a); }},
        {"Ring current", "nearest ring atom",    "Å",    [](const QtFrame& f, std::size_t a){ return f.nearestRingAtom(a); }},

        // Bond anisotropy
        {"Bond anisotropy", "mc_shielding T0",   "Å⁻³",  [](const QtFrame& f, std::size_t a){ return f.mcShielding(a).T0; }},
        {"Bond anisotropy", "mc_shielding |T2|", "Å⁻³",  [](const QtFrame& f, std::size_t a){ return T2Mag(f.mcShielding(a)); }},
        {"Bond anisotropy", "co_sum",            "",     [](const QtFrame& f, std::size_t a){ return f.mcCOSum(a); }},
        {"Bond anisotropy", "nearest C=O dist",  "Å",    [](const QtFrame& f, std::size_t a){ return f.mcNearestCODist(a); }},

        // Quadrupole / Dispersion
        {"Quadrupole",  "pq_shielding T0",       "",     [](const QtFrame& f, std::size_t a){ return f.pqShielding(a).T0; }},
        {"Quadrupole",  "pq_shielding |T2|",     "",     [](const QtFrame& f, std::size_t a){ return T2Mag(f.pqShielding(a)); }},
        {"Dispersion",  "disp_shielding T0",     "",     [](const QtFrame& f, std::size_t a){ return f.dispShielding(a).T0; }},
        {"Dispersion",  "disp_shielding |T2|",   "",     [](const QtFrame& f, std::size_t a){ return T2Mag(f.dispShielding(a)); }},

        // Electrostatics
        {"Electrostatics", "coulomb_shielding T0", "ppm",  [](const QtFrame& f, std::size_t a){ return f.coulombShielding(a).T0; }},
        {"Electrostatics", "|E| Coulomb",          "V/Å",  [](const QtFrame& f, std::size_t a){ return f.coulombETotal(a).norm(); }},
        {"Electrostatics", "|E| magnitude field",  "V/Å",  [](const QtFrame& f, std::size_t a){ return f.coulombEMagnitude(a); }},
        {"Electrostatics", "APBS EFG T0",          "V/Å²", [](const QtFrame& f, std::size_t a){ return f.apbsEfg(a).T0; }},
        {"Electrostatics", "|E| APBS",             "V/Å",  [](const QtFrame& f, std::size_t a){ return f.apbsEfield(a).norm(); }},
        {"Electrostatics", "AIMNet2 shielding T0", "ppm",  [](const QtFrame& f, std::size_t a){ return f.aimnet2Shielding(a).T0; }},

        // H-bond
        {"H-bond", "hbond_shielding T0",   "",    [](const QtFrame& f, std::size_t a){ return f.hbondShielding(a).T0; }},
        {"H-bond", "nearest dist",         "Å",   [](const QtFrame& f, std::size_t a){ return f.hbondNearestDist(a); }},
        {"H-bond", "count ≤ 3.5 Å",        "",    [](const QtFrame& f, std::size_t a){ return double(f.hbondCount35A(a)); }},

        // SASA
        {"SASA", "sasa",                   "Å²",  [](const QtFrame& f, std::size_t a){ return f.sasa(a); }},

        // Water environment
        {"Water", "|E| water",             "V/Å", [](const QtFrame& f, std::size_t a){ return f.waterEfield(a).norm(); }},
        {"Water", "n first shell",         "",    [](const QtFrame& f, std::size_t a){ return double(f.waterNFirst(a)); }},
        {"Water", "n second shell",        "",    [](const QtFrame& f, std::size_t a){ return double(f.waterNSecond(a)); }},
        {"Water", "half-shell asymmetry",  "",    [](const QtFrame& f, std::size_t a){ return f.waterHalfShellAsymmetry(a); }},
        {"Water", "dipole cos",            "",    [](const QtFrame& f, std::size_t a){ return f.waterDipoleCos(a); }},

        // Charges
        {"Charges", "AIMNet2 charge",      "e",   [](const QtFrame& f, std::size_t a){ return f.aimnet2Charge(a); }},
        {"Charges", "EEQ charge",          "e",   [](const QtFrame& f, std::size_t a){ return f.eeqCharge(a); }},
        {"Charges", "EEQ coord number",    "",    [](const QtFrame& f, std::size_t a){ return f.eeqCoordinationNumber(a); }},
    };
    return kDescs;
}

const Desc& descAt(int i) { return descs()[i]; }
int         descCount()   { return static_cast<int>(descs().size()); }

QString comboLabel(const Desc& d) {
    QString s = QStringLiteral("%1: %2").arg(QString::fromLatin1(d.group),
                                              QString::fromLatin1(d.label));
    if (d.unit && d.unit[0] != '\0') {
        s += QStringLiteral("  (%1)").arg(QString::fromLatin1(d.unit));
    }
    return s;
}

}  // namespace


QtAtomTimeSeriesDock::QtAtomTimeSeriesDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Time Series"), parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomTimeSeriesDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    auto* container = new QWidget(this);
    auto* vbox      = new QVBoxLayout(container);
    vbox->setContentsMargins(4, 4, 4, 4);
    vbox->setSpacing(4);

    // Scalar picker — grouped labels in one combo. Future work: a two-
    // tier combo or a tree if this gets unwieldy.
    scalarCombo_ = new QComboBox(container);
    for (int i = 0; i < descCount(); ++i) {
        scalarCombo_->addItem(comboLabel(descAt(i)), i);
    }
    auto* top = new QHBoxLayout;
    top->addWidget(new QLabel(QStringLiteral("Scalar:"), container));
    top->addWidget(scalarCombo_, 1);
    vbox->addLayout(top);

    // Chart setup.
    chart_       = new QChart();
    chart_->legend()->hide();
    chart_->setBackgroundRoundness(0);
    chart_->setMargins(QMargins(4, 4, 4, 4));

    dataSeries_ = new QLineSeries(chart_);
    dataSeries_->setName(QStringLiteral("data"));
    chart_->addSeries(dataSeries_);

    cursorSeries_ = new QLineSeries(chart_);
    cursorSeries_->setName(QStringLiteral("cursor"));
    QPen cursorPen(Qt::red);
    cursorPen.setWidthF(1.5);
    cursorPen.setStyle(Qt::DashLine);
    cursorSeries_->setPen(cursorPen);
    chart_->addSeries(cursorSeries_);

    xAxis_ = new QValueAxis(chart_);
    xAxis_->setTitleText(QStringLiteral("frame"));
    xAxis_->setLabelFormat(QStringLiteral("%d"));
    chart_->addAxis(xAxis_, Qt::AlignBottom);

    yAxis_ = new QValueAxis(chart_);
    yAxis_->setTitleText(QStringLiteral("value"));
    chart_->addAxis(yAxis_, Qt::AlignLeft);

    dataSeries_->attachAxis(xAxis_);
    dataSeries_->attachAxis(yAxis_);
    cursorSeries_->attachAxis(xAxis_);
    cursorSeries_->attachAxis(yAxis_);

    chartView_ = new QChartView(chart_, container);
    chartView_->setRenderHint(QPainter::Antialiasing);
    vbox->addWidget(chartView_, 1);

    valueLabel_ = new QLabel(QStringLiteral("—"), container);
    valueLabel_->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    vbox->addWidget(valueLabel_);

    setWidget(container);

    ACONNECT(scalarCombo_.data(), qOverload<int>(&QComboBox::currentIndexChanged),
             this, &QtAtomTimeSeriesDock::onScalarChanged);
}

void QtAtomTimeSeriesDock::setContext(const model::QtProtein*      protein,
                                       const model::QtConformation* conformation) {
    protein_      = protein;
    conformation_ = conformation;
    if (xAxis_ && conformation_) {
        xAxis_->setRange(0, std::max<int>(0, conformation_->frameCount() - 1));
    }
}

void QtAtomTimeSeriesDock::setPickedAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    hasSelection_ = true;
    atomIdx_      = atomIdx;
    rebuildSeries();
    updateCursor();
}

void QtAtomTimeSeriesDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    updateCursor();
}

void QtAtomTimeSeriesDock::clearSelection() {
    ASSERT_THREAD(this);
    hasSelection_ = false;
    dataSeries_->clear();
    cursorSeries_->clear();
    if (valueLabel_) valueLabel_->setText(QStringLiteral("—"));
}

void QtAtomTimeSeriesDock::onScalarChanged(int /*comboIndex*/) {
    if (hasSelection_) rebuildSeries();
    updateCursor();
}

void QtAtomTimeSeriesDock::rebuildSeries() {
    if (!dataSeries_ || !yAxis_ || !xAxis_) return;
    if (!protein_ || !conformation_ || !hasSelection_) return;
    if (atomIdx_ >= protein_->atomCount()) return;

    const int comboIdx = scalarCombo_ ? scalarCombo_->currentIndex() : 0;
    if (comboIdx < 0 || comboIdx >= descCount()) return;
    const Desc& d = descAt(comboIdx);

    const int T = static_cast<int>(conformation_->frameCount());
    QList<QPointF> pts;
    pts.reserve(T);

    double yMin = +std::numeric_limits<double>::infinity();
    double yMax = -std::numeric_limits<double>::infinity();
    for (int t = 0; t < T; ++t) {
        const auto& f = conformation_->frame(static_cast<std::size_t>(t));
        const double v = d.extract(f, atomIdx_);
        if (std::isfinite(v)) {
            pts.append(QPointF(t, v));
            if (v < yMin) yMin = v;
            if (v > yMax) yMax = v;
        }
    }
    dataSeries_->replace(pts);

    if (!std::isfinite(yMin) || !std::isfinite(yMax)) {
        yMin = 0.0; yMax = 1.0;
    }
    if (yMax - yMin < 1e-9) {
        // Flat series — pad the axis so the line is still visible.
        const double centre = yMin;
        yMin = centre - 1.0;
        yMax = centre + 1.0;
    } else {
        const double pad = 0.05 * (yMax - yMin);
        yMin -= pad; yMax += pad;
    }
    yAxis_->setRange(yMin, yMax);
    xAxis_->setRange(0, std::max(0, T - 1));

    yAxis_->setTitleText(QString::fromLatin1(d.unit ? d.unit : "value"));
    chart_->setTitle(QStringLiteral("%1  —  atom %2")
                       .arg(QString::fromLatin1(d.label))
                       .arg(atomIdx_));

    qCDebug(cTs).noquote()
        << "rebuild | atom=" << atomIdx_
        << "| scalar=" << d.label
        << "| n=" << pts.size()
        << "| y=[" << yMin << "," << yMax << "]";
}

void QtAtomTimeSeriesDock::updateCursor() {
    if (!cursorSeries_ || !yAxis_) return;

    const int T = conformation_ ? static_cast<int>(conformation_->frameCount()) : 0;
    const int t = std::clamp(frame_, 0, std::max(0, T - 1));

    const double yMin = yAxis_->min();
    const double yMax = yAxis_->max();
    QList<QPointF> cur = { QPointF(t, yMin), QPointF(t, yMax) };
    cursorSeries_->replace(cur);

    // Current-value readout — read at the cursor frame for the selected
    // scalar, if a selection is active.
    if (hasSelection_ && conformation_ && scalarCombo_ && valueLabel_) {
        const int comboIdx = scalarCombo_->currentIndex();
        if (comboIdx >= 0 && comboIdx < descCount()) {
            const Desc& d = descAt(comboIdx);
            const auto& f = conformation_->frame(static_cast<std::size_t>(t));
            const double v = d.extract(f, atomIdx_);
            valueLabel_->setText(QStringLiteral("frame %1: %2 %3")
                .arg(t)
                .arg(v, 0, 'g', 6)
                .arg(QString::fromLatin1(d.unit ? d.unit : "")));
        }
    }
}

}  // namespace h5reader::app
