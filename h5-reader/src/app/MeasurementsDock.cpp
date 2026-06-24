#include "MeasurementsDock.h"

#include "../diagnostics/ObjectCensus.h"

#include "../model/AtomSelection.h"         // GeometryKind, NameForGeometryKind
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"  // GeometryMeasurement, Measure
#include "../model/QtProtein.h"
#include "../model/QtResidueNames.h"

#include <QChar>
#include <QFont>
#include <QLabel>
#include <QLatin1Char>
#include <QString>
#include <QStringList>
#include <QVBoxLayout>
#include <QWidget>

namespace h5reader::app {

namespace {

// Non-ASCII glyphs via code points so the source stays ASCII (house rule).
QString degreeGlyph() { return QString(QChar(char16_t(0x00B0))); }    // degree sign
QString angstromGlyph() { return QString(QChar(char16_t(0x00C5))); }  // A with ring above

QString atomLabel(const model::QtProtein& protein, std::size_t atom) {
    if (atom >= protein.atomCount())
        return QStringLiteral("atom %1").arg(static_cast<qulonglong>(atom));
    const QString name = protein.atomNames(atom).iupac;
    const auto& a = protein.atom(atom);
    if (a.residueIndex < 0 || static_cast<std::size_t>(a.residueIndex) >= protein.residueCount())
        return name;
    const auto& res = protein.residue(static_cast<std::size_t>(a.residueIndex));
    const QString resName = QString::fromLatin1(model::IupacResidue3LetterFor(res.aminoAcid));
    const QString chain = res.address.chainId.isEmpty()
                              ? QString()
                              : res.address.chainId + QLatin1Char(':');
    return QStringLiteral("%1%2%3:%4")
        .arg(chain, resName)
        .arg(res.address.residueNumber)
        .arg(name);
}

}  // namespace

MeasurementsDock::MeasurementsDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Measurements"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("MeasurementsDock"));

    auto* body = new QWidget(this);
    auto* layout = new QVBoxLayout(body);
    layout->setContentsMargins(10, 8, 10, 8);
    layout->setSpacing(4);

    kindLabel_ = new QLabel(body);
    kindLabel_->setWordWrap(true);
    QFont kf = kindLabel_->font();
    kf.setBold(true);
    kindLabel_->setFont(kf);

    valueLabel_ = new QLabel(body);
    QFont vf = valueLabel_->font();
    vf.setPointSize(vf.pointSize() + 8);
    vf.setBold(true);
    valueLabel_->setFont(vf);
    valueLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);

    atomsLabel_ = new QLabel(body);
    atomsLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    atomsLabel_->setWordWrap(true);

    layout->addWidget(kindLabel_);
    layout->addWidget(valueLabel_);
    layout->addWidget(atomsLabel_);
    layout->addStretch(1);

    setWidget(body);
    clear();
}

void MeasurementsDock::setContext(const model::QtProtein* protein,
                                  const model::Conformation* conf) {
    protein_ = protein;
    conf_ = conf;
    if (!protein_ || !conf_) {
        atoms_.clear();
        frame_ = 0;
    }
    recompute();
}

void MeasurementsDock::setAtoms(const std::vector<std::size_t>& atoms) {
    atoms_ = atoms;
    recompute();
}

void MeasurementsDock::setFrame(int frame) {
    frame_ = frame < 0 ? 0 : frame;
    recompute();
}

void MeasurementsDock::clear() {
    atoms_.clear();
    recompute();
}

void MeasurementsDock::recompute() {
    if (!protein_ || !conf_ || atoms_.size() < 2 || atoms_.size() > 4) {
        kindLabel_->setText(
            QStringLiteral("Select 2-4 atoms to measure a distance, angle, or dihedral."));
        valueLabel_->clear();
        atomsLabel_->clear();
        return;
    }

    const model::GeometryMeasurement m =
        model::Measure(*conf_, static_cast<std::size_t>(frame_), atoms_);

    kindLabel_->setText(QString::fromLatin1(model::NameForGeometryKind(m.kind)));

    if (!m.valid) {
        valueLabel_->setText(QStringLiteral("--"));
    } else if (m.kind == model::GeometryKind::Distance) {
        valueLabel_->setText(
            QStringLiteral("%1 %2").arg(m.value, 0, 'f', 3).arg(angstromGlyph()));
    } else {
        valueLabel_->setText(QStringLiteral("%1%2").arg(m.value, 0, 'f', 1).arg(degreeGlyph()));
    }

    QStringList names;
    for (std::size_t i = 0; i < atoms_.size(); ++i)
        names << QStringLiteral("%1. %2").arg(i + 1).arg(atomLabel(*protein_, atoms_[i]));
    atomsLabel_->setText(names.join(QLatin1Char('\n')));
}

}  // namespace h5reader::app
