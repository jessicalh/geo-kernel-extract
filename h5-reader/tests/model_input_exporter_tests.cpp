#include "io/ModelInputExporter.h"
#include "io/QtFieldCatalog.gen.h"
#include "model/QtAtom.h"
#include "model/QtAtomNames.h"
#include "model/QtBond.h"
#include "model/QtResidue.h"
#include "model/QtResidueNames.h"
#include "model/QtRing.h"
#include "model/QtRingMembership.h"
#include "model/QtTopology.h"

#define private public
#include "model/QtProtein.h"
#undef private

#include "model/QtConformationSnapshot.h"
#include "model/SingleConformation.h"

#include <QTemporaryDir>
#include <QtTest>

#include <memory>
#include <vector>

namespace {

std::unique_ptr<h5reader::model::QtProtein> oneAtomProtein() {
    auto protein = std::make_unique<h5reader::model::QtProtein>();
    protein->atoms_.resize(1);
    protein->topology_ = std::make_unique<h5reader::model::QtTopology>(1,
                                                                       std::vector<h5reader::model::QtBond>{},
                                                                       std::vector<std::unique_ptr<h5reader::model::QtRing>>{},
                                                                       std::vector<h5reader::model::QtRingMembership>{},
                                                                       0,
                                                                       0);
    return protein;
}

}  // namespace

class ModelInputExporterTests final : public QObject {
    Q_OBJECT

private slots:
    void missingStaticPositionsAreRejected() {
        auto protein = oneAtomProtein();
        auto snapshot = std::make_shared<h5reader::model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        h5reader::model::SingleConformation conformation(protein.get(), snapshot);
        QTemporaryDir output;
        QVERIFY(output.isValid());

        h5reader::io::ModelInputExporter exporter(protein.get(), &conformation, output.path());
        QString error;
        QVERIFY(!exporter.prepare(&error));
        QCOMPARE(error, QStringLiteral("required field is absent: pos"));
    }

    void wrongStaticPositionShapeIsRejected() {
        auto protein = oneAtomProtein();
        auto snapshot = std::make_shared<h5reader::model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        h5reader::model::NpyColumn& positions = snapshot->mutableColumn(h5reader::io::FieldKind::Pos);
        positions.present = true;
        positions.rows = 1;
        positions.cols = 2;
        positions.data = {1.0, 2.0};
        h5reader::model::SingleConformation conformation(protein.get(), snapshot);
        QTemporaryDir output;
        QVERIFY(output.isValid());

        h5reader::io::ModelInputExporter exporter(protein.get(), &conformation, output.path());
        QString error;
        QVERIFY(!exporter.prepare(&error));
        QCOMPARE(error, QStringLiteral("pos has shape (1,2), expected (1,3)"));
    }
};

QTEST_GUILESS_MAIN(ModelInputExporterTests)
#include "model_input_exporter_tests.moc"
