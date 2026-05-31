#include "ChargeStore.h"

#include "../model/QtAtomNames.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"

#include <QFile>
#include <QRegularExpression>
#include <QTextStream>

namespace h5reader::rediscover {

namespace {

QString stripped(QString line) {
    const int semi = line.indexOf(QLatin1Char(';'));
    if (semi >= 0) line = line.left(semi);
    return line.trimmed();
}

bool splitTrailingNumber(const QString& s, QString& base, int& number) {
    int pos = s.size();
    while (pos > 0 && s[pos - 1].isDigit()) --pos;
    if (pos == s.size()) return false;
    bool ok = false;
    number = s.mid(pos).toInt(&ok);
    if (!ok) return false;
    base = s.left(pos);
    return true;
}

bool atomNameMatches(const model::QtAtom& atom, const model::QtAtomNames& names,
                     const QString& topolName) {
    if (names.amber == topolName || names.iupac == topolName || names.bmrb == topolName)
        return true;

    // GROMACS/AMBER topologies often number a two-proton methylene as HB1/HB2
    // while the NMR/IUPAC projection uses HB2/HB3. This is an alias at the
    // projection boundary, not a positional pick: the atom order is still 1:1
    // and the residue + typed atom identity remain fixed.
    QString topBase, modelBase;
    int topNum = 0, modelNum = 0;
    if (!splitTrailingNumber(topolName, topBase, topNum)) {
        for (const QString& modelName : {names.amber, names.iupac, names.bmrb})
            if (modelName == topolName + QStringLiteral("1")) return true;
        return false;
    }
    if (atom.element != model::Element::H) return false;
    for (const QString& modelName : {names.amber, names.iupac, names.bmrb}) {
        if (!splitTrailingNumber(modelName, modelBase, modelNum)) continue;
        if (topBase == modelBase && modelNum == topNum + 1) return true;
        if (topBase == modelBase && modelNum == 10 + topNum) return true;
        // Glycine alpha protons are the classic Markley/IUPAC inversion:
        // topology HA1/HA2 may project to model HA3/HA2.
        if (atom.locant == model::Locant::Alpha && topBase == QStringLiteral("HA")
            && modelBase == QStringLiteral("HA") && topNum == 1 && modelNum == 3)
            return true;
    }
    return false;
}

}  // namespace

bool LoadFf14sbChargesFromTopol(const QString& topolPath, model::QtProtein& protein,
                                QString* err_out) {
    QFile file(topolPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (err_out) *err_out = QStringLiteral("cannot open topol.top for FF14SB charges: %1").arg(topolPath);
        return false;
    }

    QTextStream in(&file);
    bool inAtoms = false;
    std::size_t atomRow = 0;
    int lineNo = 0;
    while (!in.atEnd()) {
        ++lineNo;
        const QString line = stripped(in.readLine());
        if (line.isEmpty()) continue;
        if (line.startsWith(QLatin1Char('['))) {
            const QString section = line.mid(1, line.indexOf(QLatin1Char(']')) - 1).trimmed();
            if (section == QStringLiteral("atoms")) {
                inAtoms = true;
                continue;
            }
            if (inAtoms) break;
            continue;
        }
        if (!inAtoms) continue;

        const QStringList toks = line.split(QRegularExpression(QStringLiteral("\\s+")), Qt::SkipEmptyParts);
        if (toks.size() < 8) {
            if (err_out)
                *err_out = QStringLiteral("malformed [ atoms ] row at %1:%2").arg(topolPath).arg(lineNo);
            return false;
        }
        if (atomRow >= protein.atomCount()) {
            if (err_out)
                *err_out = QStringLiteral("topol.top has more [ atoms ] rows than protein atoms (%1)")
                               .arg(protein.atomCount());
            return false;
        }

        bool resOk = false, chargeOk = false;
        const int resnr = toks[2].toInt(&resOk);
        const QString atomName = toks[4];
        const double charge = toks[6].toDouble(&chargeOk);
        if (!resOk || !chargeOk) {
            if (err_out)
                *err_out = QStringLiteral("cannot parse resnr/charge at %1:%2").arg(topolPath).arg(lineNo);
            return false;
        }

        const model::QtAtom& atom = protein.atom(atomRow);
        if (atom.residueIndex < 0 || static_cast<std::size_t>(atom.residueIndex) >= protein.residueCount()) {
            if (err_out) *err_out = QStringLiteral("protein atom %1 has invalid residue index").arg(atomRow);
            return false;
        }
        const model::QtResidue& residue = protein.residue(static_cast<std::size_t>(atom.residueIndex));
        const model::QtAtomNames& names = protein.atomNames(atomRow);
        const bool nameMatches = atomNameMatches(atom, names, atomName);
        if (residue.address.residueNumber != resnr || !nameMatches) {
            if (err_out) {
                *err_out = QStringLiteral("topol.top [ atoms ] mismatch at protein atom %1: "
                                          "topol=(resnr %2 atom %3), model=(resnr %4 amber %5 iupac %6 bmrb %7)")
                               .arg(atomRow)
                               .arg(resnr)
                               .arg(atomName)
                               .arg(residue.address.residueNumber)
                               .arg(names.amber)
                               .arg(names.iupac)
                               .arg(names.bmrb);
            }
            return false;
        }

        protein.setAtomPartialCharge(atomRow, charge);
        ++atomRow;
    }

    if (!inAtoms) {
        if (err_out) *err_out = QStringLiteral("topol.top has no inline [ atoms ] section: %1").arg(topolPath);
        return false;
    }
    if (atomRow != protein.atomCount()) {
        if (err_out)
            *err_out = QStringLiteral("topol.top [ atoms ] row count %1 != protein atom count %2")
                           .arg(atomRow)
                           .arg(protein.atomCount());
        return false;
    }
    return true;
}

}  // namespace h5reader::rediscover
