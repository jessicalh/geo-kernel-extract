#include "SignalDictionary.h"

#include <algorithm>

namespace h5reader::model {

std::vector<SignalSpec> CoreSignalDictionary() {
    return {
        {SignalKey{SignalFamily::Geometry, QStringLiteral("distance"), QString()},
         QStringLiteral("Distance"), QString::fromUtf8("Å"),
         SignalAnchorKind::AtomTuple, SignalValueKind::Distance, true},
        {SignalKey{SignalFamily::Geometry, QStringLiteral("angle"), QString()},
         QStringLiteral("Angle"), QString::fromUtf8("°"),
         SignalAnchorKind::AtomTuple, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::Geometry, QStringLiteral("dihedral"), QString()},
         QStringLiteral("Dihedral"), QString::fromUtf8("°"),
         SignalAnchorKind::AtomTuple, SignalValueKind::Angle, true},

        {SignalKey{SignalFamily::DftShielding, QStringLiteral("sigma"), QStringLiteral("total.T0")},
         QStringLiteral("DFT sigma total"), QStringLiteral("ppm"),
         SignalAnchorKind::Atom, SignalValueKind::Scalar, true},

        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("phi"), QString()},
         QStringLiteral("phi"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("psi"), QString()},
         QStringLiteral("psi"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("omega"), QString()},
         QStringLiteral("omega"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("chi1"), QString()},
         QStringLiteral("chi1"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("chi2"), QString()},
         QStringLiteral("chi2"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("chi3"), QString()},
         QStringLiteral("chi3"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
        {SignalKey{SignalFamily::ResidueDihedral, QStringLiteral("chi4"), QString()},
         QStringLiteral("chi4"), QString::fromUtf8("°"),
         SignalAnchorKind::Residue, SignalValueKind::Angle, true},
    };
}

std::optional<SignalSpec> FindSignalSpec(const SignalKey& key) {
    const auto dictionary = CoreSignalDictionary();
    const auto it = std::find_if(dictionary.begin(), dictionary.end(), [&key](const SignalSpec& spec) {
        return spec.key == key;
    });
    if (it == dictionary.end())
        return std::nullopt;
    return *it;
}

}  // namespace h5reader::model
