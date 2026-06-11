// ResidentIndexes — day-one immutable indexes paired with RunData.

#pragma once

#include "RingGeometryCache.h"
#include "SpatialIndexSet.h"
#include "TemporalIndex.h"
#include "TypedAtomIndex.h"

#include <QHash>
#include <QString>
#include <QStringList>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace h5reader::rediscover {

class RunData;

enum class SecondaryStructure8 : int8_t {
    H = 0,
    G = 1,
    I = 2,
    E = 3,
    B = 4,
    T = 5,
    S = 6,
    C = 7,
    Unknown = 8,
};

enum class SecondaryStructure3 : int8_t {
    Helix = 0,
    Sheet = 1,
    Coil = 2,
    Unknown = 3,
};

struct SecondaryStructureState {
    SecondaryStructure8 ss8 = SecondaryStructure8::Unknown;
    SecondaryStructure3 ss3 = SecondaryStructure3::Unknown;
    bool observed = false;
};

inline SecondaryStructure3 Ss3ForSs8(SecondaryStructure8 ss8) {
    switch (ss8) {
    case SecondaryStructure8::H:
    case SecondaryStructure8::G:
    case SecondaryStructure8::I:
        return SecondaryStructure3::Helix;
    case SecondaryStructure8::E:
    case SecondaryStructure8::B:
        return SecondaryStructure3::Sheet;
    case SecondaryStructure8::T:
    case SecondaryStructure8::S:
    case SecondaryStructure8::C:
        return SecondaryStructure3::Coil;
    case SecondaryStructure8::Unknown:
        return SecondaryStructure3::Unknown;
    }
    return SecondaryStructure3::Unknown;
}

inline SecondaryStructureState DecodeSs8(const std::array<double, 8>& oneHot) {
    double sum = 0.0;
    double best = -std::numeric_limits<double>::infinity();
    int bestIdx = -1;
    for (int i = 0; i < 8; ++i) {
        const double v = oneHot[static_cast<std::size_t>(i)];
        sum += std::abs(v);
        if (v > best) {
            best = v;
            bestIdx = i;
        }
    }
    if (!(sum > 0.0) || bestIdx < 0 || !(best > 0.0))
        return {};
    const auto ss8 = static_cast<SecondaryStructure8>(bestIdx);
    return SecondaryStructureState{ss8, Ss3ForSs8(ss8), true};
}

class AtomNameIndex {
public:
    void reset(std::size_t atomCount) { perAtom_.assign(atomCount, QString()); byName_.clear(); }
    void add(std::size_t atom, const QString& name) {
        if (atom >= perAtom_.size()) return;
        const QString key = name.isEmpty() ? QStringLiteral("<unknown>") : name;
        perAtom_[atom] = key;
        byName_[key].push_back(static_cast<int32_t>(atom));
    }
    const QString& nameForAtom(std::size_t atom) const {
        static const QString empty;
        return atom < perAtom_.size() ? perAtom_[atom] : empty;
    }
    bool containsAtom(const QString& name, std::size_t atom) const {
        return atom < perAtom_.size() && perAtom_[atom] == name;
    }
    const std::vector<int32_t>& atomsForName(const QString& name) const {
        static const std::vector<int32_t> empty;
        auto it = byName_.find(name);
        return it == byName_.end() ? empty : it.value();
    }
    QStringList names() const { return byName_.keys(); }
    std::size_t atomCount() const { return perAtom_.size(); }

private:
    std::vector<QString> perAtom_;
    QHash<QString, std::vector<int32_t>> byName_;
};

class ChemicalCategoryIndex {
public:
    void reset(std::size_t atomCount) { perAtom_.assign(atomCount, QString()); byCategory_.clear(); }
    void add(std::size_t atom, const QString& category) {
        if (atom >= perAtom_.size()) return;
        const QString key = category.isEmpty() ? QStringLiteral("unknown") : category;
        perAtom_[atom] = key;
        byCategory_[key].push_back(static_cast<int32_t>(atom));
    }
    const QString& categoryForAtom(std::size_t atom) const {
        static const QString empty;
        return atom < perAtom_.size() ? perAtom_[atom] : empty;
    }
    bool containsAtom(const QString& category, std::size_t atom) const {
        return atom < perAtom_.size() && perAtom_[atom] == category;
    }
    const std::vector<int32_t>& atomsForCategory(const QString& category) const {
        static const std::vector<int32_t> empty;
        auto it = byCategory_.find(category);
        return it == byCategory_.end() ? empty : it.value();
    }
    QStringList categories() const { return byCategory_.keys(); }

private:
    std::vector<QString> perAtom_;
    QHash<QString, std::vector<int32_t>> byCategory_;
};

class SecondaryStructureIndex {
public:
    void reset(std::size_t atomCount, std::size_t frameCount) {
        atomCount_ = atomCount;
        frameCount_ = frameCount;
        states_.assign(atomCount_ * frameCount_, {});
    }
    void set(std::size_t atom, std::size_t frame, SecondaryStructureState state) {
        if (atom >= atomCount_ || frame >= frameCount_) return;
        states_[frame * atomCount_ + atom] = state;
    }
    SecondaryStructureState state(std::size_t atom, std::size_t frame) const {
        if (atom >= atomCount_ || frame >= frameCount_) return {};
        return states_[frame * atomCount_ + atom];
    }
    std::size_t atomCount() const { return atomCount_; }
    std::size_t frameCount() const { return frameCount_; }

private:
    std::size_t atomCount_ = 0;
    std::size_t frameCount_ = 0;
    std::vector<SecondaryStructureState> states_;
};

enum class DihedralKind : int8_t {
    Phi = 0,
    Psi = 1,
    Omega = 2,
    Chi1 = 3,
    Chi2 = 4,
    Chi3 = 5,
    Chi4 = 6,
    Count = 7,
};

enum class BackboneAngleSignPolicy : int8_t {
    Unresolved = 0,
    NpyIsIupac = 1,
    NpyIsNegatedIupac = 2,
};

struct DihedralState {
    double radians = std::numeric_limits<double>::quiet_NaN();
    int8_t fixed_bin = -1;
    bool present = false;
};

inline double WrapRadians(double radians) {
    constexpr double pi = 3.14159265358979323846264338327950288;
    while (radians <= -pi) radians += 2.0 * pi;
    while (radians > pi) radians -= 2.0 * pi;
    return radians;
}

inline double AngularDistance(double a, double b) {
    return std::abs(WrapRadians(a - b));
}

inline int8_t FixedDihedralBin(double radians) {
    if (!std::isfinite(radians)) return -1;
    constexpr double pi = 3.14159265358979323846264338327950288;
    constexpr double third = pi / 3.0;
    const double x = WrapRadians(radians);
    if (x < -third) return 0;
    if (x <= third) return 1;
    return 2;
}

inline BackboneAngleSignPolicy ChooseBackboneSignPolicy(double straightError,
                                                        double negatedError,
                                                        std::size_t samples) {
    if (samples == 0 || !std::isfinite(straightError) || !std::isfinite(negatedError))
        return BackboneAngleSignPolicy::Unresolved;
    constexpr double margin = 1.0e-6;
    if (negatedError + margin < straightError)
        return BackboneAngleSignPolicy::NpyIsNegatedIupac;
    if (straightError + margin < negatedError)
        return BackboneAngleSignPolicy::NpyIsIupac;
    return BackboneAngleSignPolicy::Unresolved;
}

inline double ApplyBackboneSignPolicy(double npyAngle, BackboneAngleSignPolicy policy) {
    if (policy == BackboneAngleSignPolicy::NpyIsNegatedIupac) return WrapRadians(-npyAngle);
    if (policy == BackboneAngleSignPolicy::NpyIsIupac) return WrapRadians(npyAngle);
    return std::numeric_limits<double>::quiet_NaN();
}

class DihedralIndex {
public:
    void reset(std::size_t atomCount, std::size_t frameCount) {
        atomCount_ = atomCount;
        frameCount_ = frameCount;
        for (auto& states : states_) states.assign(atomCount_ * frameCount_, {});
    }
    void set(DihedralKind kind, std::size_t atom, std::size_t frame, DihedralState state) {
        if (kind == DihedralKind::Count || atom >= atomCount_ || frame >= frameCount_) return;
        states_[static_cast<std::size_t>(kind)][frame * atomCount_ + atom] = state;
    }
    DihedralState state(DihedralKind kind, std::size_t atom, std::size_t frame) const {
        if (kind == DihedralKind::Count || atom >= atomCount_ || frame >= frameCount_) return {};
        return states_[static_cast<std::size_t>(kind)][frame * atomCount_ + atom];
    }
    void setBackboneSignPolicy(BackboneAngleSignPolicy policy) { backboneSignPolicy_ = policy; }
    BackboneAngleSignPolicy backboneSignPolicy() const { return backboneSignPolicy_; }
    std::size_t atomCount() const { return atomCount_; }
    std::size_t frameCount() const { return frameCount_; }

private:
    std::size_t atomCount_ = 0;
    std::size_t frameCount_ = 0;
    BackboneAngleSignPolicy backboneSignPolicy_ = BackboneAngleSignPolicy::Unresolved;
    std::array<std::vector<DihedralState>, static_cast<std::size_t>(DihedralKind::Count)> states_;
};

struct ResidentIndexes {
    TypedAtomIndex typedAtoms;
    AtomNameIndex iupacNames;
    ChemicalCategoryIndex chemicalCategories;
    SecondaryStructureIndex secondaryStructure;
    DihedralIndex dihedrals;
    SpatialIndexSet spatial;
    RingGeometryCache ringGeometry;
    TemporalIndex temporal;
};

ResidentIndexes BuildResidentIndexes(const RunData& run);

}  // namespace h5reader::rediscover
