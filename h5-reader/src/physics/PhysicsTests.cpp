#include "physics/BuckinghamEfield.h"
#include "physics/ClassicalSourceMath.h"
#include "physics/CircularRingCurrent.h"
#include "physics/EfgFeature.h"
#include "physics/LiteratureAccessors.h"
#include "physics/LocalFrameBasis.h"
#include "physics/McConnellLiteratureKernel.h"
#include "physics/RingCurrentScalars.h"
#include "physics/SphericalBasis.h"
#include "constants/LiteratureConstants.h"

#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

using h5reader::model::BondCategory;
using h5reader::model::Element;
using h5reader::model::Mat3;
using h5reader::model::RingTypeIndex;
using h5reader::model::Vec3;

int failures = 0;

void fail(const std::string& message) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
}

bool sameBits(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

bool endsWith(std::string_view text, std::string_view suffix) {
    return text.size() >= suffix.size()
           && text.substr(text.size() - suffix.size()) == suffix;
}

std::string readText(const std::string& path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

double producerDefault(const std::string& source, std::string_view key) {
    const std::regex pattern("add\\(\"" + std::string(key)
                             + "\"\\s*,\\s*([-+0-9.eE]+)\\s*,");
    std::smatch match;
    if (!std::regex_search(source, match, pattern)) {
        fail("producer default key not found: " + std::string(key));
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::stod(match[1].str());
}

void expectSameBits(double actual, double expected, const std::string& label) {
    if (!sameBits(actual, expected)) {
        std::ostringstream msg;
        msg << label << " bit mismatch actual=" << actual << " expected=" << expected;
        fail(msg.str());
    }
}

void testL0ProductionFreezeIdentity() {
    const std::string source = readText(std::string(H5READER_REPO_ROOT) + "/src/CalculatorConfig.cpp");
    if (source.empty()) {
        fail("could not read src/CalculatorConfig.cpp");
        return;
    }

    struct LobeKey {
        RingTypeIndex type;
        const char* key;
    };
    const std::array<LobeKey, 8> keys{{
        {RingTypeIndex::PheBenzene, "phe_benzene_jb_lobe_offset"},
        {RingTypeIndex::TyrPhenol, "tyr_phenol_jb_lobe_offset"},
        {RingTypeIndex::TrpBenzene, "trp_benzene_jb_lobe_offset"},
        {RingTypeIndex::TrpPyrrole, "trp_pyrrole_jb_lobe_offset"},
        {RingTypeIndex::TrpPerimeter, "trp_indole_perimeter_jb_lobe_offset"},
        {RingTypeIndex::HisImidazole, "his_imidazole_jb_lobe_offset"},
        {RingTypeIndex::HidImidazole, "hid_imidazole_jb_lobe_offset"},
        {RingTypeIndex::HieImidazole, "hie_imidazole_jb_lobe_offset"},
    }};

    for (const LobeKey& item : keys) {
        expectSameBits(h5reader::physics::JohnsonBoveyLobeOffset(item.type).value,
                       producerDefault(source, item.key),
                       item.key);
    }
    expectSameBits(h5reader::physics::JohnsonBoveyLobeOffset(RingTypeIndex::ProPyrrolidine).value,
                   0.0,
                   "pro_pyrrolidine_saturated_jb_lobe_offset");
}

void testL1ClosedFormGoldens() {
    const double amideShift =
        h5reader::physics::BuckinghamShiftPpm(Element::H, "ALA", "HN", "backbone_amide_h", 1.0);
    if (!near(amideShift, -3.66))
        fail("amide 1HN Buckingham unit-field shift should be -3.66 ppm");

    std::array<double, h5reader::model::kAromaticRingTypeCount> unitT0{};
    unitT0[static_cast<std::size_t>(RingTypeIndex::PheBenzene)] = 0.25;
    const double benzeneShift = h5reader::physics::RingPerTypeT0Ppm(unitT0.data(), unitT0.size());
    if (!near(benzeneShift, -3.0))
        fail("benzene ring unit-current T0 converter should produce -3.0 ppm");

    const double mcShift =
        h5reader::physics::McConnellProducerT0ToPpm(BondCategory::PeptideCO, 0.5);
    const double expectedMc = -nmr::constants::kMcConnellMolarPrefactor.value
                              * nmr::constants::kMcConnellPeptideCO.value
                              * 0.5;
    if (!near(mcShift, expectedMc))
        fail("McConnell packed T0 to ppm converter mismatch");
}

void testL3LiteratureMagnitudeBounds() {
    if (!(nmr::constants::kCoulombKeVA.value > 10.0
          && nmr::constants::kCoulombKeVA.value < 20.0))
        fail("Coulomb ke outside expected molecular range");

    for (const auto& c : nmr::constants::kJohnsonBoveyLobeOffsetByType) {
        if (!(c.value >= 0.0 && c.value <= 1.0))
            fail(std::string(c.key) + " lobe offset outside expected range");
        if (std::string_view(c.key).find("_correction") != std::string_view::npos)
            fail(std::string(c.key) + " production key was suffixed as correction");
    }

    for (const auto& c : nmr::constants::kRingIntensityByType) {
        if (!(c.value <= 0.0 && c.value >= -25.0))
            fail(std::string(c.key) + " ring intensity outside expected range");
        if (!endsWith(c.key, "_correction"))
            fail(std::string(c.key) + " correction key missing suffix");
    }

    for (const auto& c : nmr::constants::kCorrectionLiteratureConstants) {
        const std::string_view key(c.key);
        const bool requiresSuffix =
            key.rfind("buckingham.", 0) == 0 || key.rfind("sigma0.", 0) == 0
            || key.rfind("mcconnell.", 0) == 0 || key.rfind("larsen.", 0) == 0
            || key.rfind("ring.intensity.", 0) == 0;
        if (requiresSuffix && !endsWith(key, "_correction"))
            fail(std::string(c.key) + " correction key missing suffix");
    }

    if (!(h5reader::physics::BuckinghamA(Element::H, "ALA", "HN", "backbone_amide_h").value > 0.0))
        fail("amide 1HN Buckingham A should be positive");
    if (!(h5reader::physics::BuckinghamA(Element::O, "SER", "OG", "hydroxyl_oxygen").value < 0.0))
        fail("hydroxyl oxygen Buckingham A should be negative");
    if (!(h5reader::physics::McConnellDeltaChi(BondCategory::PeptideCN).value < 0.0))
        fail("peptide C-N McConnell delta chi should be negative");
}

void testFrameAndTensorSmoke() {
    const auto frame = h5reader::physics::BuildHNFrame(Vec3(0.0, 0.0, 0.0),
                                                       Vec3(0.0, 0.0, 1.0),
                                                       Vec3(1.0, 0.0, 0.0),
                                                       Vec3(1.0, 0.0, 0.0),
                                                       false);
    if (!frame.is_valid) {
        fail("HN frame should build for orthogonal synthetic geometry");
        return;
    }

    const auto projection = h5reader::physics::ProjectBuckinghamEfield(frame, Vec3(0.0, 0.0, 2.0));
    if (!projection.present || !near(projection.e_parallel, 2.0))
        fail("Buckingham E-field projection should preserve z-axis field");

    Mat3 m = Mat3::Zero();
    m(0, 0) = 1.0;
    m(1, 1) = 2.0;
    m(2, 2) = 3.0;
    const auto st = h5reader::physics::DecomposeLibrary(m);
    const Mat3 reconstructed = h5reader::physics::ReconstructLibraryT2Matrix(st.T0, st.T2);
    if (!near((m - reconstructed).norm(), 0.0))
        fail("library T2 reconstruction should invert symmetric decomposition");
}

void testCircularRingCurrent() {
    struct ParameterCase {
        RingTypeIndex type;
        int protonationVariant;
        double radiusA;
        double currentNanoamperePerTesla;
    };
    const std::array<ParameterCase, 7> parameterCases{{
        {RingTypeIndex::PheBenzene, -1, 1.39, -15.203781575116498},
        {RingTypeIndex::TyrPhenol, -1, 1.39, -12.24375330385488},
        {RingTypeIndex::TrpBenzene, -1, 1.39, -15.876515273130503},
        {RingTypeIndex::TrpPyrrole, -1, 1.182, -14.261954397896893},
        {RingTypeIndex::HisImidazole, 2, 1.182, -17.087435929555713},
        {RingTypeIndex::HidImidazole, -1, 1.182, -16.81834245035011},
        {RingTypeIndex::HieImidazole, -1, 1.182, -16.81834245035011},
    }};
    for (const ParameterCase& expected : parameterCases) {
        const auto actual = h5reader::physics::CandidateACircularParameters(expected.type, expected.protonationVariant);
        if (!actual) {
            fail("Candidate-A parameter row unexpectedly absent");
            continue;
        }
        if (!near(actual->radiusA, expected.radiusA) || !near(actual->lobeOffsetA, 0.64)
            || !near(actual->currentNanoamperePerTesla, expected.currentNanoamperePerTesla)) {
            fail("Candidate-A parameter row mismatch");
        }
    }

    const auto trp6 = h5reader::physics::CandidateACircularParameters(RingTypeIndex::TrpBenzene);
    const auto trp5 = h5reader::physics::CandidateACircularParameters(RingTypeIndex::TrpPyrrole);
    if (!trp6 || !trp5) {
        fail("Candidate-A TRP parameters should be available");
        return;
    }
    if (h5reader::physics::CandidateACircularParameters(RingTypeIndex::TrpPerimeter)) {
        fail("Candidate-A must omit the diagnostic TRP perimeter");
    }
    if (h5reader::physics::CandidateACircularParameters(RingTypeIndex::ProPyrrolidine))
        fail("Candidate-A must omit saturated proline rings");
    if (h5reader::physics::CandidateACircularParameters(RingTypeIndex::HisImidazole, -1))
        fail("Candidate-A must require an explicit protonated-His variant");

    const Vec3 center = Vec3::Zero();
    const Vec3 normal(0.0, 0.0, 1.0);
    const h5reader::physics::CircularRingPlane plane{center, normal, 0.0};
    const auto axis = h5reader::physics::EvaluateCircularShielding(Vec3(0.0, 0.0, 3.0), plane, *trp6);
    const auto offAxis = h5reader::physics::EvaluateCircularShielding(Vec3(1.2, 0.0, 2.7), plane, *trp6);
    if (!axis || !offAxis) {
        fail("Candidate-A canonical points should be finite");
        return;
    }
    if (!near(axis->T0, 2.1064547488716219, 1e-10))
        fail("Candidate-A axial golden mismatch");
    if (!near(offAxis->T0, 1.8395779259231924, 1e-10))
        fail("Candidate-A off-axis golden mismatch");

    const h5reader::physics::CircularRingPlane rotatedPlane{center, Vec3(1.0, 0.0, 0.0), 0.0};
    const auto rotated = h5reader::physics::EvaluateCircularShielding(Vec3(2.7, 1.2, 0.0), rotatedPlane, *trp6);
    if (!rotated || !near(rotated->T0, offAxis->T0, 1e-10))
        fail("Candidate-A field should rotate with the ring frame");
}

}  // namespace

int main() {
    testL0ProductionFreezeIdentity();
    testL1ClosedFormGoldens();
    testL3LiteratureMagnitudeBounds();
    testFrameAndTensorSmoke();
    testCircularRingCurrent();

    if (failures != 0) {
        std::cerr << failures << " physics test failure(s)\n";
        return 1;
    }
    std::cout << "physics tests passed\n";
    return 0;
}
