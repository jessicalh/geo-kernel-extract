#include "OrcaShieldingParser.h"

#include "../model/SphericalDecomposition.h"

#include <cctype>
#include <sstream>
#include <string>

namespace h5reader::io {

namespace {

using model::Mat3;
using model::Vec3;

model::Element elementFromSymbol(const std::string& s) {
    if (s == "H") return model::Element::H;
    if (s == "C") return model::Element::C;
    if (s == "N") return model::Element::N;
    if (s == "O") return model::Element::O;
    if (s == "S") return model::Element::S;
    return model::Element::Unknown;
}

// Scan forward for the next line containing `label`, then read the following
// three lines as a 3x3 matrix (three whitespace-separated doubles each).
// Returns false if EOF or a new "Nucleus" block is hit before the label.
bool seekMatrix(std::istream& in, const char* label, Mat3& m) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.find(label) != std::string::npos) {
            for (int r = 0; r < 3; ++r) {
                if (!std::getline(in, line)) return false;
                std::istringstream ls(line);
                if (!(ls >> m(r, 0) >> m(r, 1) >> m(r, 2))) return false;
            }
            return true;
        }
        if (line.find("Nucleus") != std::string::npos) return false;
    }
    return false;
}

}  // namespace

model::DftShieldingFrame ParseOrcaNmrShielding(std::istream& in) {
    model::DftShieldingFrame frame;
    std::string             line;
    bool                    inSection = false;
    // The ORCA-input geometry block (precedes CHEMICAL SHIELDINGS): the
    // orientation the shielding tensors are computed in. Captured additively so
    // the extractor can frame-check the DFT tensors against the H5 positions.
    std::vector<Vec3> coords;

    while (std::getline(in, line)) {
        if (!inSection) {
            // First "CARTESIAN COORDINATES (ANGSTROEM)" block (skip the dashed
            // rule line, then read "<symbol> x y z" rows until one fails to parse;
            // ignore the later (A.U.) block — units differ).
            if (coords.empty()
                && line.find("CARTESIAN COORDINATES (ANGSTROEM)") != std::string::npos) {
                std::getline(in, line);  // dashed rule
                while (std::getline(in, line)) {
                    std::istringstream ls(line);
                    std::string sym;
                    double x = 0, y = 0, z = 0;
                    if (ls >> sym >> x >> y >> z) coords.emplace_back(x, y, z);
                    else break;
                }
                continue;
            }
            if (line.find("CHEMICAL SHIELDINGS") != std::string::npos)
                inSection = true;
            continue;
        }

        const auto np = line.find("Nucleus");
        if (np == std::string::npos)
            continue;

        // Token after "Nucleus" is index+element, e.g. "0N", "4C", "10H".
        std::istringstream ls(line.substr(np + 7));
        std::string        tok;
        if (!(ls >> tok))
            continue;
        std::size_t k = 0;
        while (k < tok.size() && std::isdigit(static_cast<unsigned char>(tok[k])))
            ++k;
        if (k == 0 || k >= tok.size())
            continue;  // malformed (no digits, or no trailing element)
        const int         idx = std::stoi(tok.substr(0, k));
        const std::string el  = tok.substr(k);

        Mat3 dia, para, total;
        if (!seekMatrix(in, "Diamagnetic", dia))             break;
        if (!seekMatrix(in, "Paramagnetic", para))           break;
        if (!seekMatrix(in, "Total shielding tensor", total)) break;

        model::DftAtomShielding a;
        a.element = elementFromSymbol(el);
        a.total   = model::DecomposeShielding(total);
        a.dia     = model::DecomposeShielding(dia);
        a.para    = model::DecomposeShielding(para);

        // Additive (DESIGN.md): keep the raw 3x3 the parser already built,
        // which DecomposeShielding above would otherwise discard. The
        // rediscover substrate re-decomposes these in the library T2 order.
        a.total_raw = total;
        a.dia_raw   = dia;
        a.para_raw  = para;

        if (idx >= 0) {
            if (static_cast<int>(frame.atoms.size()) <= idx)
                frame.atoms.resize(idx + 1);
            frame.atoms[idx] = a;
        }
    }

    // Attach the ORCA-input position to each atom (nucleus index == coordinate
    // order == topology order). Atoms beyond the parsed coords keep zero.
    for (std::size_t i = 0; i < frame.atoms.size() && i < coords.size(); ++i)
        frame.atoms[i].orca_coord = coords[i];

    frame.valid = !frame.atoms.empty();
    return frame;
}

}  // namespace h5reader::io
