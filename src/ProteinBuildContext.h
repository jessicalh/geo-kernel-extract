#pragma once
//
// ProteinBuildContext: how this protein instance was built.
// Immutable after construction. Records provenance.
//

#include <string>
#include <vector>
#include <memory>
#include <cmath>

namespace nmr {

class ProteinBuildContext {
public:
    std::string pdb_source;
    std::string deposition_date;
    std::string organism;
    double crystal_resolution = std::nan("");
    std::string protonation_tool;
    double protonation_pH = std::nan("");
    std::string force_field;
    std::string prmtop_path;          // AMBER prmtop that produced charges (provenance)
    std::string tleap_script_path;    // tleap input script used (provenance)
    std::vector<std::string> stripped;
    std::vector<std::string> assumptions;

    std::unique_ptr<ProteinBuildContext> Clone() const {
        auto c = std::make_unique<ProteinBuildContext>();
        c->pdb_source = pdb_source;
        c->deposition_date = deposition_date;
        c->organism = organism;
        c->crystal_resolution = crystal_resolution;
        c->protonation_tool = protonation_tool;
        c->protonation_pH = protonation_pH;
        c->force_field = force_field;
        c->prmtop_path = prmtop_path;
        c->tleap_script_path = tleap_script_path;
        c->stripped = stripped;
        c->assumptions = assumptions;
        return c;
    }

    std::string Describe() const {
        std::string desc = "Source: " + pdb_source;
        if (!protonation_tool.empty())
            desc += ", protonation: " + protonation_tool;
        if (!force_field.empty())
            desc += ", ff: " + force_field;
        return desc;
    }
};

}  // namespace nmr
