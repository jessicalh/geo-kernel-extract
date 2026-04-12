#include "GromacsFinalResult.h"
#include "OperationLog.h"
#include "Atom.h"
#include "Residue.h"

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

namespace nmr {

GromacsFinalResult::GromacsFinalResult(GromacsProtein& gp)
    : gp_(gp)
{}


bool GromacsFinalResult::Finalize(
        const std::string& output_dir,
        const std::string& temp_dir) {

    OperationLog::Scope scope("GromacsFinalResult::Finalize", gp_.protein_id());

    const auto& paths = gp_.frame_paths();

    OperationLog::Info(LogCalcOther, "GromacsFinalResult::Finalize",
        gp_.protein_id() + ": " +
        std::to_string(paths.size()) + " frames processed, " +
        std::to_string(gp_.protein().AtomCount()) + " atoms, " +
        "conformation 0 in Protein vector, " +
        std::to_string(gp_.protein().ConformationCount()) +
        " conformations in Protein");

    // Write per-atom trajectory catalog
    if (gp_.AtomCount() > 0) {
        std::string catalog_path = output_dir + "/atom_catalog.csv";
        std::ofstream csv(catalog_path);
        csv << "atom_idx,element,resname,resnum,atom_name,"
            << "rmsf,"
            << "water_n_first_mean,water_n_first_std,"
            << "water_n_first_min,water_n_first_min_frame,"
            << "water_n_first_max,water_n_first_max_frame,"
            << "water_emag_mean,water_emag_std,"
            << "water_emag_min,water_emag_min_frame,"
            << "water_emag_max,water_emag_max_frame,"
            << "sasa_mean,sasa_std,"
            << "sasa_min,sasa_min_frame,"
            << "sasa_max,sasa_max_frame,"
            << "half_shell_mean,half_shell_std,"
            << "dipole_cos_mean,dipole_cos_std,"
            << "nearest_ion_mean,"
            << "n_frames_dry,n_frames_exposed,n_frames_total\n";

        const auto& protein = gp_.protein();
        for (size_t ai = 0; ai < gp_.AtomCount(); ++ai) {
            const auto& ga = gp_.AtomAt(ai);
            const auto& atom = protein.AtomAt(ai);
            const auto& residue = protein.ResidueAt(atom.residue_index);

            const char* elem = "?";
            switch (atom.element) {
                case Element::H: elem = "H"; break;
                case Element::C: elem = "C"; break;
                case Element::N: elem = "N"; break;
                case Element::O: elem = "O"; break;
                case Element::S: elem = "S"; break;
                default: break;
            }

            csv << ai << ","
                << elem << ","
                << GetAminoAcidType(residue.type).three_letter_code << ","
                << residue.sequence_number << ","
                << atom.pdb_atom_name << ","
                << ga.RMSF() << ","
                << ga.water_n_first.mean << "," << ga.water_n_first.Std() << ","
                << ga.water_n_first.min_val << "," << ga.water_n_first.min_frame << ","
                << ga.water_n_first.max_val << "," << ga.water_n_first.max_frame << ","
                << ga.water_emag.mean << "," << ga.water_emag.Std() << ","
                << ga.water_emag.min_val << "," << ga.water_emag.min_frame << ","
                << ga.water_emag.max_val << "," << ga.water_emag.max_frame << ","
                << ga.sasa.mean << "," << ga.sasa.Std() << ","
                << ga.sasa.min_val << "," << ga.sasa.min_frame << ","
                << ga.sasa.max_val << "," << ga.sasa.max_frame << ","
                << ga.half_shell.mean << "," << ga.half_shell.Std() << ","
                << ga.dipole_cos.mean << "," << ga.dipole_cos.Std() << ","
                << ga.nearest_ion_dist.mean << ","
                << ga.n_frames_dry << "," << ga.n_frames_exposed << ","
                << ga.water_n_first.count << "\n";
        }

        OperationLog::Info(LogCalcOther, "GromacsFinalResult",
            "wrote atom_catalog.csv: " + std::to_string(gp_.AtomCount()) +
            " atoms x " + std::to_string(gp_.AtomAt(0).water_n_first.count) +
            " frames");
    }

    // TODO: write {protein_id}.h5 from conformation 0 + topology
    // TODO: select winners, move winner dirs from temp to output
    // TODO: clean up temp

    return true;
}

}  // namespace nmr
