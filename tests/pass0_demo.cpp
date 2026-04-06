// Pass 0 demonstration: the complete calling pattern.
//
// This code IS the pattern every subsequent agent copies.
// Load a protein, create conformations, attach results, query them.
//
// Usage: nmr_pass0_demo path/to/protein.pdb

#include "PdbFileReader.h"
#include "RuntimeEnvironment.h"
#include "GeometryResult.h"
#include "DemoResult.h"
#include <cstdio>

using namespace nmr;

int main(int argc, char* argv[]) {
    RuntimeEnvironment::Load();

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <pdb_file>\n", argv[0]);
        return 1;
    }
    const char* pdb_path = argv[1];

    // 1. Load protein
    auto r = BuildFromPdb(pdb_path);
    if (!r.Ok()) {
        fprintf(stderr, "Failed to load protein: %s\n", r.error.c_str());
        return 1;
    }
    auto& protein = *r.protein;

    printf("Loaded protein: %zu atoms, %zu residues, %zu rings, %zu bonds\n",
           protein.AtomCount(), protein.ResidueCount(),
           protein.RingCount(), protein.BondCount());

    // 2. Get the crystal conformation
    auto& conf = protein.Conformation();
    printf("Crystal conformation: %zu atoms\n", conf.AtomCount());

    // 3. Attach GeometryResult (ring/bond geometry, pre-built collections)
    conf.AttachResult(GeometryResult::Compute(conf));
    printf("GeometryResult attached\n");

    // 4. Attach DemoResult (nearest ring distance per atom)
    conf.AttachResult(DemoResult::Compute(conf));
    printf("DemoResult attached\n");

    // 5. Query: nearest ring distance for atom 42
    auto& demo = conf.Result<DemoResult>();
    printf("Atom 42: nearest ring distance = %.2f A\n",
           demo.NearestRingDistance(42));

    // 6. SphericalTensor roundtrip verification
    printf("SphericalTensor test T0 = %.6f\n", demo.TestDecomposition().T0);

    printf("\nPass 0 complete.\n");
    return 0;
}
