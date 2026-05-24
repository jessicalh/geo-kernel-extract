/*
 * C bridge for APBS — isolates FETK's Vec3/Mat3 from Eigen's.
 *
 * This header is included by both C and C++ code.
 * The implementation (apbs_bridge.c) includes APBS/FETK headers.
 * The C++ wrapper (ApbsSolver.cpp) includes Eigen headers.
 * They never meet in the same translation unit.
 *
 * In-memory path: coordinate/charge/radii arrays go directly into APBS's
 * Valist/Vpbe/Vpmg objects.  No PQR files, no input files, no temp dirs.
 *
 * Units:
 *   positions: Angstroms (double[3] per atom)
 *   charges: elementary charges (double per atom)
 *   radii: Angstroms (double per atom)
 *   potential grid: kT/e (double per grid point)
 *   E-field (derived by caller): kT/(e*Å)
 *   EFG (derived by caller): kT/(e*Å²)
 *
 * Self-potential note:
 *   The returned grid includes each atom's own Coulomb potential.
 *   At atom positions, ∇²φ_self = -(q/ε)δ(r-r_i), which the grid
 *   discretizes into a large finite Laplacian.  When computing the EFG
 *   tensor (∂²φ/∂x_i∂x_j) at atom positions, callers MUST apply a
 *   traceless projection (subtract trace/3 from diagonal) to remove
 *   this self-interaction artifact.  The external-source EFG is
 *   guaranteed traceless by Laplace's equation.
 *   See ApbsSolver::ElectricFieldGradientAtPoint for the C++ implementation.
 */

#ifndef APBS_BRIDGE_H
#define APBS_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes */
#define APBS_BRIDGE_OK      0
#define APBS_BRIDGE_ERROR  -1

/* Grid data returned from solve */
typedef struct {
    double origin[3];       /* grid origin in Angstroms */
    double spacing[3];      /* grid spacing in Angstroms */
    int    dims[3];         /* nx, ny, nz */
    double* data;           /* potential values in kT/e (caller must free) */
    int    n_points;        /* total grid points (nx*ny*nz) */
    char   error_msg[512];  /* error message if return != OK */
} ApbsGridResult;

/*
 * Solve the linearized Poisson-Boltzmann equation.
 *
 * All arrays are indexed [0..n_atoms-1].
 * Returns APBS_BRIDGE_OK on success, APBS_BRIDGE_ERROR on failure.
 * On success, result->data is allocated (caller must free with free()).
 * On failure, result->error_msg contains the reason.
 */
int apbs_solve(
    int n_atoms,
    const double* x,            /* x coordinates, Angstroms */
    const double* y,            /* y coordinates, Angstroms */
    const double* z,            /* z coordinates, Angstroms */
    const double* charges,      /* partial charges, elementary charges */
    const double* radii,        /* atomic radii, Angstroms */
    double pdie,                /* protein interior dielectric */
    double sdie,                /* solvent dielectric */
    double temperature,         /* Kelvin */
    double ionic_strength,      /* molar */
    int grid_nx, int grid_ny, int grid_nz,  /* grid dimensions */
    double fine_x, double fine_y, double fine_z,  /* fine grid lengths, Ang */
    double coarse_x, double coarse_y, double coarse_z,  /* coarse grid, Ang */
    ApbsGridResult* result
);

/* Free the grid data allocated by apbs_solve */
void apbs_free_grid(ApbsGridResult* result);

#ifdef __cplusplus
}
#endif

#endif /* APBS_BRIDGE_H */
