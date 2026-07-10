/*
 * APBS C bridge — includes APBS/FETK headers outside C++ translation units.
 * Compiled as C.
 *
 * Single solve on the caller-provided manual grid with SDH boundary
 * conditions. There is no focusing layer in this bridge.
 *
 * In-memory path: arrays → Valist → Vpbe → Vpmg → solve → extract.
 */

#include "apbs_bridge.h"

#include <apbs/apbs.h>
#include <apbs/generic/vatom.h>
#include <apbs/generic/valist.h>
#include <apbs/generic/vpbe.h>
#include <apbs/generic/mgparm.h>
#include <apbs/generic/pbeparm.h>
#include <apbs/mg/vpmg.h>
#include <apbs/mg/vpmgp.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static int compute_nlev(int dim) {
    int n = dim - 1;
    int L = 0;
    while (n > 1 && (n % 2) == 0) { n /= 2; L++; }
    return L;
}

void apbs_free_grid(ApbsGridResult* result) {
    if (result && result->data) { free(result->data); result->data = NULL; }
}

static Valist* build_atom_list(int n_atoms,
                                const double* x, const double* y, const double* z,
                                const double* charges, const double* radii) {
    Valist* alist = Valist_ctor();
    if (!alist) return NULL;
    alist->atoms = (Vatom*)Vmem_malloc(alist->vmem, n_atoms, sizeof(Vatom));
    if (!alist->atoms) { Valist_dtor(&alist); return NULL; }
    alist->number = n_atoms;
    for (int i = 0; i < n_atoms; i++) {
        Vatom_ctor2(&alist->atoms[i]);
        double pos[3] = { x[i], y[i], z[i] };
        Vatom_setPosition(&alist->atoms[i], pos);
        Vatom_setCharge(&alist->atoms[i], charges[i]);
        Vatom_setRadius(&alist->atoms[i], radii[i]);
        Vatom_setAtomID(&alist->atoms[i], i);
        Vatom_setPartID(&alist->atoms[i], 0);
        Vatom_setEpsilon(&alist->atoms[i], 0.0);
    }
    Valist_getStatistics(alist);
    return alist;
}

int apbs_solve(
    int n_atoms,
    const double* x, const double* y, const double* z,
    const double* charges, const double* radii,
    const ApbsSolveParams* params,
    ApbsGridResult* result)
{
    if (!result) return APBS_BRIDGE_ERROR;
    memset(result, 0, sizeof(ApbsGridResult));
    if (n_atoms <= 0 || !x || !y || !z || !charges || !radii) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "invalid atom arrays or n_atoms=%d", n_atoms);
        return APBS_BRIDGE_ERROR;
    }
    if (!params) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "ApbsSolveParams is null");
        return APBS_BRIDGE_ERROR;
    }
    if (params->grid_nx <= 1 || params->grid_ny <= 1 ||
        params->grid_nz <= 1 ||
        !isfinite(params->grid_len_x) || params->grid_len_x <= 0.0 ||
        !isfinite(params->grid_len_y) || params->grid_len_y <= 0.0 ||
        !isfinite(params->grid_len_z) || params->grid_len_z <= 0.0) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "invalid manual grid dimensions or lengths");
        return APBS_BRIDGE_ERROR;
    }
    if (!isfinite(params->pdie) || params->pdie <= 0.0 ||
        !isfinite(params->sdie) || params->sdie <= 0.0 ||
        !isfinite(params->temperature) || params->temperature <= 0.0) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "invalid dielectric or temperature");
        return APBS_BRIDGE_ERROR;
    }
    if (params->mobile_ion_count < 0 ||
        (params->mobile_ion_count > 0 &&
         (!params->mobile_ion_conc_M || !params->mobile_ion_radius_A ||
          !params->mobile_ion_charge_e))) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "invalid mobile-ion count or arrays");
        return APBS_BRIDGE_ERROR;
    }
    for (int ion = 0; ion < params->mobile_ion_count; ++ion) {
        if (!isfinite(params->mobile_ion_conc_M[ion]) ||
            params->mobile_ion_conc_M[ion] < 0.0 ||
            !isfinite(params->mobile_ion_radius_A[ion]) ||
            params->mobile_ion_radius_A[ion] <= 0.0 ||
            !isfinite(params->mobile_ion_charge_e[ion])) {
            snprintf(result->error_msg, sizeof(result->error_msg),
                     "invalid mobile-ion entry %d", ion);
            return APBS_BRIDGE_ERROR;
        }
    }
    for (int atom = 0; atom < n_atoms; ++atom) {
        if (!isfinite(x[atom]) || !isfinite(y[atom]) ||
            !isfinite(z[atom]) || !isfinite(charges[atom]) ||
            !isfinite(radii[atom]) || radii[atom] <= 0.0) {
            snprintf(result->error_msg, sizeof(result->error_msg),
                     "invalid atom entry %d", atom);
            return APBS_BRIDGE_ERROR;
        }
    }

    Valist* alist = build_atom_list(n_atoms, x, y, z, charges, radii);
    if (!alist) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Failed to build atom list (%d atoms)", n_atoms);
        return APBS_BRIDGE_ERROR;
    }

    Vpbe* pbe = Vpbe_ctor(alist, params->mobile_ion_count,
        (double*)params->mobile_ion_conc_M,
        (double*)params->mobile_ion_radius_A,
        (double*)params->mobile_ion_charge_e,
        params->temperature, params->pdie, params->sdie,
        1.4, 0, 10.0, 0.0, 0.0, 1.0, 0.0);
    if (!pbe) {
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpbe_ctor failed");
        Valist_dtor(&alist); return APBS_BRIDGE_ERROR;
    }

    double cx = Valist_getCenterX(alist);
    double cy = Valist_getCenterY(alist);
    double cz = Valist_getCenterZ(alist);

    /* Use the caller-provided manual grid for the single SDH-boundary solve. */
    MGparm* mgparm = MGparm_ctor(MCT_MANUAL);
    if (!mgparm) {
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "MGparm_ctor failed");
        return APBS_BRIDGE_ERROR;
    }

    const int grid_nx = params->grid_nx;
    const int grid_ny = params->grid_ny;
    const int grid_nz = params->grid_nz;
    const double grid_len_x = params->grid_len_x;
    const double grid_len_y = params->grid_len_y;
    const double grid_len_z = params->grid_len_z;

    mgparm->dime[0] = grid_nx; mgparm->dime[1] = grid_ny; mgparm->dime[2] = grid_nz;
    mgparm->setdime = 1;
    mgparm->glen[0] = grid_len_x; mgparm->glen[1] = grid_len_y; mgparm->glen[2] = grid_len_z;
    mgparm->setglen = 1;
    mgparm->cmeth = MCM_POINT;
    mgparm->center[0] = cx; mgparm->center[1] = cy; mgparm->center[2] = cz;
    mgparm->setgcent = 1;
    mgparm->chgm = VCM_BSPL2;
    mgparm->setchgm = 1;

    Vpmgp* pmgp = Vpmgp_ctor(mgparm);
    if (!pmgp) {
        MGparm_dtor(&mgparm); Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpmgp_ctor failed");
        return APBS_BRIDGE_ERROR;
    }

    pmgp->nx = grid_nx; pmgp->ny = grid_ny; pmgp->nz = grid_nz;
    pmgp->hx = grid_len_x / (grid_nx - 1);
    pmgp->hy = grid_len_y / (grid_ny - 1);
    pmgp->hzed = grid_len_z / (grid_nz - 1);
    pmgp->xcent = cx; pmgp->ycent = cy; pmgp->zcent = cz;
    pmgp->xlen = grid_len_x; pmgp->ylen = grid_len_y; pmgp->zlen = grid_len_z;
    pmgp->xmin = cx - grid_len_x/2.0; pmgp->ymin = cy - grid_len_y/2.0; pmgp->zmin = cz - grid_len_z/2.0;
    pmgp->xmax = cx + grid_len_x/2.0; pmgp->ymax = cy + grid_len_y/2.0; pmgp->zmax = cz + grid_len_z/2.0;

    int nlx = compute_nlev(grid_nx), nly = compute_nlev(grid_ny), nlz = compute_nlev(grid_nz);
    pmgp->nlev = nlx < nly ? (nlx < nlz ? nlx : nlz) : (nly < nlz ? nly : nlz);
    pmgp->nonlin = 0;
    pmgp->ipkey = -1;
    pmgp->bcfl = BCFL_SDH;
    pmgp->meth = 2;

    Vpmg* pmg = Vpmg_ctor(pmgp, pbe, 0, VNULL, mgparm, PCE_NO);
    if (!pmg) {
        Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm); Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpmg_ctor failed");
        return APBS_BRIDGE_ERROR;
    }

    if (!Vpmg_fillco(pmg, VSM_MOL, 0.3, VCM_BSPL2,
                      0,VNULL, 0,VNULL, 0,VNULL, 0,VNULL, 0,VNULL, 0,VNULL)) {
        Vpmg_dtor(&pmg); Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm);
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpmg_fillco failed");
        return APBS_BRIDGE_ERROR;
    }

    if (!Vpmg_solve(pmg)) {
        Vpmg_dtor(&pmg); Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm);
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpmg_solve failed");
        return APBS_BRIDGE_ERROR;
    }

    int nx = pmgp->nx, ny = pmgp->ny, nz = pmgp->nz;
    int n_points = nx * ny * nz;
    result->dims[0] = nx; result->dims[1] = ny; result->dims[2] = nz;
    result->n_points = n_points;
    result->spacing[0] = pmgp->hx; result->spacing[1] = pmgp->hy; result->spacing[2] = pmgp->hzed;
    result->origin[0] = pmgp->xcent - (nx-1)*pmgp->hx/2.0;
    result->origin[1] = pmgp->ycent - (ny-1)*pmgp->hy/2.0;
    result->origin[2] = pmgp->zcent - (nz-1)*pmgp->hzed/2.0;

    result->data = (double*)malloc(n_points * sizeof(double));
    if (!result->data) {
        Vpmg_dtor(&pmg); Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm);
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "malloc failed");
        return APBS_BRIDGE_ERROR;
    }

    if (!Vpmg_fillArray(pmg, result->data, VDT_POT, 0.0, PBE_LPBE, NULL)) {
        free(result->data); result->data = NULL;
        Vpmg_dtor(&pmg); Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm);
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "fillArray failed");
        return APBS_BRIDGE_ERROR;
    }

    Vpmg_dtor(&pmg); Vpmgp_dtor(&pmgp); MGparm_dtor(&mgparm);
    Vpbe_dtor(&pbe); Valist_dtor(&alist);
    return APBS_BRIDGE_OK;
}
