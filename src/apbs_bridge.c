/*
 * APBS C bridge — includes APBS/FETK headers that conflict with Eigen.
 * Compiled as C.
 *
 * Single solve on the COARSE grid with SDH boundary conditions.
 * The coarse grid (extent + 70Å) puts the boundary far enough from
 * the molecule for SDH to be accurate.
 *
 * TODO: implement two-level focusing (coarse→fine) to match APBS mg-auto.
 * Vpmg_ctor(focusFlag=1) crashes in the Debian APBS 3.4.1 build.
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
    double pdie, double sdie, double temperature, double ionic_strength,
    int grid_nx, int grid_ny, int grid_nz,
    double fine_x, double fine_y, double fine_z,
    double coarse_x, double coarse_y, double coarse_z,
    ApbsGridResult* result)
{
    memset(result, 0, sizeof(ApbsGridResult));
    if (n_atoms <= 0) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "n_atoms must be > 0, got %d", n_atoms);
        return APBS_BRIDGE_ERROR;
    }

    Valist* alist = build_atom_list(n_atoms, x, y, z, charges, radii);
    if (!alist) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Failed to build atom list (%d atoms)", n_atoms);
        return APBS_BRIDGE_ERROR;
    }

    double ionConc[2]  = { ionic_strength, ionic_strength };
    double ionRadii[2] = { 0.95, 1.81 };
    double ionQ[2]     = { 1.0, -1.0 };

    Vpbe* pbe = Vpbe_ctor(alist, 2, ionConc, ionRadii, ionQ,
        temperature, pdie, sdie, 1.4, 0, 10.0, 0.0, 0.0, 1.0, 0.0);
    if (!pbe) {
        snprintf(result->error_msg, sizeof(result->error_msg), "Vpbe_ctor failed");
        Valist_dtor(&alist); return APBS_BRIDGE_ERROR;
    }

    double cx = Valist_getCenterX(alist);
    double cy = Valist_getCenterY(alist);
    double cz = Valist_getCenterZ(alist);

    /* Use COARSE grid extent for single solve — SDH boundary is accurate
     * at ~35Å from the molecule (half of the 70Å padding). */
    MGparm* mgparm = MGparm_ctor(MCT_MANUAL);
    if (!mgparm) {
        Vpbe_dtor(&pbe); Valist_dtor(&alist);
        snprintf(result->error_msg, sizeof(result->error_msg), "MGparm_ctor failed");
        return APBS_BRIDGE_ERROR;
    }

    mgparm->dime[0] = grid_nx; mgparm->dime[1] = grid_ny; mgparm->dime[2] = grid_nz;
    mgparm->setdime = 1;
    mgparm->glen[0] = coarse_x; mgparm->glen[1] = coarse_y; mgparm->glen[2] = coarse_z;
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
    pmgp->hx = coarse_x / (grid_nx - 1);
    pmgp->hy = coarse_y / (grid_ny - 1);
    pmgp->hzed = coarse_z / (grid_nz - 1);
    pmgp->xcent = cx; pmgp->ycent = cy; pmgp->zcent = cz;
    pmgp->xlen = coarse_x; pmgp->ylen = coarse_y; pmgp->zlen = coarse_z;
    pmgp->xmin = cx - coarse_x/2.0; pmgp->ymin = cy - coarse_y/2.0; pmgp->zmin = cz - coarse_z/2.0;
    pmgp->xmax = cx + coarse_x/2.0; pmgp->ymax = cy + coarse_y/2.0; pmgp->zmax = cz + coarse_z/2.0;

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
