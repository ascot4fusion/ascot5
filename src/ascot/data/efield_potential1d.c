/**
 * Implements efield_potential1d.h.
 */
#include "efield_potential1d.h"
#include "bfield.h"
#include "defines.h"
#include "efield.h"
#include "utils/interp.h"
#include <stdio.h>
#include <stdlib.h>

int EfieldPotential1D_init(
    EfieldPotential1D *efield, size_t nrho, real rholim[2], real dvdrho[nrho])
{
    int err = 0;
    real *c = (real *)malloc(nrho * sizeof(real));
    for (size_t i = 0; i < nrho; i++)
    {
        c[i] = dvdrho[i];
    }
    Linear1D_init(&efield->dvdrho, nrho, NATURALBC, rholim, c);
    return err;
}

void EfieldPotential1D_free(EfieldPotential1D *efield)
{
    free(efield->dvdrho.c);
}

void EfieldPotential1D_offload(EfieldPotential1D *efield)
{
    SUPPRESS_UNUSED_WARNING(efield);
    GPU_MAP_TO_DEVICE(
        efield->dvdrho,
        efield->dvdrho.c [0:efield->dvrho.v * efield->dvrho.n_x]);
}

err_t EfieldPotential1D_eval_e(
    real e[3], real r, real phi, real z, real t, EfieldPotential1D *efield,
    Bfield *bfield)
{
    err_t err = 0, err2;
    real psi_dpsi[4], rho_drho[4];
    err = Bfield_eval_psi_dpsi(psi_dpsi, r, phi, z, t, bfield);
    err2 = Bfield_eval_rho_drho(rho_drho, psi_dpsi, bfield);
    err = err ? err : err2;
    rho_drho[2] = rho_drho[2] / r;

    real dvdrho;
    int interperr = Linear1D_eval_f(&dvdrho, &efield->dvdrho, rho_drho[0]);
    e[0] = -dvdrho * rho_drho[1];
    e[1] = -dvdrho * rho_drho[2];
    e[2] = -dvdrho * rho_drho[3];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_EFIELD_POTENTIAL1D_C);
    return err;
}
