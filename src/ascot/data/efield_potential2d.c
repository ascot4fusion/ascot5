/**
 * Implements "efield_potential2d.h".
 */
#include "efield_potential2d.h"
#include "defines.h"
#include "efield.h"
#include "parallel.h"
#include "utils/interp.h"
#include <stddef.h>

int EfieldPotential2D_init(
    EfieldPotential2D *efield, size_t nr, size_t nz, real rlim[2], real zlim[2],
    real vpot[nr * nz])
{
    int err = 0;
    err += Spline2D_init(
        &efield->potential, nr, nz, NATURALBC, NATURALBC, rlim, zlim, vpot);
    return err;
}

void EfieldPotential2D_free(EfieldPotential2D *efield)
{
    free(efield->potential.c);
}

void EfieldPotential2D_offload(EfieldPotential2D *efield){
    SUPPRESS_UNUSED_WARNING(efield);
    GPU_MAP_TO_DEVICE(
    efield -> potential.c
    [0:efield->potential.nx * data->efield->potential.ny * NSIZE_COMP2D])
}

err_t EfieldPotential2D_eval_e(
    real e[3], real r, real z, EfieldPotential2D *efield)
{
    real vdv[6];
    int interperr = Spline2D_eval_f_df(vdv, &efield->potential, r, z);
    e[0] = -vdv[1];
    e[1] = 0;
    e[2] = -vdv[2];

    err_t err = 0;
    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_EFIELD_POTENTIAL2D_C);
    return err;
}
