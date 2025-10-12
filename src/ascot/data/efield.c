/**
 * Implements efield.h.
 */
#include "efield.h"
#include "bfield.h"
#include "defines.h"
#include "efield_cartesian.h"
#include "efield_potential1d.h"
#include <stdio.h>

void Efield_free(Efield *efield)
{
    switch (efield->type)
    {
    case EFIELD_CARTESIAN:
        EfieldCartesian_free(efield->cartesian);
        break;
    case EFIELD_POTENTIAL1D:
        EfieldPotential1D_free(efield->potential1d);
        break;
    }
}

void Efield_offload(Efield *efield)
{
    switch (efield->type)
    {
    case EFIELD_CARTESIAN:
        EfieldCartesian_offload(efield->cartesian);
        break;
    case EFIELD_POTENTIAL1D:
        EfieldPotential1D_offload(efield->potential1d);
        break;
    }
}

err_t Efield_eval_e(
    real e[3], real r, real phi, real z, real t, Efield *efield, Bfield *bfield)
{
    err_t err = 0;
    switch (efield->type)
    {
    case EFIELD_CARTESIAN:
        err = EfieldCartesian_eval_e(e, phi, efield->cartesian);
        break;
    case EFIELD_POTENTIAL1D:
        err = EfieldPotential1D_eval_e(
            e, r, phi, z, t, efield->potential1d, bfield);
        break;
    default:
        /* If electric field is not given it evaluates to zero. */
        e[0] = 0;
        e[1] = 0;
        e[2] = 0;
        break;
    }

    e[0] = err ? 0.0 : e[0];
    e[1] = err ? 0.0 : e[1];
    e[2] = err ? 0.0 : e[2];

    return err;
}
