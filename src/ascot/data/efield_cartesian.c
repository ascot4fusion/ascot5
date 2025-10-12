/**
 * Implements efield_cartesian.h.
 */
#include "efield_cartesian.h"
#include "bfield.h"
#include "defines.h"
#include "efield.h"
#include "utils/mathlib.h"
#include <stdlib.h>

int EfieldCartesian_init(EfieldCartesian *efield, real exyz[3])
{
    efield->exyz[0] = exyz[0];
    efield->exyz[1] = exyz[1];
    efield->exyz[2] = exyz[2];
    return 0;
}

void EfieldCartesian_free(EfieldCartesian *efield)
{
    // No resources were allocated
    (void)efield;
}

void EfieldCartesian_offload(EfieldCartesian *efield)
{
    SUPPRESS_UNUSED_WARNING(efield);
    GPU_MAP_TO_DEVICE(data->exyz [0:3])
}

err_t EfieldCartesian_eval_e(real e[3], real phi, EfieldCartesian *efield)
{
    math_vec_xyz2rpz(efield->exyz, e, phi);
    return 0;
}
