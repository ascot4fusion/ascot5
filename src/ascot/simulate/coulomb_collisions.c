/**
 * @file mccc.c
 * @brief Interface for using mccc package within ascot5
 */
#include "coulomb_collisions.h"
#include <math.h>
#include <stdlib.h>

void mccc_init(
    mccc_data *mdata, int include_energy, int include_pitch, int include_gcdiff)
{
    mdata->include_energy = include_energy;
    mdata->include_pitch = include_pitch;
    mdata->include_gcdiff = include_gcdiff;

    mdata->usetabulated = 0;
}
