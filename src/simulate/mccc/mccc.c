/**
 * @file mccc.c
 * @brief Interface for using mccc package within ascot5
 */
#include <stdlib.h>
#include <math.h>
#include "mccc.h"

/**
 * @brief Set collision operator data.
 *
 * @param mdata pointer to collision operator data struct
 * @param include_energy can collisions change marker energy, either 0 or 1
 * @param include_pitch  can collisions change marker pitch, either 0 or 1
 * @param include_gcdiff can collisions change GC position, either 0 or 1
 */
void mccc_init(mccc_data* mdata, int include_energy, int include_pitch,
               int include_gcdiff) {
    mdata->include_energy = include_energy;
    mdata->include_pitch  = include_pitch;
    mdata->include_gcdiff = include_gcdiff;

    mdata->usetabulated = 0;
}