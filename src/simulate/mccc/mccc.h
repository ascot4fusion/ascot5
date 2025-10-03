/**
 * @file mccc.h
 * @brief Header file for mccc package
 */
#ifndef MCCC_H
#define MCCC_H

#include "B_field.h"
#include "ascot5.h"
#include "mccc_wiener.h"
#include "particle.h"
#include "plasma.h"
#include "random.h"

/**
 * @brief Defines minimum energy boundary condition
 *
 * This times local electron temperature is minimum energy boundary. If guiding
 * center energy goes below this, it is mirrored to prevent collision
 * coefficients from diverging.
 */
#define MCCC_CUTOFF 0.1

void mccc_init(
    mccc_data *mdata, int include_energy, int include_pitch,
    int include_gcdiff);
void mccc_fo_euler(
    particle_simd_fo *p, real *h, plasma_data *pdata, mccc_data *mdata,
    real *rnd);
void mccc_gc_euler(
    particle_simd_gc *p, real *h, B_field_data *Bdata, plasma_data *pdata,
    mccc_data *mdata, real *rnd);
void mccc_gc_milstein(
    particle_simd_gc *p, real *hin, real *hout, real tol, mccc_wienarr *w,
    B_field_data *Bdata, plasma_data *pdata, mccc_data *mdata, real *rnd);

#endif
