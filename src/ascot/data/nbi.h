/**
 * @file nbi.h
 *
 */
#ifndef NBI_H
#define NBI_H

#include "defines.h"
#include "marker.h"
#include "nbi_data.h"
#include "utils/random.h"

int nbi_init(
    Nbi *nbi, size_t ninj, size_t *id, int *anum, int *znum, real *mass, real *power,
    real *efrac, real *energy, real *div_h, real *div_v, real *div_halo_v,
    real *div_halo_h, real *div_halo_frac, size_t *nbeamlet, real *beamlet_xyz);
void nbi_free(Nbi *nbi);
void nbi_inject(real *xyz, real *vxyz, Nbi *inj, random_data *rng);

#endif
