/**
 * @file nbi.h
 * @brief Header file for nbi.c
 */
#ifndef NBI_H
#define NBI_H

#include "ascot5.h"
#include "B_field.h"
#include "plasma.h"
#include "random.h"

typedef struct {
    int id;
    int n_beamlet;
    real* beamlet_x;
    real* beamlet_y;
    real* beamlet_z;
    real* beamlet_dirx;
    real* beamlet_diry;
    real* beamlet_dirz;
    real power;
    real energy;
    real efrac[3];
    real divergence[3];
    int anum;
    int znum;
} nbi_injector;

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, real* anum, real* znum, random_data* rng);
void nbi_ionize(real* xyz, real* vxyz, int anum, int znum, B_field_data* Bdata, plasma_data* plsdata, random_data* rng);

#endif
