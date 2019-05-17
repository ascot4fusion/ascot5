/**
 * @file nbi.h
 * @brief Header file for nbi.c
 */
#ifndef NBI_H
#define NBI_H

#include "ascot5.h"
#include "B_field.h"
#include "particle.h"
#include "plasma.h"
#include "random.h"
#include "wall.h"

#define NBI_MAX_DISTANCE 100

typedef struct {
    int id;
    int n_beamlet;
    real* beamlet_x;
    real* beamlet_y;
    real* beamlet_z;
    real* beamlet_dx;
    real* beamlet_dy;
    real* beamlet_dz;
    real power;
    real energy;
    real efrac[3];
    real div_h;
    real div_v;
    real div_halo_frac;
    real div_halo_h;
    real div_halo_v;
    int anum;
    int znum;
} nbi_injector;

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, real* anum, real* znum, random_data* rng);
void nbi_ionize(real* xyz, real* vxyz, int* shinethrough, int anum, int znum,
                B_field_data* Bdata, plasma_data* plsdata, wall_data* walldata,
                random_data* rng);
void nbi_generate(int nprt, particle_state* p, int* shinethrough,
                  nbi_injector* n, B_field_data* Bdata, plasma_data* plsdata,
                  wall_data* walldata, random_data* rng);

#endif
