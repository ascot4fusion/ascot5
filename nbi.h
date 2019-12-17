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

/**
 * @brief Structure for describing an NBI injector
 */
typedef struct {
    int id;             /**< Integer identifier for this injector */
    int n_beamlet;      /**< Number of beamlets in this injector */
    real* beamlet_x;    /**< X coordinates of beamlets [m] */
    real* beamlet_y;    /**< Y coordinates of beamlets [m] */
    real* beamlet_z;    /**< Z coordinates of beamlets [m] */
    real* beamlet_dx;   /**< X components of beamlet unit direction vectors */
    real* beamlet_dy;   /**< Y components of beamlet unit direction vectors */
    real* beamlet_dz;   /**< Z components of beamlet unit direction vectors */
    real power;         /**< Total power injected [W] */
    real energy;        /**< Full energy of injected particles [J] */
    real efrac[3];      /**< Fractions of full, half and one-third energy */
    real div_h;         /**< Vertical divergence [radians] */
    real div_v;         /**< Horizontal divergence [radians] */
    real div_halo_frac; /**< Fraction of power in the halo */
    real div_halo_h;    /**< Horizontal divergence of the halo [radians] */
    real div_halo_v;    /**< Vertical divergence of the halo [radians] */
    int anum;           /**< Mass number of injected species */
    int znum;           /**< Charge number of injected species */
    real mass;          /**< Mass of injected species */
} nbi_injector;

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, int* anum, int* znum, real* mass, random_data* rng);
void nbi_ionize(real* xyz, real* vxyz, int* shinethrough, int anum, int znum,
                B_field_data* Bdata, plasma_data* plsdata, wall_data* walldata,
                random_data* rng);
void nbi_generate(int nprt, particle* p, nbi_injector* n,
                  B_field_data* Bdata, plasma_data* plsdata,
                  wall_data* walldata, random_data* rng);

#endif
