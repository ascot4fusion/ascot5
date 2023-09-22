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
#define MAX_NBI_INJ 10

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

/**
 * @brief NBI parameters consisting of a bundle of injectors.
 */
typedef struct {
    int ninj; /**< number of injectors */
    int n_beamlet[MAX_NBI_INJ]; /**< number of beamlets in a given injector */
    int id[MAX_NBI_INJ];
    real mass[MAX_NBI_INJ];    /**< plasma species masses [kg]          */
    int anum[MAX_NBI_INJ];     /**< ion species atomic number           */
    int znum[MAX_NBI_INJ];     /**< ion species charge number           */
    real power[MAX_NBI_INJ];
    real energy[MAX_NBI_INJ];
    real efrac[MAX_NBI_INJ*3];
    real div_h[MAX_NBI_INJ];
    real div_v[MAX_NBI_INJ];
    real div_halo_frac[MAX_NBI_INJ];
    real div_halo_h[MAX_NBI_INJ];
    real div_halo_v[MAX_NBI_INJ];
    int offload_array_length;  /**< number of elements in offload_array */
} nbi_offload_data;

/**
 * @brief NBI data on target.
 */
typedef struct {
    int ninj;          /**< number of injectors */
    nbi_injector inj[MAX_NBI_INJ]; /**< array of injectors */
} nbi_data;

int nbi_init_offload(nbi_offload_data* offload_data, real** offload_array);

void nbi_init(nbi_data* nbi, nbi_offload_data* offload_data,
              real* offload_array);

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, int* anum, int* znum, real* mass, random_data* rng);
void nbi_ionize(real* xyz, real* vxyz, real time, int* shinethrough, int anum, int znum,
                B_field_data* Bdata, plasma_data* plsdata, wall_data* walldata,
                random_data* rng);
void nbi_generate(int nprt, real t0, real t1, particle* p, nbi_injector* n,
                  B_field_data* Bdata, plasma_data* plsdata,
                  wall_data* walldata, random_data* rng);

#endif
