/**
 * @file nbi.h
 * @brief Header file for nbi.c
 */
#ifndef NBI_H
#define NBI_H

#include "defines.h"
#include "particle.h"
#include "random.h"

/**
 * @brief Structure for describing an NBI injector
 */
typedef struct
{
    int id;             /**< integer identifier for this injector           */
    int n_beamlet;      /**< number of beamlets in this injector            */
    real *beamlet_x;    /**< x coordinates of beamlets [m]                  */
    real *beamlet_y;    /**< y coordinates of beamlets [m]                  */
    real *beamlet_z;    /**< z coordinates of beamlets [m]                  */
    real *beamlet_dx;   /**< x components of beamlet unit direction vectors */
    real *beamlet_dy;   /**< y components of beamlet unit direction vectors */
    real *beamlet_dz;   /**< z components of beamlet unit direction vectors */
    real power;         /**< this injector's power injected [W]             */
    real energy;        /**< full energy of injected particles [J]          */
    real efrac[3];      /**< fractions of full, 1/2 and 1/3 energy          */
    real div_h;         /**< vertical divergence [rad]                      */
    real div_v;         /**< horizontal divergence [rad]                    */
    real div_halo_frac; /**< fraction of power in the halo                  */
    real div_halo_h;    /**< horizontal divergence of the halo [rad]        */
    real div_halo_v;    /**< vertical divergence of the halo [rad]          */
    int anum;           /**< mass number of injected species                */
    int znum;           /**< charge number of injected species              */
    real mass;          /**< mass of injected species [kg]                  */
} nbi_injector;

/**
 * @brief NBI data consisting of `ninj` injectors.
 */
typedef struct
{
    int ninj;          /**< number of injectors */
    nbi_injector *inj; /**< array of injectors  */
} Nbi;

int nbi_init(
    Nbi *nbi, int ninj, int *id, int *anum, int *znum, real *mass, real *power,
    real *efrac, real *energy, real *div_h, real *div_v, real *div_halo_v,
    real *div_halo_h, real *div_halo_frac, int *nbeamlet, real *beamlet_xyz);
void nbi_free(Nbi *nbi);
void nbi_inject(real *xyz, real *vxyz, nbi_injector *inj, random_data *rng);

#endif
