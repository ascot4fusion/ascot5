/**
 * @file nbi.h
 * @brief Header file for nbi.c
 */
#ifndef NBI_H
#define NBI_H

#include "ascot5.h"
#include "particle.h"
#include "random.h"

/**
 * @brief Structure for describing an NBI injector
 */
typedef struct {
    int id;             /**< integer identifier for this injector           */
    int n_beamlet;      /**< number of beamlets in this injector            */
    real* beamlet_x;    /**< x coordinates of beamlets [m]                  */
    real* beamlet_y;    /**< y coordinates of beamlets [m]                  */
    real* beamlet_z;    /**< z coordinates of beamlets [m]                  */
    real* beamlet_dx;   /**< x components of beamlet unit direction vectors */
    real* beamlet_dy;   /**< y components of beamlet unit direction vectors */
    real* beamlet_dz;   /**< z components of beamlet unit direction vectors */
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
 * @brief NBI parameters consisting of a bundle of injectors.
 */
typedef struct {
    int ninj; /**< number of injectors */
    int id[NBI_MAX_INJ];          /**< integer identifier for the injectors   */
    int n_beamlet[NBI_MAX_INJ];   /**< number of beamlets in a given injector */
    real power[NBI_MAX_INJ];      /**< injected power [W]                     */
    real energy[NBI_MAX_INJ];     /**< full energy of injected particles [J]  */
    real efrac[NBI_MAX_INJ*3];    /**< fractions of full, 1/2 and 1/3 energy  */
    real div_h[NBI_MAX_INJ];      /**< vertical divergence [rad]              */
    real div_v[NBI_MAX_INJ];      /**< horizontal divergence [rad]            */
    real div_halo_frac[NBI_MAX_INJ]; /**< fraction of power in the halo       */
    real div_halo_h[NBI_MAX_INJ]; /**< horizontal divergence of the halo [rad]*/
    real div_halo_v[NBI_MAX_INJ]; /**< vertical divergence of the halo [rad]  */
    int anum[NBI_MAX_INJ];        /**< mass number of injected species        */
    int znum[NBI_MAX_INJ];        /**< charge number of injected species      */
    real mass[NBI_MAX_INJ];       /**< mass of injected species [kg]          */
    int offload_array_length;     /**< number of elements in offload_array    */
} nbi_offload_data;

/**
 * @brief NBI data on target.
 */
typedef struct {
    int ninj;                      /**< number of injectors */
    nbi_injector inj[NBI_MAX_INJ]; /**< array of injectors  */
} nbi_data;

int nbi_init_offload(nbi_offload_data* offload_data, real** offload_array);
void nbi_init(nbi_data* nbi, nbi_offload_data* offload_data,
              real* offload_array);
void nbi_free_offload(nbi_offload_data* offload_data, real** offload_array);
void nbi_inject(real* xyz, real* vxyz, nbi_injector* inj, random_data* rng);

#endif
