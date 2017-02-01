/**
 * @file particle.h
 * @brief Header file for particle.c
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include "ascot5.h"
#include "B_field.h"

/** @brief Struct representing a single physical particle in cylindrical
 *         coordinates.
 *
 * This is primarily used for particle initialization, this will be converted
 * into SIMD structs in appropriate coordinates for the simulation.
 */
typedef struct {
    real r;        /**< r coordinate */
    real phi;      /**< phi coordinate */
    real z;        /**< z coordinate */
    real rdot;     /**< r velocity */
    real phidot;   /**< phi velocity */
    real zdot;     /**< z velocity */
    real mass;     /**< mass */
    real charge;   /**< charge */
    real weight;   /**< test particle weight */
    real time;     /**< particle simulation time */
    integer id;        /**< arbitrary id for the particle */
    integer running;  /**< 1 if the particle has hit the wall */
    integer endcond;   /**< particle end condition */
    integer walltile;  /**< id of walltile if particle hit the wall */
} particle;

/** @brief Struct representing a single guiding center in cylindrical
 *         coordinates
 */
typedef struct {
    real r;        /**< r coordinate */
    real phi;      /**< phi coordinate */
    real z;        /**< z coordinate */
    real vpar;     /**< parallel velocity */
    real mu;       /**< magnetic moment */
    real theta;    /**< gyroangle */
    real mass;     /**< mass */
    real charge;   /**< charge */
    real weight;   /**< test particle weight */
    real time;     /**< particle simulation time */
    integer id;        /**< arbitrary id for the particle */
    integer running;  /**< 1 if the particle has hit the wall */
    integer endcond;   /**< particle end condition */
    integer walltile;  /**< id of walltile if particle hit the wall */
} particle_gc;

/** @brief Struct representing a group of NSIMD particles in full orbit
 *         cylindrical coordinates.
 *
 * Convention: whenever one particle property is changed, it is the
 * responsibility of the function performing the change to make sure all
 * parameters are consistent (for example, update magnetic field when position
 * changes.
 */
typedef struct {
    real r[NSIMD] __memalign__;        /**< r coordinate */
    real phi[NSIMD] __memalign__;      /**< phi coordinate */
    real z[NSIMD] __memalign__;        /**< z coordinate */
    real rdot[NSIMD] __memalign__;     /**< r velocity */
    real phidot[NSIMD] __memalign__;   /**< phi velocity */
    real zdot[NSIMD] __memalign__;     /**< z velocity */
    real mass[NSIMD] __memalign__;     /**< mass */
    real charge[NSIMD] __memalign__;   /**< charge */
    real weight[NSIMD] __memalign__;   /**< test particle weight */
    real time[NSIMD] __memalign__;     /**< particle simulation time */
    integer id[NSIMD] __memalign__;       /**< arbitrary id for the particle */
    integer running[NSIMD] __memalign__;
    integer endcond[NSIMD] __memalign__;  /**< particle end condition */
    integer walltile[NSIMD] __memalign__; /**< id of walltile if particle hit
                                               the wall */
    real B_r[NSIMD] __memalign__;      /**< magnetic field r component at
                                            particle position */
    real B_phi[NSIMD] __memalign__;    /**< magnetic field phi component at
                                            particle position */
    real B_z[NSIMD] __memalign__;      /**< magnetic field z component at
                                            particle position */
    real prev_r[NSIMD] __memalign__;      /**< previous r coordinate */
    real prev_phi[NSIMD] __memalign__;    /**< previous phi coordinate */
    real prev_z[NSIMD] __memalign__;      /**< previous z coordinate */
    integer index[NSIMD] __memalign__;
} particle_simd_fo;

/** @brief Struct representing a group of NSIMD particles in guiding center 
 *         polar coordinates.
 
 * Convention: whenever one particle property is changed, it is the
 * responsibility of the function performing the change to make sure all
 * parameters are consistent (for example, update magnetic field when position
 * changes.
 */
typedef struct {
    real r[NSIMD] __memalign__;        /**< r coordinate */
    real phi[NSIMD] __memalign__;      /**< phi coordinate */
    real z[NSIMD] __memalign__;        /**< z coordinate */
    real vpar[NSIMD] __memalign__;     /**< parallel velocity */
    real mu[NSIMD] __memalign__;       /**< magnetic moment */
    real theta[NSIMD] __memalign__;    /**< gyroangle */
    real mass[NSIMD] __memalign__;     /**< mass */
    real charge[NSIMD] __memalign__;   /**< charge */
    real weight[NSIMD] __memalign__;   /**< test particle weight */
    real time[NSIMD] __memalign__;     /**< particle simulation time */
    integer id[NSIMD] __memalign__;       /**< arbitrary id for the particle */
    integer running[NSIMD] __memalign__; /**< 1 if the particle has hit the
                                               wall */
    integer endcond[NSIMD] __memalign__;  /**< particle end condition */
    integer walltile[NSIMD] __memalign__; /**< id of walltile if particle hit
                                               the wall */
    real B_r[NSIMD] __memalign__;      /**< magnetic field r component at
                                            particle position */
    real B_phi[NSIMD] __memalign__;    /**< magnetic field phi component at
                                            particle position */
    real B_z[NSIMD] __memalign__;      /**< magnetic field z component at
                                            particle position */
    real prev_r[NSIMD] __memalign__;      /**< previous r coordinate */
    real prev_phi[NSIMD] __memalign__;    /**< previous phi coordinate */
    real prev_z[NSIMD] __memalign__;      /**< previous z coordinate */
    integer index[NSIMD] __memalign__;
} particle_simd_gc;

#pragma omp declare target
void particle_to_gc(particle* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata);
void particle_to_gc_dummy(particle_simd_gc* p_gc, int j);
void gc_to_particle(particle_simd_gc* p_gc, int j, particle* p);
#pragma omp end declare target

#endif
