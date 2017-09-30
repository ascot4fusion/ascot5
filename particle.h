/**
 * @file particle.h
 * @brief Header file for particle.c
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include "ascot5.h"
#include "B_field.h"
#include "E_field.h"
#include "error.h"

/** @brief Struct the physical state of the particle. Includes self-consistent
 *         full-orbit and guiding-center coordinates.
 *
 * Convention: whenever one particle property is changed, it is the
 * responsibility of the function performing the change to make sure all
 * parameters are consistent (for example, update magnetic field when position
 * changes.
 *
 */
typedef struct {
    real r;           /**< guiding center r coordinate */
    real phi;         /**< guiding center phi coordinate */
    real z;           /**< guiding center z coordinate */
    real vpar;        /**< parallel velocity */
    real mu;          /**< magnetic moment */
    real theta;       /**< gyroangle */
    real rprt;        /**< particle r coordinate */
    real phiprt;      /**< particle phi coordinate */
    real zprt;        /**< particle z coordinate */
    real rdot;        /**< dr/dt */
    real phidot;      /**< dphi/dt */
    real zdot;        /**< dz/dt */
    real mass;        /**< mass */
    real charge;      /**< charge */
    real weight;      /**< test particle weight */
    real time;        /**< particle simulation time */
    real cputime;     /**< test particle cputime */
    real rho;         /**< test particle rho coordinate */
    real pol;         /**< test particle poloidal coordinate */
    integer id;       /**< arbitrary id for the particle */
    integer endcond;  /**< particle end condition */
    integer walltile; /**< id of walltile if particle hit
                                               the wall */
    real B_r;         /**< magnetic field r component at
                                            particle position */
    real B_phi;       /**< magnetic field phi component at
                                            particle position */
    real B_z;         /**< magnetic field z component at
                                            particle position */
    real B_r_dr;      /**< gradient of B_r with respect to r*/
    real B_phi_dr;    /**< gradient of B_phi with respect to r*/
    real B_z_dr;      /**< gradient of B_z with respect to r*/
    real B_r_dphi;    /**< gradient of B_r with respect to phi*/
    real B_phi_dphi;  /**< gradient of B_phi with respect to phi*/
    real B_z_dphi;    /**< gradient of B_z with respect to phi*/
    real B_r_dz;      /**< gradient of B_r with respect to z*/
    real B_phi_dz;    /**< gradient of B_phi with respect to z*/
    real B_z_dz;      /**< gradient of B_z with respect to z*/
    
    a5err err;        /**< error flag */
} particle_state;

/** @brief Struct representing a single physical particle in cylindrical
 *         coordinates.
 *
 * This is primarily used for particle initialization, this will be converted
 * into SIMD structs in appropriate coordinates for the simulation.
 *
 * @todo rdot,phidot, and zdot to v_r, v_phi and v_z? (Only affects phidot)
 */
typedef struct {
    real r;        /**< r coordinate */
    real phi;      /**< phi coordinate */
    real z;        /**< z coordinate */
    real v_r;     /**< r velocity */
    real v_phi;   /**< phi velocity */
    real v_z;     /**< z velocity */
    real mass;     /**< mass */
    real charge;   /**< charge */
    real weight;   /**< test particle weight */
    real time;     /**< particle simulation time */
    integer id;        /**< arbitrary id for the particle */
} particle;

/** @brief Struct representing a single guiding center in cylindrical
 *         coordinates
 */
typedef struct {
    real r;        /**< r coordinate */
    real phi;      /**< phi coordinate */
    real z;        /**< z coordinate */
    real energy;   /**< kinetic energy */
    real pitch;    /**< pitch angle */
    real theta;    /**< gyroangle */
    real mass;     /**< mass */
    real charge;   /**< charge */
    real weight;   /**< test particle weight */
    real time;     /**< particle simulation time */
    integer id;        /**< arbitrary id for the particle */
} particle_gc;

/** @brief Struct representing a single magnetic field line in cylindrical
 *         coordinates
 */
typedef struct {
    real r;        /**< r coordinate */
    real phi;      /**< phi coordinate */
    real z;        /**< z coordinate */
    real pitch;    /**< pitch angle */
    real weight;   /**< test particle weight */
    real time;     /**< particle simulation time */
    integer id;    /**< arbitrary id for the particle */
} particle_ml;

typedef struct {
    int n;
    particle_state** p;
    int next;
} particle_queue;

typedef enum input_particle_type {
    input_particle_type_p, input_particle_type_gc, input_particle_type_ml, input_particle_type_s
} input_particle_type;

/** @brief Struct representing either a particle or guiding center, including
 *         its type.
 */
typedef struct {
    input_particle_type type;
    union {
        particle p;
        particle_gc p_gc;
	particle_ml p_ml;
	particle_state p_s;
    };
} input_particle;

/** @brief Struct representing a group of NSIMD particles in full orbit
 *         cylindrical coordinates.
 *
 * Convention: whenever one particle property is changed, it is the
 * responsibility of the function performing the change to make sure all
 * parameters are consistent (for example, update magnetic field when position
 * changes.
 *
 */
typedef struct {
    real r[NSIMD] __memalign__;        /**< r coordinate */
    real phi[NSIMD] __memalign__;      /**< phi coordinate */
    real z[NSIMD] __memalign__;        /**< z coordinate */
    real rdot[NSIMD] __memalign__;     /**< dr/dt */
    real phidot[NSIMD] __memalign__;   /**< dphi/dt */
    real zdot[NSIMD] __memalign__;     /**< dz/dt */
    real mass[NSIMD] __memalign__;     /**< mass */
    real charge[NSIMD] __memalign__;   /**< charge */
    real weight[NSIMD] __memalign__;   /**< test particle weight */
    real time[NSIMD] __memalign__;     /**< particle simulation time */
    real cputime[NSIMD] __memalign__;  /**< particle cpu time */
    real rho[NSIMD] __memalign__;      /**< particle rho coordinate */
    real pol[NSIMD] __memalign__;      /**< particle poloidal coordinate */
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

    real B_r_dr[NSIMD] __memalign__;      /**< gradient of B_r with respect to r*/
    real B_phi_dr[NSIMD] __memalign__;    /**< gradient of B_phi with respect to r*/
    real B_z_dr[NSIMD] __memalign__;      /**< gradient of B_z with respect to r*/
    real B_r_dphi[NSIMD] __memalign__;    /**< gradient of B_r with respect to phi*/
    real B_phi_dphi[NSIMD] __memalign__;  /**< gradient of B_phi with respect to phi*/
    real B_z_dphi[NSIMD] __memalign__;    /**< gradient of B_z with respect to phi*/
    real B_r_dz[NSIMD] __memalign__;      /**< gradient of B_r with respect to z*/
    real B_phi_dz[NSIMD] __memalign__;    /**< gradient of B_phi with respect to z*/
    real B_z_dz[NSIMD] __memalign__;      /**< gradient of B_z with respect to z*/

    a5err err[NSIMD] __memalign__;        /**< error flag */

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
    real cputime[NSIMD] __memalign__;  /**< particle cpu time */
    real rho[NSIMD] __memalign__;      /**< particle rho coordinate */
    real pol[NSIMD] __memalign__;      /**< particle poloidal coordinate */
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
    real B_r_dr[NSIMD] __memalign__;      /**< gradient of B_r with respect to r*/
    real B_phi_dr[NSIMD] __memalign__;    /**< gradient of B_phi with respect to r*/
    real B_z_dr[NSIMD] __memalign__;      /**< gradient of B_z with respect to r*/
    real B_r_dphi[NSIMD] __memalign__;    /**< gradient of B_r with respect to phi*/
    real B_phi_dphi[NSIMD] __memalign__;  /**< gradient of B_phi with respect to phi*/
    real B_z_dphi[NSIMD] __memalign__;    /**< gradient of B_z with respect to phi*/
    real B_r_dz[NSIMD] __memalign__;      /**< gradient of B_r with respect to z*/
    real B_phi_dz[NSIMD] __memalign__;    /**< gradient of B_phi with respect to z*/
    real B_z_dz[NSIMD] __memalign__;      /**< gradient of B_z with respect to z*/

    a5err err[NSIMD] __memalign__;        /**< error flag */

    integer index[NSIMD] __memalign__;
} particle_simd_gc;

/** @brief Struct representing a group of NSIMD magnetic field lines in 
 *         cylinder coordinates.
 *
 */
typedef struct {
    real r[NSIMD] __memalign__;        /**< r coordinate */
    real phi[NSIMD] __memalign__;      /**< phi coordinate */
    real z[NSIMD] __memalign__;        /**< z coordinate */
    real pitch[NSIMD] __memalign__;     /**< pitc = +1 */
    real weight[NSIMD] __memalign__;   /**< test particle weight */
    real time[NSIMD] __memalign__;     /**< field line simulation "time" i.e. distance/c */
    real cputime[NSIMD] __memalign__;  /**< particle cpu time */
    real rho[NSIMD] __memalign__;      /**< particle rho coordinate */
    real pol[NSIMD] __memalign__;      /**< particle poloidal coordinate */
    integer id[NSIMD] __memalign__;       /**< arbitrary id for the field line */
    integer running[NSIMD] __memalign__; /**< 1 if the field line has hit the
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
  
    real B_r_dr[NSIMD] __memalign__;      /**< gradient of B_r with respect to r*/
    real B_phi_dr[NSIMD] __memalign__;    /**< gradient of B_phi with respect to r*/
    real B_z_dr[NSIMD] __memalign__;      /**< gradient of B_z with respect to r*/
    real B_r_dphi[NSIMD] __memalign__;    /**< gradient of B_r with respect to phi*/
    real B_phi_dphi[NSIMD] __memalign__;  /**< gradient of B_phi with respect to phi*/
    real B_z_dphi[NSIMD] __memalign__;    /**< gradient of B_z with respect to phi*/
    real B_r_dz[NSIMD] __memalign__;      /**< gradient of B_r with respect to z*/
    real B_phi_dz[NSIMD] __memalign__;    /**< gradient of B_phi with respect to z*/
    real B_z_dz[NSIMD] __memalign__;      /**< gradient of B_z with respect to z*/

    a5err err[NSIMD] __memalign__;        /**< error flag */

    integer index[NSIMD] __memalign__;
} particle_simd_ml;




#pragma omp declare target
void particle_to_fo_dummy(particle_simd_fo* p_fo, int j);
void particle_to_gc_dummy(particle_simd_gc* p_gc, int j);
void particle_to_ml_dummy(particle_simd_ml* p_ml, int j);

int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle);
int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle);
int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      B_field_data* Bdata, int* cycle);

void particle_input_to_state(input_particle* p, particle_state* ps, B_field_data* Bdata);

#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo, int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc, int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml, int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p, B_field_data* Bdata);
#pragma omp declare simd uniform(p_fo,Bdata)
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc, B_field_data* Bdata);
#pragma omp end declare target

#endif
