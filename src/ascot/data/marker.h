/**
 * @file marker.h
 * Marker simulation structs and functions to transform between those.
 */
#ifndef MARKER_H
#define MARKER_H

#include "bfield.h"
#include "defines.h"
#include "efield.h"
#include <stdint.h>

/**
 * Marker queue from which markers are pulled and returned for the simulation.
 *
 * Each time a marker has finished simulation, a new marker is chosen from this
 * queue and the old marker's data is updated. Markers are never removed from
 * the queue but an index is kept to mark where the next not yet simulated
 * marker is found. Markers are represented by State struct when they
 * are stored in the queue.
 *
 * Note: The queue can and is accessed by several threads, so make sure each
 * access is thread-safe.
 */
typedef struct
{
    size_t n;                 /**< Total number of markers in this queue.     */
    State **p;                /**< Markers in this queue.                     */
    volatile size_t next;     /**< Index for the next unsimulated marker.     */
    volatile size_t finished; /**< How many markers have finished simulation. */
} MarkerQueue;

/**
 * Vector of simulated particle markers.
 *
 * This struct contains all physical and simulation parameters necessary for
 * a gyro-orbit simulation, and additional parameters to e.g. avoid unnecessary
 * magnetic field evaluations.
 *
 * A function changing any of these parameters is responsible for making sure
 * all fields remain consistent, e.g. if position changes then the magnetic
 * field should be updated at the same time.
 *
 * All parameters are arrays of length NSIMD (CPU) or <number of markers in this
 * process> (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers (to meet NSIMD quota) and for those ``id`` and
 * ``running`` are set to zero.
 */
typedef struct
{
    real *r;           /**< Particle R coordinate [m].                        */
    real *phi;         /**< Particle phi coordinate [rad].                    */
    real *z;           /**< Particle z coordinate [m].                        */
    real *p_r;         /**< Momentum r coordinate [kg m/s].                   */
    real *p_phi;       /**< Momentum phi coordinate [kg m/s].                 */
    real *p_z;         /**< Momentum z coordinate [kg m/s].                   */
    real *mass;        /**< Mass [kg].                                        */
    real *charge;      /**< Charge [C].                                       */
    real *time;        /**< Marker simulation time [s].                       */
    int *znum;         /**< Particle atomic number.                           */
    int *anum;         /**< Particle mass number.                             */
    real *B_r;         /**< BR at marker position [T].                        */
    real *B_phi;       /**< Bphi at marker position [T].                      */
    real *B_z;         /**< Bz at marker position [T].                        */
    real *B_r_dr;      /**< dB_R/dR at marker position [T/m].                 */
    real *B_phi_dr;    /**< dB_phi/dR at marker position [T/m].               */
    real *B_z_dr;      /**< dB_z/dR at marker position [T/m].                 */
    real *B_r_dphi;    /**< dB_R/dphi at marker position [T/rad].             */
    real *B_phi_dphi;  /**< dB_phi/dphi at marker position [T/rad].           */
    real *B_z_dphi;    /**< dB_z/dphi at marker position [T/rad].             */
    real *B_r_dz;      /**< dB_R/dz at marker position [T/m].                 */
    real *B_phi_dz;    /**< dB_phi/dz at marker position [T/m].               */
    real *B_z_dz;      /**< dB_z/dz at marker position [T/m].                 */
    int *bounces;      /**< Number of times pitch sign changed.               */
    real *weight;      /**< Marker weight.                                    */
    real *cputime;     /**< Marker wall time [s].                             */
    real *rho;         /**< Marker rho coordinate [1].                        */
    real *theta;       /**< Marker cumulative poloidal coordinate [rad].      */
    size_t *id;        /**< Unique ID for the marker.                         */
    endcond_t *endcond;/**< Marker end condition.                             */
    size_t *walltile;  /**< Index (>0) of walltile if marker has hit the wall.*/
    real *mileage;     /**< Duration the marker has been simulated [s].       */
    int *running;      /**< Is the marker currently simulated.                */
    err_t *err;        /**< Error flag, zero if no error.                     */
    size_t *index;     /**< This marker's index at the marker queue.          */
    size_t n_mrk;      /**< How many markers this struct contains.            */
} MarkerGyroOrbit;

/**
 * Vector of simulated guiding center markers.
 *
 * This struct contains all physical and simulation parameters necessary for
 * a guiding-center simulation, and additional parameters to e.g. avoid
 * unnecessary magnetic field evaluations.
 *
 * A function changing any of these parameters is responsible for making sure
 * all fields remain consistent, e.g. if position changes then the magnetic
 * field should be updated at the same time.
 *
 * All parameters are arrays of length NSIMD (CPU) or <number of markers in this
 * process> (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers (to meet NSIMD quota) and for those ``id`` and
 * ``running`` are set to zero.
 */
typedef struct
{
    real *r;           /**< Guiding center R coordinate [m].                  */
    real *phi;         /**< Guiding center phi coordinate [rad].              */
    real *z;           /**< Guiding center z coordinate [m].                  */
    real *ppar;        /**< Parallel momentum [kg m/s].                       */
    real *mu;          /**< Magnetic moment [J/T] .                           */
    real *zeta;        /**< Gyroangle [rad].                                  */
    real *mass;        /**< Mass [kg].                                        */
    real *charge;      /**< Charge [C].                                       */
    real *time;        /**< Marker simulation time [s].                       */
    int *znum;         /**< Particle atomic number.                           */
    int *anum;         /**< Particle mass number.                             */
    real *B_r;         /**< BR at marker position [T].                        */
    real *B_phi;       /**< Bphi at marker position [T].                      */
    real *B_z;         /**< Bz at marker position [T].                        */
    real *B_r_dr;      /**< dB_R/dR at marker position [T/m].                 */
    real *B_phi_dr;    /**< dB_phi/dR at marker position [T/m].               */
    real *B_z_dr;      /**< dB_z/dR at marker position [T/m].                 */
    real *B_r_dphi;    /**< dB_R/dphi at marker position [T/rad].             */
    real *B_phi_dphi;  /**< dB_phi/dphi at marker position [T/rad].           */
    real *B_z_dphi;    /**< dB_z/dphi at marker position [T/rad].             */
    real *B_r_dz;      /**< dB_R/dz at marker position [T/m].                 */
    real *B_phi_dz;    /**< dB_phi/dz at marker position [T/m].               */
    real *B_z_dz;      /**< dB_z/dz at marker position [T/m].                 */
    int *bounces;      /**< Number of times pitch sign changed.               */
    real *weight;      /**< Marker weight.                                    */
    real *cputime;     /**< Marker wall time [s].                             */
    real *rho;         /**< Marker rho coordinate [1].                        */
    real *theta;       /**< Marker cumulative poloidal coordinate [rad].      */
    size_t *id;        /**< Unique ID for the marker.                         */
    endcond_t *endcond;/**< Marker end condition.                             */
    size_t *walltile;  /**< Index (>0) of walltile if marker has hit the wall.*/
    real *mileage;     /**< Duration the marker has been simulated [s].       */
    int *running;      /**< Is the marker currently simulated.                */
    err_t *err;        /**< Error flag, zero if no error.                     */
    size_t *index;     /**< This marker's index at the marker queue.          */
    size_t n_mrk;      /**< How many markers this struct contains.            */
} MarkerGuidingCenter;

/**
 * Vector of simulated field line markers.
 *
 * This struct contains all physical and simulation parameters necessary for
 * a field-line simulation, and additional parameters to e.g. avoid unnecessary
 * magnetic field evaluations.
 *
 * A function changing any of these parameters is responsible for making sure
 * all fields remain consistent, e.g. if position changes then the magnetic
 * field should be updated at the same time.
 *
 * All parameters are arrays of length NSIMD (CPU) or <number of markers in this
 * process> (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers (to meet NSIMD quota) and for those ``id`` and
 * ``running`` are set to zero.
 */
typedef struct
{
    real *r;           /**< Field line R coordinate [m].                      */
    real *phi;         /**< Field line phi coordinate [rad].                  */
    real *z;           /**< Field line z coordinate [m].                      */
    int *pitch;        /**< Is the marker traced along or opposite to B-field.*/
    real *time;        /**< The fixed time instant this field line exists [s].*/
    real *B_r;         /**< BR at marker position [T].                        */
    real *B_phi;       /**< Bphi at marker position [T].                      */
    real *B_z;         /**< Bz at marker position [T].                        */
    real *B_r_dr;      /**< dB_R/dR at marker position [T/m].                 */
    real *B_phi_dr;    /**< dB_phi/dR at marker position [T/m].               */
    real *B_z_dr;      /**< dB_z/dR at marker position [T/m].                 */
    real *B_r_dphi;    /**< dB_R/dphi at marker position [T/rad].             */
    real *B_phi_dphi;  /**< dB_phi/dphi at marker position [T/rad].           */
    real *B_z_dphi;    /**< dB_z/dphi at marker position [T/rad].             */
    real *B_r_dz;      /**< dB_R/dz at marker position [T/m].                 */
    real *B_phi_dz;    /**< dB_phi/dz at marker position [T/m].               */
    real *B_z_dz;      /**< dB_z/dz at marker position [T/m].                 */
    int *bounces;      /**< Number of times pitch sign changed.               */
    real *cputime;     /**< Marker wall time [s].                             */
    real *rho;         /**< Marker rho coordinate [1].                        */
    real *theta;       /**< Marker cumulative poloidal coordinate [rad].      */
    size_t *id;        /**< Unique ID for the marker.                         */
    endcond_t *endcond;/**< Marker end condition.                             */
    size_t *walltile;  /**< Index (>0) of walltile if marker has hit the wall.*/
    real *mileage;     /**< Distance the marker has been traced [m].          */
    int *running;      /**< Is the marker currently simulated.                */
    err_t *err;        /**< Error flag, zero if no error.                     */
    size_t *index;     /**< This marker's index at the marker queue.          */
    size_t n_mrk;      /**< How many markers this struct contains.            */
} MarkerFieldLine;

void marker_allocate_go(MarkerGyroOrbit *mrk, size_t nmrk);

void marker_allocate_gc(MarkerGuidingCenter *mrk, size_t nmrk);

void marker_allocate_fl(MarkerFieldLine *mrk, size_t nmrk);

void marker_to_go_dummy(MarkerGyroOrbit *mrk, size_t index);

void marker_to_gc_dummy(MarkerGuidingCenter *mrk, size_t index);

void marker_to_go_dummy(MarkerGyroOrbit *mrk, size_t index);

void marker_to_gc_dummy(MarkerGuidingCenter *mrk, size_t index);

void marker_to_fl_dummy(MarkerFieldLine *mrk, size_t index);

size_t marker_cycle_go(
    MarkerQueue *queue, MarkerGyroOrbit *p, Bfield *bfield, size_t *cycle);

size_t marker_cycle_gc(
    MarkerQueue *queue, MarkerGuidingCenter *p, Bfield *bfield, size_t *cycle);

size_t marker_cycle_fl(
    MarkerQueue *queue, MarkerFieldLine *p, Bfield *bfield, size_t *cycle);

void marker_offload_go(MarkerGyroOrbit *p);

void marker_onload_go(MarkerGyroOrbit *p);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
err_t state_to_go(
    State *p, size_t i, MarkerGyroOrbit *p_fo, size_t index, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
void marker_go_to_state(
    MarkerGyroOrbit *p_fo, size_t index, State *p, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
err_t state_to_gc(
    State *p, size_t i, MarkerGuidingCenter *p_gc, size_t index, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
void marker_gc_to_state(
    MarkerGuidingCenter *p_gc, size_t index, State *p, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
err_t state_to_fl(
    State *p, size_t i, MarkerFieldLine *p_ml, size_t index, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(bfield)
void marker_fl_to_state(MarkerFieldLine *p_ml, size_t index, State *p);

DECLARE_TARGET_SIMD_UNIFORM(p_fo, bfield)
int marker_go_to_gc(
    MarkerGyroOrbit *p_fo, size_t index, MarkerGuidingCenter *p_gc,
    Bfield *bfield);

GPU_DECLARE_TARGET_SIMD
void marker_copy_go(MarkerGyroOrbit *p1, size_t i, MarkerGyroOrbit *p2, size_t j);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD
void marker_copy_gc(
    MarkerGuidingCenter *p1, size_t i, MarkerGuidingCenter *p2, size_t j);

DECLARE_TARGET_SIMD
void marker_copy_fl(MarkerFieldLine *p1, size_t i, MarkerFieldLine *p2, size_t j);

#endif
