/**
 * @file marker.h
 * Marker simulation structs and functions to transform between those.
 */
#ifndef MARKER_H
#define MARKER_H

#include "bfield.h"
#include "datatypes.h"
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
 * All parameters are arrays of length NSIMD (CPU) or "number of markers in this
 * process" (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers and for those ``id`` and ``running`` are set to
 * zero.
 */
typedef struct
{
    size_t size;        /**< How many markers this struct contains.           */
    int znum;           /**< Marker species' atomic number.                   */
    int anum;           /**< Marker species' mass number.                     */
    real mass;          /**< Marker species' mass [kg].                       */
    int *charge;        /**< Charge state [e].                                */
    real *r;            /**< Particle R coordinate [m].                       */
    real *phi;          /**< Particle phi coordinate [rad].                   */
    real *z;            /**< Particle z coordinate [m].                       */
    real *p_r;          /**< Momentum r coordinate [kg m/s].                  */
    real *p_phi;        /**< Momentum phi coordinate [kg m/s].                */
    real *p_z;          /**< Momentum z coordinate [kg m/s].                  */
    real *time;         /**< Marker simulation time [s].                      */
    real *br;           /**< BR at marker position [T].                       */
    real *bphi;         /**< Bphi at marker position [T].                     */
    real *bz;           /**< Bz at marker position [T].                       */
    real *dbrdr;        /**< dBR/dR at marker position [T/m].                 */
    real *dbphidr;      /**< dBphi/dR at marker position [T/m].               */
    real *dbzdr;        /**< dBz/dR at marker position [T/m].                 */
    real *dbrdphi;      /**< dBR/dphi at marker position [T/rad].             */
    real *dbphidphi;    /**< dBphi/dphi at marker position [T/rad].           */
    real *dbzdphi;      /**< dBz/dphi at marker position [T/rad].             */
    real *dbrdz;        /**< dBR/dz at marker position [T/m].                 */
    real *dbphidz;      /**< dBphi/dz at marker position [T/m].               */
    real *dbzdz;        /**< dBz/dz at marker position [T/m].                 */
    int *bounces;       /**< Number of times pitch sign changed.              */
    real *weight;       /**< Marker weight.                                   */
    real *cputime;      /**< Marker wall time [s].                            */
    real *rho;          /**< Marker rho coordinate [1].                       */
    real *theta;        /**< Marker cumulative poloidal coordinate [rad].     */
    size_t *id;         /**< Unique ID for the marker.                        */
    endcond_t *endcond; /**< Marker end condition.                            */
    size_t *walltile;   /**< Index (>0) of tile when marker hits the wall.    */
    real *mileage;      /**< Duration the marker has been simulated [s].      */
    int *running;       /**< Is the marker currently simulated.               */
    err_t *err;         /**< Error flag, zero if no error.                    */
    size_t *index;      /**< This marker's index at the marker queue.         */
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
 * All parameters are arrays of length NSIMD (CPU) or "number of markers in this
 * process" (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers and for those ``id`` and ``running`` are set to
 * zero.
 */
typedef struct
{
    size_t size;        /**< How many markers this struct contains.           */
    int znum;           /**< Marker species' atomic number.                   */
    int anum;           /**< Marker species' mass number.                     */
    real mass;          /**< Marker species' mass [kg].                       */
    int *charge;        /**< Charge state [e].                                */
    real *r;            /**< Guiding center R coordinate [m].                 */
    real *phi;          /**< Guiding center phi coordinate [rad].             */
    real *z;            /**< Guiding center z coordinate [m].                 */
    real *ppar;         /**< Parallel momentum [kg m/s].                      */
    real *mu;           /**< Magnetic moment [J/T] .                          */
    real *zeta;         /**< Gyroangle [rad].                                 */
    real *time;         /**< Marker simulation time [s].                      */
    real *br;           /**< BR at marker position [T].                       */
    real *bphi;         /**< Bphi at marker position [T].                     */
    real *bz;           /**< Bz at marker position [T].                       */
    real *dbrdr;        /**< dBR/dR at marker position [T/m].                 */
    real *dbphidr;      /**< dBphi/dR at marker position [T/m].               */
    real *dbzdr;        /**< dBz/dR at marker position [T/m].                 */
    real *dbrdphi;      /**< dBR/dphi at marker position [T/rad].             */
    real *dbphidphi;    /**< dBphi/dphi at marker position [T/rad].           */
    real *dbzdphi;      /**< dBz/dphi at marker position [T/rad].             */
    real *dbrdz;        /**< dBR/dz at marker position [T/m].                 */
    real *dbphidz;      /**< dBphi/dz at marker position [T/m].               */
    real *dbzdz;        /**< dBz/dz at marker position [T/m].                 */
    int *bounces;       /**< Number of times pitch sign changed.              */
    real *weight;       /**< Marker weight.                                   */
    real *cputime;      /**< Marker wall time [s].                            */
    real *rho;          /**< Marker rho coordinate [1].                       */
    real *theta;        /**< Marker cumulative poloidal coordinate [rad].     */
    size_t *id;         /**< Unique ID for the marker.                        */
    endcond_t *endcond; /**< Marker end condition.                            */
    size_t *walltile;   /**< Index (>0) of tile when marker hits the wall.    */
    real *mileage;      /**< Duration the marker has been simulated [s].      */
    int *running;       /**< Is the marker currently simulated.               */
    err_t *err;         /**< Error flag, zero if no error.                    */
    size_t *index;      /**< This marker's index at the marker queue.         */
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
 * All parameters are arrays of length NSIMD (CPU) or "number of markers in this
 * process" (GPU) to faciliate vector operations that are simultaneously
 * performed for all markers in this struct. Dummy markers are used to fill
 * possible missing markers and for those ``id`` and ``running`` are set to
 * zero.
 */
typedef struct
{
    size_t size;        /**< How many markers this struct contains.           */
    real *r;            /**< Field line R coordinate [m].                     */
    real *phi;          /**< Field line phi coordinate [rad].                 */
    real *z;            /**< Field line z coordinate [m].                     */
    int *pitch;         /**< Is the marker traced along or opposite to B-field.*/
    real *time;         /**< The fixed time instant this field line exists [s].*/
    real *br;           /**< BR at marker position [T].                       */
    real *bphi;         /**< Bphi at marker position [T].                     */
    real *bz;           /**< Bz at marker position [T].                       */
    real *dbrdr;        /**< dBR/dR at marker position [T/m].                 */
    real *dbphidr;      /**< dBphi/dR at marker position [T/m].               */
    real *dbzdr;        /**< dBz/dR at marker position [T/m].                 */
    real *dbrdphi;      /**< dBR/dphi at marker position [T/rad].             */
    real *dbphidphi;    /**< dBphi/dphi at marker position [T/rad].           */
    real *dbzdphi;      /**< dBz/dphi at marker position [T/rad].             */
    real *dbrdz;        /**< dBR/dz at marker position [T/m].                 */
    real *dbphidz;      /**< dBphi/dz at marker position [T/m].               */
    real *dbzdz;        /**< dBz/dz at marker position [T/m].                 */
    real *cputime;      /**< Marker wall time [s].                            */
    real *rho;          /**< Marker rho coordinate [1].                       */
    real *theta;        /**< Marker cumulative poloidal coordinate [rad].     */
    real *mileage;      /**< Distance the marker has been traced [m].         */
    size_t *id;         /**< Unique ID for the marker.                        */
    size_t *walltile;   /**< Index (>0) of tile when marker hits the wall.    */
    size_t *index;      /**< This marker's index at the marker queue.         */
    endcond_t *endcond; /**< Marker end condition.                            */
    err_t *err;         /**< Error flag, zero if no error.                    */
    int *running;       /**< Is the marker currently simulated.               */
} MarkerFieldLine;

/**
 * Allocate field line marker simulation vector.
 *
 * @param mrk Allocated marker vector.
 * @param vector_size Number of markers in the vector.
 *
 * @return Zero on success.
 */
int MarkerGyroOrbit_allocate(MarkerGyroOrbit *mrk, size_t vector_size);

/**
 * Deallocate field line marker simulation vector.
 *
 * @param mrk Deallocated marker vector.
 */
void MarkerGyroOrbit_deallocate(MarkerGyroOrbit *mrk);

/**
 * Copy field line marker data from host to GPU.
 *
 * @param mrk Marker vector.
 */
void MarkerGyroOrbit_offload(MarkerGyroOrbit *mrk);

/**
 * Copy field line marker data from GPU to host.
 *
 * @param mrk Marker vector.
 */
void MarkerGyroOrbit_onload(MarkerGyroOrbit *mrk);

GPU_DECLARE_TARGET_SIMD_UNIFORM(copy, original)
/**
 * Copy field line marker from one simulation vector to another.
 *
 * @param copy The simulation vector copy.
 * @param original The simulation vector to copy.
 * @param index Index of the marker in the simulation vector (same in both
 *        structures).
 */
void MarkerGyroOrbit_copy(
    MarkerGyroOrbit *copy, MarkerGyroOrbit *original, size_t index);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(mrk, queue, bfield)
/**
 * Retrieve field line marker from a queue of marker states.
 *
 * The marker position in the queue is stored so that the marker can be returned
 * in the same position. In case the conversion from state to field line marker
 * failed, the error field of the state is updated and this function returns
 * non-zero value.
 *
 * @param mrk Marker simulation vector.
 * @param queue The queue of marker states.
 * @param mrk_index The index of the marker in the simulation vector.
 * @param queue_index The index of the retrieved marker in the queue.
 * @param bfield Magnetic field data.
 *
 * @return Zero on success.
 */
int MarkerGyroOrbit_from_queue(
    MarkerGyroOrbit *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(queue, mrk)
/**
 * Return field line marker back to queue and convert it to a state.
 *
 * The marker is returned on the same position it was retrieved from.
 *
 * @param queue The queue of marker states.
 * @param mrk Marker simulation vector.
 * @param index The index of the marker in the simulation vector.
 * @param bfield Magnetic field data.
 */
void MarkerGyroOrbit_to_queue(
    MarkerQueue *queue, MarkerGyroOrbit *mrk, size_t index, Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

/**
 * Allocate field line marker simulation vector.
 *
 * @param mrk Allocated marker vector.
 * @param vector_size Number of markers in the vector.
 *
 * @return Zero on success.
 */
int MarkerGuidingCenter_allocate(MarkerGuidingCenter *mrk, size_t vector_size);

/**
 * Deallocate field line marker simulation vector.
 *
 * @param mrk Deallocated marker vector.
 */
void MarkerGuidingCenter_deallocate(MarkerGuidingCenter *mrk);

/**
 * Copy field line marker data from host to GPU.
 *
 * @param mrk Marker vector.
 */
void MarkerGuidingCenter_offload(MarkerGuidingCenter *mrk);

/**
 * Copy field line marker data from GPU to host.
 *
 * @param mrk Marker vector.
 */
void MarkerGuidingCenter_onload(MarkerGuidingCenter *mrk);

GPU_DECLARE_TARGET_SIMD_UNIFORM(copy, original)
/**
 * Copy field line marker from one simulation vector to another.
 *
 * @param copy The simulation vector copy.
 * @param original The simulation vector to copy.
 * @param index Index of the marker in the simulation vector (same in both
 *        structures).
 */
void MarkerGuidingCenter_copy(
    MarkerGuidingCenter *copy, MarkerGuidingCenter *original, size_t index);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(mrk, queue, bfield)
/**
 * Retrieve field line marker from a queue of marker states.
 *
 * The marker position in the queue is stored so that the marker can be returned
 * in the same position. In case the conversion from state to field line marker
 * failed, the error field of the state is updated and this function returns
 * non-zero value.
 *
 * @param mrk Marker simulation vector.
 * @param queue The queue of marker states.
 * @param mrk_index The index of the marker in the simulation vector.
 * @param queue_index The index of the retrieved marker in the queue.
 * @param bfield Magnetic field data.
 *
 * @return Zero on success.
 */
int MarkerGuidingCenter_from_queue(
    MarkerGuidingCenter *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(queue, mrk)
/**
 * Return field line marker back to queue and convert it to a state.
 *
 * The marker is returned on the same position it was retrieved from.
 *
 * @param queue The queue of marker states.
 * @param mrk Marker simulation vector.
 * @param index The index of the marker in the simulation vector.
 * @param bfield Magnetic field data.
 */
void MarkerGuidingCenter_to_queue(
    MarkerQueue *queue, MarkerGuidingCenter *mrk, size_t index, Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

/**
 * Allocate field line marker simulation vector.
 *
 * @param mrk Allocated marker vector.
 * @param vector_size Number of markers in the vector.
 *
 * @return Zero on success.
 */
int MarkerFieldLine_allocate(MarkerFieldLine *mrk, size_t vector_size);

/**
 * Deallocate field line marker simulation vector.
 *
 * @param mrk Deallocated marker vector.
 */
void MarkerFieldLine_deallocate(MarkerFieldLine *mrk);

/**
 * Copy field line marker data from host to GPU.
 *
 * @param mrk Marker vector.
 */
void MarkerFieldLine_offload(MarkerFieldLine *mrk);

/**
 * Copy field line marker data from GPU to host.
 *
 * @param mrk Marker vector.
 */
void MarkerFieldLine_onload(MarkerFieldLine *mrk);

GPU_DECLARE_TARGET_SIMD_UNIFORM(copy, original)
/**
 * Copy field line marker from one simulation vector to another.
 *
 * @param copy The simulation vector copy.
 * @param original The simulation vector to copy.
 * @param index Index of the marker in the simulation vector (same in both
 *        structures).
 */
void MarkerFieldLine_copy(
    MarkerFieldLine *copy, MarkerFieldLine *original, size_t index);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(mrk, queue, bfield)
/**
 * Retrieve field line marker from a queue of marker states.
 *
 * The marker position in the queue is stored so that the marker can be returned
 * in the same position. In case the conversion from state to field line marker
 * failed, the error field of the state is updated and this function returns
 * non-zero value.
 *
 * @param mrk Marker simulation vector.
 * @param queue The queue of marker states.
 * @param mrk_index The index of the marker in the simulation vector.
 * @param queue_index The index of the retrieved marker in the queue.
 * @param bfield Magnetic field data.
 *
 * @return Zero on success.
 */
int MarkerFieldLine_from_queue(
    MarkerFieldLine *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(queue, mrk)
/**
 * Return field line marker back to queue and convert it to a state.
 *
 * The marker is returned on the same position it was retrieved from.
 *
 * @param queue The queue of marker states.
 * @param mrk Marker simulation vector.
 * @param index The index of the marker in the simulation vector.
 */
void MarkerFieldLine_to_queue(
    MarkerQueue *queue, MarkerFieldLine *mrk, size_t index);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(mrk_gc, mrk_go, bfield)
/**
 * Convert gyro-orbit marker to guiding-center marker.
 *
 * The conversion is done using proper guiding-center transformation.
 *
 * @param mrk_gc Vector of guiding-center markers where the result is stored.
 * @param mrk_go Vector of gyro-orbit markers used as input.
 * @param index Index of the marker in the vectors.
 * @param bfield Magnetic field data.
 *
 * @return Zero on success.
 */
int MarkerGyroOrbit_to_MarkerGuidingCenter(
    MarkerGuidingCenter *mrk_gc, const MarkerGyroOrbit *mrk_go, size_t index,
    Bfield *bfield);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(q)
/**
 * Iterate through marker array and return index of the finished marker.
 *
 * This function iterates through the markers that are currently being simulated
 * and stops when it finds the first marker that is finished.
 *
 * @param next_in_queue Pointer to the index of the next marker in the queue.
 * @param q Marker queue.
 * @param nmrk Total number of markers in the queue.
 * @param start Starting index.
 * @param ids Array of marker ids.
 * @param running Array of marker running flags.
 * @return Index of the first marker that is finished.
 */
size_t MarkerQueue_cycle(
    size_t *next_in_queue, MarkerQueue *q, size_t nmrk, size_t start,
    size_t ids[nmrk], int running[nmrk]);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif
