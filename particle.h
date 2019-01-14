/**
 * @file particle.h
 * @brief Header file for particle.c
 *
 * The relationship between the seven different marker structs is:
 *
 *    particle >--+            +--> particle_simd_fo
 *                |            |
 * particle_gc >--particle_state--> particle_simd_gc
 *                |            |
 * particle_ml >--+            +--> particle_simd_ml
 *
 * which is explained in particle.c. This file contains the definitions
 * of these structs as well as input_particle and particle_queue structs. Former
 * is a wrapper for particle, particle_gc, particle_ml, and particle_state
 * structs while latter is a queue from where markers are chosen when simulation
 * begins and updated when simulation ends.
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include "ascot5.h"
#include "B_field.h"
#include "E_field.h"
#include "error.h"

/**
 * @brief General representation of a marker
 *
 * This struct is a self-consistent representation of a marker which can be
 * constructed from any input and from which any simulation marker struct can
 * be constructed. This structure is not intended to be used during the
 * simulation, but before and after for storing marker data, and making it easy
 * to switch between different marker types and construct markers from inputs.
 *
 * Whenever marker properties are changed, it is the responsibility of the
 * function performing the change to make sure all parameters remain consistent.
 * For example, magnetic field must be updated when marker position changes.
 */
typedef struct {
    real r;           /**< Guiding center R coordinate [m]                 */
    real phi;         /**< Guiding center phi coordinate [rad]             */
    real z;           /**< Guiding center z coordinate [m]                 */
    real vpar;        /**< Parallel velocity [m/s]                         */
    real mu;          /**< Magnetic moment [J/T]                           */
    real theta;       /**< Gyroangle [rad]                                 */
    real rprt;        /**< Particle R coordinate [m]                       */
    real phiprt;      /**< Particle phi coordinate [phi]                   */
    real zprt;        /**< Particle z coordinate [m]                       */
    real rdot;        /**< dr/dt [m/s]                                     */
    real phidot;      /**< dphi/dt [rad/s]                                 */
    real zdot;        /**< dz/dt [m/s]                                     */
    real mass;        /**< Mass [kg]                                       */
    real charge;      /**< Charge [e]                                      */
    real weight;      /**< Marker weight                                   */
    real time;        /**< Marker simulation time [s]                      */
    real cputime;     /**< Marker wall-clock time [s]                      */
    real rho;         /**< Marker rho coordinate                           */
    real pol;         /**< Marker poloidal coordinate [rad]                */
    integer id;       /**< Arbitrary but unique ID for the marker          */
    integer endcond;  /**< Marker end condition                            */
    integer walltile; /**< ID of walltile if marker has hit the wall       */
    real B_r;         /**< Magnetic field R component at (r, phi, z) [T]   */
    real B_phi;       /**< Magnetic field phi component at (r, phi, z) [T] */
    real B_z;         /**< Magnetic field z component at (r, phi, z) [T]   */
    real B_r_dr;      /**< dB_R/dR at (r, phi, z) [T/m]                    */
    real B_phi_dr;    /**< dB_phi/dR at (r, phi, z) [T/m]                  */
    real B_z_dr;      /**< dB_z/dR at (r, phi, z) [T/m]                    */
    real B_r_dphi;    /**< dB_R/dphi at (r, phi, z) [T/m]                  */
    real B_phi_dphi;  /**< dB_phi/dphi at (r, phi, z) [T/m]                */
    real B_z_dphi;    /**< dB_z/dphi at (r, phi, z) [T/m]                  */
    real B_r_dz;      /**< dB_R/dz at (r, phi, z) [T/m]                    */
    real B_phi_dz;    /**< dB_phi/dz at (r, phi, z) [T/m]                  */
    real B_z_dz;      /**< dB_z/dz at (r, phi, z) [T/m]                    */

    a5err err;        /**< error flag */
} particle_state;

/**
 * @brief Particle input
 *
 * When particle marker data is read, this struct is created and filled. This
 * struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                    */
    real phi;    /**< phi coordinate [rad]                */
    real z;      /**< z coordinate [m]                    */
    real v_r;    /**< Velocity R-component [m/s]          */
    real v_phi;  /**< Velocity phi-component [m/s]        */
    real v_z;    /**< Velocity z-component [m/s]          */
    real mass;   /**< Mass [kg]                           */
    real charge; /**< Charge [e]                          */
    real weight; /**< Particle marker weight              */
    real time;   /**< Particle marker simulation time [s] */
    integer id;  /**< Unique ID for the particle marker   */
} particle;

/**
 * @brief Guiding center input
 *
 * When guiding center marker data is read, this struct is created and filled.
 * This struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                          */
    real phi;    /**< phi coordinate [rad]                      */
    real z;      /**< z coordinate [m]                          */
    real energy; /**< Kinetic energy [J]                        */
    real pitch;  /**< Pitch                                     */
    real theta;  /**< Gyroangle [rad]                           */
    real mass;   /**< Mass [kg]                                 */
    real charge; /**< Charge [e]                                */
    real weight; /**< Guiding center marker weight              */
    real time;   /**< Guiding center marker simulation time [s] */
    integer id;  /**< Unique ID for the guiding center marker   */
} particle_gc;

/**
 * @brief Field line input
 *
 * When field line marker data is read, this struct is created and filled. This
 * struct is then converted to a particle_state struct.
 */
typedef struct {
    real r;      /**< R coordinate [m]                      */
    real phi;    /**< phi coordinate [rad]                  */
    real z;      /**< z coordinate [m]                      */
    real pitch;  /**< Direction                             */
    real weight; /**< Field line marker weight              */
    real time;   /**< Field line marker simulation time [s] */
    integer id;  /**< Unique ID for the field line marker   */
} particle_ml;

/**
 * @brief Marker queue
 *
 * Each time a marker has finished simulation, a new marker is chosen from this
 * queue and the old marker's data is updated. Markers are never removed from
 * the queue but an index is kept to mark where the next not yet simulated
 * marker is found. Markers are represented by particle_state struct when they
 * are stored in the queue.
 *
 * Note: The queue can and is accessed by several threads, so make sure each
 * access is thread-safe.
 */
typedef struct {
    int n;                 /**< Total number of markers in this queue        */
    particle_state** p;    /**< Pointer to an array storing pointers to all
                                markers within this queue.                   */
    volatile int next;     /**< Index where next unsimulated marker is found */
    volatile int finished; /**< Number of markers who have finished
                                simulation                                   */
} particle_queue;

/**
 * @brief Marker types enum
 *
 * Used to indicate what marker type is stored in input_particle wrapper.
 */
typedef enum input_particle_type {
    input_particle_type_p,  /**< Type corresponding to particle struct       */
    input_particle_type_gc, /**< Type corresponding to particle_gc struct    */
    input_particle_type_ml, /**< Type corresponding to particle_ml struct    */
    input_particle_type_s   /**< Type corresponding to particle_state struct */
} input_particle_type;

/**
 * @brief Wrapper for marker structs
 *
 * This struct wraps particle_state struct and all input structs. Reason is
 * because input data can have several marker types and with this wrapper only
 * a single array is required to store them. The same array can be used when
 * the input markers are turned into marker states.
 *
 * Only a single type is stored here indicated by the "type" field. Storing a
 * a new marker struct removes the old marker struct.
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

/**
 * @brief Struct representing NSIMD particle markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_fo_fixed.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;       /**< Particle R coordinate [m]          */
    real phi[NSIMD] __memalign__;     /**< Particle phi coordinate [phi]      */
    real z[NSIMD] __memalign__;       /**< Particle z coordinate [m]          */
    real rdot[NSIMD] __memalign__;    /**< dr/dt [m/s]                        */
    real phidot[NSIMD] __memalign__;  /**< dphi/dt [rad/s]                    */
    real zdot[NSIMD] __memalign__;    /**< dz/dt [m/s]                        */
    real mass[NSIMD] __memalign__;    /**< Mass [kg]                          */
    real charge[NSIMD] __memalign__;  /**< Charge [e]                         */
    real time[NSIMD] __memalign__;    /**< Marker simulation time [s]         */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real pol[NSIMD] __memalign__;     /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
} particle_simd_fo;

/**
 * @brief Struct representing NSIMD guiding center markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_gc_fixed.c or simulate_gc_adaptive.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;      /**< Guiding center R coordinate [m]     */
    real phi[NSIMD] __memalign__;    /**< Guiding center phi coordinate [phi] */
    real z[NSIMD] __memalign__;      /**< Guiding center z coordinate [m]     */
    real vpar[NSIMD] __memalign__;   /**< Parallel velocity [m/s]             */
    real mu[NSIMD] __memalign__;     /**< Magnetic moment [J/T]               */
    real theta[NSIMD] __memalign__;  /**< Gyroangle [rad]                     */
    real mass[NSIMD] __memalign__;   /**< Mass [kg]                           */
    real charge[NSIMD] __memalign__; /**< Charge [e]                          */
    real time[NSIMD] __memalign__;   /**< Marker simulation time [s]          */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real pol[NSIMD] __memalign__;     /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
} particle_simd_gc;

/**
 * @brief Struct representing NSIMD field line markers
 *
 * This struct is used in simulation when the simulation loop in
 * simulate_ml_adaptive.c is used.
 *
 * It contains physical and simulation parameters necessary for the simulation
 * If a function makes changes to any of these parameters, it is that function's
 * responsibility to make sure all fields remain consistent, i.e., if position
 * changes then the magnetic field should be updated. Each field is a memory
 * aligned array with length NSIMD, so this struct represents NSIMD markers
 * (they can be dummy or markers whose simulation has been terminated) and so it
 * can be used within SIMD loops.
 *
 * The fields are aligned to 64 bit with __memalign__ (see ascot5.h).
 */
typedef struct {
    /* Physical coordinates and parameters */
    real r[NSIMD] __memalign__;     /**< Field line R coordinate [m]          */
    real phi[NSIMD] __memalign__;   /**< Field line phi coordinate [phi]      */
    real z[NSIMD] __memalign__;     /**< Field line z coordinate [m]          */
    real pitch[NSIMD] __memalign__; /**< Field line direction: along (1) or
                                         against (-1) magnetic field vector   */
    real time[NSIMD] __memalign__;  /**< Field line simulation "time" i.e.
                                         (distance / speed of light) [m]      */

    /* Magnetic field data */
    real B_r[NSIMD] __memalign__;        /**< Magnetic field R component at
                                              marker position [T]             */
    real B_phi[NSIMD] __memalign__;      /**< Magnetic field phi component at
                                              marker position [T]             */
    real B_z[NSIMD] __memalign__;        /**< Magnetic field z component at
                                              marker position [T]             */

    real B_r_dr[NSIMD] __memalign__;     /**< dB_R/dR at marker position [T/m]     */
    real B_phi_dr[NSIMD] __memalign__;   /**< dB_phi/dR at marker position [T/m]   */
    real B_z_dr[NSIMD] __memalign__;     /**< dB_z/dR at marker position [T/m]     */
    real B_r_dphi[NSIMD] __memalign__;   /**< dB_R/dphi at marker position [T/m]   */
    real B_phi_dphi[NSIMD] __memalign__; /**< dB_phi/dphi at marker position [T/m] */
    real B_z_dphi[NSIMD] __memalign__;   /**< dB_z/dphi at marker position [T/m]   */
    real B_r_dz[NSIMD] __memalign__;     /**< dB_R/dz at marker position [T/m]     */
    real B_phi_dz[NSIMD] __memalign__;   /**< dB_phi/dz at marker position [T/m]   */
    real B_z_dz[NSIMD] __memalign__;     /**< dB_z/dz at marker position [T/m]     */

    /* Quantities used in diagnostics */
    real weight[NSIMD] __memalign__;  /**< Marker weight                      */
    real cputime[NSIMD] __memalign__; /**< Marker wall-clock time [s]         */
    real rho[NSIMD] __memalign__;     /**< Marker rho coordinate              */
    real pol[NSIMD] __memalign__;     /**< Marker poloidal coordinate [rad]   */

    integer id[NSIMD] __memalign__;       /**< Unique ID for the marker       */
    integer endcond[NSIMD] __memalign__;  /**< Marker end condition           */
    integer walltile[NSIMD] __memalign__; /**< ID of walltile if marker has
                                               hit the wall                   */

    /* Meta data */
    integer running[NSIMD] __memalign__; /**< Indicates whether this marker is
                                              currently simulated (1) or not  */
    a5err err[NSIMD] __memalign__;       /**< Error flag, zero if no error    */
    integer index[NSIMD] __memalign__;   /**< Marker index at marker queue    */
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

void particle_input_to_state(input_particle* p, particle_state* ps,
                             B_field_data* Bdata);

#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo,
                           int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p,
                          B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc,
                           int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p,
                          B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml,
                           int j, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p,
                          B_field_data* Bdata);
#pragma omp declare simd uniform(p_fo,Bdata)
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc,
                      B_field_data* Bdata);
#pragma omp end declare target

#endif
