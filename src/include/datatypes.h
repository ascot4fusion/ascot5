/**
 * @file datatypes.h
 * Simulation data.
 */
#ifndef DATATYPES_H
#define DATATYPES_H

#include "atomic_data.h"
#include "bfield_data.h"
#include "boozer_data.h"
#include "diag_data.h"
#include "efield_data.h"
#include "mhd_data.h"
#include "nbi_data.h"
#include "neutral_data.h"
#include "options.h"
#include "plasma_data.h"
#include "rfof_data.h"
#include "wall_data.h"
#include <stddef.h>
#include <stdint.h>

/**
 * Marker end condition bit masks.
 *
 * These bit masks are used to mark specific end condition as being active.
 */
typedef enum ENDCOND
{
    ENDCOND_TLIM = (1u << 0),   /**< Fixed time limit or mileage reached.     */
    ENDCOND_EMIN = (1u << 1),   /**< Marker energy below minimum value.       */
    ENDCOND_THERM = (1u << 2),  /**< Marker energy below local thermal value. */
    ENDCOND_WALL = (1u << 3),   /**< Marker has hit wall.                     */
    ENDCOND_RHOMAX = (1u << 4), /**< Minimum rho surface crossed.             */
    ENDCOND_RHOMIN = (1u << 5), /**< Maximum rho surface crossed.             */
    ENDCOND_POLMAX = (1u << 6), /**< Poloidal orbit limit reached.            */
    ENDCOND_TORMAX = (1u << 7), /**< Toroidal orbit limit reached.            */
    ENDCOND_CPUMAX = (1u << 8), /**< Wall time exceeded.                      */
    ENDCOND_HYBRID = (1u << 9), /**< Finished GC in hybrid mode.              */
    ENDCOND_NEUTR = (1u << 10), /**< Marker has neutralized.                  */
    ENDCOND_IONIZ = (1u << 11), /**< Marker has ionized.                      */
} ENDCOND;

/**
 * Simulation end condition type.
 *
 * The error flag is a 32 bit integer where each bit represents a different
 * end condition. Multiple end conditions can activate simultaneously.
 */
typedef unsigned int endcond_t;

/**
 * Enum type for indicating type of error.
 *
 * Assign unique value for each type just in case. Do not use zero! Please use
 * running numbering and put the latest entry last.
 */
typedef enum error_type
{
    ERR_INPUT_EVALUATION = 1,  /**< Failure when evaluating input data      */
    ERR_UNKNOWN_INPUT = 2,     /**< Input data type not regonizable         */
    ERR_INPUT_UNPHYSICAL = 3,  /**< Input evaluation result is unphysical   */
    ERR_MARKER_UNPHYSICAL = 4, /**< Some of marker fields are unphysical    */
    ERR_INVALID_TIMESTEP = 5,  /**< Time step is zero, NaN or too small     */
    ERR_WIENER_ARRAY = 6,      /**< Wiener array is full or inconsistent    */
    ERR_INTEGRATION = 7,       /**< Integrating marker coordinates yield
                                    unphysical results                      */
    ERR_ATOMIC_EVALUATION = 8  /**< Failure when evaluating atomic reaction */
} error_type;

/**
 * Available reactions.
 */
typedef enum Reaction
{
    DT_He4n = 1,
    DHe3_He4p = 2,
    DD_Tp = 3,
    DD_He3n = 4,
} Reaction;

/**
 * General simulation data.
 *
 * This structure holds references to all data required to simulate markers
 * except the markers themselves. Options, input data, and diagnostics are all
 * accessable from here.
 *
 * For interfaces, the struct is stored directly as the interfaces themselves
 * store the pointer to the data.
 */
typedef struct
{
    Mhd mhd;                  /**< MHD data interface.                        */
    Wall wall;                /**< Wall data interface.                       */
    Rfof rfof;                /**< Rfof data read via RFOF library.           */
    Bfield bfield;            /**< Magnetic field data interface.             */
    Efield efield;            /**< Electric field data interface.             */
    Plasma plasma;            /**< Plasma data interface.                     */
    Boozer *boozer;           /**< Boozer data.                               */
    Atomic *atomic;           /**< Atomic data.                               */
    Neutral neutral;          /**< Neutral data interface.                    */
    Options *options;         /**< Simulation options.                        */
    Diagnostics *diagnostics; /**< Diagnostics data interface.                */
    void *random_data;        /**< Random number generator.                   */
    void *mccc_data;          /**< Tabulated special functions and collision
                                   operator parameters.                       */
} Simulation;

/**
 * Data for calculating fusion source.
 */
typedef struct
{
    int type1;          /**< Distribution type (1:beam, 2:thermal).           */
    int type2;          /**< Distribution type (1:beam, 2:thermal).           */
    int thermal1;       /**< Thermal species index for reactant 1.            */
    int thermal2;       /**< Thermal species index for reactant 2.            */
    histogram *beam1;   /**< Distribution data for reactant 1.                */
    histogram *beam2;   /**< Distribution data for reactant 2.                */
    double *r;          /**< Radial coordinate at the grid center [m].        */
    double *phi;        /**< Toroidal coordinate at the grid center [rad].    */
    double *z;          /**< Axial coordinate at the grid center [m].         */
    double *vol;        /**< Grid cell volume [m^3].                          */
    size_t volshape[3]; /**< Dimensions of r, phi, z, and volume.             */
    Reaction reaction;  /**< The fusion reaction that is modelled.            */
    float mult;         /**< Multiplication factor which is 0.5 if species is
                             interacting with itself, 1.0 otherwise.          */
} FusionSource;

/**
 * General representation of a marker.
 *
 * This structure holds all data that is used to represent a marker. Both
 * particles and guiding centers are represented by this structure.
 */
typedef struct
{
    real r;            /**< Guiding center R coordinate [m].                  */
    real z;            /**< Guiding center z coordinate [m].                  */
    real phi;          /**< Guiding center phi coordinate [rad].              */
    real rprt;         /**< Particle R coordinate [m].                        */
    real zprt;         /**< Particle z coordinate [m].                        */
    real phiprt;       /**< Particle phi coordinate [phi].                    */
    real zeta;         /**< Guiding center gyroangle [rad].                   */
    real ekin;         /**< Guiding center kinetic energy [J].                */
    real pitch;        /**< Guiding center pitch [1].                         */
    real pr;           /**< Particle momentum r component [kg m/s].           */
    real pz;           /**< Particle momentum z component [kg m/s].           */
    real pphi;         /**< Particle momentum phi component [kg m/s].         */
    real mass;         /**< Mass [kg].                                        */
    real charge;       /**< Charge [C].                                       */
    int anum;          /**< Atomic mass number of marker species.             */
    int znum;          /**< Charge number of marker species.                  */
    real time;         /**< Marker simulation time [s].                       */
    real theta;        /**< Marker poloidal coordinate [rad].                 */
    real weight;       /**< Marker weight.                                    */
    real mileage;      /**< Duration this marker has been simulated [s].      */
    real cputime;      /**< Marker wall-clock time [s].                       */
    size_t id;         /**< Arbitrary but unique ID for the marker.           */
    endcond_t endcond; /**< Marker end condition.                             */
    size_t walltile;   /**< ID of walltile if marker has hit the wall.        */
    err_t err;         /**< Error flag.                                       */
} State;

#endif
