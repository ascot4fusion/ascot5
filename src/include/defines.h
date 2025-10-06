/**
 * @file ascot5.h
 * Main header file for ASCOT5.
 *
 * This header file defines general types and structs used widely by the program
 * as well as any features dependent on compiler defined macros.
 */
#ifndef ASCOT5_H
#define ASCOT5_H

#include <omp.h>
#include <time.h>

/** Stringify arguments */
#define stringify(...) #__VA_ARGS__

/** Wrapper to stringify for converting compiler macros to strings */
#define str_macro(c) stringify(c)

/** Macro to suppress unused variable warnings. */
#ifndef GPU
#define SUPPRESS_UNUSED_WARNING(x) (void)(x)
#else
#define SUPPRESS_UNUSED_WARNING(x)
#endif

/** Stringify argument and convert it to a directive */
#define str_pragma(c) _Pragma(stringify(c))

/**
 * This is used to tell the compiler that we want a variable aligned to
 * 64 bits for Xeon Phi; may not be always necessary */
#define __memalign__ __attribute__((aligned(64)))

/** We use a custom type real to represent floating point numbers; precision
 * can be defined compile-time. */
#if defined SINGLEPRECISION
typedef int integer; /**< Single precision integer */
typedef float real;  /**< Singe precision float    */
#else
typedef long integer; /**< Double precision integer */
typedef double real;  /**< Double precision float   */
#endif

/** Number of particles simulated simultaneously in a particle group
 *         operations */
#ifndef NSIMD
#define NSIMD 16
#endif

/** Maximum number of plasma species */
#define MAX_SPECIES 8

/** Maximum number of MHD modes */
#define MHD_MODES_MAX_NUM 512

/** Maximum number of atomic reactions */
#define MAX_ATOMIC 32

/** Maximum distance BBNBI traces markers in seconds */
#define NBI_MAX_DISTANCE 1e-3

/** Maximum number of NBI injectors */
#define NBI_MAX_INJ 32

/** Maximum number of Wiener processes stored (effectively number
 *  of time step reductions) */
#define WIENERSLOTS 20

/** Determine whether to use geometric method of Box-Muller to
 *  to generate normal random numbers */
#define A5_CCOL_USE_GEOBM 1

/** If adaptive time step falls below this value, produce an error */
#define A5_EXTREMELY_SMALL_TIMESTEP 1e-12

/** How often progress is being written (s) in the stdout file */
#define A5_PRINTPROGRESSINTERVAL 20

/** Wall time */
#define A5_WTIME omp_get_wtime()

/** Choose whether magnetic field spline interpolation is done
 *  with an explicit (1) or compact (0) way */
#define INTERP_SPL_EXPL 0

/** Choose whether to use tabulated values for collision coefficients */
#define A5_CCOL_USE_TABULATED 0

/** Default depth of octree struct */
#define WALL_OCTREE_DEPTH 7

/** Input filename where RFOF parameters are stored */
#define RFOF_CODEPARAM_XML "rfof_codeparam.xml"

/**
 * @brief Parameters and data required to evaluate Coulomb collisions
 */
typedef struct {
    int usetabulated;   /**< Use tabulated values for special functions    */
    int include_energy; /**< Let collisions change energy                  */
    int include_pitch;  /**< Let collisions change pitch                   */
    int include_gcdiff; /**< Let collisions change guiding center position */
} mccc_data;

#endif
