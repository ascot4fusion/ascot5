/**
 * @file defines.h
 * Common definitions and macros.
 *
 * This header file defines general types and structs used widely by the program
 * as well as any features dependent on compiler defined macros.
 */
#ifndef DEFINES_H
#define DEFINES_H

#include <omp.h>
#include <time.h>

/**
 * Stringify arguments.
 */
#define stringify(...) #__VA_ARGS__

/**
 * Stringify argument and convert it to a directive.
 */
#define str_pragma(c) _Pragma(stringify(c))

/**
 * Align variable to 64 bits.
 *
 * This might not be needed anymore.
 */
#define __memalign__ __attribute__((aligned(64)))

/**
 * Suppress warnings from unused variables arising when the code is not compiled
 * for the GPU.
 *
 * This is mostly used in offload functions as those are only relevant for the
 * GPU.
 */
#ifdef GPU
#define SUPPRESS_UNUSED_WARNING(x)
#else
#define SUPPRESS_UNUSED_WARNING(x) (void)(x)
#endif

/**
 * Real number that is either in double or single precision.
 */
#ifdef USE_SINGLE_PRECISION
typedef float real;
#else
typedef double real;
#endif

/**
 * Number of particles simulated simultaneously in CPU by a single thread when
 * vector instructions are enabled.
 *
 * If the code was compiled in double precision and the CPU supports AVX512
 * instructions, this number should be 8.
 */
#ifndef NSIMD
#define NSIMD 8
#endif

/**
 * Maximum number of plasma species.
 */
#define MAX_SPECIES 16

/**
 * Maximum number of MHD modes.
 */
#define MHD_MODES_MAX_NUM 512

/**
 * Maximum number of atomic reactions.
 */
#define MAX_ATOMIC 32

/**
 * Maximum distance BBNBI traces markers in seconds.
 */
#define NBI_MAX_DISTANCE 1e-3

/**
 * Maximum number of NBI injectors.
 */
#define NBI_MAX_INJ 32

/**
 * Maximum number of Wiener processes stored (effectively number of time step
 * reductions).
 */
#define WIENERSLOTS 20

/**
 * Determine whether to use geometric method of Box-Muller to to generate normal
 * random numbers.
 */
#define A5_CCOL_USE_GEOBM 1

/**
 * If adaptive time step falls below this value, produce an error.
 */
#define A5_EXTREMELY_SMALL_TIMESTEP 1e-12

/**
 * Get current wall time.
 */
#define A5_WTIME omp_get_wtime()

/**
 * Spcifies whether to use tabulated values for collision coefficients.
 */
#define A5_CCOL_USE_TABULATED 0

/**
 * Default depth of octree struct.
 */
#define WALL_OCTREE_DEPTH 7

/**
 * Input filename where RFOF parameters are stored.
 */
#define RFOF_CODEPARAM_XML "rfof_codeparam.xml"

/**
 * Parameters and data required to evaluate Coulomb collisions.
 */
typedef struct
{
    int usetabulated;   /**< Use tabulated values for special functions    */
    int include_energy; /**< Let collisions change energy                  */
    int include_pitch;  /**< Let collisions change pitch                   */
    int include_gcdiff; /**< Let collisions change guiding center position */
} mccc_data;

/**
 * Simulation error flag type.
 *
 * The error flag is a 32 bit integer that contains following information:
 *
 * - 0 - 7 bits store error type.
 * - 8 - 17 bits store the file name where error was raised.
 * - 18 - 31 bits store line number where error was raised.
 */
typedef unsigned int err_t;

/**
 * Generate an error.
 *
 * @param type Type of the error.
 * @param file The current file.
 */
#define ERROR_RAISE(type, file)                                                \
    (((err_t)(type) & 0xFFu) | (((err_t)(file) & 0x3FFu) << 8) |               \
     (((err_t)(__LINE__) & 0x3FFFu) << 18))

/**
 * Generate an error if condition is met and no previous error exists.
 *
 * @param err Existing error flag.
 *        A new error is raised only if this is zero.
 * @param condition Raise error if this condition is true.
 * @param type Type of the error.
 * @param file The current file.
 */
#define ERROR_CHECK(err, condition, type, file_id)                             \
    ((err) ? (err) : ((condition) ? ERROR_RAISE(type, file_id) : 0))

/**
 * Error types.
 */
typedef enum
{
    ERR_IMAGINARY_RHO = 1,
    ERR_INTERPOLATED_OUTSIDE_RANGE,
    ERR_TOO_MANY_TIME_STEP_REDUCTIONS,
} Error_type;

/**
 * Enum type for indicating which file error originated from.
 *
 * This list is used to generate error codes during the simulation. For
 * convenience all source files are listed here even though not all of them
 * raise errors.
 *
 * To make this list maintainable, the list is sorted alphabetically (based on
 * the full path). To list the source files in order:
 * "find . -type f -name '*.c' | sort".
 */
typedef enum
{
    COREAPI_INTERPOLATE_C,
    COREAPI_SOLVE_C,
    COREAPI_TOOLS_C,
    DATA_ATOMIC_C,
    DATA_BFIELD_C,
    DATA_BFIELD_ANALYTICAL_C,
    DATA_BFIELD_CARTESIAN_C,
    DATA_BFIELD_SPLINE2D_C,
    DATA_BFIELD_SPLINE3D_C,
    DATA_BFIELD_STELLARATOR_C,
    DATA_BOOZER_C,
    DATA_DIAG_C,
    DATA_DIAG_ORB_C,
    DATA_DIAG_TRANSCOEF_C,
    DATA_DIST_5D_C,
    DATA_DIST_6D_C,
    DATA_DIST_COM_C,
    DATA_DIST_RHO5D_C,
    DATA_DIST_RHO6D_C,
    DATA_EFIELD_C,
    DATA_EFIELD_CARTESIAN_C,
    DATA_EFIELD_POTENTIAL1D_C,
    DATA_HIST_C,
    DATA_MARKER_C,
    DATA_MHD_C,
    DATA_MHD_DYNAMIC_C,
    DATA_MHD_STATIONARY_C,
    DATA_NBI_C,
    DATA_NEUTRAL_C,
    DATA_NEUTRAL_ARBITRARY_C,
    DATA_NEUTRAL_RADIAL_C,
    DATA_PLASMA_C,
    DATA_PLASMA_DYNAMIC1D_C,
    DATA_PLASMA_LINEAR1D_C,
    DATA_RFOF_C,
    DATA_WALL_C,
    DATA_WALL_CONTOUR2D_C,
    DATA_WALL_TRIANGULAR3D_C,
    ENDCOND_C,
    SIMULATE_ATOMIC_REACTIONS_C,
    SIMULATE_CCOLL_GC_EULER_C,
    SIMULATE_CCOLL_GC_MILSTEIN_C,
    SIMULATE_CCOLL_GO_EULER_C,
    SIMULATE_CCOLL_WIENER_C,
    SIMULATE_COULOMB_COLLISIONS_C,
    SIMULATE_FIELD_LINE_ADAPTIVE_C,
    SIMULATE_FUSION_SOURCE_C,
    SIMULATE_GUIDING_CENTER_ADAPTIVE_C,
    SIMULATE_GUIDING_CENTER_FIXED_C,
    SIMULATE_GYRO_ORBIT_FIXED_C,
    SIMULATE_NBI_SOURCE_C,
    SIMULATE_ORBIT_FL_CASHKARP_C,
    SIMULATE_ORBIT_GC_CASHKARP_C,
    SIMULATE_ORBIT_GC_RK4_C,
    SIMULATE_ORBIT_GO_VPA_C,
    UTILS_BOSCHHAL_CE_C,
    UTILS_GCTRANSFORM_C,
    UTILS_INTERP1DCOMP_C,
    UTILS_INTERP2DCOMP_C,
    UTILS_INTERP3DCOMP_C,
    UTILS_LININT_C,
    UTILS_LIST_C,
    UTILS_MATH_C,
    UTILS_OCTREE_C,
    UTILS_RANDOM_C,
    UTILS_SPLINECOMP_C,
    UTILS_SUZUKI_C,
    C_FILE_COUNT
} Source_file;

/**
 * List of source files.
 *
 * The list is in the same order as the ``Source_file`` enum so this can be used
 * to decode the enum values to strings.
 */
static const char *file_names[C_FILE_COUNT] = {
    "./coreapi/interpolate.c",
    "./coreapi/solve.c",
    "./coreapi/tools.c",
    "./data/atomic.c",
    "./data/bfield.c",
    "./data/bfield_analytical.c",
    "./data/bfield_cartesian.c",
    "./data/bfield_spline2d.c",
    "./data/bfield_spline3d.c",
    "./data/bfield_stellarator.c",
    "./data/boozer.c",
    "./data/diag.c",
    "./data/diag_orb.c",
    "./data/diag_transcoef.c",
    "./data/dist_5D.c",
    "./data/dist_6D.c",
    "./data/dist_com.c",
    "./data/dist_rho5D.c",
    "./data/dist_rho6D.c",
    "./data/efield.c",
    "./data/efield_cartesian.c",
    "./data/efield_potential1d.c",
    "./data/hist.c",
    "./data/marker.c",
    "./data/mhd.c",
    "./data/mhd_dynamic.c",
    "./data/mhd_stationary.c",
    "./data/nbi.c",
    "./data/neutral.c",
    "./data/neutral_arbitrary.c",
    "./data/neutral_radial.c",
    "./data/plasma.c",
    "./data/plasma_dynamic1d.c",
    "./data/plasma_linear1d.c",
    "./data/rfof.c",
    "./data/wall.c",
    "./data/wall_contour2d.c",
    "./data/wall_triangular3d.c",
    "./endcond.c",
    "./simulate/atomic_reactions.c",
    "./simulate/ccoll_gc_euler.c",
    "./simulate/ccoll_gc_milstein.c",
    "./simulate/ccoll_go_euler.c",
    "./simulate/ccoll_wiener.c",
    "./simulate/coulomb_collisions.c",
    "./simulate/field_line_adaptive.c",
    "./simulate/fusion_source.c",
    "./simulate/guiding_center_adaptive.c",
    "./simulate/guiding_center_fixed.c",
    "./simulate/gyro_orbit_fixed.c",
    "./simulate/nbi_source.c",
    "./simulate/orbit_fl_cashkarp.c",
    "./simulate/orbit_gc_cashkarp.c",
    "./simulate/orbit_gc_rk4.c",
    "./simulate/orbit_go_vpa.c",
    "./utils/boschhale.c",
    "./utils/gctransform.c",
    "./utils/interp1Dcomp.c",
    "./utils/interp2Dcomp.c",
    "./utils/interp3Dcomp.c",
    "./utils/linint.c",
    "./utils/list.c",
    "./utils/math.c",
    "./utils/octree.c",
    "./utils/random.c",
    "./utils/splinecomp.c",
    "./utils/suzuki.c"};

#endif
