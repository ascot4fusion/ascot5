/**
 * @mainpage ASCOT5
 *
 * @section intro Introduction
 *
 * ASCOT5 is a test-particle Monte Carlo orbit-following code for solving the
 * Fokker-Planck equation in toroidal plasmas.
 *
 * @section compilation Compilation
 *
 * Currently only Intel's compiler version
 * 14.0 and up support all the OpenMP features used, but gcc can be used to
 * compile some test programs and the full code with limited features.
 *
 * Before compilation the necessary compiler modules need to be loaded.
 * The command on CSC Bull and Helios is "module load intel/<version>".
 * Available versions can be checked with "module avail".
 *
 * Due to the use of compiler flags to define code functionality, the code
 * should always be fully recompiled by calling <make clean> first if the
 * parameters are changed.
 *
 * Syntax for compilation:
 *
 * make clean && make <program_name> <parameters>
 *
 * Available parameters:
 *
 *  - NSIMD=n   number of particles in a group (default 16); these are
 *              processed simultaneously by each thread
 *  - CC=...    compiler (default icc)
 *  - TARGET=1  Offload computation to Xeon Phi accelerator
 *  - VERBOSE=n print increasing amounts of progress information:
 *                0: No information except bare essentials
 *                1: Standard information; everything happening outside
 *                   simulation loops is printed
 *                2: Extensive information; a record of simulation progress
 *                   is written process-specific *.stdout files
 *  - MPI=1     enable MPI
 *
 *  Available programs:
 *
 *   - ascot5_main : The stand-alone ASCOT5 main program
 *
 *  Example:
 *
 *  make clean && make ascot5_main CC=icc NSIMD=16
 */

/**
 * @file ascot5.h
 * @brief Main header file for ASCOT5
 *
 * This header file defines general types and structs used widely by the program
 * as well as any features dependent on compiler defined macros.
 */
#ifndef ASCOT5_H
#define ASCOT5_H

#include <omp.h>
#include <time.h>

/** This is used to tell the compiler that we want a variable aligned to
 * 64 bits for Xeon Phi; may not be always necessary */
#define __memalign__ __attribute__((aligned(64)))

/** We use a custom type real to represent floating point numbers; precision
 * can be defined compile-time. */
#if defined SINGLEPRECISION
typedef int integer;
typedef float real;
#else
typedef long integer;
typedef double real;
#endif

/** @brief Number of particles simulated simultaneously in a particle group
 *         operations */
#ifndef NSIMD
#define NSIMD 16
#endif

/** @brief Maximum number of plasma species */
#define MAX_SPECIES 8

/** @brief Maximum number of Wiener processes stored (effectively number
 *  of time step reductions) */
#define WIENERSLOTS 20

/** @brief Determine whether to use geometric method of Box-Muller to
 *  to generate normal random numbers */
#define A5_CCOL_USE_GEOBM 1

/** @brief Turn off energy component of Coulomb collisions */
#define A5_CCOL_NOENERGY 0

/** @brief Turn off pitch component of Coulomb collisions */
#define A5_CCOL_NOPITCH  0

/** @brief Turn off spatial component in gc picture (classical diffusion)
    of Coulomb collisions */
#define A5_CCOL_NOGCDIFF 0

/** @brief If adaptive time step falls below this value, produce an error */
#define A5_EXTREMELY_SMALL_TIMESTEP 1e-12

/** @brief How often progress is being written (s) in the stdout file */
#define A5_PRINTPROGRESSINTERVAL 60

/** @brief Wall time */
#define A5_WTIME omp_get_wtime()

/** @brief Choose whether magnetic field spline interpolation is done
 *  with an explicit (1) or compact (0) way */
#define INTERP_SPL_EXPL 0

/** @brief Choose whether to use tabulated values for collision coefficients */
#define A5_CCOL_USE_TABULATED 0

#endif
