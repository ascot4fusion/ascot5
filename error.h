/**
 * @file error.h
 * @brief Error module for ASCOT5. 
 *
 * Contains error codes, their interpretation, 
 * and functions to raise errors.
 */
#ifndef ERROR_H
#define ERROR_H

#include <math.h>

typedef unsigned long int a5err;

static const unsigned int ERRMOD_REJECTED = 1;/** @brief Error arised when input marker was first transformed into a state in particle.c. */
static const unsigned int ERRMOD_STATE    = 2;/** @brief Error arised when transforming between simd and state structures in particle.c. */
static const unsigned int ERRMOD_ORBSTEP  = 5;/** @brief Error originates from step_gc_cashkarp.c, step_gc_rk4.c or step_fo_vpa.c. */
static const unsigned int ERRMOD_COLLSTEP = 6;/** @brief Error originates from mccc.c or associated modules. */

static const unsigned int ERR_UNPHYSICAL_B     = 1;/** @brief Magnetic field is zero or NaN */
static const unsigned int ERR_OUTSIDE_BFIELD   = 2;/** @brief Tried to evaluate magnetic field at outside the grid */
static const unsigned int ERR_UNPHYSICAL_PSI   = 3;/** @brief Psi is zero or NaN */
static const unsigned int ERR_OUTSIDE_PSIFIELD = 4;/** @brief Tried to evaluate psi at outside the grid */
static const unsigned int ERR_UNPHYSICAL_E     = 5;/** @brief Electric field is zero or NaN */
static const unsigned int ERR_OUTSIDE_EFIELD   = 6;/** @brief Tried to evaluate electric field at outside the grid */
static const unsigned int ERR_UNPHYSICAL_TE    = 7;/** @brief Electron temperature is non-positive or NaN */
static const unsigned int ERR_UNPHYSICAL_TI    = 8;/** @brief Ion temperature is non-positive or NaN */
static const unsigned int ERR_UNPHYSICAL_TN    = 9;/** @brief Neutral temperature is negative or NaN */
static const unsigned int ERR_UNPHYSICAL_NE    = 10;/** @brief Electron density is non-positive or NaN */
static const unsigned int ERR_UNPHYSICAL_NI    = 11;/** @brief Ion density is non-positive or NaN */
static const unsigned int ERR_UNPHYSICAL_NN    = 12;/** @brief Neutral density is negative or NaN */
static const unsigned int ERR_OUTSIDE_PLASMA   = 13;/** @brief Tried to evaluate plasma profiles at outside the grid */
static const unsigned int ERR_UNPHYSICAL_FO    = 14;/** @brief R <= 0 or vtot > c or at least one fo parameter is NaN */ 
static const unsigned int ERR_UNPHYSICAL_GC    = 15;/** @brief R <= 0 or vpar > c or mu < 0 or at least one gc parameter is NaN */ 
static const unsigned int ERR_UNPHYSICAL_ML    = 16;/** @brief R <= 0 or pitch == 0 or at least one ml parameter is NaN */ 
static const unsigned int ERR_EXTREMELY_SMALL_TIMESTEP = 17;/** @brief Time step smaller than the A5_EXTREMELY_SMALL_TIMESTEP defined in ascot5.h */
static const unsigned int ERR_TIMESTEP_ZERO_OR_NAN = 18;/** @brief Time step is zero or NaN. */
static const unsigned int ERR_CCOEF_EVAL_FAIL  = 19;/** @brief Collision coefficient evaluation failed */
static const unsigned int ERR_CCOL_EVAL_FAIL   = 20;/** @brief Collision evaluation failed */
static const unsigned int ERR_WIENARR_EXCEEDED = 21;/** @brief Exceeded wiener array capacity */
static const unsigned int ERR_WIENARR_NOPROC   = 22;/** @brief Cannot clean wiener processes for t < t0 since W(t0) does not exists */

#pragma omp declare target

/**
 * @brief Raise a new error
 * 
 * @param errtype type of error as defined in error.h
 * @param linenumber line where this function was called
 */
#pragma omp declare simd
static inline a5err error_raise(unsigned int errtype, int linenumber) {
    a5err err = (a5err) errtype;
    err += (a5err)(linenumber*256);
    return err;
}

/**
 * @brief Define from which module error appears
 * 
 * @param a5err a raised error 
 * @param modnumber module identifier as defined in error.h
 */
#pragma omp declare simd
static inline a5err error_module(a5err err, unsigned int modnumber) {
    return err + (a5err)(modnumber*256*1024);
}

static inline int error_getmod(a5err err) {
    return (int)(floor(err / (256*1024)));
}

static inline int error_getline(a5err err) {
    err -= (a5err)(floor(err / (256*1024))*1024*256);
    return  (int)(floor(err / (256)));
}

static inline int error_getmsg(a5err err) {
    return (int) ( (err % (256*1024)) % 256 );
}

#pragma omp end declare target
#endif
