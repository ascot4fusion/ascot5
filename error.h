/**
 * @file error.h
 * @brief Error module for ASCOT5.
 *
 * Contains error codes, their interpretation, and a function to raise errors.
 */
#ifndef ERROR_H
#define ERROR_H

#include <math.h>
#include <stdio.h>

/**
 * @brief Simulation error flag
 */
typedef unsigned long int a5err;

/**
 * @brief Enum type for indicating which file error originated from.
 *
 * Assign unique value for each type just in case. Do not use zero! Please use
 * running numbering and put the latest entry last.
 */
typedef enum error_file {
    EF_MCCC_WIENER       =   1, /**< Error is from mccc_wiener.c              */
    EF_MCCC_PUSH         =   2, /**< Error is from mccc_push.c                */
    EF_MCCC_COEFS        =   3, /**< Error is from mccc_coefs.c               */
    EF_MCCC              =   4, /**< Error is from mccc.c                     */
    EF_STEP_FO_VPA       =   5, /**< Error is from step_fo_vpa.c              */
    EF_STEP_GC_CASHKARP  =   6, /**< Error is from step_gc_cashkarp.c         */
    EF_STEP_GC_RK4       =   7, /**< Error is from step_gc_rk4.c              */
    EF_N0_3D             =   8, /**< Error is from N0_3D.c                    */
    EF_N0_ST             =   9, /**< Error is from N0_ST.c                    */
    EF_B_3DS             =  10, /**< Error is from B_3DS.c                    */
    EF_B_2DS             =  11, /**< Error is from B_2DS.c                    */
    EF_B_STS             =  12, /**< Error is from B_STS.c                    */
    EF_B_GS              =  13, /**< Error is from B_GS.c                     */
    EF_PLASMA_1DS        =  14, /**< Error is from plasma_1DS.c               */
    EF_PLASMA            =  15, /**< Error is from plasma.c                   */
    EF_E_FIELD           =  16, /**< Error is from E_field.c                  */
    EF_NEUTRAL           =  17, /**< Error is from neutral.c                  */
    EF_E_1DS             =  18, /**< Error is from E_1DS.c                    */
    EF_B_FIELD           =  19, /**< Error is from B_field.c                  */
    EF_PARTICLE          =  20  /**< Error is from particle.c                 */
}error_file;

/**
 * @brief Enum type for indicating type of error.
 *
 * Assign unique value for each type just in case. Do not use zero! Please use
 * running numbering and put the latest entry last.
 */
typedef enum error_type {
    ERR_INTERPOLATION_FAILED = 1 /**< Failure when interpolating input data   */
}error_type;


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
static const unsigned int ERR_OUTSIDE_N0DATA   = 23;/** @brief Tried to evaluate neutral density outside the grid */
static const unsigned int ERR_OUTSIDE_AXISGRID = 24;/** @brief Tried to evaluate neutral density outside the grid */
static const unsigned int ERR_UNKNOWN_INPUT    = 25;/** @brief Unknown input */
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
    //err -= (a5err)(floor(err / (256*1024))*1024*256);
    return  (int)(floor(err / (256)));
}

static inline int error_getmsg(a5err err) {
    return (int) ( (err % (256*1024)) % 256 );
}

/**
 * @brief Raise a new error
 *
 * This is a SIMD function.
 *
 * @param type type of error
 * @param line line number this function is called from; use macro __LINE__ here
 * @param file file this function is called from
 *
 * @return error containing info on error type and line and file error happened
 */
#pragma omp declare simd
static inline a5err error_throw(error_type type, int line, error_file file) {
    a5err err = 0;
    err += (a5err)(type);
    err += (a5err)(line*256);
    err += (a5err)(file*256*1024);
    return err;
}
#pragma omp end declare target

void error_parse(a5err err, int* msg, int* line, int* file);

void error_parse2str(a5err err, char* msg, char* line, char* file);

#endif
