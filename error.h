/**
 * @file error.h
 * @brief Error module for ASCOT5. 
 *
 * Contains error codes, their interpretation, 
 * and functions to raise errors.
 */
#ifndef ERROR_H
#define ERROR_H

typedef unsigned long int a5err;

static const unsigned int ERRMOD_ORBSTEP  = 1;
static const unsigned int ERRMOD_COLLSTEP = 2;
static const unsigned int ERRMOD_DIAGORB  = 3;
static const unsigned int ERRMOD_DIAGDIST = 4;

static const unsigned int ERR_UNPHYSICAL_B     = 1;
static const unsigned int ERR_OUTSIDE_BFIELD   = 2;
static const unsigned int ERR_UNPHYSICAL_PSI   = 3;
static const unsigned int ERR_OUTSIDE_PSIFIELD = 4;
static const unsigned int ERR_UNPHYSICAL_E     = 5;
static const unsigned int ERR_OUTSIDE_EFIELD   = 6;
static const unsigned int ERR_UNPHYSICAL_TE    = 7;
static const unsigned int ERR_UNPHYSICAL_TI    = 8;
static const unsigned int ERR_UNPHYSICAL_TN    = 9;
static const unsigned int ERR_UNPHYSICAL_NE    = 10;
static const unsigned int ERR_UNPHYSICAL_NI    = 11;
static const unsigned int ERR_UNPHYSICAL_NN    = 12;
static const unsigned int ERR_OUTSIDE_PLASMA   = 13;
static const unsigned int ERR_FO_ZERO_OR_NAN   = 14;
static const unsigned int ERR_GC_ZERO_OR_NAN   = 15;
static const unsigned int ERR_ML_ZERO_OR_NAN   = 16;
static const unsigned int ERR_EXTREMELY_SMALL_TIMESTEP = 17;
static const unsigned int ERR_TIMESTEP_ZERO_OR_NAN = 18;

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

#pragma omp end declare target
#endif
