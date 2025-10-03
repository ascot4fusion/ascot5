/**
 * @file error.h
 * Error module for ASCOT5.
 *
 * This modules handles marker simulation time errors, i.e., errors that arise
 * when a marker is being simulated and which do not lead to termination of the
 * whole program, just to the termination of that marker.
 *
 * Error flag is stored in particle_state struct and care should be taken that
 * no error flag overwriting occurs by ceasing the simulation of that marker.
 * Error flag is represented by a5err which is a 64 bit integer. Zero means
 * no error. When error occurs, an error is raised with error_raise() which
 * uses a5err to store the following data as follows:
 *
 * - 0 - 7 bits store error message
 * - 8 - 17 bits store line number where error was raised
 * - 18 - 49 bits store the file name where error was raised
 * - 50 - 63 bits are left empty.
 *
 * Whenever adding a new error type, include that type to error_type enum. These
 * types should be generic as line and file are enough to locate what exactly
 * went wrong. Then add description of that error type to error_parse2str().
 *
 * Whenever adding a new file where error can arise, add that file to error_file
 * enum and write the file name to error_parse2str().
 */
#ifndef ERROR_H
#define ERROR_H

#include "offload.h"
#include <math.h>
#include <stdio.h>

/**
 * Simulation error flag
 */
typedef unsigned long int a5err;

/**
 * Enum type for indicating which file error originated from.
 *
 * Assign unique value for each type just in case. Do not use zero! Please use
 * running numbering and put the latest entry last.
 */
typedef enum error_file
{
    EF_MCCC_WIENER = 1,      /**< Error is from mccc_wiener.c              */
    EF_MCCC_PUSH = 2,        /**< Error is from mccc_push.c                */
    EF_MCCC_COEFS = 3,       /**< Error is from mccc_coefs.c               */
    EF_MCCC = 4,             /**< Error is from mccc.c                     */
    EF_STEP_FO_VPA = 5,      /**< Error is from step_fo_vpa.c              */
    EF_STEP_GC_CASHKARP = 6, /**< Error is from step_gc_cashkarp.c         */
    EF_STEP_GC_RK4 = 7,      /**< Error is from step_gc_rk4.c              */
    EF_N0_1D = 8,            /**< Error is from N0_1D.c                    */
    EF_N0_3D = 9,            /**< Error is from N0_3D.c                    */
    EF_N0_ST = 10,           /**< Error is from N0_ST.c                    */
    EF_B_3DS = 11,           /**< Error is from B_3DS.c                    */
    EF_B_2DS = 12,           /**< Error is from B_2DS.c                    */
    EF_B_STS = 13,           /**< Error is from B_STS.c                    */
    EF_B_GS = 14,            /**< Error is from B_GS.c                     */
    EF_PLASMA_1D = 15,       /**< Error is from plasma_1DS.c               */
    EF_PLASMA_1DS = 16,      /**< Error is from plasma_1DS.c               */
    EF_PLASMA = 17,          /**< Error is from plasma.c                   */
    EF_E_FIELD = 18,         /**< Error is from E_field.c                  */
    EF_NEUTRAL = 19,         /**< Error is from neutral.c                  */
    EF_E_1DS = 20,           /**< Error is from E_1DS.c                    */
    EF_B_FIELD = 21,         /**< Error is from B_field.c                  */
    EF_PARTICLE = 22,        /**< Error is from particle.c                 */
    EF_BOOZER = 23,          /**< Error is from boozer.c                   */
    EF_MHD = 24,             /**< Error is from mhd.c                      */
    EF_ATOMIC = 25,          /**< Error is from atomic.c                   */
    EF_ASIGMA = 26,          /**< Error is from asigma.c                   */
    EF_ASIGMA_LOC = 27,      /**< Error is from asigma_loc.c               */
    EF_SUZUKI = 28           /**< Error is from suzuki.c                   */
} error_file;

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
 * Raise a new error
 *
 * This is a SIMD function.
 *
 * @param type Type of error.
 * @param line Line number this function is called from
 *        Use macro __LINE__ here.
 * @param file File this function is called from.
 *
 * @return Error containing info on error type and line and file error happened
 */
DECLARE_TARGET_SIMD
static inline a5err error_raise(error_type type, int line, error_file file)
{
    a5err err = 0;
    err += (a5err)(type);
    err += (a5err)(line * 256);
    err += (a5err)(file * 256 * 1024);
    return err;
}

/**
 * Retrieve stored data from the error flag.
 *
 * The values of type and file corresponds to the ones that were in error_type
 * and error_file enums that were used when error was thrown with error_raise().
 *
 * This function is host only.
 *
 * @param err The error flag.
 * @param msg Pointer for storing type of the error.
 * @param line Pointer for storing the line where error originated from.
 * @param file Pointer for storing the file where error originated from.
 */
void error_parse(a5err err, int *msg, int *line, int *file);

/**
 * Convert error flag in string format.
 *
 * This function is host only.
 *
 * @param err The error flag.
 * @param msg Pointer for error message.
 * @param line Pointer for storing the line where error originated from.
 * @param file Pointer for storing the file where error originated from.
 */
void error_parse2str(a5err err, char *msg, char *line, char *file);

#endif
