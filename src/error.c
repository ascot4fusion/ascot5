/**
 * @file error.c
 * @brief Error module for ASCOT5.
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
#include <stdlib.h>
#include <math.h>
#include "error.h"

/**
 * @brief Retrieve stored data from the error flag
 *
 * The values of type and file corresponds to the ones that were in error_type
 * and error_file enums that were used when error was thrown with error_raise().
 *
 * This function is host only.
 *
 * @param err the error flag
 * @param msg pointer for storing type of the error
 * @param line pointer for storing the line where error originated from
 * @param file pointer for storing the file where error originated from
 */
void error_parse(a5err err, int* msg, int* line, int* file) {
    *file = (int)( floor( err / ( 256 * 1024 ) ) );
    err -= *file * ( 256 * 1024 );
    *line = (int)( floor( err / 256 ) );
    err -= *line * ( 256 );
    *msg = err;
}

/**
 * @brief Convert error flag in string format
 *
 * This function is host only.
 *
 * @param err the error flag
 * @param msg pointer for error message
 * @param line pointer for storing the line where error originated from
 * @param file pointer for storing the file where error originated from
 */
void error_parse2str(a5err err, char* msg, char* line, char* file) {
    int msg_i, line_i, file_i;
    error_parse(err, &msg_i, &line_i, &file_i);

    /* Find the error message */
    switch(msg_i) {

        case ERR_INPUT_EVALUATION:
            sprintf(msg, "Input evaluation failed "
                    "(marker could be outside input data grid)");
            break;

        case ERR_UNKNOWN_INPUT:
            sprintf(msg, "Input was not regonized "
                    "(offload or target data could be uninitialized)");
            break;

        case ERR_INPUT_UNPHYSICAL:
            sprintf(msg, "Input evaluation yields unphysical results "
                    "(something could be wrong with the input)");
            break;

        case ERR_MARKER_UNPHYSICAL:
            sprintf(msg, "One or more of marker's fields are unphysical or "
                    "inconsistent (marker input could be corrupted)");
            break;

        case ERR_INVALID_TIMESTEP:
            sprintf(msg, "Time step is zero, NaN, or smaller than MIN_ALLOWED"
                    "_TIME_STEP in ascot5.h (time step limits could be too "
                    "conservative)");
            break;

        case ERR_WIENER_ARRAY:
            sprintf(msg, "Wiener array is full of rejected steps or  "
                    "limits could be too conservative or initial step too "
                    "large)");
            break;

        case ERR_ATOMIC_EVALUATION:
            sprintf(msg, "Atomic reaction evaluation failed "
                    "(atomic reactions for the current charge state might"
                    "not have been implemented)");
            break;

        default:
            sprintf(msg, "Unknown error");
            break;
    }

    /* Line number to string */
    sprintf(line, "%4d", line_i);

    /* Find the file where error originated from */
    switch(file_i) {

        case EF_MCCC_WIENER:
            sprintf(file, "mccc_wiener.c");
            break;

        case EF_MCCC_PUSH:
            sprintf(file, "mccc_push.c");
            break;

        case EF_MCCC_COEFS:
            sprintf(file, "mccc_coefs.c");
            break;

        case EF_MCCC:
            sprintf(file, "mccc.c");
            break;

        case EF_STEP_FO_VPA:
            sprintf(file, "step_fo_vpa.c");
            break;

        case EF_STEP_GC_RK4:
            sprintf(file, "step_gc_rk4.c");
            break;

        case EF_STEP_GC_CASHKARP:
            sprintf(file, "step_gc_cashkarp.c");
            break;

        case EF_PLASMA:
            sprintf(file, "plasma.c");
            break;

        case EF_PLASMA_1D:
            sprintf(file, "plasma_1D.c");
            break;

        case EF_PLASMA_2D:
            sprintf(file, "plasma_2D.c");
            break;

        case EF_PLASMA_1DS:
            sprintf(file, "plasma_1DS.c");
            break;

        case EF_E_FIELD:
            sprintf(file, "E_field.c");
            break;

        case EF_E_1DS:
            sprintf(file, "E_1DS.c");
            break;

        case EF_NEUTRAL:
            sprintf(file, "neutral.c");
            break;

        case EF_N0_1D:
            sprintf(file, "N0_1D.c");
            break;

        case EF_N0_3D:
            sprintf(file, "N0_3D.c");
            break;

        case EF_N0_ST:
            sprintf(file, "N0_ST.c");
            break;

        case EF_B_FIELD:
            sprintf(file, "B_field.c");
            break;

        case EF_B_GS:
            sprintf(file, "B_GS.c");
            break;

        case EF_B_STS:
            sprintf(file, "B_STS.c");
            break;

        case EF_B_2DS:
            sprintf(file, "B_2DS.c");
            break;

        case EF_B_3DS:
            sprintf(file, "B_3DS.c");
            break;

        case EF_PARTICLE:
            sprintf(file, "particle.c");
            break;

        case EF_BOOZER:
            sprintf(file, "boozer.c");
            break;

        case EF_MHD:
            sprintf(file, "mhd.c");
            break;

        case EF_ATOMIC:
            sprintf(file, "atomic.c");
            break;

        case EF_ASIGMA:
            sprintf(file, "asigma.c");
            break;

        case EF_ASIGMA_LOC:
            sprintf(file, "asigma_loc.c");
            break;

        default:
            sprintf(file, "unknown file");
            break;
    }
}
