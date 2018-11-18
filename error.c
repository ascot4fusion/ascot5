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

        case ERR_INTERPOLATION_FAILED:
            sprintf(msg, "Interpolation failed");
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

        default:
            sprintf(file, "unknown file");
            break;
    }
}
