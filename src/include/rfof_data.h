/**
 * @file rfof_data.h
 * RFOF simulation data.
 */
#ifndef RFOF_DATA_H
#define RFOF_DATA_H

#include "defines.h"

/**
 * Reusable struct for storing marker specific data during the simulation loop.
 *
 * The data in this struct is altered during the simulation. Only pointers are
 * stored as the actual data is stored on the Fortran side.
 */
typedef struct rfof_marker
{
    void *p[NSIMD]; /**< The marker struct in a format required by librfof */
    void *history_array[NSIMD]; /**< Stores values of the resonance function for
                                     estimating the next time-step */
    void *diag_array[NSIMD]; /**< C equivalents of Fortran diagnostics pointers
                                  which are required but unused at the moment */
    int nrow[NSIMD]; /**< Number of rows in an resonance history matrix    */
    int ncol[NSIMD]; /**< Number of columns in an resonance history matrix */

} rfof_marker;

/** RFOF simulation input data.
 *
 * Immutable input data shared between all markers. The actual data is stored
 * in the Fortran side and this struct only stores the pointers.
 */
typedef struct
{
    void *rfof_input_params; /**< Pointer to rfof_input_param struct on
                                  the fortran side.                 */
    void *rfglobal;          /**< Wave field; same for all markers. */
} Rfof;

#endif
