/**
 * @file hdf5_asigma.c
 * @brief Module for reading atomic reaction cross-section data from HDF5 file
 *
 * Atomic data must be read by calling hdf5_asigma_init_offload()
 * contained in this module. This module contains routines to read atomic
 * reactions data from an ASCOT5 HDF5 file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../asigma.h"
#include "hdf5_helpers.h"
#include "hdf5_asigma.h"

/**
 * @brief Read atomic data from HDF5 file
 *
 * This function reads atomic cross-section (sigma) data with given qid
 * while also initializing offload data and allocating and filling
 * offload array. The file is opened and closed outside this function.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 *
 * @todo
 * - Include AMNS possibility and other alternative asigma types
 */
int hdf5_asigma_init_offload(hid_t f,
                             asigma_offload_data* offload_data,
                             real** offload_array, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    /* Initialize if data was read succesfully */
    if(!err) {
        err = asigma_init_offload(offload_data, offload_array);
    }

    return err;
}
