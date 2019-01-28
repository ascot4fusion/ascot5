/**
 * @file hdf5_boozer.c
 * @brief Module for reading Boozer data from HDF5 file
 *
 */
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../boozer.h"
#include "hdf5_helpers.h"
#include "hdf5_boozer.h"

/**
 * @brief Initialize Boozer offload data from HDF5 file
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is initialized here
 * @param offload_array pointer to offload array which is allocated and
 *                      initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 *
 * @todo Konsta will write this.
 */
int hdf5_boozer_init_offload(hid_t f, boozer_offload_data* offload_data,
                             real** offload_array, char* qid) {
    return 0;
}
