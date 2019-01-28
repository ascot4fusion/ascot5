/**
 * @file hdf5_mhd.c
 * @brief Module for reading MHD data from HDF5 file
 *
 */
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../mhd.h"
#include "hdf5_helpers.h"
#include "hdf5_mhd.h"

/**
 * @brief Initialize MHD offload data from HDF5 file
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
int hdf5_mhd_init_offload(hid_t f, mhd_offload_data* offload_data,
                          real** offload_array, char* qid) {
    return 0;
}
