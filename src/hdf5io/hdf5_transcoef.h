/**
 * @file hdf5_transcoef.h
 * @brief Header file for hdf5_transcoef.c
 */
#ifndef HDF5_TRANSCOEF_H
#define HDF5_TRANSCOEF_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../diag/diag_transcoef.h"

int hdf5_transcoef_write(hid_t f, char* run,
                         diag_transcoef_offload_data* diag,
                         real* coefarr);

#endif
