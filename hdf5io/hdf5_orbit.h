/**
 * @file hdf5_orbits.h
 * @brief Header file for hdf5_orbits.c
 */
#ifndef HDF5_ORBITS_H
#define HDF5_ORBITS_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../diag/diag_orb.h"

int hdf5_orbits_write(hid_t f, char* qid, diag_orb_offload_data* diag,
                      real* orbits);

#endif
