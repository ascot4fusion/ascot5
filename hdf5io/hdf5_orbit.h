/**
 * @file hdf5_orbit.h
 * @brief Header file for hdf5_orbit.c
 */
#ifndef HDF5_ORBIT_H
#define HDF5_ORBIT_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../diag/diag_orb.h"

int hdf5_orbit_write(hid_t f, char* qid, diag_orb_offload_data* diag,
                     real* orbits);

#endif
