/**
 * @brief hdf5_dist.h
 * @brief Header file for hdf5_dist.c
 */
#ifndef HDF5_DIST_H
#define HDF5_DIST_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../diag/dist_5D.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_rho5D.h"
#include "../diag/dist_rho6D.h"
#include "../diag/dist_com.h"

int hdf5_dist_write_5D(hid_t f, char* path, dist_5D_data* dist);
int hdf5_dist_write_6D(hid_t f, char* path, dist_6D_data* dist);
int hdf5_dist_write_rho5D(hid_t f, char* path, dist_rho5D_data* dist);
int hdf5_dist_write_rho6D(hid_t f, char* path, dist_rho6D_data* dist);
int hdf5_dist_write_COM(hid_t f, char* path, dist_COM_data* dist);

#endif
