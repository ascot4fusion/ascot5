/**
 * @file hdf5_bfield.h
 * @brief Header file for hdf5_bfield.c
 */
#ifndef HDF5_BFIELD_H
#define HDF5_BFIELD_H
#include "ascot5.h"
#include "B_2D.h"
#include "B_3D.h"
#include "B_GS.h"
#include "B_ST.h"

void hdf5_bfield_init_offload_ST(B_ST_offload_data* offload_data, real** offload_array);

#endif
