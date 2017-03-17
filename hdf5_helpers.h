/**
 * @file hdf5_helpers.h
 * @brief Header file for hdf5_helpers.h
 */
#ifndef HDF5_HELPERS_H
#define HDF5_HELPERS_H

#include <hdf5.h>

hid_t hdf5_create(const char* filename);
hid_t hdf5_open(const char* filename);
herr_t hdf5_close(hid_t file_id);
hid_t hdf5_create_group(hid_t loc, const char* path);

#endif
