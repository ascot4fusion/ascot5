/**
 * @file hdf5_helpers.h
 * @brief Header file for hdf5_helpers.h
 */
#ifndef HDF5_HELPERS_H
#define HDF5_HELPERS_H

#include <hdf5.h>
#include "../ascot5.h"

void hdf5_init(void);
hid_t hdf5_create(const char* filename);
hid_t hdf5_open(const char* filename);
hid_t hdf5_open_ro(const char* filename);
herr_t hdf5_close(hid_t file_id);
hid_t hdf5_create_group(hid_t loc, const char* path);
herr_t hdf5_find_group(hid_t loc, const char* path);
char* hdf5_generate_qid_path(const char* original, char* qid, char* path);
char* hdf5_gen_path(const char* original, char* qid, char* path);
int hdf5_read_double(const char* var, real* ptr, hid_t file, char* qid,
                     const char* errfile, int errline);
int hdf5_read_short(const char* var, short* ptr, hid_t file, char* qid,
                    const char* errfile, int errline);
int hdf5_read_int(const char* var, int* ptr, hid_t file, char* qid,
                  const char* errfile, int errline);
int hdf5_read_long(const char* var, long* ptr, hid_t file, char* qid,
                   const char* errfile, int errline);
herr_t  hdf5_write_string_attribute(hid_t loc, const char* path,
                                    const char* attrname,  const char* string);
herr_t hdf5_write_extendible_dataset_double(hid_t group,
                                            const char* datasetname,
                                            int length, double* data);
herr_t hdf5_write_extendible_dataset_long(hid_t group,
                                          const char* datasetname,
                                          int length, long* data);
herr_t hdf5_write_extendible_dataset_int(hid_t group,
                                         const char* datasetname,
                                         int length, int* data);

#endif
