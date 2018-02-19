#ifndef HDF5_DIST_H
#define HDF5_DIST_H

#include "../ascot5.h"
#include "../dist_5D.h"
#include "../dist_6D.h"

void hdf5_dist_write_5D(dist_5D_offload_data* dist, real* hist, char* filename,
                        char* qid);
void hdf5_dist_write_6D(dist_6D_offload_data* dist, real* hist, char* filename,
                        char* qid);

#endif
