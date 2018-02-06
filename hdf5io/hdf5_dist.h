
#ifndef HDF5_DIST_H
#define HDF5_DIST_H

#include "../ascot5.h"
#include "../distributions.h"

void hdf5_dist_write_rzvv(dist_rzvv_offload_data* dist, real* hist,
			  char* filename, char* qid);

#endif
