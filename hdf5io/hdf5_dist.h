#ifndef HDF5_DIST_H
#define HDF5_DIST_H

#include "../ascot5.h"
#include "../diag/dist_5D.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_rho5D.h"
#include "../diag/dist_rho6D.h"

void hdf5_dist_write_5D(dist_5D_offload_data* dist, real* hist, char* filename,
                        char* qid);
void hdf5_dist_write_6D(dist_6D_offload_data* dist, real* hist, char* filename,
                        char* qid);
void hdf5_dist_write_rho5D(dist_rho5D_offload_data* dist, real* hist, char* filename,
                        char* qid);
void hdf5_dist_write_rho6D(dist_rho6D_offload_data* dist, real* hist, char* filename,
                        char* qid);
                        
#endif
