#ifndef HDF5_DIAG_H
#define HDF5_DIAG_H

#include "../ascot5.h"
#include "../simulate.h"

void hdf5_diag_write(sim_offload_data* sim, real* diag_offload_array, char* out, char* qid);

#endif
