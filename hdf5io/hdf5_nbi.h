#ifndef HDF5_NBI_H
#define HDF5_NBI_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../nbi.h"

int hdf5_nbi_read(hid_t f, int* n_inj, nbi_injector** inj);

#endif
