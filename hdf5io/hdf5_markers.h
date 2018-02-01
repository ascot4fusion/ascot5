#ifndef HDF5_MARKERS_H
#define HDF5_MARKERS_H

#include <hdf5.h>
#include "../particle.h"

void hdf5_markers_init(hid_t f, int *n, input_particle** p);

void hdf5_markers_init_particle(hid_t f, int n, input_particle* p, char* qid);

void hdf5_markers_init_guiding_center(hid_t f, int n, input_particle* p, char* qid);

void hdf5_markers_init_field_line(hid_t f, int n, input_particle* p, char* qid);

#endif
