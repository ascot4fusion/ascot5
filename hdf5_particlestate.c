/**
 * @file hdf5_particlestate.c
 * @brief HDF5 particle state IO
 */
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "ascot5.h"
#include "hdf5_helpers.h"
#include "particle.h"

/**
   @brief Writes the particle state to an ASCOT4 HDF5 file.
*/
int hdf5_particlestate_write(hid_t file, char *state, int n, particle* p) {
    hid_t state_group = hdf5_create_group(file, state);
    if(state_group < 0) {
        return 0;
    }

    hsize_t dims[1];
    dims[0] = n;
    real* data = (real*) malloc(n * sizeof(real));

    int i;
    for(i = 0; i < n; i++) {
        data[i] = p[i].r;
    }
    H5LTmake_dataset(state_group, "Rprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "Rprt", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].phi;
    }
    H5LTmake_dataset(state_group, "phiprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "phiprt", "unit", "deg");

    for(i = 0; i < n; i++) {
        data[i] = p[i].z;
    }
    H5LTmake_dataset(state_group, "zprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "zprt", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].rdot;
    }
    H5LTmake_dataset(state_group, "vR", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vR", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].phidot;
    }
    H5LTmake_dataset(state_group, "vphi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vphi", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].zdot;
    }
    H5LTmake_dataset(state_group, "vz", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vz", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].weight;
    }
    H5LTmake_dataset(state_group, "weight", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "weight", "unit", "1/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].time;
    }
    H5LTmake_dataset(state_group, "time", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "time", "unit", "s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].id;
    }
    H5LTmake_dataset(state_group, "id", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "id", "unit", "");

    for(i = 0; i < n; i++) {
        data[i] = p[i].endcond;
    }
    H5LTmake_dataset(state_group, "endCond", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "endCond", "unit", "");

    for(i = 0; i < n; i++) {
        if(p[i].running)
            data[i] = -1;
        else
            data[i] = 1;
    }
    H5LTmake_dataset(state_group, "wallTile", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "wallTile", "unit", "");

    for(i = 0; i < n; i++) {
        data[i] = -1;
    }
    H5LTmake_dataset(state_group, "R", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "R", "unit", "m");

    H5Gclose(state_group);
}
