/**
 * @file hdf5_orbit.c
 * @brief Module for writing orbits diagnostics data to a HDF5 file
 */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_helpers.h"
#include "../simulate.h"
#include "../particle.h"
#include "../ascot5.h"
#include "../diag/diag_orb.h"
#include "../consts.h"
#include "hdf5_orbit.h"

void hdf5_orbit_writeset(hid_t group,  const char* name, const char* unit,
                         int type, int arraylength, real confac,
                         integer* mask, integer size, real* data);

/**
 * @brief Write orbit diagnostics data to a HDF5 file
 *
 * @param f hdf5 file
 * @param path path to group which is created here and where the data is stored
 * @param data orbit diagnostics offload data
 * @param orbits array storing the orbit data (format same as in init() in
 *        diag_orb.c)
 *
 * @return zero on success
 */
int hdf5_orbit_write(hid_t f, char* path, diag_orb_offload_data* data,
                     real* orbits) {
    hid_t group = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(group < 0) {
        return 1;
    }

    int arraylength = data->Nmrk*data->Npnt;
    int datasize = 0;
    integer* mask = malloc(arraylength*sizeof(real));

    for(integer i=0; i < arraylength; i++) {
        mask[i] = orbits[i] > 0;
        if(mask[i]) {
            datasize++;
        }
    }

    int dtypef64 = 0;
    int dtypei32 = 1;
    int dtypei16 = 2;
    if(data->record_mode == simulate_mode_fo) {

        hdf5_orbit_writeset(group,  "ids", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*0]);
        hdf5_orbit_writeset(group,  "mileage", "s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*1]);
        hdf5_orbit_writeset(group,  "r", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*2]);
        hdf5_orbit_writeset(group,  "phi", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*3]);
        hdf5_orbit_writeset(group,  "z", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*4]);
        hdf5_orbit_writeset(group,  "pr", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*5]);
        hdf5_orbit_writeset(group,  "pphi", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*6]);
        hdf5_orbit_writeset(group,  "pz", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*7]);
        hdf5_orbit_writeset(group,  "weight", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*8]);
        hdf5_orbit_writeset(group,  "charge", "e", dtypei16, arraylength,
                            1.0/CONST_E,
                            mask, datasize, &orbits[arraylength*9]);
        hdf5_orbit_writeset(group,  "rho", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*10]);
        hdf5_orbit_writeset(group,  "theta", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*11]);
        hdf5_orbit_writeset(group,  "br", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*12]);
        hdf5_orbit_writeset(group,  "bphi", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*13]);
        hdf5_orbit_writeset(group,  "bz", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*14]);
        hdf5_orbit_writeset(group,  "simmode", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*15]);

    }

    if(data->record_mode == simulate_mode_gc) {
        hdf5_orbit_writeset(group,  "ids", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*0]);
        hdf5_orbit_writeset(group,  "mileage", "s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*1]);
        hdf5_orbit_writeset(group,  "r", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*2]);
        hdf5_orbit_writeset(group,  "phi", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*3]);
        hdf5_orbit_writeset(group,  "z", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*4]);
        hdf5_orbit_writeset(group,  "ppar", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*5]);
        hdf5_orbit_writeset(group,  "mu", "eV/T", dtypef64, arraylength,
                            1.0/CONST_E,
                            mask, datasize, &orbits[arraylength*6]);
        hdf5_orbit_writeset(group,  "zeta", "rad", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*7]);
        hdf5_orbit_writeset(group,  "weight", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*8]);
        hdf5_orbit_writeset(group,  "charge", "e", dtypei16, arraylength,
                            1.0/CONST_E,
                            mask, datasize, &orbits[arraylength*9]);
        hdf5_orbit_writeset(group,  "rho", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*10]);
        hdf5_orbit_writeset(group,  "theta", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*11]);
        hdf5_orbit_writeset(group,  "br", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*12]);
        hdf5_orbit_writeset(group,  "bphi", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*13]);
        hdf5_orbit_writeset(group,  "bz", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*14]);
        hdf5_orbit_writeset(group,  "simmode", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*15]);
    }

    if(data->record_mode == simulate_mode_ml) {
        hdf5_orbit_writeset(group,  "ids", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*0]);
        hdf5_orbit_writeset(group,  "mileage", "s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*1]);
        hdf5_orbit_writeset(group,  "r", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*2]);
        hdf5_orbit_writeset(group,  "phi", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*3]);
        hdf5_orbit_writeset(group,  "z", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*4]);
        hdf5_orbit_writeset(group,  "rho", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*5]);
        hdf5_orbit_writeset(group,  "theta", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*6]);
        hdf5_orbit_writeset(group,  "br", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*7]);
        hdf5_orbit_writeset(group,  "bphi", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*8]);
        hdf5_orbit_writeset(group,  "bz", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*9]);
        hdf5_orbit_writeset(group,  "simmode", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*10]);
    }

    if(data->record_mode == simulate_mode_hybrid) {
        hdf5_orbit_writeset(group,  "ids", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*0]);
        hdf5_orbit_writeset(group,  "mileage", "s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*1]);
        hdf5_orbit_writeset(group,  "r", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*2]);
        hdf5_orbit_writeset(group,  "phi", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*3]);
        hdf5_orbit_writeset(group,  "z", "m", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*4]);
        hdf5_orbit_writeset(group,  "pr", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*5]);
        hdf5_orbit_writeset(group,  "pphi", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*6]);
        hdf5_orbit_writeset(group,  "pz", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*7]);
        hdf5_orbit_writeset(group,  "ppar", "kg*m/s", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*8]);
        hdf5_orbit_writeset(group,  "mu", "eV/T", dtypef64, arraylength,
                            1.0/CONST_E,
                            mask, datasize, &orbits[arraylength*9]);
        hdf5_orbit_writeset(group,  "zeta", "rad", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*10]);
        hdf5_orbit_writeset(group,  "weight", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*11]);
        hdf5_orbit_writeset(group,  "charge", "e", dtypei16, arraylength,
                            1.0/CONST_E,
                            mask, datasize, &orbits[arraylength*12]);
        hdf5_orbit_writeset(group,  "rho", "1", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*13]);
        hdf5_orbit_writeset(group,  "theta", "deg", dtypef64, arraylength,
                            180.0/CONST_PI,
                            mask, datasize, &orbits[arraylength*14]);
        hdf5_orbit_writeset(group,  "br", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*15]);
        hdf5_orbit_writeset(group,  "bphi", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*16]);
        hdf5_orbit_writeset(group,  "bz", "T", dtypef64, arraylength, 1,
                            mask, datasize, &orbits[arraylength*17]);
        hdf5_orbit_writeset(group,  "simmode", "1", dtypei32, arraylength, 1,
                            mask, datasize, &orbits[arraylength*18]);

    }

    if(data->mode == DIAG_ORB_POINCARE) {
        hdf5_orbit_writeset(group,  "pncrid", "1", dtypei16, arraylength, 1,
                            mask, datasize, &orbits[arraylength * data->Nfld]);
        hdf5_orbit_writeset(group,  "pncrdi", "1", dtypei16, arraylength, 1,
                            mask, datasize, &orbits[arraylength * (data->Nfld+1)]);
    }

    free(mask);
    H5Gclose (group);

    return 0;
}

/**
 * @brief Helper function for writing orbit diagnostic datasets
 *
 * @param group HDF5 group where dataset is written
 * @param name name the data is wrote with
 * @param unit unit the the data is wrote with
 * @param type is data double (0), int (1), or integer (2)
 * @param arraylength length of the orbit data array
 * @param confac unit conversion factor orbits is multiplied before writing
 * @param mask flag array indicating which elements in orbits array contain data
 * @param size number of data elements
 * @param orbits orbit data array
 */
void hdf5_orbit_writeset(hid_t group,  const char* name, const char* unit,
                         int type, int arraylength, real confac,
                         integer* mask, integer size, real* orbits) {
    if(type == 0) {
        real* data = (real*) malloc(size * sizeof(real));
        integer j = 0;
        for(integer i = 0; i < arraylength; i++) {
            if(mask[i]) {
                data[j] = confac*orbits[i];
                j++;
            }
        }
        hdf5_write_extendible_dataset_double(group, name, size, data);
        H5LTset_attribute_string(group, name, "unit", unit);
        free(data);
    }
    else if(type == 1) {
        int* data = (int*) malloc(size * sizeof(int));
        integer j = 0;
        for(integer i = 0; i < arraylength; i++) {
            if(mask[i]) {
                data[j] = (int)(confac*orbits[i]);
                j++;
            }
        }
        hdf5_write_extendible_dataset_int(group, name, size, data);
        H5LTset_attribute_string(group, name, "unit", unit);
        free(data);
    }
    else if(type == 2) {
        integer* data = (integer*) malloc(size * sizeof(integer));
        integer j = 0;
        for(integer i = 0; i < arraylength; i++) {
            if(mask[i]) {
                data[j] = (integer)(confac*orbits[i]);
                j++;
            }
        }
        hdf5_write_extendible_dataset_long(group, name, size, data);
        H5LTset_attribute_string(group, name, "unit", unit);
        free(data);
    }
}
