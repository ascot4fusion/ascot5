/**
 * @file hdf5_particlestate.c
 * @brief Module for writing marker state to HDF5 file
 */
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../particle.h"
#include "hdf5_helpers.h"

/**
 * @brief Writes marker state to an ASCOT5 HDF5 file.
 *
 * Markers are written in /results/run-XXXXXXXXXX/state/ group where X's are
 * the QID.
 *
 * All fields in particle_state struct are written. Fields are written as they
 * are except some unit conversions. Each field is written as n length array
 * (the order is same between all fields) and each dataset includes string
 * attribute that indicates the unit of the stored quantity.
 *
 * @param f output HDF5 file id
 * @param qid QID of the run group
 * @param state name of the state
 * @param n number of markers in state array
 * @param p array holding marker states
 *
 * @return Zero on success
*/
int hdf5_particlestate_write(hid_t f, char* qid, char* state, integer n,
                             particle_state* p) {

    char path[256];
    hdf5_gen_path("/results/run-XXXXXXXXXX/", qid, path);
    strcat(path, state);

    hid_t state_group = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT);
    if(state_group < 0) {
        return 1;
    }

    real* data = (real*) malloc(n * sizeof(real));

    /* Particle coordinates */
    integer i;
    for(i = 0; i < n; i++) {
        data[i] = p[i].rprt;
    }
    hdf5_write_extendible_dataset_double(state_group, "Rprt", n, data);
    H5LTset_attribute_string(state_group, "Rprt", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].phiprt*(180/CONST_PI);
    }
    hdf5_write_extendible_dataset_double(state_group, "phiprt", n, data);
    H5LTset_attribute_string(state_group, "phiprt", "unit", "deg");

    for(i = 0; i < n; i++) {
        data[i] = p[i].zprt;
    }
    hdf5_write_extendible_dataset_double(state_group, "zprt", n, data);
    H5LTset_attribute_string(state_group, "zprt", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].rdot;
    }
    hdf5_write_extendible_dataset_double(state_group, "vR", n, data);
    H5LTset_attribute_string(state_group, "vR", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].phidot * p[i].rprt;
    }
    hdf5_write_extendible_dataset_double(state_group, "vphi", n, data);
    H5LTset_attribute_string(state_group, "vphi", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].zdot;
    }
    hdf5_write_extendible_dataset_double(state_group, "vz", n, data);
    H5LTset_attribute_string(state_group, "vz", "unit", "m/s");

    /* Guiding center coordinates */
    for(i = 0; i < n; i++) {
        data[i] = p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "R", n, data);
    H5LTset_attribute_string(state_group, "R", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].phi*(180/CONST_PI);
    }
    hdf5_write_extendible_dataset_double(state_group, "phi", n, data);
    H5LTset_attribute_string(state_group, "phi", "unit", "deg");

    for(i = 0; i < n; i++) {
        data[i] = p[i].z;
    }
    hdf5_write_extendible_dataset_double(state_group, "z", n, data);
    H5LTset_attribute_string(state_group, "z", "unit", "m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].vpar;
    }
    hdf5_write_extendible_dataset_double(state_group, "vpar", n, data);
    H5LTset_attribute_string(state_group, "vpar", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].mu/CONST_E;
    }
    hdf5_write_extendible_dataset_double(state_group, "mu", n, data);
    H5LTset_attribute_string(state_group, "mu", "unit", "eV/T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].theta;
    }
    hdf5_write_extendible_dataset_double(state_group, "theta", n, data);
    H5LTset_attribute_string(state_group, "theta", "unit", "rad");

    /* Common */
    for(i = 0; i < n; i++) {
        data[i] = p[i].weight;
    }
    hdf5_write_extendible_dataset_double(state_group, "weight", n, data);
    H5LTset_attribute_string(state_group, "weight", "unit", "markers/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].time;
    }
    hdf5_write_extendible_dataset_double(state_group, "time", n, data);
    H5LTset_attribute_string(state_group, "time", "unit", "s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].cputime;
    }
    hdf5_write_extendible_dataset_double(state_group, "cputime", n, data);
    H5LTset_attribute_string(state_group, "cputime", "unit", "s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].rho;
    }
    hdf5_write_extendible_dataset_double(state_group, "rho", n, data);
    H5LTset_attribute_string(state_group, "rho", "unit", "1");

    for(i = 0; i < n; i++) {
        data[i] = p[i].pol*(180/CONST_PI);
    }
    hdf5_write_extendible_dataset_double(state_group, "pol", n, data);
    H5LTset_attribute_string(state_group, "pol", "unit", "deg");

    for(i = 0; i < n; i++) {
        data[i] = p[i].mass/CONST_U;
    }
    hdf5_write_extendible_dataset_double(state_group, "mass", n, data);
    H5LTset_attribute_string(state_group, "mass", "unit", "amu");

    /* Magnetic field */
    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_R", n, data);
    H5LTset_attribute_string(state_group, "B_R", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_phi", n, data);
    H5LTset_attribute_string(state_group, "B_phi", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_z", n, data);
    H5LTset_attribute_string(state_group, "B_z", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_R_dR", n, data);
    H5LTset_attribute_string(state_group, "B_R_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_R_dphi", n, data);
    H5LTset_attribute_string(state_group, "B_R_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_R_dz", n, data);
    H5LTset_attribute_string(state_group, "B_R_dz", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_phi_dR", n, data);
    H5LTset_attribute_string(state_group, "B_phi_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_phi_dphi", n, data);
    H5LTset_attribute_string(state_group, "B_phi_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_phi_dz", n, data);
    H5LTset_attribute_string(state_group, "B_phi_dz", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_z_dR", n, data);
    H5LTset_attribute_string(state_group, "B_z_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_z_dphi", n, data);
    H5LTset_attribute_string(state_group, "B_z_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "B_z_dz", n, data);
    H5LTset_attribute_string(state_group, "B_z_dz", "unit", "T/m");

    free(data);

    /* Integer quantities */
    integer* intdata = (integer*) malloc(n * sizeof(integer));

    for(i = 0; i < n; i++) {
        intdata[i] = p[i].id;
    }
    hdf5_write_extendible_dataset_long(state_group, "id", n, intdata);
    H5LTset_attribute_string(state_group, "id", "unit", "1");

    for(i = 0; i < n; i++) {
        intdata[i] = p[i].endcond;
    }
    hdf5_write_extendible_dataset_long(state_group, "endcond", n, intdata);
    H5LTset_attribute_string(state_group, "endcond", "unit", "1");

    for(i = 0; i < n; i++) {
        intdata[i] = p[i].walltile;
    }
    hdf5_write_extendible_dataset_long(state_group, "walltile", n, intdata);
    H5LTset_attribute_string(state_group, "walltile", "unit", "1");

    free(intdata);

    int* intdata32 = (int*) malloc(n * sizeof(int));
    for(i = 0; i < n; i++) {
        intdata32[i] = (int)(p[i].charge/CONST_E);
    }
    hdf5_write_extendible_dataset_int(state_group, "charge", n, intdata32);
    H5LTset_attribute_string(state_group, "charge", "unit", "e");

    /* Error data */
    int d1, d2;// Dummies
    for(i = 0; i < n; i++) {
        error_parse(p[i].err, &(intdata32[i]), &d1, &d2);
    }
    hdf5_write_extendible_dataset_int(state_group, "errormsg", n, intdata32);
    H5LTset_attribute_string(state_group, "errormsg", "unit", "1");

    for(i = 0; i < n; i++) {
        error_parse(p[i].err, &d1, &(intdata32[i]), &d2);
    }
    hdf5_write_extendible_dataset_int(state_group, "errorline", n, intdata32);
    H5LTset_attribute_string(state_group, "errorline", "unit", "1");

    for(i = 0; i < n; i++) {
        error_parse(p[i].err, &d1, &d2, &(intdata32[i]));
    }
    hdf5_write_extendible_dataset_int(state_group, "errormod", n, intdata32);
    H5LTset_attribute_string(state_group, "errormod", "unit", "1");

    free(intdata32);

    H5Gclose(state_group);

    return 0;
}
