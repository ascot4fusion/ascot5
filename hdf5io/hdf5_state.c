/**
 * @file hdf5_state.c
 * @brief Module for writing marker state to HDF5 file
 */
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "../particle.h"
#include "hdf5_helpers.h"
#include "hdf5_state.h"

/**
 * @brief Writes marker state to an ASCOT5 HDF5 file.
 *
 * Markers are written in /results/run_XXXXXXXXXX/state/ group where X's are
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
int hdf5_state_write(hid_t f, char* qid, char* state, integer n,
                     particle_state* p) {

    char path[256];
    hdf5_gen_path("/results/run_XXXXXXXXXX/", qid, path);
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
    hdf5_write_extendible_dataset_double(state_group, "rprt", n, data);
    H5LTset_attribute_string(state_group, "rprt", "unit", "m");

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
        real pnorm =  sqrt(p[i].p_r*p[i].p_r + p[i].p_phi*p[i].p_phi + p[i].p_z*p[i].p_z);
        real gamma = physlib_gamma_pnorm(p[i].mass, pnorm);
        data[i] = p[i].p_r / (gamma*p[i].mass);
    }
    hdf5_write_extendible_dataset_double(state_group, "vr", n, data);
    H5LTset_attribute_string(state_group, "vr", "unit", "m/s");

    for(i = 0; i < n; i++) {
        real pnorm =  sqrt(p[i].p_r*p[i].p_r + p[i].p_phi*p[i].p_phi + p[i].p_z*p[i].p_z);
        real gamma = physlib_gamma_pnorm(p[i].mass, pnorm);
        data[i] = p[i].p_phi / (gamma*p[i].mass);
    }
    hdf5_write_extendible_dataset_double(state_group, "vphi", n, data);
    H5LTset_attribute_string(state_group, "vphi", "unit", "m/s");

    for(i = 0; i < n; i++) {
        real pnorm = sqrt(p[i].p_r*p[i].p_r + p[i].p_phi*p[i].p_phi + p[i].p_z*p[i].p_z);
        real gamma = physlib_gamma_pnorm(p[i].mass, pnorm);
        data[i] = p[i].p_z / (gamma*p[i].mass);
    }
    hdf5_write_extendible_dataset_double(state_group, "vz", n, data);
    H5LTset_attribute_string(state_group, "vz", "unit", "m/s");

    /* Guiding center coordinates */
    for(i = 0; i < n; i++) {
        data[i] = p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "r", n, data);
    H5LTset_attribute_string(state_group, "r", "unit", "m");

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
        real Bnorm = sqrt(
            p[i].B_r * p[i].B_r + p[i].B_phi * p[i].B_phi + p[i].B_z * p[i].B_z);
        real gamma = physlib_gamma_ppar(p[i].mass, p[i].mu, p[i].ppar, Bnorm);
        data[i] = p[i].ppar  / (gamma*p[i].mass);
    }
    hdf5_write_extendible_dataset_double(state_group, "vpar", n, data);
    H5LTset_attribute_string(state_group, "vpar", "unit", "m/s");

    for(i = 0; i < n; i++) {
        data[i] = p[i].mu/CONST_E;
    }
    hdf5_write_extendible_dataset_double(state_group, "mu", n, data);
    H5LTset_attribute_string(state_group, "mu", "unit", "eV/T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].zeta;
    }
    hdf5_write_extendible_dataset_double(state_group, "zeta", n, data);
    H5LTset_attribute_string(state_group, "zeta", "unit", "rad");

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
        data[i] = p[i].theta*(180/CONST_PI);
    }
    hdf5_write_extendible_dataset_double(state_group, "theta", n, data);
    H5LTset_attribute_string(state_group, "theta", "unit", "deg");

    for(i = 0; i < n; i++) {
        data[i] = p[i].mass/CONST_U;
    }
    hdf5_write_extendible_dataset_double(state_group, "mass", n, data);
    H5LTset_attribute_string(state_group, "mass", "unit", "amu");

    /* Magnetic field */
    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r;
    }
    hdf5_write_extendible_dataset_double(state_group, "br", n, data);
    H5LTset_attribute_string(state_group, "br", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi;
    }
    hdf5_write_extendible_dataset_double(state_group, "bphi", n, data);
    H5LTset_attribute_string(state_group, "bphi", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z;
    }
    hdf5_write_extendible_dataset_double(state_group, "bz", n, data);
    H5LTset_attribute_string(state_group, "bz", "unit", "T");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "brdr", n, data);
    H5LTset_attribute_string(state_group, "brdr", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "brdphi", n, data);
    H5LTset_attribute_string(state_group, "brdphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_r_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "brdz", n, data);
    H5LTset_attribute_string(state_group, "brdz", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "bphidr", n, data);
    H5LTset_attribute_string(state_group, "bphidr", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "bphidphi", n, data);
    H5LTset_attribute_string(state_group, "bphidphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_phi_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "bphidz", n, data);
    H5LTset_attribute_string(state_group, "bphidz", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dr;
    }
    hdf5_write_extendible_dataset_double(state_group, "bzdr", n, data);
    H5LTset_attribute_string(state_group, "bzdr", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dphi/p[i].r;
    }
    hdf5_write_extendible_dataset_double(state_group, "bzdphi", n, data);
    H5LTset_attribute_string(state_group, "bzdphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
        data[i] = p[i].B_z_dz;
    }
    hdf5_write_extendible_dataset_double(state_group, "bzdz", n, data);
    H5LTset_attribute_string(state_group, "bzdz", "unit", "T/m");

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
        intdata32[i] = (int)round(p[i].charge/CONST_E);
    }
    hdf5_write_extendible_dataset_int(state_group, "charge", n, intdata32);
    H5LTset_attribute_string(state_group, "charge", "unit", "e");

    for(i = 0; i < n; i++) {
        intdata32[i] = p[i].anum;
    }
    hdf5_write_extendible_dataset_int(state_group, "anum", n, intdata32);
    H5LTset_attribute_string(state_group, "anum", "unit", "1");

    for(i = 0; i < n; i++) {
        intdata32[i] = p[i].znum;
    }
    hdf5_write_extendible_dataset_int(state_group, "znum", n, intdata32);
    H5LTset_attribute_string(state_group, "znum", "unit", "1");

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
