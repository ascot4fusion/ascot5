/**
 * @file hdf5_particlestate.c
 * @brief HDF5 particle state IO
 */
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../consts.h"
#include "hdf5_helpers.h"
#include "../particle.h"

/**
   @brief Writes the particle state to an ASCOT5 HDF5 file.
*/
int hdf5_particlestate_write(char* fn, char *state, int n, input_particle* p) {
    hid_t file = hdf5_open(fn);

    hid_t state_group = hdf5_create_group(file, state);
    if(state_group < 0) {
        return 0;
    }

    hsize_t dims[1];
    dims[0] = n;

    int i;
    for(i = 0; i < n; i++) {
	if(p[i].type != input_particle_type_s) {
	    return -1;
	}
    }

    real* data = (real*) malloc(n * sizeof(real));

    /* Particle coordinates */
    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.rprt;
    }
    H5LTmake_dataset(state_group, "Rprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "Rprt", "unit", "m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.phiprt*(180/CONST_PI);
    }
    H5LTmake_dataset(state_group, "phiprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "phiprt", "unit", "deg");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.zprt;
    }
    H5LTmake_dataset(state_group, "zprt", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "zprt", "unit", "m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.rdot;
    }
    H5LTmake_dataset(state_group, "vR", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vR", "unit", "m/s");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.phidot * p[i].p_s.rprt;
    }
    H5LTmake_dataset(state_group, "vphi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vphi", "unit", "m/s");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.zdot;
    }
    H5LTmake_dataset(state_group, "vz", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vz", "unit", "m/s");

    /* Guiding center coordinates */
    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.r;
    }
    H5LTmake_dataset(state_group, "R", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "R", "unit", "m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.phi*(180/CONST_PI);
    }
    H5LTmake_dataset(state_group, "phi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "phi", "unit", "deg");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.z;
    }
    H5LTmake_dataset(state_group, "z", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "z", "unit", "m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.vpar;
    }
    H5LTmake_dataset(state_group, "vpar", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "vpar", "unit", "m/s");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.mu/CONST_E;
    }
    H5LTmake_dataset(state_group, "mu", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "mu", "unit", "eV/T");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.theta;
    }
    H5LTmake_dataset(state_group, "theta", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "theta", "unit", "rad");

    /* Common */
    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.weight;
    }
    H5LTmake_dataset(state_group, "weight", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "weight", "unit", "markers/s");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.time;
    }
    H5LTmake_dataset(state_group, "time", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "time", "unit", "s");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.cputime;
    }
    H5LTmake_dataset(state_group, "cputime", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "cputime", "unit", "s");

   for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.rho;
    }
    H5LTmake_dataset(state_group, "rho", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "rho", "unit", "1");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.pol*(180/CONST_PI);
    }
    H5LTmake_dataset(state_group, "pol", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "pol", "unit", "deg");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.mass/CONST_U;
    }
    H5LTmake_dataset(state_group, "mass", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "mass", "unit", "amu");

    /* Magnetic field */
    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_r;
    }
    H5LTmake_dataset(state_group, "B_R", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_R", "unit", "T");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_phi;
    }
    H5LTmake_dataset(state_group, "B_phi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_phi", "unit", "T");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_z;
    }
    H5LTmake_dataset(state_group, "B_z", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_z", "unit", "T");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_r_dr;
    }
    H5LTmake_dataset(state_group, "B_R_dR", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_R_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_r_dphi/p[i].p_s.r;
    }
    H5LTmake_dataset(state_group, "B_R_dphi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_R_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_r_dz;
    }
    H5LTmake_dataset(state_group, "B_R_dz", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_R_dz", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_phi_dr;
    }
    H5LTmake_dataset(state_group, "B_phi_dR", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_phi_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_phi_dphi/p[i].p_s.r;
    }
    H5LTmake_dataset(state_group, "B_phi_dphi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_phi_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_phi_dz;
    }
    H5LTmake_dataset(state_group, "B_phi_dz", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_phi_dz", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_z_dr;
    }
    H5LTmake_dataset(state_group, "B_z_dR", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_z_dR", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_z_dphi/p[i].p_s.r;
    }
    H5LTmake_dataset(state_group, "B_z_dphi", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_z_dphi", "unit", "T/m");

    for(i = 0; i < n; i++) {
	data[i] = p[i].p_s.B_z_dz;
    }
    H5LTmake_dataset(state_group, "B_z_dz", 1, dims, H5T_IEEE_F64LE, data);
    H5LTset_attribute_string(state_group, "B_z_dz", "unit", "T/m");

    free(data);

    /* Integer quantities */
    integer* intdata = (integer*) malloc(n * sizeof(integer));

    for(i = 0; i < n; i++) {
	intdata[i] = p[i].p_s.id;
    }
    H5LTmake_dataset(state_group, "id", 1, dims, H5T_STD_I64LE, intdata);
    H5LTset_attribute_string(state_group, "id", "unit", "1");

    for(i = 0; i < n; i++) {
	intdata[i] = (int)(p[i].p_s.charge/CONST_E);
    }
    H5LTmake_dataset(state_group, "charge", 1, dims, H5T_STD_I32LE, intdata);
    H5LTset_attribute_string(state_group, "charge", "unit", "e");

    for(i = 0; i < n; i++) {
	intdata[i] = p[i].p_s.endcond;
    }
    H5LTmake_dataset(state_group, "endCond", 1, dims, H5T_STD_I64LE, intdata);
    H5LTset_attribute_string(state_group, "endCond", "unit", "1");

    for(i = 0; i < n; i++) {
	intdata[i] = p[i].p_s.walltile;
    }
    H5LTmake_dataset(state_group, "wallTile", 1, dims, H5T_STD_I64LE, intdata);
    H5LTset_attribute_string(state_group, "wallTile", "unit", "1");

    free(intdata);

    H5Gclose(state_group);
    hdf5_close(file);
    return 1;
}
