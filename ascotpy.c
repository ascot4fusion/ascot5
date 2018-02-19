/**
 * @file ascotpy.c
 * @brief Python interface
 */
#include <stdio.h>
#include <string.h>

#include "ascot5.h"
#include "simulate.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "wall.h"
#include "offload.h"

#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_input.h"
#include "hdf5io/hdf5_bfield.h"
#include "hdf5io/hdf5_efield.h"
#include "hdf5io/hdf5_plasma.h"
#include "hdf5io/hdf5_wall.h"

static sim_offload_data sim_offload;
static sim_data sim;

static real* B_offload_array;
static real* E_offload_array;
static real* plasma_offload_array;
static real* wall_offload_array;


int ascotpy_init_bfield(char* fn) {
    hdf5_init();
    hid_t f = hdf5_open(fn);
    int err = hdf5_bfield_init_offload(f, &sim_offload.B_offload_data, &B_offload_array);
    
    B_field_init(&sim.B_data, &sim_offload.B_offload_data, B_offload_array);
    B_field_free_offload(&sim_offload.B_offload_data, &B_offload_array);
    
    err = hdf5_close(f);
    return err;
}

int ascotpy_free_bfield() {

    return 0;
}

int ascotpy_bfield_eval_B(int Neval, real* R, real* phi, real* z, real* BR, real* Bphi, real* Bz) {
    real B[3];

    for(int k = 0; k < Neval; k++) {
	B_field_eval_B(B, R[k], phi[k], z[k], &sim.B_data);
	BR[k]   = B[0];
	Bphi[k] = B[1];
	Bz[k]   = B[2];
    }
    
    return 0;
}

int ascotpy_init_efield(char* fn) {
    hdf5_init();
    hid_t f = hdf5_open(fn);
    int err = hdf5_efield_init_offload(f, &sim_offload.E_offload_data, &E_offload_array);
    
    E_field_init(&sim.E_data, &sim_offload.E_offload_data, E_offload_array);
    E_field_free_offload(&sim_offload.E_offload_data, &E_offload_array);
    
    err = hdf5_close(f);
    return err;
}

int ascotpy_free_efield() {

    return 0;
}

int ascotpy_init_plasma(char* fn) {
    hdf5_init();
    hid_t f = hdf5_open(fn);
    int err = hdf5_plasma_init_offload(f, &sim_offload.plasma_offload_data, &plasma_offload_array);
    
    plasma_init(&sim.plasma_data, &sim_offload.plasma_offload_data, plasma_offload_array);
    plasma_free_offload(&sim_offload.plasma_offload_data, &plasma_offload_array);
    
    err = hdf5_close(f);
    return err;
}

int ascotpy_free_plasma() {

    return 0;
}

int ascotpy_init_wall(char* fn) {
    hdf5_init();
    hid_t f = hdf5_open(fn);
    int err = hdf5_wall_init_offload(f, &sim_offload.wall_offload_data, &wall_offload_array);
    
    wall_init(&sim.wall_data, &sim_offload.wall_offload_data, wall_offload_array);
    wall_free_offload(&sim_offload.wall_offload_data, &wall_offload_array);
    
    err = hdf5_close(f);
    return err;
}

int ascotpy_free_wall() {

    return 0;
}
