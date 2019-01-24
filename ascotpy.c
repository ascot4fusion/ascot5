/**
 * @file ascotpy.c
 * @brief C side of the interactive Python interface.
 *
 * Functions in this file can be called from Python directly via ascotpy
 * module.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf5.h>

#include "ascot5.h"
#include "simulate.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "wall.h"
#include "neutral.h"

#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_bfield.h"
#include "hdf5io/hdf5_efield.h"
#include "hdf5io/hdf5_plasma.h"
#include "hdf5io/hdf5_wall.h"

/** Simulation offload struct for holding offload data structs. */
static sim_offload_data sim_offload;

/** Simulation struct for holding data structs. */
static sim_data sim;

static real* Bdata;      /**< Magnetic field data (i.e. offload array) */
//static real* Edata;      /**< Electric field data (i.e. offload array) */
//static real* plasmadata; /**< Plasma data (i.e. offload array)         */
//static real* walldata;   /**< Wall data (i.e. offload array)           */

/**
 * @brief Initialize input data.
 *
 * This is similar to what is found in hdf5_interface, except we want to control
 * what data fields are initialized. It is user's responsibility not to
 * initialize same data twice.
 *
 * @param fn name of the HDF5 file.
 * @param bfield flag for initializing magnetic field.
 * @param efield flag for initializing electric field.
 * @param plasma flag for initializing plasma data.
 * @param wall flag for initializing wall data.
 * @param neutral flag for initializing neutral data.
 *
 * @return zero if initialization succeeded.
 */
int ascotpy_init(char* fn, int bfield, int efield, int plasma, int wall,
                 int neutral) {
    hdf5_init();
    hid_t f = hdf5_open(fn);
    if(f < 0) {
        /* Failed to open file. */
        return 1;
    }

    char qid[11];

    /* Initialize magnetic field data if requested. */
    if(bfield) {
        if( hdf5_find_group(f, "/bfield/") ) {
            /* Data does not exist. */
            return 1;
        }
        if( hdf5_get_active_qid(f, "/bfield/", qid) ) {
            /* No active qid found. */
            return 1;
        }
        if( hdf5_bfield_init_offload(f, &sim_offload.B_offload_data, &Bdata,
                                     qid) ) {
            /* Initialization failed. */
            return 1;
        }
        B_field_init(&sim.B_data, &sim_offload.B_offload_data, Bdata);
    }

    if(hdf5_close(f)) {
        return 1;
    }
    return 0;
}

/**
 * @brief Free input data.
 *
 * @param bfield flag for initializing magnetic field.
 * @param efield flag for initializing electric field.
 * @param plasma flag for initializing plasma data.
 * @param wall flag for initializing wall data.
 * @param neutral flag for initializing neutral data.
 */
int ascotpy_free(int bfield, int efield, int plasma, int wall, int neutral) {
    if(bfield) {
        free(Bdata);
    }
    return 0;
}

/**
 * @brief Evaluate magnetic field vector and derivatives at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param BR
 * @param Bphi
 * @param Bz
 * @param BR_dR
 * @param BR_dphi
 * @param BR_dz
 * @param Bphi_dR
 * @param Bphi_dphi
 * @param Bphi_dz
 * @param Bz_dR
 * @param Bz_dphi
 * @param Bz_dz
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_bfield_eval_B_dB(int Neval, real* R, real* phi, real* z,
                             real* BR, real* Bphi, real* Bz,
                             real* BR_dR, real* BR_dphi, real* BR_dz,
                             real* Bphi_dR, real* Bphi_dphi, real* Bphi_dz,
                             real* Bz_dR, real* Bz_dphi, real* Bz_dz) {
    real B[3];

    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_B(B, R[k], phi[k], z[k], &sim.B_data) ) {
            return 1;
        }
        BR[k]   = B[0];
        Bphi[k] = B[1];
        Bz[k]   = B[2];
    }
    return 0;
}
