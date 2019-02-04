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
#include "hdf5io/hdf5_neutral.h"

/** Simulation offload struct for holding offload data structs. */
static sim_offload_data sim_offload;

/** Simulation struct for holding data structs. */
static sim_data sim;

static real* Bdata;       /**< Magnetic field data (i.e. offload array) */
static real* Edata;       /**< Electric field data (i.e. offload array) */
static real* plasmadata;  /**< Plasma data (i.e. offload array)         */
static real* walldata;    /**< Wall data (i.e. offload array)           */
static real* neutraldata; /**< Neutral data (i.e. offload array)        */

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
        return 1;
    }

    char qid[11];

    /* Initialize magnetic field data if requested. */
    if(bfield) {
        if( hdf5_find_group(f, "/bfield/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/bfield/", qid) ) {
            return 1;
        }
        if( hdf5_bfield_init_offload(f, &sim_offload.B_offload_data,
                &Bdata, qid) ) {
            return 1;
        }
        B_field_init(&sim.B_data, &sim_offload.B_offload_data, Bdata);
    }

    /* Initialize electric field data if requested. */
    if(efield) {
        if( hdf5_find_group(f, "/efield/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/efield/", qid) ) {
            return 1;
        }
        if( hdf5_efield_init_offload(f, &sim_offload.E_offload_data,
                &Edata, qid) ) {
            return 1;
        }
        E_field_init(&sim.E_data, &sim_offload.E_offload_data, Edata);
    }

    /* Initialize plasma data if requested. */
    if(plasma) {
        if( hdf5_find_group(f, "/plasma/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/plasma/", qid) ) {
            return 1;
        }
        if( hdf5_plasma_init_offload(f, &sim_offload.plasma_offload_data,
                &plasmadata, qid) ) {
            return 1;
        }
        plasma_init(&sim.plasma_data, &sim_offload.plasma_offload_data,
                plasmadata);
    }

    /* Initialize wall data if requested. */
    if(wall) {
        if( hdf5_find_group(f, "/wall/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/wall/", qid) ) {
            return 1;
        }
        if( hdf5_wall_init_offload(f, &sim_offload.wall_offload_data,
                &walldata, qid) ) {
            return 1;
        }
        wall_init(&sim.wall_data, &sim_offload.wall_offload_data,
                walldata);
    }

    /* Initialize neutral data if requested. */
    if(neutral) {
        if( hdf5_find_group(f, "/neutral/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/neutral/", qid) ) {
            return 1;
        }
        if( hdf5_neutral_init_offload(f, &sim_offload.neutral_offload_data,
                &neutraldata, qid) ) {
            return 1;
        }
        neutral_init(&sim.neutral_data, &sim_offload.neutral_offload_data,
                neutraldata);
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
    if(efield) {
        free(Edata);
    }
    if(plasma) {
        free(plasmadata);
    }
    if(wall) {
        free(walldata);
    }
    if(neutral) {
        free(neutraldata);
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
 * @param t time coordinates of the evaluation points [s].
 * @param BR output array [T].
 * @param Bphi output array [T].
 * @param Bz output array [T].
 * @param BR_dR output array [T].
 * @param BR_dphi output array [T].
 * @param BR_dz output array [T].
 * @param Bphi_dR output array [T].
 * @param Bphi_dphi output array [T].
 * @param Bphi_dz output array [T].
 * @param Bz_dR output array [T].
 * @param Bz_dphi output array [T].
 * @param Bz_dz output array [T].
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_B_field_eval_B_dB(int Neval, real* R, real* phi, real* z, real* t,
                              real* BR, real* Bphi, real* Bz,
                              real* BR_dR, real* BR_dphi, real* BR_dz,
                              real* Bphi_dR, real* Bphi_dphi, real* Bphi_dz,
                              real* Bz_dR, real* Bz_dphi, real* Bz_dz) {
    real B[12];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_B_dB(B, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            return 1;
        }
        BR[k]        = B[0];
        Bphi[k]      = B[4];
        Bz[k]        = B[8];
        BR_dR[k]     = B[1];
        BR_dphi[k]   = B[2];
        BR_dz[k]     = B[3];
        Bphi_dR[k]   = B[5];
        Bphi_dphi[k] = B[6];
        Bphi_dz[k]   = B[7];
        Bz_dR[k]     = B[9];
        Bz_dphi[k]   = B[10];
        Bz_dz[k]     = B[11];
    }
    return 0;
}

/**
 * @brief Evaluate poloidal flux at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param psi output array.
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_B_field_eval_psi(int Neval, real* R, real* phi, real* z, real* t,
                             real* psi) {
    real psival[1];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psival, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            return 1;
        }
        psi[k] = psival[0];
    }
    return 0;
}

/**
 * @brief Evaluate normalized poloidal flux at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param psi output array.
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_B_field_eval_rho(int Neval, real* R, real* phi, real* z, real* t,
                             real* rho) {
    real rhoval[1];
    real psival[1];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psival, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            return 1;
        }
        if( B_field_eval_rho(rhoval, psival[0], &sim.B_data) ) {
            return 1;
        }
        rho[k] = rhoval[0];
    }
    return 0;
}

/**
 * @brief Get magnetic axis at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param phi phi coordinates of the evaluation points [rad].
 * @param psi output array.
 */
void ascotpy_B_field_get_axis(int Neval, real* phi, real* Raxis, real* zaxis) {

    for(int k = 0; k < Neval; k++) {
        Raxis[k] = B_field_get_axis_r(&sim.B_data, phi[k]);
        zaxis[k] = B_field_get_axis_z(&sim.B_data, phi[k]);
    }
}

/**
 * @brief Evaluate electric field vector at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param ER output array [V/m].
 * @param Ephi output array [V/m].
 * @param Ez output array [V/m].
 */
int ascotpy_E_field_eval_E(int Neval, real* R, real* phi, real* z, real* t,
                           real* ER, real* Ephi, real* Ez) {

    real E[3];
    for(int k = 0; k < Neval; k++) {
        if( E_field_eval_E(E, R[k], phi[k], z[k], t[k],
                           &sim.E_data, &sim.B_data) ) {
            return 1;
        }
        ER[k]   = E[0];
        Ephi[k] = E[1];
        Ez[k]   = E[2];
    }
    return 0;
}

/**
 * @brief Get number of plasma species.
 *
 * @return number of plasma species.
 */
int ascotpy_plasma_get_n_species() {
    return plasma_get_n_species(&sim.plasma_data);
}

/**
 * @brief Get mass and charge of all plasma species.
 *
 * @param mass output array [kg].
 * @param charge output array [C].
 */
void ascotpy_plasma_get_species_mass_and_charge(real* mass, real* charge) {

    int n_species = plasma_get_n_species(&sim.plasma_data);
    real* m = plasma_get_species_mass(&sim.plasma_data);
    real* q = plasma_get_species_charge(&sim.plasma_data);
    for(int i=0; i<n_species; i++) {
        mass[i]   = m[i];
        charge[i] = q[i];
    }
}

/**
 * @brief Evaluate plasma density and temperature at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param dens output array [m^-3].
 * @param temp output array [eV].
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_plasma_eval_background(int Neval, real* R, real* phi, real* z,
                                   real* t, real* dens, real* temp) {

    int n_species = plasma_get_n_species(&sim.plasma_data);
    real psi[1];
    real rho[1];
    real n[n_species];
    real T[n_species];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psi, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            return 1;
        }
        if( B_field_eval_rho(rho, psi[0], &sim.B_data) ) {
            return 1;
        }
        if( plasma_eval_densandtemp(rho[0], t[k], &sim.plasma_data, n, T) ) {
            return 1;
        }
        for(int i=0; i<n_species; i++) {
            dens[k + i*Neval] = n[i];
            temp[k + i*Neval] = T[i];
        }
    }
    return 0;
}

/**
 * @brief Evaluate neutral density at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param dens output array [m^-3].
 *
 * @return zero if evaluation succeeded.
 */
int ascotpy_neutral_eval_density(int Neval, real* R, real* phi, real* z,
                                 real* t, real* dens) {

    real n0[1];
    for(int k = 0; k < Neval; k++) {
        if( neutral_eval_n0(n0, R[k], phi[k], z[k], t[k], &sim.neutral_data) ) {
            return 1;
        }
        dens[k] = n0[0];
    }
    return 1;
}
