/**
 * @file libascot.c
 * @brief Library of Ascot5 functions for external use.
 *
 * Functions in this file allows to evaluate input data and quantities using
 * the same methods as is used in actual simulation.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>

#include "ascot5.h"
#include "simulate.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "wall.h"
#include "neutral.h"
#include "boozer.h"
#include "mhd.h"

#include "simulate/mccc/mccc.h"

#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_bfield.h"
#include "hdf5io/hdf5_efield.h"
#include "hdf5io/hdf5_plasma.h"
#include "hdf5io/hdf5_wall.h"
#include "hdf5io/hdf5_neutral.h"
#include "hdf5io/hdf5_boozer.h"
#include "hdf5io/hdf5_mhd.h"

/** Simulation offload struct for holding offload data structs. */
static sim_offload_data sim_offload;

/** Simulation struct for holding data structs. */
static sim_data sim;

static real* Bdata;       /**< Magnetic field data (i.e. offload array) */
static real* Edata;       /**< Electric field data (i.e. offload array) */
static real* plasmadata;  /**< Plasma data (i.e. offload array)         */
static real* walldata;    /**< Wall data (i.e. offload array)           */
static real* neutraldata; /**< Neutral data (i.e. offload array)        */
static real* boozerdata;  /**< Boozer data (i.e. offload array)         */
static real* mhddata;     /**< MHD data (i.e. offload array)            */

/**
 * @brief Initialize input data.
 *
 * This is similar to what is found in hdf5_interface, except we want to control
 * what data fields are initialized. It is user's responsibility not to
 * initialize same data twice.
 *
 * @param fn name of the HDF5 file.
 * @param bfield  flag for initializing magnetic field.
 * @param efield  flag for initializing electric field.
 * @param plasma  flag for initializing plasma data.
 * @param wall    flag for initializing wall data.
 * @param neutral flag for initializing neutral data.
 * @param boozer  flag for initializing boozer data.
 * @param mhd     flag for initializing mhd data.
 *
 * @return zero if initialization succeeded.
 */
int libascot_init(char* fn, int bfield, int efield, int plasma, int wall,
                  int neutral, int boozer, int mhd) {
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

    /* Initialize boozer data if requested. */
    if(boozer) {
        if( hdf5_find_group(f, "/boozer/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/boozer/", qid) ) {
            return 1;
        }
        if( hdf5_boozer_init_offload(f, &sim_offload.boozer_offload_data,
                &boozerdata, qid) ) {
            return 1;
        }
        boozer_init(&sim.boozer_data, &sim_offload.boozer_offload_data,
                    boozerdata);
    }

    /* Initialize mhd data if requested. */
    if(mhd) {
        if( hdf5_find_group(f, "/mhd/") ) {
            return 1;
        }
        if( hdf5_get_active_qid(f, "/mhd/", qid) ) {
            return 1;
        }
        if( hdf5_mhd_init_offload(f, &sim_offload.mhd_offload_data,
                &mhddata, qid) ) {
            return 1;
        }
        mhd_init(&sim.mhd_data, &sim_offload.mhd_offload_data,
                 mhddata);
    }

    if(hdf5_close(f)) {
        return 1;
    }
    return 0;
}

/**
 * @brief Free input data.
 *
 * @param bfield  flag for initializing magnetic field.
 * @param efield  flag for initializing electric field.
 * @param plasma  flag for initializing plasma data.
 * @param wall    flag for initializing wall data.
 * @param neutral flag for initializing neutral data.
 * @param boozer  flag for initializing boozer data.
 * @param mhd     flag for initializing mhd data.
 */
int libascot_free(int bfield, int efield, int plasma, int wall, int neutral,
                  int boozer, int mhd) {
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
    if(boozer) {
        free(boozerdata);
    }
    if(mhd) {
        free(mhddata);
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
 */
void libascot_B_field_eval_B_dB(int Neval, real* R, real* phi, real* z, real* t,
                                real* BR, real* Bphi, real* Bz,
                                real* BR_dR, real* BR_dphi, real* BR_dz,
                                real* Bphi_dR, real* Bphi_dphi, real* Bphi_dz,
                                real* Bz_dR, real* Bz_dphi, real* Bz_dz) {
    real B[12];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_B_dB(B, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
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
}

/**
 * @brief Evaluate normalized poloidal flux at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param rho output array for the normalized poloidal flux.
 * @param psi output array for the poloidal flux.
 */
void libascot_B_field_eval_rho(int Neval, real* R, real* phi, real* z, real* t,
                               real* rho, real* psi) {
    real rhoval[1];
    real psival[1];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psival, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(rhoval, psival[0], &sim.B_data) ) {
            continue;
        }
        rho[k] = rhoval[0];
        psi[k] = psival[0];
    }
}

/**
 * @brief Get magnetic axis at given coordinates.
 *
 * @param Neval number of evaluation points.
 * @param phi phi coordinates of the evaluation points [rad].
 * @param Raxis output array for axis R coordinates.
 * @param zaxis output array for axis z coordinates.
 */
void libascot_B_field_get_axis(int Neval, real* phi, real* Raxis, real* zaxis) {

    for(int k = 0; k < Neval; k++) {
        Raxis[k] = B_field_get_axis_r(&sim.B_data, phi[k]);
        zaxis[k] = B_field_get_axis_z(&sim.B_data, phi[k]);
    }
}

/**
 * @brief Get Rz coordinates for uniform grid of rho values.
 *
 * Purpose of this function is to form a connection between rho and Rz
 * coordinates.
 *
 * This function starts from magnetic axis (at given phi), and takes small steps
 * in the direction of theta until the rho evaluated at that point exceeds
 * rhomin. Rz coordinates of that point are recorded, and steps are continued
 * until rhomax is reached and coordinates of that point are recorded.
 * Interval [(Rmin, zmin), (Rmax,zmax)] is then divided to nrho uniform grid
 * points and rho is evaluated at those points.
 *
 * @param nrho number of evaluation points.
 * @param minrho minimum rho value.
 * @param maxrho maximum rho value.
 * @param theta poloidal angle [rad].
 * @param phi toroidal angle [rad].
 * @param t time coordinate [s].
 * @param r output array for R coordinates [m].
 * @param z output array for R coordinates [m].
 * @param rho output array for rho values.
 */
void libascot_B_field_eval_rhovals(int nrho, real minrho, real maxrho,
                                   real theta, real phi, real t,
                                   real* r, real* z, real* rho) {
    /* Maximum number of steps and step length [m] */
    int NSTEP = 500;
    real step = 0.01;

    real Raxis = B_field_get_axis_r(&sim.B_data, phi);
    real zaxis = B_field_get_axis_z(&sim.B_data, phi);

    real psival, rho0;

    /* Start evaluation from axis */
    real R1 = Raxis;
    real z1 = zaxis;
    if(B_field_eval_psi(&psival, R1, phi, z1, t, &sim.B_data)) {
        return;
    }
    if( B_field_eval_rho(&rho0, psival, &sim.B_data) ) {
        return;
    }

    int iter = 0;
    while(rho0 < minrho && iter < NSTEP) {
        iter++;

        R1 = Raxis + iter*step*cos(theta);
        z1 = zaxis + iter*step*sin(theta);
        if( B_field_eval_psi(&psival, R1, phi, z1, t, &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(&rho0, psival, &sim.B_data) ) {
            continue;
        }
    }
    if(iter > 0) {
        /* Value stored in rho0 is > minrho so we take previous R,z coordinates,
         * since it is the last point where rho0 < minrho. */
        R1 = Raxis + (iter-1)*step*cos(theta);
        z1 = zaxis + (iter-1)*step*sin(theta);
    }

    if(iter == NSTEP) {
        /* Maximum iterations reached. Abort. */
        return;
    }

    real R2 = R1;
    real z2 = z1;

    iter = 0;
    while(rho0 < maxrho && iter < NSTEP) {
        iter++;

        R2 = R1 + iter*step*cos(theta);
        z2 = z1 + iter*step*sin(theta);
        if( B_field_eval_psi(&psival, R2, phi, z2, t, &sim.B_data) ) {
            break;
        }
        if( B_field_eval_rho(&rho0, psival, &sim.B_data) ) {
            break;
        }
    }

    /* With limits determined, we can make a simple linspace from (R1, z1) to
     * (R2, z2) and evaluate rho at those points.                             */
    real rstep = (R2 - R1) / (nrho-1);
    real zstep = (z2 - z1) / (nrho-1);
    for(int i=0; i<nrho; i++) {
        r[i] = R1 + rstep*i;
        z[i] = z1 + zstep*i;
        if( B_field_eval_psi(&psival, r[i], phi, z[i], t, &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(&rho0, psival, &sim.B_data) ) {
            continue;
        }
        rho[i] = rho0;
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
void libascot_E_field_eval_E(int Neval, real* R, real* phi, real* z, real* t,
                             real* ER, real* Ephi, real* Ez) {

    real E[3];
    for(int k = 0; k < Neval; k++) {
        if( E_field_eval_E(E, R[k], phi[k], z[k], t[k],
                           &sim.E_data, &sim.B_data) ) {
            continue;
        }
        ER[k]   = E[0];
        Ephi[k] = E[1];
        Ez[k]   = E[2];
    }
}

/**
 * @brief Get number of plasma species.
 *
 * @return number of plasma species.
 */
int libascot_plasma_get_n_species() {
    return plasma_get_n_species(&sim.plasma_data);
}

/**
 * @brief Get mass and charge of all plasma species.
 *
 * @param mass output array [kg].
 * @param charge output array [C].
 */
void libascot_plasma_get_species_mass_and_charge(real* mass, real* charge) {

    int n_species = plasma_get_n_species(&sim.plasma_data);
    const real* m = plasma_get_species_mass(&sim.plasma_data);
    const real* q = plasma_get_species_charge(&sim.plasma_data);
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
void libascot_plasma_eval_background(int Neval, real* R, real* phi, real* z,
                                     real* t, real* dens, real* temp) {

    int n_species = plasma_get_n_species(&sim.plasma_data);
    real psi[1];
    real rho[1];
    real n[MAX_SPECIES];
    real T[MAX_SPECIES];

    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psi, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(rho, psi[0], &sim.B_data) ) {
            continue;
        }
        if( plasma_eval_densandtemp(n, T, rho[0], R[k], phi[k], z[k], t[k],
                                    &sim.plasma_data) ) {
            continue;
        }
        for(int i=0; i<n_species; i++) {
            dens[k + i*Neval] = n[i];
            temp[k + i*Neval] = T[i];
        }
    }
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
void libascot_neutral_eval_density(int Neval, real* R, real* phi, real* z,
                                   real* t, real* dens) {

    real n0[1];
    for(int k = 0; k < Neval; k++) {
        if( neutral_eval_n0(n0, R[k], phi[k], z[k], t[k], &sim.neutral_data) ) {
            continue;
        }
        dens[k] = n0[0];
    }
}

/**
 * @brief Evaluate boozer coordinates and derivatives.
 */
void libascot_boozer_eval_psithetazeta(int Neval, real* r, real* phi, real* z,
                                       real* t, real* psi, real* theta,
                                       real* zeta, real* dpsidr, real* dpsidphi,
                                       real* dpsidz, real* dthetadr,
                                       real* dthetadphi, real* dthetadz,
                                       real* dzetadr, real* dzetadphi,
                                       real* dzetadz) {
    real psithetazeta[12];
    int isinside;
    for(int k = 0; k < Neval; k++) {
        if( boozer_eval_psithetazeta(psithetazeta, &isinside, r[k], phi[k],
                                     z[k], &sim.boozer_data) ) {
            continue;
        }
        if(!isinside) {
            continue;
        }
        psi[k]        = psithetazeta[0];
        theta[k]      = psithetazeta[4];
        zeta[k]       = psithetazeta[8];
        dpsidr[k]     = psithetazeta[1];
        dpsidphi[k]   = psithetazeta[2];
        dpsidz[k]     = psithetazeta[3];
        dthetadr[k]   = psithetazeta[5];
        dthetadphi[k] = psithetazeta[6];
        dthetadz[k]   = psithetazeta[7];
        dzetadr[k]    = psithetazeta[9];
        dzetadphi[k]  = psithetazeta[10];
        dzetadz[k]    = psithetazeta[11];
    }
}

/**
 * @brief Evaluate MHD perturbation EM-field components
 */
void libascot_mhd_eval_perturbation(int Neval, real* R, real* phi, real* z,
                                    real* t, real* mhd_br, real* mhd_bphi,
                                    real* mhd_bz, real* mhd_er, real* mhd_ephi,
                                    real* mhd_ez) {
    real pert_field[6];
    for(int k = 0; k < Neval; k++) {
        if( mhd_perturbations(pert_field, R[k], phi[k], z[k], t[k],
                              &sim.boozer_data, &sim.mhd_data, &sim.B_data) ) {
            continue;
        }
        mhd_br[k]   = pert_field[0];
        mhd_bphi[k] = pert_field[1];
        mhd_bz[k]   = pert_field[2];
        mhd_er[k]   = pert_field[3];
        mhd_ephi[k] = pert_field[4];
        mhd_ez[k]   = pert_field[5];
    }
}

/**
 * @brief Evaluate collision coefficients.
 */
int libascot_eval_collcoefs(int Neval, real* va, real R, real phi, real z,
                            real t, real ma, real qa, real* F, real* Dpara,
                            real* Dperp, real* K, real* nu, real* Q, real* dQ,
                            real* dDpara, real* clog, real* mu0, real* mu1,
                            real* dmu0) {


    return mccc_eval_coefs(ma, qa, R, phi, z, t, va, Neval,
                           &sim.plasma_data, &sim.B_data,
                           F, Dpara, Dperp, K, nu, Q, dQ, dDpara, clog,
                           mu0, mu1, dmu0);
}
