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
#include "math.h"
#include "simulate.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "wall.h"
#include "neutral.h"
#include "boozer.h"
#include "mhd.h"
#include "consts.h"

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


/**
 * @brief Evaluate magnetic field vector and derivatives at given coordinates.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
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
void libascot_B_field_eval_B_dB(
    sim_offload_data* sim_offload_data, real* B_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* BR, real* Bphi, real* Bz,
    real* BR_dR, real* BR_dphi, real* BR_dz, real* Bphi_dR, real* Bphi_dphi,
    real* Bphi_dz, real* Bz_dR, real* Bz_dphi, real* Bz_dz) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    real B[15];
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
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param rho output array for the normalized poloidal flux.
 * @param psi output array for the poloidal flux.
 */
void libascot_B_field_eval_rho(
    sim_offload_data* sim_offload_data, real* B_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* rho, real* psi) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    real rhoval[2];
    real psival[1];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psival, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }
        psi[k] = psival[0];
        if( B_field_eval_rho(rhoval, psival[0], &sim.B_data) ) {
            continue;
        }
        rho[k] = rhoval[0];
    }
}

/**
 * @brief Get magnetic axis at given coordinates.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param Neval number of evaluation points.
 * @param phi phi coordinates of the evaluation points [rad].
 * @param Raxis output array for axis R coordinates.
 * @param zaxis output array for axis z coordinates.
 */
void libascot_B_field_get_axis(
    sim_offload_data* sim_offload_data, real* B_offload_array, int Neval,
    real* phi, real* Raxis, real* zaxis) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    real axisrz[2];
    for(int k = 0; k < Neval; k++) {
        if( B_field_get_axis_rz(axisrz, &sim.B_data, phi[k]) ) {
            continue;
        }
        Raxis[k] = axisrz[0];
        zaxis[k] = axisrz[1];
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
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
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
void libascot_B_field_eval_rhovals(
    sim_offload_data* sim_offload_data, real* B_offload_array, int nrho,
    real minrho, real maxrho, real theta, real phi, real t, real* r, real* z,
    real* rho) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    /* Maximum number of steps and step length [m] */
    int NSTEP = 50000;
    real step = 0.0001;

    real axisrz[2];
    if( B_field_get_axis_rz(axisrz, &sim.B_data, phi) ) {
        return;
    }
    real Raxis = axisrz[0];
    real zaxis = axisrz[1];

    real psival, rho0[2];

    /* Start evaluation from axis */
    real R1 = Raxis;
    real z1 = zaxis;
    if(B_field_eval_psi(&psival, R1, phi, z1, t, &sim.B_data)) {
        return;
    }
    if( B_field_eval_rho(rho0, psival, &sim.B_data) ) {
        return;
    }

    int iter = 0;
    while(rho0[0] < minrho && iter < NSTEP) {
        iter++;

        R1 = Raxis + iter*step*cos(theta);
        z1 = zaxis + iter*step*sin(theta);
        if( B_field_eval_psi(&psival, R1, phi, z1, t, &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(rho0, psival, &sim.B_data) ) {
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
    while(rho0[0] < maxrho && iter < NSTEP) {
        iter++;

        R2 = R1 + iter*step*cos(theta);
        z2 = z1 + iter*step*sin(theta);
        if( B_field_eval_psi(&psival, R2, phi, z2, t, &sim.B_data) ) {
            break;
        }
        if( B_field_eval_rho(rho0, psival, &sim.B_data) ) {
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
        if( B_field_eval_rho(rho0, psival, &sim.B_data) ) {
            continue;
        }
        rho[i] = rho0[0];
    }
}

/**
 * @brief Evaluate electric field vector at given coordinates.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param E_offload_array initialized electric field offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param ER output array [V/m].
 * @param Ephi output array [V/m].
 * @param Ez output array [V/m].
 */
void libascot_E_field_eval_E(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* E_offload_array, int Neval, real* R, real* phi, real* z, real* t,
    real* ER, real* Ephi, real* Ez) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    E_field_init(&sim.E_data, &sim_offload_data->E_offload_data,
                 E_offload_array);
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
 * @param sim_offload_data initialized simulation offload data struct
 * @param plasma_offload_array initialized plasma offload data
 *
 * @return number of plasma species.
 */
int libascot_plasma_get_n_species(
    sim_offload_data* sim_offload_data, real* plasma_offload_array) {

    sim_data sim;
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
    return plasma_get_n_species(&sim.plasma_data);
}

/**
 * @brief Get mass and charge of all plasma species.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param plasma_offload_array initialized plasma offload data
 * @param mass output array [kg].
 * @param charge output array [C].
 */
void libascot_plasma_get_species_mass_and_charge(
    sim_offload_data* sim_offload_data, real* plasma_offload_array, real* mass,
    real* charge) {

    sim_data sim;
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
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
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized mangetic field offload data
 * @param plasma_offload_array initialized plasma offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param dens output array [m^-3].
 * @param temp output array [eV].
 */
void libascot_plasma_eval_background(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* plasma_offload_array, int Neval, real* R, real* phi, real* z, real* t,
    real* dens, real* temp) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
    int n_species = plasma_get_n_species(&sim.plasma_data);
    real psi[1];
    real rho[2];
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
            temp[k + i*Neval] = T[i]/CONST_E;
        }
    }
}

/**
 * @brief Evaluate neutral density at given coordinates.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param neutral_offload_array initialized neutral offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param dens output array [m^-3].
 */
void libascot_neutral_eval_density(
    sim_offload_data* sim_offload_data, real* neutral_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* dens) {

    sim_data sim;
    neutral_init(&sim.neutral_data, &sim_offload_data->neutral_offload_data,
                 neutral_offload_array);
    real psi[1];
    real rho[2];
    real n0[1];
    for(int k = 0; k < Neval; k++) {
        if( B_field_eval_psi(psi, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(rho, psi[0], &sim.B_data) ) {
            continue;
        }
        if( neutral_eval_n0(n0, rho[0], R[k], phi[k], z[k], t[k],
                            &sim.neutral_data) ) {
            continue;
        }
        dens[k] = n0[0];
    }
}

/**
 * @brief Evaluate boozer coordinates and derivatives.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param boozer_offload_array initialized boozer offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param psi output array
 * @param theta output array
 * @param zeta output array
 * @param dpsidr output array
 * @param dpsidphi output array
 * @param dpsidz output array
 * @param dthetadr output array
 * @param dthetadphi output array
 * @param dthetadz output array
 * @param dzetadr output array
 * @param dzetadphi output array
 * @param dzetadz output array
 * @param rho output array
 */
void libascot_boozer_eval_psithetazeta(
    sim_offload_data* sim_offload_data, real* boozer_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* psi, real* theta, real* zeta,
    real* dpsidr, real* dpsidphi, real* dpsidz, real* dthetadr,
    real* dthetadphi, real* dthetadz, real* dzetadr, real* dzetadphi,
    real* dzetadz, real* rho) {

    sim_data sim;
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    real psithetazeta[12];
    real rhoval[2];
    int isinside;
    for(int k = 0; k < Neval; k++) {
        if( boozer_eval_psithetazeta(psithetazeta, &isinside, R[k], phi[k],
                                     z[k], &sim.B_data, &sim.boozer_data) ) {
            continue;
        }
        if(!isinside) {
            continue;
        }
        if( B_field_eval_rho(rhoval, psithetazeta[0], &sim.B_data) ) {
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
        rho[k]        = rhoval[0];
    }
}

/**
 * @brief Evaluate boozer coordinates related quantities.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param boozer_offload_array initialized boozer offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param qprof array for storing the (flux averaged) safety factor.
 * @param jac array for storing the coordinate Jacobian.
 * @param jacb2 array for storing the coordinate Jacobian multiplied with B^2.
 */
void libascot_boozer_eval_fun(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* boozer_offload_array, int Neval, real* R, real* phi, real* z, real* t,
    real* qprof, real* jac, real* jacB2) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    real psithetazeta[12];
    real B[15];
    int isinside;
    for(int k = 0; k < Neval; k++) {
        if( boozer_eval_psithetazeta(psithetazeta, &isinside, R[k], phi[k],
                                     z[k], &sim.B_data, &sim.boozer_data) ) {
            continue;
        }
        if(!isinside) {
            continue;
        }
        if( B_field_eval_B_dB(B, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }

        real bvec[]      = {B[0], B[4], B[8]};
        real gradpsi[]   = {psithetazeta[1],
                            psithetazeta[2]/R[k],
                            psithetazeta[3]};
        real gradtheta[] = {psithetazeta[5],
                            psithetazeta[6]/R[k],
                            psithetazeta[7]};
        real gradzeta[]  = {psithetazeta[9],
                            psithetazeta[10]/R[k],
                            psithetazeta[11]};

        real veca[3], vecb[3];

        math_cross(gradpsi, gradzeta, veca);
        math_cross(gradpsi, gradtheta, vecb);
        qprof[k] = (veca[1] - bvec[1]) / vecb[1];

        math_cross(gradtheta, gradzeta, veca);
        jac[k]   = 1/math_dot(veca, gradpsi);
        jacB2[k] = jac[k]*math_norm(bvec)*math_norm(bvec);
    }
}

/**
 * @brief Evaluate MHD perturbation potentials
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param boozer_offload_array initialized boozer offload data
 * @param mhd_offload_array initialized MHD offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param alpha output array
 * @param dadr output array
 * @param dadphi output array
 * @param dadz output array
 * @param dadt output array
 * @param Phi output array
 * @param dPhidr output array
 * @param dPhidphi output array
 * @param dPhidz output array
 * @param dPhidt output array
 */
void libascot_mhd_eval(
    sim_offload_data* sim_offload_data, real* boozer_offload_array,
    real* mhd_offload_array, int Neval, real* R, real* phi, real* z, real* t,
    real* alpha, real* dadr, real* dadphi, real* dadz, real* dadt, real* Phi,
    real* dPhidr, real* dPhidphi, real* dPhidz, real* dPhidt) {

    sim_data sim;
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);

    real mhd_dmhd[10];
    for(int k = 0; k < Neval; k++) {
        if( mhd_eval(mhd_dmhd, R[k], phi[k], z[k], t[k],
                     &sim.boozer_data, &sim.mhd_data, &sim.B_data) ) {
            continue;
        }
        alpha[k]    = mhd_dmhd[0];
        dadr[k]     = mhd_dmhd[2];
        dadphi[k]   = mhd_dmhd[3];
        dadz[k]     = mhd_dmhd[4];
        dadt[k]     = mhd_dmhd[1];
        Phi[k]      = mhd_dmhd[5];
        dPhidr[k]   = mhd_dmhd[7];
        dPhidphi[k] = mhd_dmhd[8];
        dPhidz[k]   = mhd_dmhd[9];
        dPhidt[k]   = mhd_dmhd[6];
    }
}

/**
 * @brief Evaluate MHD perturbation EM-field components
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param boozer_offload_array initialized boozer offload data
 * @param mhd_offload_array initialized MHD offload data
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param mhd_br output array
 * @param mhd_bphi output array
 * @param mhd_bz output array
 * @param mhd_er output array
 * @param mhd_ephi output array
 * @param mhd_ez output array
 * @param mhd_phi output array
 */
void libascot_mhd_eval_perturbation(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* boozer_offload_array, real* mhd_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* mhd_br, real* mhd_bphi,
    real* mhd_bz, real* mhd_er, real* mhd_ephi, real* mhd_ez, real* mhd_phi) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);
    real pert_field[7];
    int onlypert = 1;
    for(int k = 0; k < Neval; k++) {
        if( mhd_perturbations(pert_field, R[k], phi[k], z[k], t[k], onlypert,
                              &sim.boozer_data, &sim.mhd_data, &sim.B_data) ) {
            continue;
        }
        mhd_br[k]   = pert_field[0];
        mhd_bphi[k] = pert_field[1];
        mhd_bz[k]   = pert_field[2];
        mhd_er[k]   = pert_field[3];
        mhd_ephi[k] = pert_field[4];
        mhd_ez[k]   = pert_field[5];
        mhd_phi[k]  = pert_field[6];
    }
}

/**
 * @brief Evaluate collision coefficients.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param plasma_offload_array initialized plasma offload data
 * @param Neval number of evaluation points.
 * @param va test particle velocities at the evaluation point [m/s].
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param F output array
 * @param Dpara output array
 * @param Dperp output array
 * @param K output array
 * @param nu output array
 * @param Q output array
 * @param dQ output array
 * @param dDpara output array
 * @param clog output array
 * @param mu0 output array
 * @param mu1 output array
 * @param dmu0 output array
 */
int libascot_eval_collcoefs(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* plasma_offload_array, int Neval, real* va, real R, real phi, real z,
    real t, real ma, real qa, real* F, real* Dpara, real* Dperp, real* K,
    real* nu, real* Q, real* dQ, real* dDpara, real* clog, real* mu0, real* mu1,
    real* dmu0) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
    return mccc_eval_coefs(
        ma, qa, R, phi, z, t, va, Neval, &sim.plasma_data, &sim.B_data,
        F, Dpara, Dperp, K, nu, Q, dQ, dDpara, clog, mu0, mu1, dmu0);
}
