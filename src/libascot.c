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
#include "gitver.h"
#include "math.h"
#include "simulate.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "wall.h"
#include "neutral.h"
#include "boozer.h"
#include "mhd.h"
#include "asigma.h"
#include "consts.h"
#include "physlib.h"

#include "simulate/mccc/mccc_coefs.h"

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

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real B[15];
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
 * @param drhodpsi output array for normalized poloidal flux psi derivative.
 * @param psi output array for the poloidal flux [Wb].
 * @param dpsidr output array for the poloidal flux R derivative [Wb/m].
 * @param dpsidphi output array for the poloidal flux phi derivative [Wb/rad].
 * @param dpsidz output array for the poloidal flux z derivative [Wb/m].
 */
void libascot_B_field_eval_rho(
    sim_offload_data* sim_offload_data, real* B_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* rho, real* drhodpsi, real* psi,
    real* dpsidr, real* dpsidphi, real* dpsidz) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real rhoval[2], psival[4];
        if( B_field_eval_psi_dpsi(psival, R[k], phi[k], z[k], t[k],
                                  &sim.B_data) ) {
            continue;
        }
        psi[k]      = psival[0];
        dpsidr[k]   = psival[1];
        dpsidphi[k] = psival[2];
        dpsidz[k]   = psival[3];
        if( B_field_eval_rho(rhoval, psival[0], &sim.B_data) ) {
            continue;
        }
        rho[k]      = rhoval[0];
        drhodpsi[k] = rhoval[1];
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
    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real axisrz[2];
        if( B_field_get_axis_rz(axisrz, &sim.B_data, phi[k]) ) {
            continue;
        }
        Raxis[k] = axisrz[0];
        zaxis[k] = axisrz[1];
    }
}

/**
 * @brief Map (rho, theta, phi) to (R,z) coordinates.
 *
 * This function implements the Newton method. If the function fails on
 * a given position, the corresponding (R,z) values in the output arrays are
 * not altered.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param Neval number of query points.
 * @param rho the square root of the normalized poloidal flux values.
 * @param theta poloidal angles [rad].
 * @param phi toroidal angles [rad].
 * @param t time coordinate (same for all) [s].
 * @param maxiter maximum number of iterations in Newton algorithm.
 * @param tol algorithm is stopped when |rho - rho(r,z)| < tol
 * @param r output array for R coordinates [m].
 * @param z output array for z coordinates [m].
 */
void libascot_B_field_rhotheta2rz(
    sim_offload_data* sim_offload_data, real* B_offload_array, int Neval,
    real* rho, real* theta, real* phi, real t, int maxiter, real tol,
    real* r, real* z) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);

    #pragma omp parallel for
    for(int j=0; j<Neval; j++) {
        real axisrz[2];
        real rhodrho[4];
        if( B_field_get_axis_rz(axisrz, &sim.B_data, phi[j]) ) {
            continue;
        }
        if( B_field_eval_rho_drho(rhodrho, axisrz[0], phi[j], axisrz[1],
                                  &sim.B_data)) {
            continue;
        }
        if( rhodrho[0] > rho[j] ) {
            /* Due to padding, rho might not be exactly zero on the axis so we
             * return the axis position for small values of queried rho */
            r[j] = axisrz[0];
            z[j] = axisrz[1];
            continue;
        }

        real x = 1e-1;
        real rj, zj;
        real costh = cos(theta[j]);
        real sinth = sin(theta[j]);
        for(int i=0; i<maxiter; i++) {
            rj = axisrz[0] + x * costh;
            zj = axisrz[1] + x * sinth;
            if( B_field_eval_rho_drho(rhodrho, rj, phi[j], zj, &sim.B_data) ) {
                break;
            }
            if( fabs(rho[j] - rhodrho[0]) < tol ) {
                r[j] = rj;
                z[j] = zj;
                break;
            }

            real drhodx = costh * rhodrho[1] + sinth * rhodrho[3];
            x = x - (rhodrho[0] - rho[j]) / drhodx;
            if( x < 0 ) {
                /* Try again starting closer from the axis */
                x = (x + (rhodrho[0] - rho[j]) / drhodx) / 2;
            }
        }
    }
}

/**
 * @brief Find psi on axis using the gradient descent method
 *
 * Note that the psi value is not returned in case this algorithm fails.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param psi value of psi on axis if this function did not fail
 * @param rz initial (R,z) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void libascot_B_field_gradient_descent(
    sim_offload_data* sim_offload_data, real* B_offload_array, real psi[1],
    real rz[2], real step, real tol, int maxiter, int ascent) {
    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);

    if(ascent) {
        step = -1 * step;
    }

    real phi = 0.0, time = 0.0;
    real psidpsi[4], nextrz[2];
    B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, &sim.B_data);

    int iter = 0;
    while (1) {
        if( B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time,
                                  &sim.B_data) ) {
            break;
        }
        nextrz[0] = rz[0] - step * psidpsi[1];
        nextrz[1] = rz[1] - step * psidpsi[3];

        // Check convergence
        if(sqrt( (nextrz[0] - rz[0]) * (nextrz[0] - rz[0])
                + (nextrz[1] - rz[1]) * (nextrz[1] - rz[1]) ) < tol) {
            psi[0] = psidpsi[0];
            rz[0] = nextrz[0];
            rz[1] = nextrz[1];

            // Add a bit of padding
            B_field_eval_psi_dpsi(
                psidpsi, rz[0], phi, rz[1], time, &sim.B_data);
            if(ascent) {
                psi[0] = psi[0] + (tol * psidpsi[1] + tol * psidpsi[3]);
            } else {
                psi[0] = psi[0] - (tol * psidpsi[1] + tol * psidpsi[3]);
            }
            break;
        }

        rz[0] = nextrz[0];
        rz[1] = nextrz[1];
        iter++;

        if(iter == maxiter) {
            break;
        }
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

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real E[3];
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
 * @param mass mass output array [kg].
 * @param charge charge output array [C].
 * @param anum atomic mass number output array [1].
 * @param znum charge number output array [1].
 */
void libascot_plasma_get_species_mass_and_charge(
    sim_offload_data* sim_offload_data, real* plasma_offload_array, real* mass,
    real* charge, int* anum, int* znum) {

    sim_data sim;
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
    int n_species = plasma_get_n_species(&sim.plasma_data);
    const real* m = plasma_get_species_mass(&sim.plasma_data);
    const real* q = plasma_get_species_charge(&sim.plasma_data);
    const int* a  = plasma_get_species_anum(&sim.plasma_data);
    const int* z  = plasma_get_species_znum(&sim.plasma_data);
    mass[0]   = CONST_M_E;
    charge[0] = -CONST_E;
    anum[0]   = 0;
    znum[0]   = 0;
    for(int i=1; i<n_species; i++) {
        mass[i]   = m[i];
        charge[i] = q[i];
        anum[i]   = a[i-1];
        znum[i]   = z[i-1];
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

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real psi[1], rho[2], n[MAX_SPECIES], T[MAX_SPECIES];
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
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* neutral_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* dens) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    neutral_init(&sim.neutral_data, &sim_offload_data->neutral_offload_data,
                 neutral_offload_array);

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real psi[1], rho[2], n0[1];
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
 * @param B_offload_array initialized magnetic field offload data
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
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* boozer_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, real* psi, real* theta, real* zeta,
    real* dpsidr, real* dpsidphi, real* dpsidz, real* dthetadr,
    real* dthetadphi, real* dthetadz, real* dzetadr, real* dzetadphi,
    real* dzetadz, real* rho) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        int isinside;
        real psithetazeta[12], rhoval[2];
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
 * @param jacB2 array for storing the coordinate Jacobian multiplied with B^2.
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

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        int isinside;
        real psithetazeta[12], B[15];
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
        jac[k]   = -1.0 / math_dot(veca, gradpsi);
        jacB2[k] = jac[k]*math_norm(bvec)*math_norm(bvec);
    }
}

/**
 * @brief Get number of MHD modes.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param mhd_offload_array initialized mhd offload data
 *
 * @return number of MHD modes
 */
int libascot_mhd_get_n_modes(
    sim_offload_data* sim_offload_data, real* mhd_offload_array) {

    sim_data sim;
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);
    return mhd_get_n_modes(&sim.mhd_data);
}

/**
 * @brief Get MHD mode amplitude, frequency, phase, and mode numbers
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param mhd_offload_array initialized mhd offload data
 * @param nmode output array for toroidal mode number
 * @param mmode output array for poloidal mode number
 * @param amplitude output array for mode amplitude
 * @param omega output array for mode frequency
 * @param phase output array for mode phase
 */
void libascot_mhd_get_mode_specs(
    sim_offload_data* sim_offload_data, real* mhd_offload_array, int* nmode,
    int* mmode, real* amplitude, real* omega, real* phase) {

    sim_data sim;
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);
    int n_modes   = mhd_get_n_modes(&sim.mhd_data);
    const int* n  = mhd_get_nmode(&sim.mhd_data);
    const int* m  = mhd_get_mmode(&sim.mhd_data);
    const real* a = mhd_get_amplitude(&sim.mhd_data);
    const real* o = mhd_get_frequency(&sim.mhd_data);
    const real* p = mhd_get_phase(&sim.mhd_data);
    for(int i=0; i<n_modes; i++) {
        nmode[i]     = n[i];
        mmode[i]     = m[i];
        amplitude[i] = a[i];
        omega[i]     = o[i];
        phase[i]     = p[i];
    }
}

/**
 * @brief Evaluate MHD perturbation potentials
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
 * @param includemode mode index to include or MHD_INCLUDE_ALL
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
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* boozer_offload_array, real* mhd_offload_array, int Neval,
    real* R, real* phi, real* z, real* t, int includemode,
    real* alpha, real* dadr, real* dadphi, real* dadz, real* dadt, real* Phi,
    real* dPhidr, real* dPhidphi, real* dPhidz, real* dPhidt) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);

    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real mhd_dmhd[10];
        if( mhd_eval(mhd_dmhd, R[k], phi[k], z[k], t[k], includemode,
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
 * @param includemode mode index to include or MHD_INCLUDE_ALL
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
    real* R, real* phi, real* z, real* t, int includemode, real* mhd_br,
    real* mhd_bphi, real* mhd_bz, real* mhd_er, real* mhd_ephi, real* mhd_ez,
    real* mhd_phi) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    boozer_init(&sim.boozer_data, &sim_offload_data->boozer_offload_data,
                boozer_offload_array);
    mhd_init(&sim.mhd_data, &sim_offload_data->mhd_offload_data,
             mhd_offload_array);
    int onlypert = 1;
    #pragma omp parallel for
    for(int k = 0; k < Neval; k++) {
        real pert_field[7];
        if( mhd_perturbations(pert_field, R[k], phi[k], z[k], t[k], onlypert,
                              includemode, &sim.boozer_data, &sim.mhd_data,
                              &sim.B_data) ) {
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
 * @brief Evaluate collision coefficients
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param plasma_offload_array initialized plasma offload data
 * @param Neval number of evaluation points
 * @param R R coordinates of the evaluation points [m]
 * @param phi phi coordinates of the evaluation points [rad]
 * @param z z coordinates of the evaluation points [m]
 * @param t time coordinates of the evaluation points [s]
 * @param Nv number of velocity grid points
 * @param va test particle velocities at the evaluation point [m/s]
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
void libascot_eval_collcoefs(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* plasma_offload_array, int Neval, real* R, real* phi, real* z, real* t,
    int Nv, real* va, real ma, real qa, real* F, real* Dpara, real* Dperp,
    real* K, real* nu, real* Q, real* dQ, real* dDpara, real* clog,
    real* mu0, real* mu1, real* dmu0) {

    sim_data sim;
    sim.mccc_data.usetabulated = 0;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);

    /* Evaluate plasma parameters */
    int n_species  = plasma_get_n_species(&sim.plasma_data);
    const real* qb = plasma_get_species_charge(&sim.plasma_data);
    const real* mb = plasma_get_species_mass(&sim.plasma_data);

    #pragma omp parallel for
    for(int k=0; k<Neval; k++) {
        real mufun[3] = {0., 0., 0.};

        /* Evaluate rho as it is needed to evaluate plasma parameters */
        real psi, rho[2];
        if( B_field_eval_psi(&psi, R[k], phi[k], z[k], t[k], &sim.B_data) ) {
            continue;
        }
        if( B_field_eval_rho(rho, psi, &sim.B_data) ) {
            continue;
        }

        real nb[MAX_SPECIES], Tb[MAX_SPECIES];
        if( plasma_eval_densandtemp(nb, Tb, rho[0], R[k], phi[k], z[k], t[k],
                                    &sim.plasma_data) ) {
            continue;
        }

        /* Evaluate coefficients for different velocities */
        for(int iv=0; iv<Nv; iv++) {

            /* Loop through all plasma species */
            for(int ib=0; ib<n_species; ib++) {

                /* Coulomb logarithm */
                real clogab[MAX_SPECIES];
                mccc_coefs_clog(clogab, ma, qa, va[iv], n_species, mb, qb,
                                nb, Tb);

                /* Special functions */
                real vb = sqrt( 2 * Tb[ib] / mb[ib] );
                real x  = va[iv] / vb;
                mccc_coefs_mufun(mufun, x, &sim.mccc_data);

                /* Coefficients */
                real Fb      = mccc_coefs_F(ma, qa, mb[ib], qb[ib], nb[ib], vb,
                                            clogab[ib], mufun[0]);
                real Qb      = mccc_coefs_Q(ma, qa, mb[ib], qb[ib], nb[ib], vb,
                                            clogab[ib], mufun[0]);
                real dQb     = mccc_coefs_dQ(ma, qa, mb[ib], qb[ib], nb[ib], vb,
                                             clogab[ib], mufun[2]);
                real Dparab  = mccc_coefs_Dpara(ma, qa, va[iv], qb[ib], nb[ib],
                                                vb, clogab[ib], mufun[0]);
                real Dperpb  = mccc_coefs_Dperp(ma, qa, va[iv], qb[ib], nb[ib],
                                                vb, clogab[ib], mufun[1]);
                real dDparab = mccc_coefs_dDpara(ma, qa, va[iv], qb[ib], nb[ib],
                                                 vb, clogab[ib], mufun[0],
                                                 mufun[2]);
                real Kb      = mccc_coefs_K(va[iv], Dparab, dDparab, Qb);
                real nub     = mccc_coefs_nu(va[iv], Dperpb);

                /* Store requested quantities */
                int idx = ib*Nv*Neval + Nv * k + iv;
                if(mu0 != NULL)    { mu0[idx]    = mufun[0];   }
                if(mu1 != NULL)    { mu1[idx]    = mufun[1];   }
                if(dmu0 != NULL)   { dmu0[idx]   = mufun[2];   }
                if(clog != NULL)   { clog[idx]   = clogab[ib]; }
                if(F != NULL)      { F[idx]      = Fb;         }
                if(Dpara != NULL)  { Dpara[idx]  = Dparab;     }
                if(Dperp != NULL)  { Dperp[idx]  = Dperpb;     }
                if(K != NULL)      { K[idx]      = Kb;         }
                if(nu != NULL)     { nu[idx]     = nub;        }
                if(Q != NULL)      { Q[idx]      = Qb;         }
                if(dQ != NULL)     { dQ[idx]     = dQb;        }
                if(dDpara != NULL) { dDpara[idx] = dDparab;    }
            }
        }
    }
}

/**
 * @brief Evaluate cross-section for atomic reactions.
 *
 * @param sim_offload_data initialized simulation offload data struct
 * @param B_offload_array initialized magnetic field offload data
 * @param plasma_offload_array initialized magnetic field offload data
 * @param neutral_offload_array initialized magnetic field offload data
 * @param asigma_offload_array initialized magnetic field offload data
 * @param Neval number of evaluation points in (R, phi, z, t).
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param Nv number of evaluation points in velocity space.
 * @param va test particle velocities [m/s].
 * @param Aa test particle atomic mass number.
 * @param Za test particle charge number.
 * @param ma test particle mass.
 * @param ib target ion index nuber in plasma input.
 * @param reac_type reaction type
 * @param sigmav output array where evaluated values are stored [1/m^2].
 */
void libascot_eval_sigmav(
    sim_offload_data* sim_offload_data, real* B_offload_array,
    real* plasma_offload_array, real* neutral_offload_array,
    real* asigma_offload_array,
    int Neval, real* R, real* phi, real* z, real* t, int Nv, real* va,
    int Aa, int Za, real ma, int ib, int reac_type, real* sigmav) {

    sim_data sim;
    B_field_init(&sim.B_data, &sim_offload_data->B_offload_data,
                 B_offload_array);
    plasma_init(&sim.plasma_data, &sim_offload_data->plasma_offload_data,
                plasma_offload_array);
    neutral_init(&sim.neutral_data, &sim_offload_data->neutral_offload_data,
                 neutral_offload_array);
    asigma_init(&sim.asigma_data, &sim_offload_data->asigma_offload_data,
                asigma_offload_array);

    const int* Zb  = plasma_get_species_znum(&sim.plasma_data);
    const int* Ab  = plasma_get_species_anum(&sim.plasma_data);

    #pragma omp parallel for
    for (int k=0; k < Neval; k++) {
        real psi[1], rho[2], T0[1], n[MAX_SPECIES], T[MAX_SPECIES];
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
        if( neutral_eval_t0(T0, rho[0], R[k], phi[k], z[k], t[k],
                            &sim.neutral_data) ) {
            continue;
        }
        for (int j=0; j < Nv; j++) {
            real E = (physlib_gamma_vnorm(va[j]) - 1.0) * ma * CONST_C*CONST_C;
            real val;
            if( asigma_eval_sigmav(
                    &val, Za, Aa, ma, Zb[ib], Ab[ib],
                    E, T[0], T0[0], n[ib+1], reac_type, &sim.asigma_data) ) {
                continue;
            }
            sigmav[Nv*k + j] = val;
        }
    }

}
