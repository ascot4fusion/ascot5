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


#define STORE(index, val, ptr) if ((ptr)) (ptr)[(index)] = (val)
/** Store value in output array if the output array is allocated. */

/**
 * @brief Interpolate input quantities at the given coordinates.
 *
 * @param[in] bfield Magnetic field data.
 * @param[in] efield Electric field data.
 * @param[in] plasma plasma data.
 * @param[in] neutral Neutral data.
 * @param[in] boozer Boozer data.
 * @param[in] mhd MHD data.
 * @param[in] npnt Number of evaluation points.
 * @param[in] modenumber Evaluate mhd perturbation only for this mode.
 *
 *        If modenumber < 0, evaluate all modes.
 * @param[in] R R coordinates where the inputs are interpolated [m].
 * @param[in] phi phi coordinates where the inputs are interpolated [rad].
 * @param[in] z z coordinates where the inputs are interpolated [m].
 * @param[in] t time coordinates where the inputs are interpolated [s].
 * @param[out] B Magnetic field vector (Br, Bphi, Bz) [T].
 * @param[out] Bjac Magnetic field Jacobian (dBr/dr, dBr/dphi, dBr/dz, dBphi/dr,
 *             dBphi/dphi, dBphi/dz, dBz/dr, dBz/dphi, dBz/dz) [T/m].
 * @param[out] psi Poloidal flux and its derivatives (psi, dpsi/dr, dpsi/dphi,
 *             dpsi/dz) [Wb/rad, Wb/rad m].
 * @param[out] rho Square root of the normalized poloidal flux and its
 *             derivative (rho, drho/dpsi) [1, rad/Wb].
 * @param[out] E Electric field vector (Er, Ephi, Ez) [V/m].
 * @param[out] n Density of plasma species (ne, ni1, ni2, ...) [m^-3].
 * @param[out] T Temperature of plasma species (Te, Ti) [eV].
 * @param[out] n0 Density of neutral species (n1, n2, ...) [m^-3].
 * @param[out] T0 Temperature of neutral species (T1, T2, ...) [eV].
 * @param[out] theta Poloidal Boozer coordinate and its derivatives (theta,
 *             dtheta/dr, dtheta/dphi, dtheta/dz) [rad, rad/m].
 * @param[out] zeta output array [rad].
 * @param[out] alpha Magnetic perturbation eigenfunction and its derivatives
 *             (Phi, dPhi/dr, dPhi/dphi, dPhi/dz) [1]
 * @param[out] Phi Electric perturbation eigenfunction and its derivatives
 *             (Phi, dPhi/dr, dPhi/dphi, dPhi/dz) [1].
 * @param[out] mhd_br Magnetic field perturbation components due to MHD
 *             (Br, Bphi, Bz) [T].
 * @param[out] mhd_er Electric field perturbation components due to MHD
 *             (Er, Ephi, Ez) [V/m].
 * @param[out] mhd_phi Electric field perturbation potential [V/m].
 */
void libascot_interpolate(
    B_field_data* bfield, E_field_data* efield, plasma_data* plasma,
    neutral_data* neutral, boozer_data* boozer, mhd_data* mhd,
    //asigma_data* atomic,
    int npnt, int modenumber,
    real R[npnt], real phi[npnt], real z[npnt], real t[npnt],
    real B[3][npnt], real Bjac[9][npnt], real psi[4][npnt], real rho[2][npnt],
    real E[3][npnt],
    real n[][npnt], real T[2][npnt], real n0[][npnt], real T0[][npnt],
    real theta[4][npnt], real zeta[4][npnt],
    real alpha[5][npnt], real Phi[5][npnt],
    real mhd_b[3][npnt], real mhd_e[3][npnt], real mhd_phi[npnt]
    ) {

    int ONLY_PERTURBATIONS = 1;
    #pragma omp parallel for
    for(int k = 0; k < npnt; k++) {
        real Bq[15], psival[4], rhoval[2], Eq[3], ns[MAX_SPECIES],
            Ts[MAX_SPECIES], psithetazeta[12], pert_field[7], mhd_dmhd[10];
        bool psi_valid = false, rho_valid = false;
        int n_species, isinside;

        if( bfield &&
            !B_field_eval_B_dB(Bq, R[k], phi[k], z[k], t[k], bfield) ) {
            STORE(0*npnt + k, Bq[0], *B);
            STORE(1*npnt + k, Bq[4], *B);
            STORE(2*npnt + k, Bq[8], *B);
            STORE(0*npnt + k, Bq[1], *Bjac);
            STORE(1*npnt + k, Bq[2], *Bjac);
            STORE(2*npnt + k, Bq[3], *Bjac);
            STORE(3*npnt + k, Bq[5], *Bjac);
            STORE(4*npnt + k, Bq[6], *Bjac);
            STORE(5*npnt + k, Bq[7], *Bjac);
            STORE(6*npnt + k, Bq[9], *Bjac);
            STORE(7*npnt + k, Bq[10], *Bjac);
            STORE(8*npnt + k, Bq[11], *Bjac);
        }
        if( bfield &&
            !B_field_eval_psi_dpsi(psival, R[k], phi[k], z[k], t[k], bfield) ) {
            psi_valid = true;
            STORE(0*npnt + k, psival[0], *psi);
            STORE(1*npnt + k, psival[1], *psi);
            STORE(2*npnt + k, psival[2], *psi);
            STORE(3*npnt + k, psival[3], *psi);
        }
        if( bfield && psi_valid &&
            !B_field_eval_rho(rhoval, psival[0], bfield) ) {
            rho_valid = true;
            STORE(0*npnt + k, rhoval[0], *rho);
            STORE(1*npnt + k, rhoval[1], *rho);
        }
        if( efield && bfield &&
            !E_field_eval_E(Eq, R[k], phi[k], z[k], t[k], efield, bfield) ) {
            STORE(0*npnt + k, Eq[0], *E);
            STORE(1*npnt + k, Eq[1], *E);
            STORE(2*npnt + k, Eq[2], *E);
        }
        if(plasma) {
            n_species = plasma_get_n_species(plasma);
        }
        if( plasma && rho_valid &&
            !plasma_eval_densandtemp(ns, Ts, rhoval[0], R[k], phi[k], z[k], t[k],
                                     plasma) ) {
            STORE(0*npnt + k, Ts[0] / CONST_E, *T);
            STORE(1*npnt + k, Ts[1] / CONST_E, *T);
            for(int i=0; i < n_species; i++) {
                STORE(i*npnt + k, ns[i], *n);
            }
        }
        if(neutral) {
            n_species = neutral_get_n_species(neutral);
        }
        if( neutral && rho_valid &&
            !neutral_eval_n0(ns, rhoval[0], R[k], phi[k], z[k], t[k], neutral) ) {
            for(int i=0; i < n_species; i++) {
                STORE(i*npnt + k, ns[i], *n0);
            }
        }
        if( neutral && rho_valid &&
            !neutral_eval_t0(Ts, rhoval[0], R[k], phi[k], z[k], t[k], neutral) ) {
            for(int i=0; i < n_species; i++) {
                STORE(i*npnt + k, Ts[i], *T0);
            }
        }
        if( boozer && bfield &&
            !boozer_eval_psithetazeta(psithetazeta, &isinside, R[k], phi[k],
                                      z[k], bfield, boozer) ) {
            if(isinside) {
                STORE(0*npnt + k, psithetazeta[4], *theta);
                STORE(0*npnt + k, psithetazeta[8], *zeta);
                STORE(1*npnt + k, psithetazeta[5], *theta);
                STORE(2*npnt + k, psithetazeta[6], *theta);
                STORE(3*npnt + k, psithetazeta[7], *theta);
                STORE(1*npnt + k, psithetazeta[9], *zeta);
                STORE(2*npnt + k, psithetazeta[10], *zeta);
                STORE(3*npnt + k, psithetazeta[11], *zeta);
            }
        }
        if( mhd && boozer && bfield &&
            !mhd_eval(mhd_dmhd, R[k], phi[k], z[k], t[k], modenumber,
                      boozer, mhd, bfield) ) {
            STORE(0*npnt + k, mhd_dmhd[0], *alpha);
            STORE(1*npnt + k, mhd_dmhd[2], *alpha);
            STORE(2*npnt + k, mhd_dmhd[3], *alpha);
            STORE(3*npnt + k, mhd_dmhd[4], *alpha);
            STORE(4*npnt + k, mhd_dmhd[1], *alpha);
            STORE(0*npnt + k, mhd_dmhd[5], *Phi);
            STORE(1*npnt + k, mhd_dmhd[7], *Phi);
            STORE(2*npnt + k, mhd_dmhd[8], *Phi);
            STORE(3*npnt + k, mhd_dmhd[9], *Phi);
            STORE(4*npnt + k, mhd_dmhd[6], *Phi);
        }
        if( mhd && boozer && bfield &&
            !mhd_perturbations(
                pert_field, R[k], phi[k], z[k], t[k], ONLY_PERTURBATIONS,
                modenumber, boozer, mhd, bfield
            ) ) {
            STORE(0*npnt + k, pert_field[0], *mhd_b);
            STORE(1*npnt + k, pert_field[1], *mhd_b);
            STORE(2*npnt + k, pert_field[2], *mhd_b);
            STORE(0*npnt + k, pert_field[3], *mhd_e);
            STORE(1*npnt + k, pert_field[4], *mhd_e);
            STORE(2*npnt + k, pert_field[5], *mhd_e);
            STORE(0*npnt + k, pert_field[6], mhd_phi);
        }
    }
}

/**
 * @brief Get magnetic axis at given coordinates.
 *
 * @param bfield Magnetic field data
 * @param nphi Number of evaluation points.
 * @param phi phi coordinates of the evaluation points [rad].
 * @param Raxis output array for axis R coordinates.
 */
void libascot_find_axis(
    B_field_data* bfield, int nphi, real phi[nphi], real axisRz[2][nphi]) {

    #pragma omp parallel for
    for(int k = 0; k < nphi; k++) {
        real axisrz[2];
        if( !B_field_get_axis_rz(axisrz, bfield, phi[k]) ) {
            axisRz[0][k] = axisrz[0];
            axisRz[1][k] = axisrz[1];
        }
    }
}

/**
 * @brief Map (rho, theta, phi) to (R,z) coordinates.
 *
 * This function implements the Newton method. If the function fails on
 * a given position, the corresponding (R,z) values in the output arrays are
 * not altered.
 *
 * @param bfield magnetic field data
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
void libascot_rhotheta2rz(
    B_field_data* bfield, int Neval, int maxiter, real tol, real t,
    real* rho, real* theta, real* phi, real* r, real* z) {

    #pragma omp parallel for
    for(int j=0; j<Neval; j++) {
        real axisrz[2];
        real rhodrho[4];
        if( B_field_get_axis_rz(axisrz, bfield, phi[j]) ) {
            continue;
        }
        if( B_field_eval_rho_drho(rhodrho, axisrz[0], phi[j], axisrz[1],
                                  bfield)) {
            continue;
        }
        if( rhodrho[0] > rho[j] ) {
            /* Due to padding, rho might not be exactly zero on the axis so we
             * return the axis position for small values of queried rho */
            r[j] = axisrz[0];
            z[j] = axisrz[1];
            continue;
        }

        real a = 0.0, b = 5.0;
        real costh = cos(theta[j]);
        real sinth = sin(theta[j]);
        for(int i=0; i<maxiter; i++) {
	    real c = 0.5*(a + b);
	    real rj = axisrz[0] + c * costh;
            real zj = axisrz[1] + c * sinth;
	    if(rj < 0) {
	        b = c;
		continue;
	    }
	    if( B_field_eval_rho_drho(rhodrho, rj, phi[j], zj, bfield) ) {
	        b = c;
                continue;
            }
            if( fabs(rho[j] - rhodrho[0]) < tol ) {
                r[j] = rj;
                z[j] = zj;
                break;
            }
	    if( rho[j] < rhodrho[0]) {
	        b = c;
	    } else {
	        a = c;
	    }
        }
    }
}

/**
 * @brief Find psi on axis using the gradient descent method
 *
 * Note that the psi value is not returned in case this algorithm fails.
 *
 * @param bfield magnetic field data
 * @param psi value of psi on axis if this function did not fail
 * @param rz initial (R,z) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void libascot_gradient_descent(
    B_field_data* bfield, int maxiter, int ascent, real step, real tol,
    real psi[1], real rz[2]) {

    if(ascent) {
        step = -1 * step;
    }

    real phi = 0.0, time = 0.0;
    real psidpsi[4], nextrz[2];
    B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time, bfield);

    int iter = 0;
    while(1) {
        if( B_field_eval_psi_dpsi(psidpsi, rz[0], phi, rz[1], time,
                                  bfield) ) {
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
                psidpsi, rz[0], phi, rz[1], time, bfield);
            psi[0] = psi[0] + (tol * psidpsi[1] + tol * psidpsi[3]);
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
 * @brief Find one psi minimum using the gradient descent method for 3D field
 * inside a sector (phimin < phi < phimax). Not guaranteed to find the global
 * minimum inside the given sector.
 *
 * @param bfield magnetic field data
 * @param psi value of psi on axis if this function did not fail
 * @param rzphi initial (R,z,phi) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void libascot_gradient_descent_3d(
    B_field_data* bfield, int maxiter, int ascent, real phimin, real phimax,
    real step, real tol, real psi[1], real rzphi[3]) {

    if(ascent) {
        step = -1 * step;
    }

    real time = 0.0;
    real psidpsi[4], nextrzphi[3];
    B_field_eval_psi_dpsi(psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield);

    int iter = 0;
    while(1) {
        if( B_field_eval_psi_dpsi(psidpsi, rzphi[0], rzphi[2], rzphi[1], time,
                                  bfield) ) {
            break;
        }
        nextrzphi[0] = rzphi[0] - step * psidpsi[1];          // R
        nextrzphi[1] = rzphi[1] - step * psidpsi[3];          // z
        nextrzphi[2] = rzphi[2] - step/rzphi[0] * psidpsi[2]; /* phi. phidpsi[2]
                                                               is dimensionless,
                                                               must divide by R
                                                               because in
                                                               cylindrical
                                                               co-ordinates   */

        /* Check that phi remained inside the sector. If not, use the value on
           the sector boundary. */
        if (nextrzphi[2] > phimax) {nextrzphi[2] = phimax;}
        if (nextrzphi[2] < phimin) {nextrzphi[2]=phimin;}

        /* Check convergence (phi difference must be multiplied by R to get
        the arc length which has dimensions of L) */
        if(sqrt( (nextrzphi[0] - rzphi[0]) * (nextrzphi[0] - rzphi[0])
                + (nextrzphi[1] - rzphi[1]) * (nextrzphi[1] - rzphi[1])
                + rzphi[0]*(nextrzphi[2] - rzphi[2]) *
                rzphi[0]*(nextrzphi[2] - rzphi[2])) < tol){
            psi[0] = psidpsi[0];
            rzphi[0] = nextrzphi[0];
            rzphi[1] = nextrzphi[1];
            rzphi[2] = nextrzphi[2];

            // Add a bit of padding
            B_field_eval_psi_dpsi(
                psidpsi, rzphi[0], rzphi[2], rzphi[1], time, bfield);
            psi[0] = psi[0]
                + (tol * ( psidpsi[1] + psidpsi[2]/rzphi[0] + psidpsi[3] ));
            break;
        }

        rzphi[0] = nextrzphi[0];
        rzphi[1] = nextrzphi[1];
        rzphi[2] = nextrzphi[2];
        iter++;

        if(iter == maxiter) {
            break;
        }
    }
}

/**
 * @brief Evaluate collision coefficients
 *
 * @param sim_data initialized simulation data struct
 * @param Neval number of evaluation points
 * @param R R coordinates of the evaluation points [m]
 * @param phi phi coordinates of the evaluation points [rad]
 * @param z z coordinates of the evaluation points [m]
 * @param t time coordinates of the evaluation points [s]
 * @param Nv number of velocity grid points
 * @param va test particle velocities at the evaluation point [m/s]
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param coefficients The collision coefficients (F, Dpara, Dperp, K, nu, Q,
 *        dQ, dDpara, clog) [].
 * @param clog The Coulomb logarithm.
 * @param mu The special functions (mu0, mu1, dmu0)
 */
void libascot_eval_collcoefs(
    B_field_data* bfield, plasma_data* plasma, int npnt, real ma, real qa,
    real R[npnt], real phi[npnt], real z[npnt], real t[npnt], real va[npnt],
    real coefficients[6][npnt], real clog[][npnt], real mu[3][npnt]) {

    /* Evaluate plasma parameters */
    int n_species  = plasma_get_n_species(plasma);
    const real* qb = plasma_get_species_charge(plasma);
    const real* mb = plasma_get_species_mass(plasma);
    mccc_data mccc;

    #pragma omp parallel for
    for(int k=0; k<npnt; k++) {
        real mufun[3] = {0., 0., 0.};

        /* Evaluate rho as it is needed to evaluate plasma parameters */
        real psi, rho[2];
        if( B_field_eval_psi(&psi, R[k], phi[k], z[k], t[k], bfield) ) {
            continue;
        }
        if( B_field_eval_rho(rho, psi, bfield) ) {
            continue;
        }

        real nb[MAX_SPECIES], Tb[MAX_SPECIES];
        if( plasma_eval_densandtemp(nb, Tb, rho[0], R[k], phi[k], z[k], t[k],
                                    plasma) ) {
            continue;
        }

        /* Loop through all plasma species */
        for(int ib=0; ib<n_species; ib++) {

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, ma, qa, va[k], n_species, mb, qb, nb, Tb);

            /* Special functions */
            real vb = sqrt( 2 * Tb[ib] / mb[ib] );
            real x  = va[k] / vb;
            mccc_coefs_mufun(mufun, x, &mccc);

            /* Coefficients */
            real Fb = mccc_coefs_F(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real Qb = mccc_coefs_Q(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real dQb = mccc_coefs_dQ(
                ma, qa, mb[ib], qb[ib], nb[ib], vb, clogab[ib], mufun[2]);
            real Dparab = mccc_coefs_Dpara(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[0]);
            real Dperpb = mccc_coefs_Dperp(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[1]);
            real dDparab = mccc_coefs_dDpara(
                ma, qa, va[k], qb[ib], nb[ib], vb, clogab[ib], mufun[0],
                mufun[2]);
            real Kb = mccc_coefs_K(va[k], Dparab, dDparab, Qb);
            real nub = mccc_coefs_nu(va[k], Dperpb);

            int idx = ib*npnt + k;
            //STORE(n_species*npnt*0 + ib*npnt + k, mufun[0], **mu);
            /**
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
            if(dDpara != NULL) { dDpara[idx] = dDparab;    } */
        }
    }
}

/**
 * @brief Evaluate atomic reaction rate coefficient.
 *
 * @param sim_data initialized simulation data struct
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
 * @param reac_type reaction type
 * @param ratecoeff output array where evaluated values are stored [1/m^2].
 */
void libascot_eval_ratecoeff(
    B_field_data* bfield, plasma_data* plasma, neutral_data* neutral,
    asigma_data* atomic,
    int Neval, real* R, real* phi, real* z, real* t, int Nv, real* va,
    int Aa, int Za, real ma, int reac_type, real* ratecoeff) {

    const int* Zb = plasma_get_species_znum(plasma);
    const int* Ab = plasma_get_species_anum(plasma);
    int nion  = plasma_get_n_species(plasma) - 1;
    int nspec = neutral_get_n_species(neutral);

    #pragma omp parallel for
    for (int k=0; k < Neval; k++) {
        real psi[1], rho[2], T0[1], n[MAX_SPECIES], T[MAX_SPECIES],
            n0[MAX_SPECIES];
        if( B_field_eval_psi(psi, R[k], phi[k], z[k], t[k], bfield) ) {
            continue;
        }
        if( B_field_eval_rho(rho, psi[0], bfield) ) {
            continue;
        }
        if( plasma_eval_densandtemp(n, T, rho[0], R[k], phi[k], z[k], t[k],
                                    plasma) ) {
            continue;
        }
        if( neutral_eval_t0(T0, rho[0], R[k], phi[k], z[k], t[k], neutral) ) {
            continue;
        }
        if( neutral_eval_n0(n0, rho[0], R[k], phi[k], z[k], t[k], neutral) ) {
            continue;
        }
        for (int j=0; j < Nv; j++) {
            real E = (physlib_gamma_vnorm(va[j]) - 1.0) * ma * CONST_C*CONST_C;
            real val;
            switch (reac_type) {
            case sigmav_CX:
                if( asigma_eval_cx(
                        &val, Za, Aa, E, ma, nspec, Zb, Ab, T0[0], n0, atomic) ) {
                    continue;
                }
                ratecoeff[Nv*k + j] = val;
                break;
            case sigmav_BMS:
                if( asigma_eval_bms(
                        &val, Za, Aa, E, ma, nion, Zb, Ab, T[0], n, atomic) ) {
                    continue;
                }
                ratecoeff[Nv*k + j] = val * n[0];
                break;
            default:
                break;
            }
        }
    }

}

/**
 * @brief Evaluate ICRH electric field and the resonance condition.
 *
 * The evaluated electric field consists of left-hand (+) and right-hand (-)
 * circularly polarized components. The circularly polarised components are
 * notorious of their negatice effects on the mental health of their consumers.
 * In the context of this function, they have the following definitions:
 * E_plus_real = Re{E_LH ( cos( phase(E_LH) ) + i sin( phase(E_LH) ) } and
 * E_minus_real = Re{E_RH ( cos( phase(E_RH) ) - i sin( phase(E_RH) ) }.
 * E_LH, E_RH and their phases are all real. Furthermore, E_LH is called E_+ and
 * R_RH is called E_- in the context of RFOF and also more generally. Sometimes
 * the definition E+- = E_x +- iE_y. It is seen that |E+| = sqrt(E+  E+*) where
 * E+* is the complex conjugate. But because E+*=E-, one has
 * |E+| = sqrt(E+  E-) = |E-|. This is obviously not the case in this function,
 * as the magnitudes are different but this definition is also used sometimes.
 *
 * The resonance condition is given by
 *
 * omega_wave - n * omega_gyro - k_parallel * v_parallel - k_perp dot v_drift
 * = 0. Whether the Doppler shift (k_par v_par) or/and the drift term
 * (k_per v_drif) is used, is specified in the rfof_codeparam.xml. The drift
 * has not yet been implemented (Jan 2025).
 *
 * @param bfield magnetic field data.
 * @param rfof RFOF data.
 * @param Neval number of evaluation points.
 * @param R R coordinates of the evaluation points [m].
 * @param phi phi coordinates of the evaluation points [rad].
 * @param z z coordinates of the evaluation points [m].
 * @param t time coordinates of the evaluation points [s].
 * @param mass test particle mass (for computing resonance) [kg].
 * @param q test particle charge (for computing resonance) [C].
 * @param vpar test particle parallel velocity (for computing resonance) [m/s].
 * @param Eplus_real Real part of the right-handed electric field component of
 *        the wave [V/m]. (See comment above!)
 * @param Eminus_real Real part of the left-handed electric field component of
 *        the wave [V/m]. (See comment above!)
 * @param Eplus_imag Imaginary part of the right-handed electric field component
 *        of the wave [V/m]. (See comment above!)
 * @param Eminus_imag Imaginary part of the left-handed electric field
 *        component of the wave [V/m]. (See comment above!)
 * @param res_cond value of the resonance condition where zero is the resonance
 *        [1].
 */
void libascot_eval_rfof(
    B_field_data* bfield, rfof_data* rfof, int Neval,
    real* R, real* phi, real* z, real* t, real mass, real q, real vpar,
    real* Eplus_real, real* Eminus_real, real* Eplus_imag, real* Eminus_imag,
    real* res_cond) {

    #pragma omp parallel
    {
        /* The function that evaluates resonance condition takes an RFOF marker
        * as an input. However, only the R and vpar values are actually used.
        * Therefore, we initialize a dummy marker and adjust only the values of
        * R and vpar. */
        rfof_marker rfof_mrk;
        int dummy_int   = 1;
        real dummy_real = -999.0;  /*-999.0 to be on the safe side */
        rfof_set_up(&rfof_mrk, rfof);

        #pragma omp for
        for(int k = 0; k < Neval; k++) {
            real B[3];
            if( B_field_eval_B(B, R[k], phi[k], z[k], t[k], bfield) ) {
                continue;
            }
            real B_magn = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
            real gyrofreq = q * B_magn / mass;
            rfof_set_marker_manually(&rfof_mrk, &dummy_int,
                &dummy_real, &(R[k]), &dummy_real, &dummy_real, &dummy_real,
                &dummy_real, &dummy_real, &dummy_real, &dummy_real, &dummy_real,
                &dummy_real, &vpar, &dummy_real, &gyrofreq, &dummy_real,
                &dummy_real, &dummy_int, &dummy_int);

            int nharm; /* For storing return value which is not used */
            rfof_eval_resonance_function(
                &(res_cond[k]), &nharm, &rfof_mrk, rfof);

            // TODO: this should return a non-zero value for failed evaluations
            rfof_eval_rf_wave(
                &(Eplus_real[k]), &(Eminus_real[k]), &(Eplus_imag[k]),
                &(Eminus_imag[k]), R[k], z[k], rfof);
        }
        rfof_tear_down(&rfof_mrk);
    }
}
