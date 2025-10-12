/**
 * @file step_gc_rk4.c
 * @brief Guiding center integration with RK4
 **/
#include <stdio.h>
#include <math.h>
#include "defines.h"
#include "consts.h"
#include "data/bfield.h"
#include "data/efield.h"
#include "data/boozer.h"
#include "data/mhd.h"
#include "utils/mathlib.h"
#include "data/marker.h"
#include "utils/physlib.h"
#include "orbit_following.h"

/**
 * @brief Integrate a guiding center step for a struct of markers with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with RK4 simultaneously using SIMD instructions. All arrays in the
 * function are of NSIMD length so vectorization can be performed directly
 * without gather and scatter operations.
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_rk4(MarkerGuidingCenter* p, real* h, Bfield* bfield,
                 Efield* efield, int aldforce) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd aligned(h : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            err_t errflag = 0;

            real k1[6], k2[6], k3[6], k4[6];
            real tempy[6];
            real yprev[6];
            real y[6];

            real mass   = p->mass[i];
            real charge = p->charge[i];

            real B_dB[15];
            real E[3];

            real R0   = p->r[i];
            real z0   = p->z[i];
            real t0   = p->time[i];

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->ppar[i];
            yprev[4] = p->mu[i];
            yprev[5] = p->zeta[i];

            /* Magnetic field at initial position already known */
            B_dB[0]  = p->B_r[i];
            B_dB[1]  = p->B_r_dr[i];
            B_dB[2]  = p->B_r_dphi[i];
            B_dB[3]  = p->B_r_dz[i];

            B_dB[4]  = p->B_phi[i];
            B_dB[5]  = p->B_phi_dr[i];
            B_dB[6]  = p->B_phi_dphi[i];
            B_dB[7]  = p->B_phi_dz[i];

            B_dB[8]  = p->B_z[i];
            B_dB[9]  = p->B_z_dr[i];
            B_dB[10] = p->B_z_dphi[i];
            B_dB[11] = p->B_z_dz[i];

            if(!errflag) {
                errflag = Efield_eval_e(E, yprev[0], yprev[1], yprev[2],
                                         t0, efield, bfield);
            }
            if(!errflag) {
                step_gceom(k1, yprev, mass, charge, B_dB, E, aldforce);
            }


            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k1[j]/2.0;
            }


            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, efield, bfield);
            }
            if(!errflag) {
                step_gceom(k2, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k2[j]/2.0;
            }


            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, efield, bfield);
            }
            if(!errflag) {
                step_gceom(k3, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k3[j];
            }


            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], efield, bfield);}
            if(!errflag) {
                step_gceom(k4, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                y[j] = yprev[j]
                    + h[i]/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            }


            /* Test that results are physical */
            if(!errflag && y[0] <= 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }
            if(!errflag && y[4] < 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }

            /* Update gc phase space position */
            if(!errflag) {
                p->r[i]    = y[0];
                p->phi[i]  = y[1];
                p->z[i]    = y[2];
                p->ppar[i] = y[3];
                p->mu[i]   = y[4];
                p->zeta[i] = fmod(y[5],CONST_2PI);
                if(p->zeta[i]<0) {
                    p->zeta[i] = CONST_2PI + p->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, p->r[i], p->phi[i], p->z[i],
                                            t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Bfield_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }

            if(!errflag) {
                p->B_r[i]        = B_dB[0];
                p->B_r_dr[i]     = B_dB[1];
                p->B_r_dphi[i]   = B_dB[2];
                p->B_r_dz[i]     = B_dB[3];

                p->B_phi[i]      = B_dB[4];
                p->B_phi_dr[i]   = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i]   = B_dB[7];

                p->B_z[i]        = B_dB[8];
                p->B_z_dr[i]     = B_dB[9];
                p->B_z_dphi[i]   = B_dB[10];
                p->B_z_dz[i]     = B_dB[11];

                p->rho[i] = rho[0];

                /* Evaluate theta angle so that it is cumulative */
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}

/**
 * @brief Integrate a guiding center step with RK4 with MHD modes present.
 *
 * Same as previous function but with MHD present
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param boozer pointer to boozer data
 * @param mhd pointer to MHD data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_rk4_mhd(
    MarkerGuidingCenter* p, real* h, Bfield* bfield, Efield* efield,
    Boozer* boozer, Mhd* mhd, int aldforce) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd aligned(h : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            err_t errflag = 0;

            real k1[6], k2[6], k3[6], k4[6];
            real tempy[6];
            real yprev[6];
            real y[6];

            real mass   = p->mass[i];
            real charge = p->charge[i];

            real B_dB[15];
            real E[3];
            real mhd_dmhd[10];

            real R0 = p->r[i];
            real z0 = p->z[i];
            real t0 = p->time[i];

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->ppar[i];
            yprev[4] = p->mu[i];
            yprev[5] = p->zeta[i];

            /* Magnetic field at initial position already known */
            B_dB[0]  = p->B_r[i];
            B_dB[1]  = p->B_r_dr[i];
            B_dB[2]  = p->B_r_dphi[i];
            B_dB[3]  = p->B_r_dz[i];

            B_dB[4]  = p->B_phi[i];
            B_dB[5]  = p->B_phi_dr[i];
            B_dB[6]  = p->B_phi_dphi[i];
            B_dB[7]  = p->B_phi_dz[i];

            B_dB[8]  = p->B_z[i];
            B_dB[9]  = p->B_z_dr[i];
            B_dB[10] = p->B_z_dphi[i];
            B_dB[11] = p->B_z_dz[i];

            if(!errflag) {
                errflag = Efield_eval_e(E, yprev[0], yprev[1], yprev[2], t0,
                                         efield, bfield);
            }
            if(!errflag) {
                errflag = Mhd_eval(mhd_dmhd, yprev[0], yprev[1], yprev[2], t0,
                                   MHD_INCLUDE_ALL, boozer, mhd, bfield);
            }
            if(!errflag) {
                step_gceom_mhd(k1, yprev, mass, charge, B_dB, E, mhd_dmhd,
                               aldforce);
            }

            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]/2.0*k1[j];
            }

            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, efield, bfield);
            }
            if(!errflag) {
                errflag = Mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + h[i]/2.0, MHD_INCLUDE_ALL, boozer,
                                   mhd, bfield);
            }
            if(!errflag) {
                step_gceom_mhd(k2, tempy, mass, charge, B_dB, E, mhd_dmhd,
                               aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]/2.0*k2[j];
            }

            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, efield, bfield);
            }
            if(!errflag) {
                errflag = Mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + h[i]/2.0, MHD_INCLUDE_ALL, boozer,
                                   mhd, bfield);
            }
            if(!errflag) {
                step_gceom_mhd(k3, tempy, mass, charge, B_dB, E, mhd_dmhd,
                               aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k3[j];
            }

            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Efield_eval_e(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], efield, bfield);
            }
            if(!errflag) {
                errflag = Mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + h[i], MHD_INCLUDE_ALL, boozer,
                                   mhd, bfield);
            }
            if(!errflag) {
                step_gceom_mhd(k4, tempy, mass, charge, B_dB, E, mhd_dmhd,
                               aldforce);
            }
            for(int j = 0; j < 6; j++) {
                y[j] = yprev[j]
                    + h[i]/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            }

            /* Test that results are physical */
            if(!errflag && y[0] <= 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }
            else if(!errflag && fabs(y[4]) >= CONST_C) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }
            else if(!errflag && y[4] < 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }

            /* Update gc phase space position */
            if(!errflag) {
                p->r[i]    = y[0];
                p->phi[i]  = y[1];
                p->z[i]    = y[2];
                p->ppar[i] = y[3];
                p->mu[i]   = y[4];
                p->zeta[i] = fmod(y[5],CONST_2PI);
                if(p->zeta[i]<0){
                    p->zeta[i] = CONST_2PI + p->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = Bfield_eval_b_db(B_dB, p->r[i], p->phi[i], p->z[i],
                                            t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Bfield_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           t0 + h[i], bfield);
            }
            if(!errflag) {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }

            if(!errflag) {
                p->B_r[i]        = B_dB[0];
                p->B_r_dr[i]     = B_dB[1];
                p->B_r_dphi[i]   = B_dB[2];
                p->B_r_dz[i]     = B_dB[3];

                p->B_phi[i]      = B_dB[4];
                p->B_phi_dr[i]   = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i]   = B_dB[7];

                p->B_z[i]        = B_dB[8];
                p->B_z_dr[i]     = B_dB[9];
                p->B_z_dphi[i]   = B_dB[10];
                p->B_z_dz[i]     = B_dB[11];

                p->rho[i] = rho[0];

                /* Evaluate pol angle so that it is cumulative */
                real axisrz[2];
                errflag  = Bfield_eval_axis_rz(axisrz, bfield, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }

}
