/**
 * @file step_gc_cashkarp.c
 * @brief Guiding center integration with adaptive Cash Karp method
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../particle.h"
#include "../../error.h"
#include "step_gc_cashkarp.h"
#include "step_gceom.h"
#include "step_gceom_mhd.h"

/**
 * @brief Integrate a guiding center step for a struct of markers
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions.
 * All arrays in the function are of NSIMD length so vectorization can be performed
 * directly without gather and scatter operations. Informs whther time step was accepted or
 * rejected and provides a suggestion for the next time step.
 *
 * @param p marker struct that will be updated
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step. Negative sign indicates current step was rejected
 * @param tol error tolerance
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_cashkarp(particle_simd_gc* p, real* h, real* hnext, real tol,
                      B_field_data* Bdata, E_field_data* Edata, int aldforce) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd aligned(h, hnext : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real k1[6], k2[6], k3[6], k4[6], k5[6], k6[6];
            real tempy[6];
            real yprev[6];

            real mass   = p->mass[i];;
            real charge = p->charge[i];

            real B_dB[15];
            real E[3];

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
            B_dB[0] = p->B_r[i];
            B_dB[1] = p->B_r_dr[i];
            B_dB[2] = p->B_r_dphi[i];
            B_dB[3] = p->B_r_dz[i];

            B_dB[4] = p->B_phi[i];
            B_dB[5] = p->B_phi_dr[i];
            B_dB[6] = p->B_phi_dphi[i];
            B_dB[7] = p->B_phi_dz[i];

            B_dB[8] = p->B_z[i];
            B_dB[9] = p->B_z_dr[i];
            B_dB[10] = p->B_z_dphi[i];
            B_dB[11] = p->B_z_dz[i];

            if(!errflag) {
                errflag = E_field_eval_E(E, yprev[0], yprev[1], yprev[2],
                                         t0, Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k1, yprev, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                        (1.0/5) * k1[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (1.0/5)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (1.0/5)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k2, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          (3.0/40) * k1[j]
                        + (9.0/40) * k2[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (3.0/10)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/10)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k3, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          ( 3.0/10) * k1[j]
                        + (-9.0/10) * k2[j]
                        + ( 6.0/5 ) * k3[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (3.0/5)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/5)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k4, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          (-11.0/54) * k1[j]
                        + (  5.0/2 ) * k2[j]
                        + (-70.0/27) * k3[j]
                        + ( 35.0/27) * k4[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k5, tempy, mass, charge, B_dB, E, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          ( 1631.0/55296 ) * k1[j]
                        + (  175.0/512   ) * k2[j]
                        + (  575.0/13824 ) * k3[j]
                        + (44275.0/110592) * k4[j]
                        + (  253.0/4096  ) * k5[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (7.0/8)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (7.0/8)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k6, tempy, mass, charge, B_dB, E, aldforce);
            }

            /* Error estimate is a difference between RK4 and RK5 solutions. If
             * time-step is accepted, the RK5 solution will be used to advance
             * marker. */
            real rk5[6], rk4[6];
            if(!errflag) {
                real err = 0.0;
                for(int j = 0; j < 6; j++) {
                    rk5[j] = yprev[j]
                        + h[i]*(
                              ( 37.0/378 ) * k1[j]
                            + (250.0/621 ) * k3[j]
                            + (125.0/594 ) * k4[j]
                            + (512.0/1771) * k6[j] );

                    rk4[j] = yprev[j] +
                        h[i]*(
                              ( 2825.0/27648) * k1[j]
                            + (18575.0/48384) * k3[j]
                            + (13525.0/55296) * k4[j]
                            + (  277.0/14336) * k5[j]
                            + (    1.0/4    ) * k6[j] );
                    if(j==3) {
                        real yerr = fabs(rk5[j] - rk4[j]);
                        real ytol = fabs(yprev[j]) + fabs(k1[j]*h[i])
                                    + DBL_EPSILON;
                        err = fmax( err, yerr/ytol );
                    }
                    else if(j==2) {
                        real rk1[3] = {k1[0]*h[i], k1[1]*h[i], k1[2]*h[i]};
                        real yerr =
                              rk5[0] * rk5[0] + rk4[0] * rk4[0]
                            - 2 * rk5[0] * rk4[0] * cos(rk5[1] - rk4[1])
                              + ( rk5[2] - rk4[2] ) * ( rk5[2] - rk4[2] );
                        real ytol =
                              yprev[0] * yprev[0] + rk1[0] * rk1[0]
                            - 2 * yprev[0] * rk1[0] * cos(yprev[1] - rk1[1])
                            + ( yprev[2] - rk1[2] ) * ( yprev[2] - rk1[2] )
                            + DBL_EPSILON;
                        err = fmax( err, sqrt(yerr/ytol) );
                    }
                }

                err = err/tol;
                if(err <= 1){
                    /* Time step accepted */
                    hnext[i] = 0.85*h[i]*pow(err,-0.2);

                    /* Make sure we don't make a huge jump */
                    if(hnext[i] > 1.5*h[i]) {
                        hnext[i] = 1.5*h[i];
                    }
                }
                else{
                    /* Time step rejected */
                    hnext[i] = -0.85*h[i]*pow(err,-0.25);
                }
            }

            /* Test that results are physical */
            if(!errflag && fabs(hnext[i]) < A5_EXTREMELY_SMALL_TIMESTEP) {
                errflag = error_raise(
                    ERR_INVALID_TIMESTEP, __LINE__, EF_STEP_GC_CASHKARP);
            }
            else if(!errflag && rk5[0] <= 0) {
                errflag = error_raise(
                    ERR_INTEGRATION, __LINE__, EF_STEP_GC_CASHKARP);
            }
            else if(!errflag && rk5[4] < 0) {
                errflag = error_raise(
                    ERR_INTEGRATION, __LINE__, EF_STEP_GC_CASHKARP);
            }

            /* Update gc phase space position */
            if(!errflag) {
                p->r[i]     = rk5[0];
                p->phi[i]   = rk5[1];
                p->z[i]     = rk5[2];
                p->ppar[i]  = rk5[3];
                p->mu[i]    = rk5[4];
                p->zeta[i]  = fmod( rk5[5], CONST_2PI );
                if(p->zeta[i]<0) {
                    p->zeta[i] = CONST_2PI + p->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i],
                                            p->time[i] + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           p->time[i] + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
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
                errflag = B_field_get_axis_rz(axisrz, Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
                hnext[i]      = h[i];
            }
        }
    }
}

/**
 * @brief Integrate a guiding center step for a struct of markers with MHD
 *
 * Rejected step has a negative suggestion for the next time-step. The negative
 * sign is only used to indicate a rejected step and absolute value should be
 * used for the next time-step.
 *
 * @param p marker struct that will be updated
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step.
 * @param tol error tolerance
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 * @param boozer pointer to Boozer data
 * @param mhd pointer to MHD data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_cashkarp_mhd(
    particle_simd_gc* p, real* h, real* hnext, real tol, B_field_data* Bdata,
    E_field_data* Edata, boozer_data* boozer, mhd_data* mhd, int aldforce) {

    int i;
    /* Following loop will be executed simultaneously for all i */
#pragma omp simd aligned(h, hnext : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real k1[6], k2[6], k3[6], k4[6], k5[6], k6[6];
            real tempy[6];
            real yprev[6];

            real mass   = p->mass[i];;
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
            B_dB[0] = p->B_r[i];
            B_dB[1] = p->B_r_dr[i];
            B_dB[2] = p->B_r_dphi[i];
            B_dB[3] = p->B_r_dz[i];

            B_dB[4] = p->B_phi[i];
            B_dB[5] = p->B_phi_dr[i];
            B_dB[6] = p->B_phi_dphi[i];
            B_dB[7] = p->B_phi_dz[i];

            B_dB[8] = p->B_z[i];
            B_dB[9] = p->B_z_dr[i];
            B_dB[10] = p->B_z_dphi[i];
            B_dB[11] = p->B_z_dz[i];

            if(!errflag) {
                errflag = E_field_eval_E(E, yprev[0], yprev[1], yprev[2],
                                         t0, Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, yprev[0], yprev[1], yprev[2],
                                   t0, MHD_INCLUDE_ALL, boozer, mhd, Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k1, yprev, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                        (1.0/5) * k1[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (1.0/5)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (1.0/5)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + (1.0/5)*h[i], MHD_INCLUDE_ALL, boozer,
                                   mhd, Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k2, tempy, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          (3.0/40) * k1[j]
                        + (9.0/40) * k2[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (3.0/10)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/10)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + (3.0/10)*h[i], MHD_INCLUDE_ALL, boozer,
                                   mhd, Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k3, tempy, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          ( 3.0/10) * k1[j]
                        + (-9.0/10) * k2[j]
                        + ( 6.0/5 ) * k3[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (3.0/5)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/5)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + (3.0/5)*h[i], MHD_INCLUDE_ALL, boozer,
                                   mhd, Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k4, tempy, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          (-11.0/54) * k1[j]
                        + (  5.0/2 ) * k2[j]
                        + (-70.0/27) * k3[j]
                        + ( 35.0/27) * k4[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + h[i], MHD_INCLUDE_ALL, boozer, mhd,
                                   Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k5, tempy, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j]
                    + h[i]*(
                          ( 1631.0/55296 ) * k1[j]
                        + (  175.0/512   ) * k2[j]
                        + (  575.0/13824 ) * k3[j]
                        + (44275.0/110592) * k4[j]
                        + (  253.0/4096  ) * k5[j] );
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + (7.0/8)*h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + (7.0/8)*h[i], Edata, Bdata);
            }
            if(!errflag) {
                errflag = mhd_eval(mhd_dmhd, tempy[0], tempy[1], tempy[2],
                                   t0 + (7.0/8)*h[i], MHD_INCLUDE_ALL, boozer,
                                   mhd, Bdata);
            }
            if(!errflag) {
                step_gceom_mhd(
                    k6, tempy, mass, charge, B_dB, E, mhd_dmhd, aldforce);
            }

            /* Error estimate is a difference between RK4 and RK5 solutions. If
             * time-step is accepted, the RK5 solution will be used to advance
             * marker. */
            real rk5[6], rk4[6];
            if(!errflag) {
                real err = 0.0;
                for(int j = 0; j < 5; j++) {
                    rk5[j] = yprev[j]
                        + h[i]*(
                              ( 37.0/378 ) * k1[j]
                            + (250.0/621 ) * k3[j]
                            + (125.0/594 ) * k4[j]
                            + (512.0/1771) * k6[j] );

                    rk4[j] = yprev[j] +
                        h[i]*(
                              ( 2825.0/27648) * k1[j]
                            + (18575.0/48384) * k3[j]
                            + (13525.0/55296) * k4[j]
                            + (  277.0/14336) * k5[j]
                            + (    1.0/4    ) * k6[j] );
                    if(j==3) {
                        real yerr = fabs(rk5[j] - rk4[j]);
                        real ytol = fabs(yprev[j]) + fabs(k1[j]*h[i])
                                    + DBL_EPSILON;
                        err = fmax( err, yerr/ytol );
                    }
                    else if(j==2) {
                        real rk1[3] = {k1[0]*h[i], k1[1]*h[i], k1[2]*h[i]};
                        real yerr =
                              rk5[0] * rk5[0] + rk4[0] * rk4[0]
                            - 2 * rk5[0] * rk4[0] * cos(rk5[1] - rk4[1])
                              + ( rk5[2] - rk4[2] ) * ( rk5[2] - rk4[2] );
                        real ytol =
                              yprev[0] * yprev[0] + rk1[0] * rk1[0]
                            - 2 * yprev[0] * rk1[0] * cos(yprev[1] - rk1[1])
                            + ( yprev[2] - rk1[2] ) * ( yprev[2] - rk1[2] )
                            + DBL_EPSILON;
                        err = fmax( err, sqrt(yerr/ytol) );
                    }
                }

                err = err/tol;
                if(err <= 1){
                    /* Time step accepted */
                    hnext[i] = 0.85*h[i]*pow(err,-0.2);

                    /* Make sure we don't make a huge jump */
                    if(hnext[i] > 1.5*h[i]) {
                        hnext[i] = 1.5*h[i];
                    }
                }
                else{
                    /* Time step rejected */
                    hnext[i] = -0.85*h[i]*pow(err,-0.25);
                }
            }

            /* Update gc phase space position */
            if(!errflag) {
                p->r[i]     = rk5[0];
                p->phi[i]   = rk5[1];
                p->z[i]     = rk5[2];
                p->ppar[i]  = rk5[3];
                p->mu[i]    = rk5[4];
                p->zeta[i]  = fmod( rk5[5], CONST_2PI );
                if(p->zeta[i]<0) {
                    p->zeta[i] = CONST_2PI + p->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i],
                                            p->time[i] + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           p->time[i] + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
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
                errflag = B_field_get_axis_rz(axisrz, Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
                hnext[i]      = h[i];
            }
        }
    }
}
