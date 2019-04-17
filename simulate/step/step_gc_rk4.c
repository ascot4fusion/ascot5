/**
 * @file step_gc_rk4.c
 * @brief Guiding center integration with RK4
 **/
#include <stdio.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../math.h"
#include "../../particle.h"
#include "../../error.h"
#include "step_gceom.h"
#include "step_gc_rk4.h"

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
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_gc_rk4(particle_simd_gc* p, real* h, B_field_data* Bdata,
                 E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd aligned(h : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

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
            yprev[3] = p->vpar[i];
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
                errflag = E_field_eval_E(E, yprev[0], yprev[1], yprev[2],
                                         t0, Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k1, yprev, mass, charge, B_dB, E);
            }


            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k1[j]/2.0;
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k2, tempy, mass, charge, B_dB, E);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k2[j]/2.0;
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i]/2.0, Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i]/2.0, Edata, Bdata);
            }
            if(!errflag) {
                step_gceom(k3, tempy, mass, charge, B_dB, E);
            }
            for(int j = 0; j < 6; j++) {
                tempy[j] = yprev[j] + h[i]*k3[j];
            }


            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(E, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], Edata, Bdata);}
            if(!errflag) {step_gceom(k4, tempy, mass, charge, B_dB, E);
            }
            for(int j = 0; j < 6; j++) {
                y[j] = yprev[j]
                    + h[i]/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            }


            /* Test that results are physical */
            if(!errflag && y[0] <= 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }
            if(!errflag && fabs(y[4]) >= CONST_C) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }
            if(!errflag && y[4] < 0) {
                errflag = error_raise(ERR_INTEGRATION, __LINE__, EF_STEP_GC_RK4);
            }

            /* Update gc phase space position */
            if(!errflag) {
                p->r[i] = y[0];
                p->phi[i] = y[1];
                p->z[i] = y[2];
                p->vpar[i] = y[3];
                p->mu[i] = y[4];
                p->zeta[i] = fmod(y[5],CONST_2PI);
                if(p->zeta[i]<0) {
                    p->zeta[i] = CONST_2PI + p->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[1];
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           t0 + h[i], Bdata);
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
                real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
                real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axis_r) * (p->z[i]-axis_z)
                                      - (z0-axis_z) * (p->r[i]-axis_r),
                                        (R0-axis_r) * (p->r[i]-axis_r)
                                      + (z0-axis_z) * (p->z[i]-axis_z) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}
