/**
 * @file step_ml_cashkarp.c
 * @brief Field line integrator implemented with Cash-Karp method.
 **/
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../math.h"
#include "../../particle.h"
#include "step_ml_cashkarp.h"

/**
 * @brief Integrate a magnetic field line step for a struct of markers
 *
 * This function calculates a magnetic field line step for a struct of NSIMD
 * markers with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions.
 * All arrays in the function are of NSIMD length so vectorization can be
 * performed directly without gather and scatter operations. Informs whether
 * time step was accepted or rejected and provides a suggestion for the next
 * time step.
 *
 * @param p marker struct that will be integrated
 * @param t array containting time values
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step. Negative if rejected.
 * @param tol error tolerance
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_ml_cashkarp(particle_simd_ml* p, real* h, real* hnext, real tol,
                      B_field_data* Bdata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real k1[3], k2[3], k3[3], k4[3], k5[3], k6[3];
            real tempy[3];
            real yprev[3];

            real normB;

            real R0 = p->r[i];
            real z0 = p->z[i];
            real t0 = p->time[i];

            /* Direction */
            int direction = 1;
            if(p->pitch[i] < 0) {
                direction = -1;
            }

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];

            /* Magnetic field at initial position already known */
            k1[0] = p->B_r[i];
            k1[1] = p->B_phi[i];
            k1[2] = p->B_z[i];

            normB = (math_normc(k1[0], k1[1], k1[2])) * direction;
            k1[0] /= normB;
            k1[1] /= normB;
            k1[2] /= normB;
            k1[1] /= yprev[0];

            for(int j = 0; j < 3; j++) {
                tempy[j] = yprev[j]
                    + h[i] * ( (1.0/5) * k1[j] );
            }
            if(!errflag) {
                errflag = B_field_eval_B(k2, tempy[0], tempy[1], tempy[2],
                                         t0 + (1.0/5)*h[i], Bdata);
            }
            normB = (math_normc(k2[0], k2[1], k2[2])) * direction;
            k2[0] /= normB;
            k2[1] /= normB;
            k2[2] /= normB;
            k2[1] /= tempy[0];

            for(int j = 0; j < 3; j++) {
                tempy[j] = yprev[j]
                    + h[i] * (
                          (3.0/40) * k1[j]
                        + (9.0/40) * k2[j] );
            }
            if(!errflag) {
                errflag = B_field_eval_B(k3, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/10)*h[i], Bdata);
            }
            normB = (math_normc(k3[0], k3[1], k3[2])) * direction;
            k3[0] /= normB;
            k3[1] /= normB;
            k3[2] /= normB;
            k3[1] /= tempy[0];

            for(int j = 0; j < 3; j++) {
                tempy[j] = yprev[j]
                    + h[i] * (
                          ( 3.0/10) * k1[j]
                        + (-9.0/10) * k2[j]
                        + ( 6.0/5 ) * k3[j] );
            }
            if(!errflag) {
                errflag = B_field_eval_B(k4, tempy[0], tempy[1], tempy[2],
                                         t0 + (3.0/5)*h[i], Bdata);
            }
            normB = (math_normc(k4[0], k4[1], k4[2])) * direction;
            k4[0] /= normB;
            k4[1] /= normB;
            k4[2] /= normB;
            k4[1] /= tempy[0];

            for(int j = 0; j < 3; j++) {
                tempy[j] = yprev[j]
                    + h[i] * (
                          (-11.0/54) * k1[j]
                        + (  5.0/2 ) * k2[j]
                        + (-70.0/27) * k3[j]
                        + ( 35.0/27) * k4[j] );
            }
            if(!errflag) {
                errflag = B_field_eval_B(k5, tempy[0], tempy[1], tempy[2],
                                         t0 + h[i], Bdata);
            }
            normB = (math_normc(k5[0], k5[1], k5[2])) * direction;
            k5[0] /= normB;
            k5[1] /= normB;
            k5[2] /= normB;
            k5[1] /= tempy[0];

            for(int j = 0; j < 3; j++) {
                tempy[j] = yprev[j]
                    + h[i] * (
                          ( 1631.0/55296 ) * k1[j]
                        + (  175.0/512   ) * k2[j]
                        + (  575.0/13824 ) * k3[j]
                        + (44275.0/110592) * k4[j]
                        + (  253.0/4096  ) * k5[j] );
            }
            if(!errflag) {
                errflag = B_field_eval_B(k6, tempy[0], tempy[1], tempy[2],
                               t0 + (7.0/8)*h[i], Bdata);
            }
            normB = (math_normc(k6[0], k6[1], k6[2])) * direction;
            k6[0] /= normB;
            k6[1] /= normB;
            k6[2] /= normB;
            k6[1] /= tempy[0];

            /* Error estimate is a difference between RK4 and RK5 solutions. If
             * time-step is accepted, the RK5 solution will be used to advance
             * marker. */
            real rk5[3];
            real err = 0.0;
            for(int j = 0; j < 3; j++) {
                rk5[j] = yprev[j]
                    + h[i] * (
                          ( 37.0/378 ) * k1[j]
                        + (250.0/621 ) * k3[j]
                        + (125.0/594 ) * k4[j]
                        + (512.0/1771) * k6[j] );

                real rk4 = yprev[j] +
                            h[i] * (
                                  ( 2825.0/27648) * k1[j]
                                + (18575.0/48384) * k3[j]
                                + (13525.0/55296) * k4[j]
                                + (  277.0/14336) * k5[j]
                                + (    1.0/4    ) * k6[j] );

                real yerr = fabs(rk5[j] - rk4);
                real ytol = fabs(yprev[j]) + fabs( k1[j] * h[i] ) + DBL_EPSILON;
                err  = fmax( err, yerr/ytol );
            }

            err = err/tol;
            if(err <= 1){
                /* Time step accepted */
                hnext[i] = 0.85*h[i]*pow(err,-0.2);

            }
            else{
                /* Time step rejected */
                hnext[i] = -0.85*h[i]*pow(err,-0.25);

            }

            p->r[i]   = rk5[0];
            p->phi[i] = rk5[1];
            p->z[i]   = rk5[2];

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[15];
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i],
                                            p->time[i] + h[i], Bdata);
            }
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


            real psi[1];
            real rho[1];
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }
            p->rho[i] = rho[0];


            /* Evaluate pol angle so that it is cumulative */
            real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
            real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
            p->pol[i] += atan2(   (R0-axis_r) * (p->z[i]-axis_z)
                                - (z0-axis_z) * (p->r[i]-axis_r),
                                  (R0-axis_r) * (p->r[i]-axis_r)
                                + (z0-axis_z) * (p->z[i]-axis_z) );

        }
    }
}
