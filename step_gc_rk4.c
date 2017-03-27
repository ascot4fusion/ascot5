/**
 * @file step_gc_rk4.c
 * @brief Calculate a guiding center step for a struct of particles with RK4
 **/
#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "step_gc_rk4.h"
#include "phys_orbit.h"
#include "B_field.h"
#include "E_field.h"
#include "math.h"
#include "particle.h"

/**
 * @brief Integrate a guiding center step for a struct of particles with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD 
 * particles with RK4 simultaneously using SIMD instructions. All arrays in the 
 * function are of NSIMD length so vectorization can be performed directly 
 * without gather and scatter operations.
 *
 * @param p particle struct that will be updated
 * @param t time
 * @param h length of time step
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_gc_rk4(particle_simd_gc* p, real t, real h, B_field_data* Bdata, E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd 
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real k1[5];
            real k2[5];
            real k3[5];
            real k4[5];
            real tempy[5];
            real yprev[5];
            real y[5];

            real mass;
            real charge;

            real B[3];
            real B_dB[12];
            real rho_drho[4];
	    real E[3];

            /* Coordinates are copied from the struct into an array to make 
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->vpar[i];
            yprev[4] = p->mu[i];
            mass = p->mass[i];
            charge = p->charge[i];

            B_field_eval_B_dB(B_dB, yprev[0], yprev[1], yprev[2], Bdata);
            B_field_eval_rho_drho(rho_drho, yprev[0], yprev[1], yprev[2], Bdata);
	    E_field_eval_E(E, rho_drho, Edata);
            phys_eomgc(k1, t, yprev, mass, charge, B_dB, E);
            int j;
            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h/2.0*k1[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            B_field_eval_rho_drho(rho_drho, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, rho_drho, Edata);
            phys_eomgc(k2, t+h/2.0, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h/2.0*k2[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            B_field_eval_rho_drho(rho_drho, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, rho_drho, Edata);
            phys_eomgc(k3, t+h/2.0, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h*k3[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            B_field_eval_rho_drho(rho_drho, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, rho_drho, Edata);
            phys_eomgc(k4, t+h, tempy, mass, charge, B_dB, E);
            for(j = 0; j < 5; j++) {
                y[j] = yprev[j]
                    + h/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            } 

            p->r[i] = y[0];
            p->phi[i] = y[1];
            p->z[i] = y[2];
            p->vpar[i] = y[3];

            p->time[i] = p->time[i] + h;

            /* Update other particle parameters to be consistent */
            B_field_eval_B(B, y[0], y[1], y[2], Bdata);
            p->B_r[i] = B[0];
            p->B_phi[i] = B[1];
            p->B_z[i] = B[2]; 

            p->prev_r[i] = yprev[0];
            p->prev_phi[i] = yprev[1];
            p->prev_z[i] = yprev[2];

        }
    }
}
