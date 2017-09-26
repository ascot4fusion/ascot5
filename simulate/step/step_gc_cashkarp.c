/**
 * @file step_gc_cashkarp.c
 * @brief Calculate a guiding center step for a struct of particles with adaptive Cash Karp method
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../particle.h"
#include "step_gc_cashkarp.h"
#include "step_gceom.h"

/**
 * @brief Integrate a guiding center step for a struct of particles with adaptive Cash Karp method
 *
 * This function calculates a guiding center step for a struct of NSIMD 
 * particles with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions. 
 * All arrays in the function are of NSIMD length so vectorization can be performed 
 * directly without gather and scatter operations. Informs whther time step was accepted or
 * rejected and provides a suggestion for the next time step.
 *
 * @param p particle struct that will be updated
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step. Sign only indicates whether current step was rejected
 * @param tol error tolerance
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_gc_cashkarp(particle_simd_gc* p, real* h, real* hnext, real tol, B_field_data* Bdata, E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd aligned(h, hnext : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real k1[6];
            real k2[6];
            real k3[6];
            real k4[6];
	    real k5[6];
	    real k6[6];
            real tempy[6];
	    real yprev[6];

            real mass;
            real charge;

            real B_dB[12];
	    real E[3];

	    real R0   = p->r[i];
	    real z0   = p->z[i];

            /* Coordinates are copied from the struct into an array to make 
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->vpar[i];
            yprev[4] = p->mu[i];
	    yprev[5] = p->theta[i];
            mass = p->mass[i];
            charge = p->charge[i];

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

	    E_field_eval_E(E, yprev[0], yprev[1], yprev[2], Edata, Bdata);
	    step_gceom(k1, yprev, mass, charge, B_dB, E);
	    int j;

	    for(j = 0; j < 6; j++) {
		tempy[j] = yprev[j] + ((1.0/5)*k1[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
	    step_gceom(k2, tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 6; j++) {
		tempy[j] = yprev[j] + ((3.0/40)*k1[j]+(9.0/40)*k2[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
	    step_gceom(k3, tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 6; j++) {
		tempy[j] = yprev[j] + ((3.0/10)*k1[j]+(-9.0/10)*k2[j]+(6.0/5)*k3[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
	    step_gceom(k4, tempy, mass, charge, B_dB, E);
	
	    for(j = 0; j < 6; j++) {
		tempy[j] = yprev[j] + ((-11.0/54)*k1[j]+(5.0/2)*k2[j]+(-70.0/27)*k3[j]+(35.0/27)*k4[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
	    step_gceom(k5, tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 6; j++) {
		tempy[j] = yprev[j] + ((1631.0/55296)*k1[j]+(175.0/512)*k2[j]+(575.0/13824)*k3[j]+(44275.0/110592)*k4[j]+(253.0/4096)*k5[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata, Bdata);
	    step_gceom(k6, tempy, mass, charge, B_dB, E);

	    real yout[6];
	    real yerr;
	    real ytol;
	    real err = 0.0;
	    for(j = 0; j < 6; j++) {
		yout[j] = yprev[j] + ( (37.0/378)*k1[j] + (250.0/621)*k3[j] + (125.0/594)*k4[j] + (512.0/1771)*k6[j] )*h[i] ;
		yerr = fabs(yprev[j] + 
			    ( (2825.0/27648)*k1[j] + (18575.0/48384)*k3[j] + (13525.0/55296)*k4[j] + (277.0/14336)*k5[j] + (1.0/4)*k6[j] )*h[i] 
			    - yout[j]);
		ytol = fabs(yprev[j]) + fabs(k1[j]*h[i]);
		err = fmax(err,yerr/ytol);
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

	    p->r[i] = yout[0];
	    p->phi[i] = yout[1];
	    p->z[i] = yout[2];
	    p->vpar[i] = yout[3];
	    p->mu[i] = yout[4];
	    p->theta[i] = fmod(yout[5],CONST_2PI);
	    if(p->theta[i]<0){p->theta[i] = CONST_2PI + p->theta[i];}
		

	    /* Evaluate magnetic field (and gradient) and rho at new position */
	    B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i], Bdata);
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
	    B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
	    B_field_eval_rho(rho, psi[0], Bdata);
	    p->rho[i] = rho[0];

	    
	    /* Evaluate pol angle so that it is cumulative */
	    real axis_r = B_field_get_axis_r(Bdata);
	    real axis_z = B_field_get_axis_z(Bdata);
	    p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r), 
				(R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );

        }
    }
}
