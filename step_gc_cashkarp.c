/**
 * @file step_gc_cashkarp.c
 * @brief Calculate a guiding center step for a struct of particles with adaptive Cash Karp method
 **/
#include <math.h>
#include "ascot5.h"
#include "step_gc_cashkarp.h"
#include "step_gc_rk4.h"
#include "B_field.h"
#include "math.h"
#include "particle.h"

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
 * @param t array containting time values
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step. Sign only indicates whether current step was rejected
 * @param tol error tolerance
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_gc_cashkarp(particle_simd_gc* p, real* t, real* h, real* hnext, real tol, B_field_data* Bdata, E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd 
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real k1[5];
            real k2[5];
            real k3[5];
            real k4[5];
	    real k5[5];
	    real k6[5];
            real tempy[5];
	    real yprev[5];

           /* Mass and charge need to have the same structure as in the previous
               arrays so that the compiler can vectorize the function call, but
               we only use the first row for actual data */
            real mass[5];
            real charge[5];

            real B[3];
            real B_dB[12];
	    real E[3];

            /* Coordinates are copied from the struct into an array to make 
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->vpar[i];
            yprev[4] = p->mu[i];
            mass[0] = p->mass[i];
            charge[0] = p->charge[i];

	    B_field_eval_B_dB(B_dB, yprev[0], yprev[1], yprev[2], Bdata);
	    E_field_eval_E(E, yprev[0], yprev[1], yprev[2], Edata);
	    ydot_gc(k1, t[0], yprev, mass, charge, B_dB, E);
	    int j;

	    for(j = 0; j < 5; j++) {
		tempy[j] = yprev[j] + ((1.0/5)*k1[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata);
	    ydot_gc(k2, t[i]+h[i]/5.0, tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 5; j++) {
		tempy[j] = yprev[j] + ((3.0/40)*k1[j]+(9.0/40)*k2[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata);
	    ydot_gc(k3, t[i]+h[i]*(3.0/10), tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 5; j++) {
		tempy[j] = yprev[j] + ((3.0/10)*k1[j]+(-9.0/10)*k2[j]+(6.0/5)*k3[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata);
	    ydot_gc(k4, t[i]+h[i]*(3.0/5), tempy, mass, charge, B_dB, E);
	
	    for(j = 0; j < 5; j++) {
		tempy[j] = yprev[j] + ((-11.0/54)*k1[j]+(5.0/2)*k2[j]+(-70.0/27)*k3[j]+(35.0/27)*k4[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata);
	    ydot_gc(k5, t[i]+h[i], tempy, mass, charge, B_dB, E);

	    for(j = 0; j < 5; j++) {
		tempy[j] = yprev[j] + ((1631.0/55296)*k1[j]+(175.0/512)*k2[j]+(575.0/13824)*k3[j]+(44275.0/110592)*k4[j]+(253.0/4096)*k5[j])*h[i];
	    }
	    B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
	    E_field_eval_E(E, tempy[0], tempy[1], tempy[2], Edata);
	    ydot_gc(k6, t[i]+h[i]*(7.0/8), tempy, mass, charge, B_dB, E);

	    real yout[5];
	    real yerr;
	    real ytol;
	    real err = 0.0;
	    for(j = 0; j < 5; j++) {
		yout[j] = yprev[j] + (37.0/378)    *k1[j] + (250.0/621)    *k3[j]
		    + (125.0/594)    *k4[j] + (512.0/1771)*k6[j];
		yerr = fabs((yprev[j] + (2825.0/27648)*k1[j] + (18575.0/48384)*k3[j]
				+ (13525.0/55296)*k4[j] + (277.0/14336)*k5[j] +      (1.0/4)*k6[j])*h[i] - yout[j]);
		ytol = fabs(yprev[j]) + fabs(k1[j]*h[i]);
		err = fmax(err,yerr/ytol);
	    }
        
	    if(err <= 1.0){
		/* Time step accepted */
	        hnext[i] = h[i]*pow(err,-0.2);
		p->r[i] = yout[0];
		p->phi[i] = yout[1];
		p->z[i] = yout[2];
		p->vpar[i] = yout[3];
		
		p->time[i] = p->time[i] + h[i];
		
		/* Update other particle parameters to be consistent */
		B_field_eval_B(B, yout[0], yout[1], yout[2], Bdata);
		p->B_r[i] = B[0];
		p->B_phi[i] = B[1];
		p->B_z[i] = B[2]; 
		
		p->prev_r[i] = yprev[0];
		p->prev_phi[i] = yprev[1];
		p->prev_z[i] = yprev[2];
	    }
	    else{
		/* Time step accepted */
	        hnext[i] = -h[i]*pow(err,-0.25);
	    }

        }
    }
}
