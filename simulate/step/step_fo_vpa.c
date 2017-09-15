/**
 * @file step_fo_lf.c
 * @brief Calculate a full orbit step for a struct of particles with leap-frog
 **/
#include <math.h>
#include <stdio.h>
#include "../../math.h"
#include "../../consts.h"
#include "step_fo_vpa.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../particle.h"

/**
 * @brief Integrate a full orbit step for a struct of particles with VPA
 *
 * The integration is performed for a struct of NSIMD particles using the 
 * volume preserving algorithm (Boris method for relativistic particles) see Zhang 2015.
 *
 * @param p particle_simd_fo struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_fo_vpa(particle_simd_fo* p, real* h, B_field_data* Bdata, E_field_data* Edata) {

  int i;
  /* Following loop will be executed simultaneously for all i */
#pragma omp simd 
  for(i = 0; i < NSIMD; i++) {
    if(p->running[i]) {

	real phi0 = p->phi[i];
	real R0   = p->r[i];
	real z0   = p->z[i];
	    
	/* Take a half step and evaluate fields at that position */
	real xhalf[3];
	xhalf[0]= p->r[i] + p->rdot[i]*h[i]/2;
	xhalf[1]= p->phi[i] + p->phidot[i]*h[i]/2;
	xhalf[2]= p->z[i] + p->zdot[i]*h[i]/2;
	    
	real Brpz[3];
	real Erpz[3];
	B_field_eval_B(Brpz, xhalf[0], xhalf[1], xhalf[2], Bdata);
	E_field_eval_E(Erpz, xhalf[0], xhalf[1], xhalf[2], Edata, Bdata);
      

	/* Electromagnetic fields to cartesian coordinates */  
	real Bxyz[3];
	real Exyz[3];

	math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
	math_vec_rpz2xyz(Erpz, Exyz, p->phi[i]);

	/* Convert velocity to cartesian coordinates */
	real vxyz[3];
	real vrpz[3] = {p->rdot[i], p->phidot[i]*p->r[i], p->zdot[i]};
	math_vec_rpz2xyz(vrpz, vxyz, p->phi[i]);

	/* Positions to cartesian coordinates */
	real posxyz[3];
	math_rpz2xyz(xhalf,posxyz);

	/* Evaluate helper variable pminus */
	real pminus[3];
	real gamma = 1 / sqrt( 1 - math_dot(vxyz,vxyz)/CONST_C2 );
	real sigma = p->charge[i]*h[i]/2;

	pminus[0] = p->mass[i]*vxyz[0]*gamma + sigma*Exyz[0];
	pminus[1] = p->mass[i]*vxyz[1]*gamma + sigma*Exyz[1];
	pminus[2] = p->mass[i]*vxyz[2]*gamma + sigma*Exyz[2];

	/* Second helper variable pplus*/
	real d = (p->charge[i]*h[i]/2) / 
	    sqrt( p->mass[i]*p->mass[i] + math_dot(pminus,pminus)/CONST_C2 );
	real d2 = d*d;

	real Bhat[9] = {       0,  Bxyz[2], -Bxyz[1], 
			-Bxyz[2],        0,  Bxyz[0], 
			 Bxyz[1], -Bxyz[0],       0};
	real Bhat2[9];
	math_matmul(Bhat, Bhat, 3, 3, 3, Bhat2);
	
	real B2 = Bxyz[0]*Bxyz[0] + Bxyz[1]*Bxyz[1] + Bxyz[2]*Bxyz[2];
	
	real A[9];
	for(int j=0; j<9; j++) {
	    A[j] = (d*Bhat[j] + d2*Bhat2[j]) * (2.0/(1+d2*B2));
	}

	// Add identity matrix to the mix
	A[0] += 1;
	A[4] += 1;
	A[8] += 1;

	real pplus[3];
	math_matmul(pminus, A, 1, 3, 3, pplus);

	/* Take the step */
	real pfinal[3];
	pfinal[0] = pplus[0] + sigma*Exyz[0];
	pfinal[1] = pplus[1] + sigma*Exyz[1];
	pfinal[2] = pplus[2] + sigma*Exyz[2];

	// gamma = sqrt(1+(p/mc)^2)
	gamma = sqrt( 1 + math_dot(pfinal,pfinal)/(p->mass[i]*p->mass[i]*CONST_C2) );

	vxyz[0] = pfinal[0]/(gamma*p->mass[i]);
	vxyz[1] = pfinal[1]/(gamma*p->mass[i]);
	vxyz[2] = pfinal[2]/(gamma*p->mass[i]);

	posxyz[0] = posxyz[0] + h[i]*vxyz[0]/2;
	posxyz[1] = posxyz[1] + h[i]*vxyz[1]/2;
	posxyz[2] = posxyz[2] + h[i]*vxyz[2]/2;

	/* Back to cylindrical coordinates */
	p->r[i] = sqrt(posxyz[0]*posxyz[0]+posxyz[1]*posxyz[1]);
	p->phi[i] = atan2(posxyz[1], posxyz[0]);
	p->z[i] = posxyz[2];
	p->rdot[i] = vxyz[0] * cos(p->phi[i]) + vxyz[1] * sin(p->phi[i]);
	p->phidot[i] = ( -vxyz[0] * sin(p->phi[i]) + vxyz[1] * cos(p->phi[i]) ) / p->r[i];
	p->zdot[i] = vxyz[2];


	/* Evaluate magnetic field (and gradient) and rho at new position */
	real BdBrpz[12];
	B_field_eval_B_dB(BdBrpz, p->r[i], p->phi[i], p->z[i], Bdata);
	p->B_r[i]        = BdBrpz[0];
	p->B_r_dr[i]     = BdBrpz[1];
	p->B_r_dphi[i]   = BdBrpz[2];
	p->B_r_dz[i]     = BdBrpz[3];

	p->B_phi[i]      = BdBrpz[4];
	p->B_phi_dr[i]   = BdBrpz[5];
	p->B_phi_dphi[i] = BdBrpz[6];
	p->B_phi_dz[i]   = BdBrpz[7];

	p->B_z[i]        = BdBrpz[8];
	p->B_z_dr[i]     = BdBrpz[9];
	p->B_z_dphi[i]   = BdBrpz[10];
	p->B_z_dz[i]     = BdBrpz[11];

	real psi[1];
	real rho[1];
	B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
	B_field_eval_rho(rho, psi[0], Bdata);
	p->rho[i] = rho[0];

	/* Evaluate phi and pol angles so that they are cumulative */
	real axis_r = B_field_get_axis_r(Bdata);
	real axis_z = B_field_get_axis_z(Bdata);
	p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r), 
			    (R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );
	real tphi = fmod(phi0 , CONST_2PI );
	if(tphi < 0){tphi = CONST_2PI+tphi;}
	tphi = fmod(p->phi[i]+CONST_2PI,CONST_2PI) -  tphi;
        
	p->phi[i] = phi0 + tphi;
    }
  }
}
