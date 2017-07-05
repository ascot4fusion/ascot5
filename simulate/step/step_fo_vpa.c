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
      real pxyz[3];
      math_rpz2xyz(xhalf,pxyz);

      /* Precompute some values that will be used repeatedly */
      real sigma = p->charge[i]*h[i]/(2*p->mass[i]*CONST_C);
      real uminus[3];
      real g = (1/sqrt(1-math_dot(vxyz,vxyz)/CONST_C2))/CONST_C;
      uminus[0] = vxyz[0]*g + sigma*Exyz[0];
      uminus[1] = vxyz[1]*g + sigma*Exyz[1];
      uminus[2] = vxyz[2]*g + sigma*Exyz[2];
      real d = sigma*CONST_C/sqrt(1+math_dot(uminus,uminus));

      real Bhat[9] = {       0, -Bxyz[2],  Bxyz[1], 
		     Bxyz[2],        0,   -Bxyz[0], 
		      -Bxyz[1], Bxyz[0],        0};
      real dd[9];
      real a = (1+pow(d,2)*math_dot(Bxyz,Bxyz));
      int j;
      for(j = 0; j < 9 ; j++) dd[j]= 2*d*Bhat[j]/a;

      real m1[9];math_matmul(dd,Bhat,3,3,3,m1);
      for(j = 0; j < 9 ; j++) m1[j] = d*m1[j] + dd[j];
      real m2[3];math_matmul(m1,uminus,3,3,1,m2);
      real uplus[3];
      uplus[0] = uminus[0]+m2[0];
      uplus[1] = uminus[1]+m2[1];
      uplus[2] = uminus[2]+m2[2];

      /* Take the step */
      vxyz[0] = uplus[0] + sigma*Exyz[0];
      vxyz[1] = uplus[1] + sigma*Exyz[1];
      vxyz[2] = uplus[2] + sigma*Exyz[2];

      g = sqrt(1+math_dot(vxyz,vxyz))/CONST_C;
      vxyz[0] = vxyz[0]/(g);
      vxyz[1] = vxyz[1]/(g);
      vxyz[2] = vxyz[2]/(g);

      pxyz[0] = pxyz[0] + h[i]*vxyz[0]/2;
      pxyz[1] = pxyz[1] + h[i]*vxyz[1]/2;
      pxyz[2] = pxyz[2] + h[i]*vxyz[2]/2;

      /* Back to cylindrical coordinates */
      p->r[i] = sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]);
      p->phi[i] = atan2(pxyz[1], pxyz[0]);
      p->z[i] = pxyz[2];
      p->rdot[i] = vxyz[0] * cos(p->phi[i]) + vxyz[1] * sin(p->phi[i]);
      p->phidot[i] = ( -vxyz[0] * sin(p->phi[i]) + vxyz[1] * cos(p->phi[i]) ) / p->r[i];
      p->zdot[i] = vxyz[2];


      /* Evaluate magnetic field (and gradient) at new position */
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
    }
  }
}
