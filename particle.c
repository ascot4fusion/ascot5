/**
 * @file particle.c
 * @brief Particle representations and helper functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "math.h"
#include "consts.h"
#include "particle.h"
#include "B_field.h"
#include "E_field.h"

/**
 * @brief Transforms particle struct into a full-orbit simulation struct
 *
 * @param p     pointer to the particle being processed
 * @param i     index
 * @param p_fo  pointer to particle_simd_fo array
 * @param j     index of the new fo struct in the SIMD array
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void particle_to_fo(particle* p, int i, particle_simd_fo* p_fo, int j,
                    B_field_data* Bdata, E_field_data* Edata){
    p_fo->r[j] = p->r;
    p_fo->phi[j] = p->phi;
    p_fo->z[j] = p->z;
    p_fo->rdot[j] = p->rdot;
    p_fo->phidot[j] = p->phidot;
    p_fo->zdot[j] = p->zdot;
    p_fo->mass[j] = p->mass;
    p_fo->charge[j] = p->charge;
    p_fo->weight[j] = p->weight;
    p_fo->time[j] = p->time;
    p_fo->id[j] = p->id; 
    p_fo->running[j] = p->running;
    p_fo->endcond[j] = p->endcond; 
    p_fo->walltile[j] = p->walltile;

    real B[3];
    B_field_eval_B(B, p->r, p->phi, p->z, Bdata);
    real E[3];
    E_field_eval_E(E, p->r, p->phi, p->z, Edata);

    p_fo->B_r[j] = B[0];					  
    p_fo->B_phi[j] = B[1];				
    p_fo->B_z[j] = B[2];				
    p_fo->E_r[j] = E[0];					  
    p_fo->E_phi[j] = E[1]; 					 
    p_fo->E_z[j] = E[2];
    p_fo->prev_r[j] = p->r;
    p_fo->prev_phi[j] = p->phi;
    p_fo->prev_z[j] = p->z;
    p_fo->index[j] = i;
}

/**
 * @brief Makes a dummy full-orbit simulation struct
 *
 * @param p_fo  pointer to particle_simd_fo array
 * @param j     index of the new fo struct in the SIMD array
 */
void particle_to_fo_dummy(particle_simd_fo* p_fo, int j){
    p_fo->r[j] = 1;
    p_fo->phi[j] = 1;
    p_fo->z[j] = 1;
    p_fo->rdot[j] = 1;
    p_fo->phidot[j] = 1;
    p_fo->zdot[j] = 1;
    p_fo->mass[j] = 1;
    p_fo->charge[j] = 1;
    p_fo->weight[j] = 0;
    p_fo->time[j] = 0;
    p_fo->id[j] = -1; 
    p_fo->running[j] = 0;
    p_fo->endcond[j] = 0; 
    p_fo->walltile[j] = 0;
    p_fo->B_r[j] = 1;					  
    p_fo->B_phi[j] = 1;				
    p_fo->B_z[j] = 1;				
    p_fo->E_r[j] = 1;					  
    p_fo->E_phi[j] = 1; 					 
    p_fo->E_z[j] = 1;
    p_fo->prev_r[j] = 1;
    p_fo->prev_phi[j] = 1;
    p_fo->prev_z[j] = 1;
    p_fo->index[j] = -1;
}


void particle_to_gc(particle* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata) {
    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, Bdata);
    real gcpos[5];
    phasespace_particle_to_guidingcenter(p->mass, p->charge, p->r, p->phi, p->z, 
					 p->rdot, p->phidot*p->r, p->zdot, B_dB, gcpos);

    p_gc->r[j] = gcpos[0];
    p_gc->phi[j] = gcpos[1];
    p_gc->z[j] =gcpos[2];
    
    real B[3];
    B_field_eval_B(B, p_gc->r[j], p_gc->phi[j], p_gc->z[j], Bdata);
    real gamma = sqrt(1+2*gcpos[4]/(p->mass*CONST_C2)+pow(gcpos[3]/(p->mass*CONST_C),2) );
    p_gc->vpar[j] = gcpos[3]/(p->mass*gamma);
    p_gc->mu[j] = gcpos[4];

    //p_gc->theta[j] = 0atan2(math_norm(ezcrossrho), math_dot(ez, rho));
    p_gc->theta[j] = 0.0; // Ill defined quantity
/*    p_gc->theta[j] = acos(math_dot(ez, rho)
                          / (math_norm(ez)*math_norm(rho))); */

    p_gc->mass[j] = p->mass;
    p_gc->charge[j] = p->charge;
    p_gc->weight[j] = p->weight;
    p_gc->time[j] = 0;
    p_gc->id[j] = p->id;
    p_gc->running[j] = p->running;
    p_gc->B_r[j] = B_dB[0];
    p_gc->B_phi[j] = B_dB[4];
    p_gc->B_z[j] = B_dB[8];
    p_gc->index[j] = i;
}

void particle_to_gc_dummy(particle_simd_gc* p_gc, int j) {
    p_gc->r[j] = 1;
    p_gc->phi[j] = 1;
    p_gc->z[j] = 1;
    p_gc->vpar[j] = 1;
    p_gc->mu[j] = 1;
    p_gc->theta[j] = 1;
    p_gc->B_r[j] = 1;
    p_gc->B_phi[j] = 1;
    p_gc->B_z[j] = 1;
    p_gc->mass[j] = 1;
    p_gc->charge[j] = 1;
    p_gc->time[j] = 0;
    p_gc->weight[j] = 0;
    p_gc->id[j] = -1;
    p_gc->running[j] = 0;
    p_gc->index[j] = -1;
}

void gc_to_particle(particle_simd_gc* p_gc, int j, particle* p) {
    real B[3];
    B[0] = p_gc->B_r[j];
    B[1] = p_gc->B_phi[j];
    B[2] = p_gc->B_z[j];
    real normB = sqrt(math_dot(B, B));
    real Bxyz[3];
    math_vec_rpz2xyz(B, Bxyz, p_gc->phi[j]);

    real ez[3];
    ez[0] = -Bxyz[2]/normB * Bxyz[0]/normB;
    ez[1] = -Bxyz[2]/normB * Bxyz[1]/normB;
    ez[2] = 1 - Bxyz[2]/normB * Bxyz[2]/normB;

    real a[3];
    a[0] = ez[0]/math_norm(ez);
    a[1] = ez[1]/math_norm(ez);
    a[2] = ez[2]/math_norm(ez);
    
    real b[3];
    math_cross(B,ez,b);
    real tmp = math_norm(b);
    b[0] = b[0]/tmp;
    b[1] = b[1]/tmp;
    b[2] = b[2]/tmp;

    real rho = fabs(sqrt(2*p_gc->mu[j]*p_gc->mass[j]/normB)/p_gc->charge[j]);
    real rhovec[3];
    rhovec[0] = rho*(sin(p_gc->theta[j]) * a[0] + cos(p_gc->theta[j]) * b[0]);
    rhovec[1] = rho*(sin(p_gc->theta[j]) * a[1] + cos(p_gc->theta[j]) * b[1]);
    rhovec[2] = rho*(sin(p_gc->theta[j]) * a[2] + cos(p_gc->theta[j]) * b[2]);

    real rhovecrpz[3];
    math_vec_xyz2rpz(rhovec,rhovecrpz,p_gc->phi[j]);

/*    p->r = p_gc->r[j] + rhovecrpz[0];
    p->phi = p_gc->phi[j] + rhovecrpz[1];
    p->z = p_gc->z[j] + rhovecrpz[2];*/

    p->r = p_gc->r[j];
    p->phi = p_gc->phi[j];
    p->z = p_gc->z[j];
    p->rdot = p_gc->mu[j];
    p->phidot = p_gc->vpar[j];
    p->zdot = 0;
    p->time = p_gc->time[j];
    p->running = p_gc->running[j];
    p->endcond = p_gc->endcond[j];
    p->walltile = p_gc->walltile[j];
}

/**
 * @brief Transforms particle to guiding center phase space 
 *
 * The transformation is done between coordinates [r,phi,z,p_r,p_phi,p_z]
 * and [R,Phi,Z,p_para,mu]. The transformation is done up to first order.
 *
 * @todo Should get momentum as an input
 *  
 */
void phasespace_particle_to_guidingcenter(real mass, real charge, real r, real phi, real z, 
					  real v_r, real v_phi, real v_z, real* B_dB, real* gcpos){
    /* Temporary re-used variables */
    int sum_index;
    real TEMP_S1;
    real TEMP_V1[3];

    /* Helper variables */
    real velocity[3] = {v_r, v_phi, v_z};
    real v_unit[3];
    math_unit(velocity, v_unit);
    real momentum = sqrt( 1/(1-math_dot(velocity,velocity)/CONST_C2) - 1 ) * mass*CONST_C;
    real gamma = sqrt( 1 + pow(momentum/(mass*CONST_C), 2) );

    real B_vec[3]  = {B_dB[0], B_dB[4], B_dB[8]};
    real B_norm = math_norm(B_vec);
    real B_unit[3];
    math_unit(B_vec, B_unit);

    real pitch = math_dot(v_unit,B_unit);

    /* Magnetic field jacobian, gradient and curl in cylindrical coordinates */
    real jacB[9] = {B_dB[1], B_dB[2]/r, B_dB[3],
		    B_dB[5], B_dB[6]/r, B_dB[7],
		    B_dB[9], B_dB[10]/r, B_dB[11]};

    real gradB[3];
    gradB[0] = B_unit[0]*jacB[0] + B_unit[0]*jacB[3] + B_unit[0]*jacB[6];
    gradB[1] = B_unit[1]*jacB[1] + B_unit[1]*jacB[4] + B_unit[1]*jacB[7];
    gradB[2] = B_unit[2]*jacB[2] + B_unit[2]*jacB[5] + B_unit[2]*jacB[8];
    math_matmul(jacB,B_unit,3,3,1,gradB);

    real curlB[3] = {jacB[7]-jacB[5], jacB[2]-jacB[6], jacB[3]+B_vec[1]/r-jacB[1]};

    real tau_B = math_dot(B_unit, curlB)/B_norm;
    real nablabhat[9];
    nablabhat[0]=(jacB[0] - gradB[0]*B_unit[0])/B_norm;
    nablabhat[1]=(jacB[1] - gradB[0]*B_unit[1])/B_norm;
    nablabhat[2]=(jacB[2] - gradB[0]*B_unit[2])/B_norm;
    nablabhat[3]=(jacB[3] - gradB[1]*B_unit[0])/B_norm;
    nablabhat[4]=(jacB[4] - gradB[1]*B_unit[1])/B_norm;
    nablabhat[5]=(jacB[5] - gradB[1]*B_unit[2])/B_norm;
    nablabhat[6]=(jacB[6] - gradB[2]*B_unit[0])/B_norm;
    nablabhat[7]=(jacB[7] - gradB[2]*B_unit[1])/B_norm;
    nablabhat[8]=(jacB[8] - gradB[2]*B_unit[2])/B_norm; 

    real kappa[3];
    math_matmul(nablabhat,B_unit,3,3,1,kappa);
    
    /* Zeroth order momentum terms */
    real p_para0 = pitch*momentum;
    real mu_0 = ( 1 - pow(pitch,2) )*pow(momentum,2)/(2*mass*B_norm);

    /* Make the spatial transformation */
    real rho[3];
    math_cross(v_unit,B_unit,rho);
    math_prod(rho,momentum/(charge*B_norm));

    real rho_unit[3];
    math_unit(rho,rho_unit);
    real rho_norm = math_norm(rho);

    gcpos[0] = r + rho[0];
    gcpos[1] = phi + atan(rho[1]/r);
    gcpos[2] = z + rho[2];

    /* First order momentum terms */
    real perphat[3];
    math_cross(rho_unit, B_unit, perphat);

    real a1ddotgradb = 0.5*(2*(rho_unit[0]*perphat[0]*nablabhat[0] + rho_unit[1]*perphat[1]*nablabhat[4]
			      +rho_unit[2]*perphat[2]*nablabhat[8]) 
			   +(rho_unit[0]*perphat[1]+rho_unit[1]*perphat[0])*(nablabhat[1]+nablabhat[3])
			   +(rho_unit[0]*perphat[2]+rho_unit[2]*perphat[0])*(nablabhat[2]+nablabhat[6])
			   +(rho_unit[1]*perphat[2]+rho_unit[2]*perphat[1])*(nablabhat[5]+nablabhat[7]));

    real p_para1 = p_para0*math_dot(rho,kappa) - ((mass*mu_0)/charge)*(tau_B + a1ddotgradb);
    TEMP_S1 = pow(p_para0,2)/(mass*B_norm);
    TEMP_V1[0] = TEMP_S1 * kappa[0] + mu_0*gradB[0]/B_norm;
    TEMP_V1[1] = TEMP_S1 * kappa[1] + mu_0*gradB[1]/B_norm;
    TEMP_V1[2] = TEMP_S1 * kappa[2] + mu_0*gradB[2]/B_norm;
    real mu_1 = -math_dot(rho, TEMP_V1) + ( (mu_0*p_para0)/(charge*B_norm) )*(tau_B + a1ddotgradb);

    /* Make the momentum transformation */
    gcpos[3] = p_para0 + p_para1;
    gcpos[4] = mu_0 + mu_1;

    /* Finally, change p_para into v_para as the gcpos[3] is defined */
    //gamma = sqrt(1+2*gcpos[4]/(mass*CONST_C2)+pow(gcpos[3]/(mass*CONST_C),2) );
    //gcpos[3] = gcpos[3]/(mass*gamma);
}

void phasespace_guidingcenter_to_particle(real R, real Phi, real Z, real v_para, real mu, 
					  real B_dB, real* gcpos){


}
