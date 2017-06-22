/**
 * @file particle.c
 * @brief Particle representations and helper functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "phys_orbit.h"
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
    p_fo->rdot[j] = p->v_r;
    p_fo->phidot[j] = p->v_phi/p->r;
    p_fo->zdot[j] = p->v_z;
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
    real rho_drho[4];
    real E[3];
    B_field_eval_rho_drho(rho_drho, p->r, p->phi, p->z, Bdata);
    /* Convert partial derivative to gradient */
    rho_drho[2] = rho_drho[2]/p->r;
    E_field_eval_E(E, rho_drho, Edata);

    p_fo->B_r[j] = B[0];					  
    p_fo->B_phi[j] = B[1];				
    p_fo->B_z[j] = B[2];
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
    p_fo->index[j] = -1;
}

void fo_to_particle(particle_simd_fo* p_fo, int j, particle* p) {
    p->r = p_fo->r[j];
    p->phi = p_fo->phi[j];
    p->z = p_fo->z[j];
    p->v_r = p_fo->rdot[j];
    p->v_phi = p_fo->phidot[j] * p_fo->r[j];
    p->v_z = p_fo->zdot[j];
    p->time = p_fo->time[j];
    p->running = p_fo->running[j];
    p->endcond = p_fo->endcond[j];
    p->walltile = p_fo->walltile[j];
}

void particle_to_gc(particle* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata) {
    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, Bdata);
    real gcpos[5];
    real gamma = phys_gammaprtv(sqrt(p->v_r*p->v_r + p->v_phi*p->v_phi + p->v_z*p->v_z));
    phys_prttogc(p->mass, p->charge, p->r, p->phi, p->z, 
		 gamma*p->mass*p->v_r, gamma*p->mass*p->v_phi, gamma*p->mass*p->v_z, B_dB, gcpos);

    p_gc->r[j] = gcpos[0];
    p_gc->phi[j] = gcpos[1];
    p_gc->z[j] =gcpos[2];
    
    real B[3];
    B_field_eval_B(B, p_gc->r[j], p_gc->phi[j], p_gc->z[j], Bdata);
    gamma = phys_gammagcp(p->mass, gcpos[3], gcpos[4]);
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

void guiding_center_to_gc(particle_gc* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata) {
    real B_dB[12];

    p_gc->r[j] = p->r;
    p_gc->phi[j] = p->phi;
    p_gc->z[j] = p->z;
    
    real B[3];
    B_field_eval_B(B, p_gc->r[j], p_gc->phi[j], p_gc->z[j], Bdata);
    p_gc->vpar[j] = p->vpar;
    p_gc->mu[j] = p->mu;

    //p_gc->theta[j] = 0atan2(math_norm(ezcrossrho), math_dot(ez, rho));
    p_gc->theta[j] = p->theta; // Ill defined quantity
/*    p_gc->theta[j] = acos(math_dot(ez, rho)
                          / (math_norm(ez)*math_norm(rho))); */

    p_gc->mass[j] = p->mass;
    p_gc->charge[j] = p->charge;
    p_gc->weight[j] = p->weight;
    p_gc->time[j] = 0;
    p_gc->id[j] = p->id;
    p_gc->running[j] = p->running;
    p_gc->B_r[j] = B[0];
    p_gc->B_phi[j] = B[1];
    p_gc->B_z[j] = B[2];
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
    p->v_r = p_gc->mu[j];
    p->v_phi = p_gc->vpar[j];
    p->v_z = 0;
    p->time = p_gc->time[j];
    p->running = p_gc->running[j];
    p->endcond = p_gc->endcond[j];
    p->walltile = p_gc->walltile[j];
}

