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
                    B_field_data* Bdata){
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

    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, Bdata);

    p_fo->B_r[j]        = B_dB[0];
    p_fo->B_r_dr[j]     = B_dB[1];
    p_fo->B_r_dphi[j]   = B_dB[2];
    p_fo->B_r_dz[j]     = B_dB[3];

    p_fo->B_phi[j]      = B_dB[4];
    p_fo->B_phi_dr[j]   = B_dB[5];
    p_fo->B_phi_dphi[j] = B_dB[6];
    p_fo->B_phi_dz[j]   = B_dB[7];

    p_fo->B_z[j]        = B_dB[8];
    p_fo->B_z_dr[j]     = B_dB[9];
    p_fo->B_z_dphi[j]   = B_dB[10];
    p_fo->B_z_dz[j]     = B_dB[11];

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
    p_gc->endcond[j] = 0; 
    p_gc->walltile[j] = 0;
    p_gc->B_r[j]        = B_dB[0];
    p_gc->B_r_dr[j]     = B_dB[1];
    p_gc->B_r_dphi[j]   = B_dB[2];
    p_gc->B_r_dz[j]     = B_dB[3];

    p_gc->B_phi[j]      = B_dB[4];
    p_gc->B_phi_dr[j]   = B_dB[5];
    p_gc->B_phi_dphi[j] = B_dB[6];
    p_gc->B_phi_dz[j]   = B_dB[7];

    p_gc->B_z[j]        = B_dB[8];
    p_gc->B_z_dr[j]     = B_dB[9];
    p_gc->B_z_dphi[j]   = B_dB[10];
    p_gc->B_z_dz[j]     = B_dB[11];
    p_gc->index[j] = i;
}

void particle_to_gc_dummy(particle_simd_gc* p_gc, int j) {
    p_gc->r[j] = 1;
    p_gc->phi[j] = 1;
    p_gc->z[j] = 1;
    p_gc->vpar[j] = 1;
    p_gc->mu[j] = 1;
    p_gc->theta[j] = 1;
    p_gc->mass[j] = 1;
    p_gc->charge[j] = 1;
    p_gc->time[j] = 0;
    p_gc->weight[j] = 0;
    p_gc->id[j] = -1;
    p_gc->running[j] = 0;
    p_gc->endcond[j] = 0; 
    p_gc->walltile[j] = 0;
    p_gc->B_r[j]        = 1;
    p_gc->B_r_dr[j]     = 1;
    p_gc->B_r_dphi[j]   = 1;
    p_gc->B_r_dz[j]     = 1;

    p_gc->B_phi[j]      = 1;
    p_gc->B_phi_dr[j]   = 1;
    p_gc->B_phi_dphi[j] = 1;
    p_gc->B_phi_dz[j]   = 1;

    p_gc->B_z[j]        = 1;
    p_gc->B_z_dr[j]     = 1;
    p_gc->B_z_dphi[j]   = 1;
    p_gc->B_z_dz[j]     = 1;
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

void particle_gc_to_gc(particle_gc* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata){
    p_gc->r[j] = p->r;
    p_gc->phi[j] = p->phi;
    p_gc->z[j] = p->z;
    p_gc->theta[j] = p->theta;
    p_gc->mass[j] = p->mass;
    p_gc->charge[j] = p->charge;
    p_gc->weight[j] = p->weight;
    p_gc->time[j] = p->time;
    p_gc->id[j] = p->id;
    p_gc->running[j] = p->running;
    p_gc->endcond[j] = p->endcond;
    p_gc->walltile[j] = p->walltile;

    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, Bdata);
    p_gc->B_r[j]        = B_dB[0];
    p_gc->B_r_dr[j]     = B_dB[1];
    p_gc->B_r_dphi[j]   = B_dB[2];
    p_gc->B_r_dz[j]     = B_dB[3];

    p_gc->B_phi[j]      = B_dB[4];
    p_gc->B_phi_dr[j]   = B_dB[5];
    p_gc->B_phi_dphi[j] = B_dB[6];
    p_gc->B_phi_dz[j]   = B_dB[7];

    p_gc->B_z[j]        = B_dB[8];
    p_gc->B_z_dr[j]     = B_dB[9];
    p_gc->B_z_dphi[j]   = B_dB[10];
    p_gc->B_z_dz[j]     = B_dB[11];

    real v_para, v_perp, magnv, magnB;
    magnv = phys_Ekintovelocity(p->mass,p->energy);
    v_para = magnv*p->pitch;
    v_perp = magnv*sqrt(1 - p->pitch*p->pitch);
    p_gc->vpar[j] = v_para;

    magnB = sqrt(B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]);
    p_gc->mu[j] = phys_mu(v_para,v_perp,p->mass,magnB);

    p_gc->index[j] = i;
}

void gc_to_particle_gc(particle_simd_gc* p_gc, int j, particle_gc* p) {
    p->r = p_gc->r[j];
    p->phi = p_gc->phi[j];
    p->z = p_gc->z[j];

    real v_perp, magnB, magnv;
    magnB = sqrt(p_gc->B_r[j]*p_gc->B_r[j]
                 + p_gc->B_phi[j]*p_gc->B_phi[j]
                 + p_gc->B_z[j]*p_gc->B_z[j]);
    v_perp = phys_mutovperp(p_gc->mass[j], p_gc->vpar[j], p_gc->mu[j], magnB);
    magnv = sqrt( p_gc->vpar[j]*p_gc->vpar[j] + v_perp*v_perp);

    p->pitch = p_gc->vpar[j]/magnv;

    real gammar;
    gammar = phys_gammagcv(p_gc->mass[j],p_gc->vpar[j],p_gc->mu[j]);
    
    p->energy = p_gc->mass[j]*CONST_C2*(gammar - 1);
    
    p->theta = p_gc->theta[j];
    p->time = p_gc->time[j];
    p->running = p_gc->running[j];
    p->endcond = p_gc->endcond[j];
    p->walltile = p_gc->walltile[j];
}


/**
 * @brief Transforms particle struct into a magnetic field line simulation struct
 *
 * @param p     pointer to the particle being processed
 * @param i     index
 * @param p_ml  pointer to particle_simd_ml array
 * @param j     index of the new ml struct in the SIMD array
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void particle_to_ml(particle* p, int i, particle_simd_ml* p_ml, int j,
                    B_field_data* Bdata){
    p_ml->r[j] = p->r;
    p_ml->phi[j] = p->phi;
    p_ml->z[j] = p->z;
    p_ml->distance[j] = p->time;
    p_ml->id[j] = p->id; 
    p_ml->running[j] = p->running;
    p_ml->endcond[j] = p->endcond; 
    p_ml->walltile[j] = p->walltile;

    real B_dB[3];
    B_field_eval_B(B_dB, p->r, p->phi, p->z, Bdata);

    p_ml->B_r[j]        = B_dB[0];
    p_ml->B_r_dr[j]     = B_dB[1];
    p_ml->B_r_dphi[j]   = B_dB[2];
    p_ml->B_r_dz[j]     = B_dB[3];

    p_ml->B_phi[j]      = B_dB[4];
    p_ml->B_phi_dr[j]   = B_dB[5];
    p_ml->B_phi_dphi[j] = B_dB[6];
    p_ml->B_phi_dz[j]   = B_dB[7];

    p_ml->B_z[j]        = B_dB[8];
    p_ml->B_z_dr[j]     = B_dB[9];
    p_ml->B_z_dphi[j]   = B_dB[10];
    p_ml->B_z_dz[j]     = B_dB[11];
    p_ml->index[j] = i;
}

/**
 * @brief Makes a dummy magnetic field line simulation struct
 *
 * @param p_fo  pointer to particle_simd_ml array
 * @param j     index of the new ml struct in the SIMD array
 */
void particle_to_ml_dummy(particle_simd_ml* p_ml, int j){
    p_ml->r[j] = 1;
    p_ml->phi[j] = 1;
    p_ml->z[j] = 1;
    p_ml->distance[j] = 0;
    p_ml->id[j] = -1; 
    p_ml->running[j] = 0;
    p_ml->endcond[j] = 0; 
    p_ml->walltile[j] = 0;  
    p_ml->B_r[j]        = 1;
    p_ml->B_r_dr[j]     = 1;
    p_ml->B_r_dphi[j]   = 1;
    p_ml->B_r_dz[j]     = 1;

    p_ml->B_phi[j]      = 1;
    p_ml->B_phi_dr[j]   = 1;
    p_ml->B_phi_dphi[j] = 1;
    p_ml->B_phi_dz[j]   = 1;

    p_ml->B_z[j]        = 1;
    p_ml->B_z_dr[j]     = 1;
    p_ml->B_z_dphi[j]   = 1;
    p_ml->B_z_dz[j]     = 1;
    p_ml->index[j] = -1;
}

void ml_to_particle(particle_simd_ml* p_ml, int j, particle* p) {
    p->r = p_ml->r[j];
    p->phi = p_ml->phi[j];
    p->z = p_ml->z[j];
    p->time = p_ml->distance[j];
    p->running = p_ml->running[j];
    p->endcond = p_ml->endcond[j];
    p->walltile = p_ml->walltile[j];
}

int particle_cycle_fo(particle_queue_fo* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle) {
    for(int i = 0; i < NSIMD; i++) {
	int i_prt;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    #pragma omp critical
	    i_prt = q->next++;
	    particle_state_to_fo(q->p[i_prt], i_prt, p,i);
	    cycle[i] = 1;
	    continue;
	}

	cycle[i] = 0;
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_fo_to_state(p, i, q->p[p->index[i]]);

            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		particle_state_to_fo(q->p[i_prt], i_prt, p, i);
		cycle[i] = 1;
            }
            else {
                p->running[i] = 0;
                p->id[i] = -1;
		cycle[i] = -1;
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
        n_running += p->running[i];
    }
    
    return n_running;
}

int particle_cycle_gc(particle_queue_gc* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle) {
    for(int i = 0; i < NSIMD; i++) {
	int i_prt;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    #pragma omp critical
	    i_prt = q->next++;
	    particle_state_to_gc(q->p[i_prt], i_prt, p,i);
	    cycle[i] = 1;
	    continue;
	}

	cycle[i] = 0;
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_gc_to_state(p, i, q->p[p->index[i]]);

            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		particle_state_to_gc(q->p[i_prt], i_prt, p, i);
		cycle[i] = 1;
            }
            else {
                p->running[i] = 0;
                p->id[i] = -1;
		cycle[i] = -1;
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
        n_running += p->running[i];
    }
    
    return n_running;
}

/**
 * @brief Transforms input to state. The integer state indicates whether we are going to simulate fo (1), gc (2) or mf(3)
 */
void particle_marker_to_state(input_particle* p, int i_prt, B_field_data* Bdata, int state) {

    if(state == 1) {
	if(p[i_prt].type == input_particle_type_p) {
	    /* Particle to fo */
	    p[i_prt].type = input_particle_type_ps;

	    real r      = p[i_prt].p.r;
	    real phi    = p[i_prt].p.phi;
	    real z      = p[i_prt].p.z;
	    real rdot   = p[i_prt].p.v_r;
	    real phidot = p[i_prt].p.v_phi/r;
	    real zdot   = p[i_prt].p.v_z;
	    real mass   = p[i_prt].p.mass;
	    real charge = p[i_prt].p.charge;
	    real weight = p[i_prt].p.weight;
	    real time   = p[i_prt].p.time;
	    int id      = p[i_prt].p.id;

	    real B_dB[12];
	    B_field_eval_B_dB(B_dB, r, phi, z, Bdata);
      
	    p[i_prt].p_s.rprt       = r;     
	    p[i_prt].p_s.phiprt     = phi;      
	    p[i_prt].p_s.zprt       = z;   
	    p[i_prt].p_s.rdot       = rdot; 
	    p[i_prt].p_s.phidot     = phidot;     
	    p[i_prt].p_s.zdot       = zdot;     
	    p[i_prt].p_s.mass       = mass;      
	    p[i_prt].p_s.charge     = charge;  
	    p[i_prt].p_s.weight     = weight;    
	    p[i_prt].p_s.time       = time;     
	    p[i_prt].p_s.id         = id;    
	    p[i_prt].p_s.endcond    = 0; 
	    p[i_prt].p_s.walltile   = 0;
	    p[i_prt].p_s.B_r        = B_dB[0];     
	    p[i_prt].p_s.B_phi      = B_dB[4];     
	    p[i_prt].p_s.B_z        = B_dB[8];      
	    p[i_prt].p_s.B_r_dr     = B_dB[1];     
	    p[i_prt].p_s.B_phi_dr   = B_dB[5];   
	    p[i_prt].p_s.B_z_dr     = B_dB[9];     
	    p[i_prt].p_s.B_r_dphi   = B_dB[2];  
	    p[i_prt].p_s.B_phi_dphi = B_dB[6];  
	    p[i_prt].p_s.B_z_dphi   = B_dB[10];   
	    p[i_prt].p_s.B_r_dz     = B_dB[3];      
	    p[i_prt].p_s.B_phi_dz   = B_dB[7];   
	    p[i_prt].p_s.B_z_dz     = B_dB[11];
	}
        else if(p[i_prt].type == input_particle_type_gc) {
	    /* Guiding center to fo */
	    p[i_prt].type = input_particle_type_ps;

	    real r      = p[i_prt].p_gc.r;
	    real phi    = p[i_prt].p_gc.phi;
	    real z      = p[i_prt].p_gc.z;
	    real pitch  = p[i_prt].p_gc.pitch;
	    real energy = p[i_prt].p_gc.energy;
	    real theta  = p[i_prt].p_gc.theta;
	    real mass   = p[i_prt].p_gc.mass;
	    real charge = p[i_prt].p_gc.charge;
	    real weight = p[i_prt].p_gc.weight;
	    real time   = p[i_prt].p_gc.time;
	    integer id  = p[i_prt].p_gc.id;

	    real B_dB[12];
	    B_field_eval_B_dB(B_dB, r, phi, z, Bdata);

	    real B_norm = sqrt(B_dB[0] * B_dB[0] + B_dB[4] * B_dB[4] + B_dB[8] * B_dB[8]);

	    real gamma = 1 + energy / (mass * CONST_C2);
	    real vgc = sqrt(1 - 1 / (gamma * gamma)) * CONST_C;
	    real vpar = pitch*vgc;
	    real mu = (1 - pitch * pitch) * vgc * vgc / (B_norm);
	    
	    real prtpos[6];
	    phys_gctoprt(mass, charge, r, phi, z, vpar, mu, theta, B_dB, prtpos);
	     
	    p[i_prt].p_s.r          = r;     
	    p[i_prt].p_s.phi        = phi;      
	    p[i_prt].p_s.z          = z;   
	    p[i_prt].p_s.mu         = mu; 
	    p[i_prt].p_s.vpar       = vpar;     
	    p[i_prt].p_s.theta      = theta;     
	    p[i_prt].p_s.mass       = mass;      
	    p[i_prt].p_s.charge     = charge;  
	    p[i_prt].p_s.weight     = weight;    
	    p[i_prt].p_s.time       = time;     
	    p[i_prt].p_s.id         = id;    
	    p[i_prt].p_s.endcond    = 0; 
	    p[i_prt].p_s.walltile   = 0;
	    p[i_prt].p_s.B_r        = B_dB[0];     
	    p[i_prt].p_s.B_phi      = B_dB[4];     
	    p[i_prt].p_s.B_z        = B_dB[8];      
	    p[i_prt].p_s.B_r_dr     = B_dB[1];     
	    p[i_prt].p_s.B_phi_dr   = B_dB[5];   
	    p[i_prt].p_s.B_z_dr     = B_dB[9];     
	    p[i_prt].p_s.B_r_dphi   = B_dB[2];  
	    p[i_prt].p_s.B_phi_dphi = B_dB[6];  
	    p[i_prt].p_s.B_z_dphi   = B_dB[10];   
	    p[i_prt].p_s.B_r_dz     = B_dB[3];      
	    p[i_prt].p_s.B_phi_dz   = B_dB[7];   
	    p[i_prt].p_s.B_z_dz     = B_dB[11];
           
	}
	else if(p[i_prt].type == input_particle_type_ml) {
	    /* Magnetic field line to fo */
	    /* NOT ACCEPTABLE TRANSFORMATION */
	    // TODO generate error
	}
    }
    else if(state == 2) {
	if(p[i_prt].type == input_particle_type_p) {
	    /* Particle to gc */
	    p[i_prt].type = input_particle_type_gcs;

	    real r      = p[i_prt].p.r;
	    real phi    = p[i_prt].p.phi;
	    real z      = p[i_prt].p.z;
	    real v_r    = p[i_prt].p.v_r;
	    real v_phi  = p[i_prt].p.v_phi;
	    real v_z    = p[i_prt].p.v_z;
	    real mass   = p[i_prt].p.mass;
	    real charge = p[i_prt].p.charge;
	    real weight = p[i_prt].p.weight;
	    real time   = p[i_prt].p.time;
	    int id      = p[i_prt].p.id;

	    real B_dB[12];
	    B_field_eval_B_dB(B_dB, r, phi, z, Bdata);
	    real gcpos[5];
	    real gamma = phys_gammaprtv(sqrt( v_r * v_r + v_phi * v_phi + v_z * v_z));
	    phys_prttogc(mass, charge, r, phi, z, 
			 gamma * mass * v_r, gamma * mass * v_phi, gamma * mass * v_z, 
			 B_dB, gcpos);



	    gamma = phys_gammagcp(mass, gcpos[3], gcpos[4]);

	    p[i_prt].p_s.r          = gcpos[0];
	    p[i_prt].p_s.phi        = gcpos[1];
	    p[i_prt].p_s.z          = gcpos[2];
	    p[i_prt].p_s.vpar       = gcpos[3]/(mass*gamma);
	    p[i_prt].p_s.mu         = gcpos[4];
	    p[i_prt].p_s.theta      = 0.0; // TODO Ill defined quantity
	    p[i_prt].p_s.mass       = mass;
	    p[i_prt].p_s.charge     = charge;
	    p[i_prt].p_s.weight     = weight;
	    p[i_prt].p_s.time       = time;
	    p[i_prt].p_s.id         = id;
	    p[i_prt].p_s.endcond    = 0; 
	    p[i_prt].p_s.walltile   = 0;
	    p[i_prt].p_s.B_r        = B_dB[0];     
	    p[i_prt].p_s.B_phi      = B_dB[4];     
	    p[i_prt].p_s.B_z        = B_dB[8];      
	    p[i_prt].p_s.B_r_dr     = B_dB[1];     
	    p[i_prt].p_s.B_phi_dr   = B_dB[5];   
	    p[i_prt].p_s.B_z_dr     = B_dB[9];     
	    p[i_prt].p_s.B_r_dphi   = B_dB[2];  
	    p[i_prt].p_s.B_phi_dphi = B_dB[6];  
	    p[i_prt].p_s.B_z_dphi   = B_dB[10];   
	    p[i_prt].p_s.B_r_dz     = B_dB[3];      
	    p[i_prt].p_s.B_phi_dz   = B_dB[7];   
	    p[i_prt].p_s.B_z_dz     = B_dB[11];
	}
	else if(p[i_prt].type == input_particle_type_gc) {
	    /* Guiding center to gc */
	    p[i_prt].type = input_particle_type_gcs;

	    real r      = p[i_prt].p_gc.r;
	    real phi    = p[i_prt].p_gc.phi;
	    real z      = p[i_prt].p_gc.z;
	    real pitch  = p[i_prt].p_gc.pitch;
	    real energy = p[i_prt].p_gc.energy;
	    real theta  = p[i_prt].p_gc.theta;
	    real mass   = p[i_prt].p_gc.mass;
	    real charge = p[i_prt].p_gc.charge;
	    real weight = p[i_prt].p_gc.weight;
	    real time   = p[i_prt].p_gc.time;
	    integer id  = p[i_prt].p_gc.id;

	    real B_dB[12];
	    B_field_eval_B_dB(B_dB, r, phi, z, Bdata);

	    real B_norm = sqrt(B_dB[0] * B_dB[0] + B_dB[4] * B_dB[4] + B_dB[8] * B_dB[8]);

	    real gamma = 1 + energy / (mass * CONST_C2);
	    real vgc = sqrt(1 - 1 / (gamma * gamma)) * CONST_C;
	    real vpar = pitch*vgc;
	    real mu = (1 - pitch * pitch) * vgc * vgc / (B_norm);

	    p[i_prt].p_s.r          = r;     
	    p[i_prt].p_s.phi        = phi;      
	    p[i_prt].p_s.z          = z;   
	    p[i_prt].p_s.mu         = mu; 
	    p[i_prt].p_s.vpar       = vpar;     
	    p[i_prt].p_s.theta      = theta;     
	    p[i_prt].p_s.mass       = mass;      
	    p[i_prt].p_s.charge     = charge;  
	    p[i_prt].p_s.weight     = weight;    
	    p[i_prt].p_s.time       = time;     
	    p[i_prt].p_s.id         = id;    
	    p[i_prt].p_s.endcond    = 0; 
	    p[i_prt].p_s.walltile   = 0;
	    p[i_prt].p_s.B_r        = B_dB[0];     
	    p[i_prt].p_s.B_phi      = B_dB[4];     
	    p[i_prt].p_s.B_z        = B_dB[8];      
	    p[i_prt].p_s.B_r_dr     = B_dB[1];     
	    p[i_prt].p_s.B_phi_dr   = B_dB[5];   
	    p[i_prt].p_s.B_z_dr     = B_dB[9];     
	    p[i_prt].p_s.B_r_dphi   = B_dB[2];  
	    p[i_prt].p_s.B_phi_dphi = B_dB[6];  
	    p[i_prt].p_s.B_z_dphi   = B_dB[10];   
	    p[i_prt].p_s.B_r_dz     = B_dB[3];      
	    p[i_prt].p_s.B_phi_dz   = B_dB[7];   
	    p[i_prt].p_s.B_z_dz     = B_dB[11];
	}
	else if(p[i_prt].type == input_particle_type_ml) {
	    /* Magnetic field line to gc */
	    /* NOT ACCEPTABLE TRANSFORMATION */
	    // TODO generate error
	}

    }
    else if(state == 3) {
        if(p[i_prt].type == input_particle_type_p) {
	    
	}
	else if(p[i_prt].type == input_particle_type_gc) {
	    
	}
	else if(p[i_prt].type == input_particle_type_ml) {

	}
    }
}

void particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo, int j) {
    p_fo->r[j] = p->rprt;
    p_fo->phi[j] = p->phiprt;
    p_fo->z[j] = p->zprt;
    p_fo->rdot[j] = p->rdot;
    p_fo->phidot[j] = p->phidot;
    p_fo->zdot[j] = p->zdot;
    p_fo->mass[j] = p->mass;
    p_fo->charge[j] = p->charge;
    p_fo->weight[j] = p->weight;
    p_fo->time[j] = p->time;
    p_fo->id[j] = p->id; 
    p_fo->endcond[j] = p->endcond;

    p_fo->running[j] = 1;
    if(p->endcond) {
	p_fo->running[j] = 0;
    }
    
    p_fo->walltile[j] = p->walltile;

    p_fo->B_r[j]        = p->B_r;
    p_fo->B_r_dr[j]     = p->B_r_dr;
    p_fo->B_r_dphi[j]   = p->B_r_dphi;
    p_fo->B_r_dz[j]     = p->B_r_dz;

    p_fo->B_phi[j]      = p->B_phi;
    p_fo->B_phi_dr[j]   = p->B_phi_dr;
    p_fo->B_phi_dphi[j] = p->B_phi_dphi;
    p_fo->B_phi_dz[j]   = p->B_phi_dz;

    p_fo->B_z[j]        = p->B_z;
    p_fo->B_z_dr[j]     = p->B_z_dr;
    p_fo->B_z_dphi[j]   = p->B_z_dphi;
    p_fo->B_z_dz[j]     = p->B_z_dz;
	        
    p_fo->index[j] = i;
}

void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p) {
    p->rprt    = p_fo->r[j];
    p->phiprt  = p_fo->phi[j];
    p->zprt    = p_fo->z[j];
    p->rdot    = p_fo->rdot[j];
    p->phidot  = p_fo->phidot[j];
    p->zdot    = p_fo->zdot[j];
    p->mass    = p_fo->mass[j];
    p->charge  = p_fo->charge[j];
    p->weight  = p_fo->weight[j];
    p->time    = p_fo->time[j];
    p->id      = p_fo->id[j]; 
    p->endcond = p_fo->endcond[j];
 
    p->walltile = p_fo->walltile[j];

    p->B_r      = p_fo->B_r[j];
    p->B_r_dr   = p_fo->B_r_dr[j];
    p->B_r_dphi = p_fo->B_r_dphi[j];
    p->B_r_dz   = p_fo->B_r_dz[j];

    p->B_phi      = p_fo->B_phi[j];
    p->B_phi_dr   = p_fo->B_phi_dr[j];
    p->B_phi_dphi = p_fo->B_phi_dphi[j];
    p->B_phi_dz   = p_fo->B_phi_dz[j];

    p->B_z      = p_fo->B_z[j];
    p->B_z_dr   = p_fo->B_z_dr[j];
    p->B_z_dphi = p_fo->B_z_dphi[j];
    p->B_z_dz   = p_fo->B_z_dz[j];
}

void particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc, int j) {
    
    p_gc->r[j]        = p->r;
    p_gc->phi[j]      = p->phi;
    p_gc->z[j]        = p->z;
    p_gc->vpar[j]     = p->vpar;
    p_gc->mu[j]       = p->mu;
    p_gc->theta[j]    = p->theta;
    p_gc->mass[j]     = p->mass;
    p_gc->charge[j]   = p->charge;
    p_gc->time[j]     = p->time;
    p_gc->weight[j]   = p->weight;
    p_gc->id[j]       = p->id;
    p_gc->endcond[j]  = p->endcond; 
    p_gc->walltile[j] = p->walltile;

    p_gc->running[j] = 1;
    if(p->endcond) {
	p_gc->running[j] = 0;
    }

    p_gc->B_r[j]        = p->B_r;
    p_gc->B_r_dr[j]     = p->B_r_dr;
    p_gc->B_r_dphi[j]   = p->B_r_dphi;
    p_gc->B_r_dz[j]     = p->B_r_dz;

    p_gc->B_phi[j]      = p->B_phi;
    p_gc->B_phi_dr[j]   = p->B_phi_dr;
    p_gc->B_phi_dphi[j] = p->B_phi_dphi;
    p_gc->B_phi_dz[j]   = p->B_phi_dz;

    p_gc->B_z[j]        = p->B_z;
    p_gc->B_z_dr[j]     = p->B_z_dr;
    p_gc->B_z_dphi[j]   = p->B_z_dphi;
    p_gc->B_z_dz[j]     = p->B_z_dz;

    p_gc->index[j] = i;
}

void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p) {

    p->r       = p_gc->r[j];
    p->phi     = p_gc->phi[j];
    p->z       = p_gc->z[j];
    p->vpar    = p_gc->vpar[j];
    p->mu      = p_gc->mu[j];
    p->theta   = p_gc->theta[j];
    p->mass    = p_gc->mass[j];
    p->charge  = p_gc->charge[j];
    p->time    = p_gc->time[j];
    p->weight  = p_gc->weight[j];
    p->id      = p_gc->id[j];
    p->endcond = p_gc->endcond[j]; 
    p->walltile = p_gc->walltile[j];

    p->B_r      = p_gc->B_r[j];
    p->B_r_dr   = p_gc->B_r_dr[j];
    p->B_r_dphi = p_gc->B_r_dphi[j];
    p->B_r_dz   = p_gc->B_r_dz[j];

    p->B_phi      = p_gc->B_phi[j];
    p->B_phi_dr   = p_gc->B_phi_dr[j];
    p->B_phi_dphi = p_gc->B_phi_dphi[j];
    p->B_phi_dz   = p_gc->B_phi_dz[j];

    p->B_z      = p_gc->B_z[j];
    p->B_z_dr   = p_gc->B_z_dr[j];
    p->B_z_dphi = p_gc->B_z_dphi[j];
    p->B_z_dz   = p_gc->B_z_dz[j];
}
