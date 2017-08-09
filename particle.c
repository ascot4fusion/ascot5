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
    p_ml->time[j] = p->time;
    p_ml->id[j] = p->id; 

    real B_dB[12];
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
    p_ml->time[j] = 0;
    p_ml->id[j] = -1; 
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
    p->time = p_ml->time[j];
}

int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle) {
    for(int i = 0; i < NSIMD; i++) {
	int i_prt;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    #pragma omp critical
	    i_prt = q->next++;
	    particle_state_to_fo(q->p[i_prt], i_prt, p, i, Bdata);
	    cycle[i] = 1;
	    continue;
	}

	cycle[i] = 0;
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_fo_to_state(p, i, q->p[p->index[i]], Bdata);

            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		particle_state_to_fo(q->p[i_prt], i_prt, p, i, Bdata);
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

int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle) {
    for(int i = 0; i < NSIMD; i++) {
	int i_prt;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    #pragma omp critical
	    i_prt = q->next++;
	    particle_state_to_gc(q->p[i_prt], i_prt, p, i, Bdata);
	    cycle[i] = 1;
	    continue;
	}

	cycle[i] = 0;
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_gc_to_state(p, i, q->p[p->index[i]], Bdata);

            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		particle_state_to_gc(q->p[i_prt], i_prt, p, i, Bdata);
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

int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      B_field_data* Bdata, int* cycle) {
    for(int i = 0; i < NSIMD; i++) {
	int i_prt;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    #pragma omp critical
	    i_prt = q->next++;
	    particle_state_to_ml(q->p[i_prt], i_prt, p, i, Bdata);
	    cycle[i] = 1;
	    continue;
	}

	cycle[i] = 0;
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_ml_to_state(p, i, q->p[p->index[i]], Bdata);

            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		particle_state_to_ml(q->p[i_prt], i_prt, p, i, Bdata);
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
 * @brief Transforms input to state
 */
void particle_input_to_state(input_particle* p, particle_state* ps, B_field_data* Bdata) {

    if(p->type == input_particle_type_p) {
	/* Particle to state */
	p->type = input_particle_type_s;

	real rprt   = p->p.r;
	real phiprt = p->p.phi;
	real zprt   = p->p.z;
	real rdot   = p->p.v_r;
	real phidot = p->p.v_phi/rprt;
	real zdot   = p->p.v_z;
	real mass   = p->p.mass;
	real charge = p->p.charge;
	real weight = p->p.weight;
	real time   = p->p.time;
	int id      = p->p.id;

	ps->rprt       = rprt;     
	ps->phiprt     = phiprt;      
	ps->zprt       = zprt;   
	ps->rdot       = rdot; 
	ps->phidot     = phidot;     
	ps->zdot       = zdot;
     
	ps->mass       = mass;      
	ps->charge     = charge;  
	ps->weight     = weight;    
	ps->time       = time; 
	ps->pol        = atan2(zprt-B_field_get_axis_r(Bdata),
					rprt-B_field_get_axis_z(Bdata));   
	ps->id         = id;    
	ps->endcond    = 0; 
	ps->walltile   = 0;
	ps->cputime    = 0;

	/* Guiding center transformation */
	real B_dB[12];
	B_field_eval_B_dB(B_dB, rprt, phiprt, zprt, Bdata);
	real gcpos[6];
	real gamma = phys_gammaprtv(sqrt( rdot * rdot + phidot*phidot * rprt*rprt + zdot * zdot));
	phys_prttogc(mass, charge, rprt, phiprt, zprt, 
		     gamma * mass * rdot, gamma * mass * phidot * rprt, gamma * mass * zdot, 
		     B_dB, gcpos);
        
	B_field_eval_B_dB(B_dB, gcpos[0], gcpos[1], gcpos[2], Bdata);
	gamma = phys_gammagcp(mass, gcpos[3], gcpos[4]);

	ps->r          = gcpos[0];
	ps->phi        = gcpos[1];
	ps->z          = gcpos[2];
	ps->vpar       = gcpos[3]/(mass*gamma);
	ps->mu         = gcpos[4];
	ps->theta      = gcpos[5];

	real psi[1];
	real rho[1];
	B_field_eval_psi(psi, gcpos[0], gcpos[1], gcpos[2], Bdata);
	B_field_eval_rho(rho, psi[0], Bdata);

	ps->rho        = rho[0];
	ps->B_r        = B_dB[0];     
	ps->B_phi      = B_dB[4];     
	ps->B_z        = B_dB[8];      
	ps->B_r_dr     = B_dB[1];     
	ps->B_phi_dr   = B_dB[5];   
	ps->B_z_dr     = B_dB[9];     
	ps->B_r_dphi   = B_dB[2];  
	ps->B_phi_dphi = B_dB[6];  
	ps->B_z_dphi   = B_dB[10];   
	ps->B_r_dz     = B_dB[3];      
	ps->B_phi_dz   = B_dB[7];   
	ps->B_z_dz     = B_dB[11];
    }
    else if(p->type == input_particle_type_gc) {
        /* Guiding center to state */
	p->type = input_particle_type_s;

	real r      = p->p_gc.r;
	real phi    = p->p_gc.phi;
	real z      = p->p_gc.z;
	real pitch  = p->p_gc.pitch;
	real energy = p->p_gc.energy;
	real theta  = p->p_gc.theta;
	real mass   = p->p_gc.mass;
	real charge = p->p_gc.charge;
	real weight = p->p_gc.weight;
	real time   = p->p_gc.time;
	integer id  = p->p_gc.id;

	real B_dB[12];
	B_field_eval_B_dB(B_dB, r, phi, z, Bdata);

	real B_norm = sqrt(B_dB[0] * B_dB[0] + B_dB[4] * B_dB[4] + B_dB[8] * B_dB[8]);

	real gamma = 1 + energy / (mass * CONST_C2);
	real vgc = sqrt(1 - 1 / (gamma * gamma)) * CONST_C;
	real vpar = pitch*vgc;
	real mu = (1 - pitch * pitch) * mass * vgc * vgc / (2 * B_norm);

	real prtpos[6];
	phys_gctoprt(mass, charge, r, phi, z, vpar, mu, theta, B_dB, prtpos);
	     
	ps->rprt       = prtpos[0];     
	ps->phiprt     = prtpos[1];      
	ps->zprt       = prtpos[2];   
	ps->rdot       = prtpos[3]/mass; 
	ps->phidot     = (prtpos[4]/mass)/ps->rprt;     
	ps->zdot       = prtpos[5]/mass;

	ps->r          = r;     
	ps->phi        = phi;      
	ps->z          = z;   
	ps->mu         = mu; 
	ps->vpar       = vpar;     
	ps->theta      = theta;
   
	ps->mass       = mass;      
	ps->charge     = charge;  
	ps->weight     = weight;    
	ps->time       = time; 
	ps->pol        = atan2(z-B_field_get_axis_z(Bdata),
					r-B_field_get_axis_r(Bdata));   
	ps->id         = id;    
	ps->endcond    = 0; 
	ps->walltile   = 0;
	ps->cputime    = 0;

	real psi[1];
	real rho[1];
	B_field_eval_psi(psi, r, phi, z, Bdata);
	B_field_eval_rho(rho, psi[0], Bdata);

	ps->rho        = rho[0];
	ps->B_r        = B_dB[0];
	ps->B_phi      = B_dB[4];     
	ps->B_z        = B_dB[8];      
	ps->B_r_dr     = B_dB[1];     
	ps->B_phi_dr   = B_dB[5];   
	ps->B_z_dr     = B_dB[9];     
	ps->B_r_dphi   = B_dB[2];  
	ps->B_phi_dphi = B_dB[6];  
	ps->B_z_dphi   = B_dB[10];   
	ps->B_r_dz     = B_dB[3];      
	ps->B_phi_dz   = B_dB[7];   
	ps->B_z_dz     = B_dB[11];
           
    }
    else if(p->type == input_particle_type_ml) {
	/* Magnetic field line to state */
	p->type = input_particle_type_s;

	real r      = p->p_ml.r;
	real phi    = p->p_ml.phi;
	real z      = p->p_ml.z;
	real pitch  = p->p_ml.pitch;
	real time   = p->p_ml.time;
	real weight = p->p_ml.weight;
	integer id  = p->p_ml.id;

	ps->rprt       = r;     
	ps->phiprt     = phi;      
	ps->zprt       = z;   
	ps->rdot       = 0; 
	ps->phidot     = 0;     
	ps->zdot       = 0;
     
	ps->mass       = 0;      
	ps->charge     = 0;  
	ps->weight     = weight;    
	ps->time       = time;     
	ps->id         = id; 
	ps->pol        = atan2(z-B_field_get_axis_z(Bdata), 
					r-B_field_get_axis_r(Bdata)); 
	ps->endcond    = 0; 
	ps->walltile   = 0;
	ps->cputime    = 0;

	real B_dB[12];
	B_field_eval_B_dB(B_dB, r, phi, z, Bdata);
	real psi[1];
	real rho[1];
	B_field_eval_psi(psi, r, phi, z, Bdata);
	B_field_eval_rho(rho, psi[0], Bdata);

	ps->r          = r;
	ps->phi        = phi;
	ps->z          = z;
	ps->vpar       = pitch >= 0;
	ps->mu         = 0;
	ps->theta      = 0;

	ps->rho        = rho[0]; 
	ps->B_r        = B_dB[0];     
	ps->B_phi      = B_dB[4];     
	ps->B_z        = B_dB[8];      
	ps->B_r_dr     = B_dB[1];     
	ps->B_phi_dr   = B_dB[5];   
	ps->B_z_dr     = B_dB[9];     
	ps->B_r_dphi   = B_dB[2];  
	ps->B_phi_dphi = B_dB[6];  
	ps->B_z_dphi   = B_dB[10];   
	ps->B_r_dz     = B_dB[3];      
	ps->B_phi_dz   = B_dB[7];   
	ps->B_z_dz     = B_dB[11];
    }
}

void particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo, int j, 
			  B_field_data* Bdata) {
    p_fo->r[j]          = p->rprt;
    p_fo->phi[j]        = p->phiprt;
    p_fo->z[j]          = p->zprt;
    p_fo->rdot[j]       = p->rdot;
    p_fo->phidot[j]     = p->phidot;
    p_fo->zdot[j]       = p->zdot;

    p_fo->mass[j]       = p->mass;
    p_fo->charge[j]     = p->charge;
    p_fo->weight[j]     = p->weight;
    p_fo->time[j]       = p->time;
    p_fo->pol[j]        = p->pol; 
    p_fo->id[j]         = p->id; 
    p_fo->endcond[j]    = p->endcond;
    p_fo->walltile[j]   = p->walltile;

    /* Magnetic field stored in state is for the gc position */
    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, Bdata);

    real psi[1];
    real rho[1];
    B_field_eval_psi(psi, p->rprt, p->phiprt, p->zprt, Bdata);
    B_field_eval_rho(rho, psi[0], Bdata);

    p_fo->rho[j]        = rho[0];

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
	        
    p_fo->running[j] = 1;
    if(p->endcond) {
	p_fo->running[j] = 0;
    }
    p_fo->cputime[j] = p->cputime;
    p_fo->index[j]   = i;
}

void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p, 
			  B_field_data* Bdata) {
    p->rprt       = p_fo->r[j];
    p->phiprt     = p_fo->phi[j];
    p->zprt       = p_fo->z[j];
    p->rdot       = p_fo->rdot[j];
    p->phidot     = p_fo->phidot[j];
    p->zdot       = p_fo->zdot[j];

    p->mass       = p_fo->mass[j];
    p->charge     = p_fo->charge[j];
    p->weight     = p_fo->weight[j];
    p->time       = p_fo->time[j];
    p->pol        = p_fo->pol[j];
    p->id         = p_fo->id[j]; 
    p->endcond    = p_fo->endcond[j];
    p->walltile   = p_fo->walltile[j];
    p->cputime    = p_fo->cputime[j];

    /* Particle to guiding center */
    real B_dB[12];
    B_dB[0]       = p_fo->B_r[j];
    B_dB[1]       = p_fo->B_r_dr[j];
    B_dB[2]       = p_fo->B_r_dphi[j];
    B_dB[3]       = p_fo->B_r_dz[j];
    B_dB[4]       = p_fo->B_phi[j];
    B_dB[5]       = p_fo->B_phi_dr[j];
    B_dB[6]       = p_fo->B_phi_dphi[j];
    B_dB[7]       = p_fo->B_phi_dz[j];
    B_dB[8]       = p_fo->B_z[j];
    B_dB[9]       = p_fo->B_z_dr[j];
    B_dB[10]      = p_fo->B_z_dphi[j];
    B_dB[11]      = p_fo->B_z_dz[j];

    real gcpos[6];
    real gamma = phys_gammaprtv(sqrt( p->rdot * p->rdot + p->phidot*p->phidot * p->rprt*p->rprt + p->zdot * p->zdot));
    phys_prttogc(p->mass, p->charge, p->rprt, p->phiprt, p->zprt, 
		 gamma * p->mass * p->rdot, gamma * p->mass * p->phidot * p->rprt, gamma * p->mass * p->zdot, 
		 B_dB, gcpos);
    gamma = phys_gammagcp(p->mass, gcpos[3], gcpos[4]);

    p->r          = gcpos[0];
    p->phi        = gcpos[1];
    p->z          = gcpos[2];
    p->vpar       = gcpos[3]/(p->mass*gamma);
    p->mu         = gcpos[4];
    p->theta      = gcpos[5];

    /* Magnetic field stored in state is for the gc position */
    B_field_eval_B_dB(B_dB, gcpos[0], gcpos[1], gcpos[2], Bdata);

    real psi[1];
    real rho[1];
    B_field_eval_psi(psi, gcpos[0], gcpos[1], gcpos[2], Bdata);
    B_field_eval_rho(rho, psi[0], Bdata);

    p->rho        = rho[0];

    p->B_r        = p_fo->B_r[j];
    p->B_r_dr     = p_fo->B_r_dr[j];
    p->B_r_dphi   = p_fo->B_r_dphi[j];
    p->B_r_dz     = p_fo->B_r_dz[j];

    p->B_phi      = p_fo->B_phi[j];
    p->B_phi_dr   = p_fo->B_phi_dr[j];
    p->B_phi_dphi = p_fo->B_phi_dphi[j];
    p->B_phi_dz   = p_fo->B_phi_dz[j];

    p->B_z        = p_fo->B_z[j];
    p->B_z_dr     = p_fo->B_z_dr[j];
    p->B_z_dphi   = p_fo->B_z_dphi[j];
    p->B_z_dz     = p_fo->B_z_dz[j];
}

void particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc, int j, 
			  B_field_data* Bdata) {
    
    p_gc->r[j]          = p->r;
    p_gc->phi[j]        = p->phi;
    p_gc->z[j]          = p->z;
    p_gc->vpar[j]       = p->vpar;
    p_gc->mu[j]         = p->mu;
    p_gc->theta[j]      = p->theta;

    p_gc->mass[j]       = p->mass;
    p_gc->charge[j]     = p->charge;
    p_gc->time[j]       = p->time;
    p_gc->weight[j]     = p->weight;
    p_gc->rho[j]        = p->rho;
    p_gc->pol[j]        = p->pol;
    p_gc->id[j]         = p->id;
    p_gc->endcond[j]    = p->endcond; 
    p_gc->walltile[j]   = p->walltile;

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

    p_gc->running[j] = 1;
    if(p->endcond) {
	p_gc->running[j] = 0;
    }
    p_gc->cputime[j] = p->cputime;
    p_gc->index[j]   = i;
}

void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p, 
			  B_field_data* Bdata) {

    p->r          = p_gc->r[j];
    p->phi        = p_gc->phi[j];
    p->z          = p_gc->z[j];
    p->vpar       = p_gc->vpar[j];
    p->mu         = p_gc->mu[j];
    p->theta      = p_gc->theta[j];

    p->mass       = p_gc->mass[j];
    p->charge     = p_gc->charge[j];
    p->time       = p_gc->time[j];
    p->weight     = p_gc->weight[j];
    p->id         = p_gc->id[j];
    p->cputime    = p_gc->cputime[j];
    p->rho        = p_gc->rho[j];
    p->pol        = p_gc->pol[j];
    p->endcond    = p_gc->endcond[j]; 
    p->walltile   = p_gc->walltile[j];

    p->B_r        = p_gc->B_r[j];
    p->B_r_dr     = p_gc->B_r_dr[j];
    p->B_r_dphi   = p_gc->B_r_dphi[j];
    p->B_r_dz     = p_gc->B_r_dz[j];

    p->B_phi      = p_gc->B_phi[j];
    p->B_phi_dr   = p_gc->B_phi_dr[j];
    p->B_phi_dphi = p_gc->B_phi_dphi[j];
    p->B_phi_dz   = p_gc->B_phi_dz[j];

    p->B_z        = p_gc->B_z[j];
    p->B_z_dr     = p_gc->B_z_dr[j];
    p->B_z_dphi   = p_gc->B_z_dphi[j];
    p->B_z_dz     = p_gc->B_z_dz[j];
    
    /* Guiding center to particle transformation */
    real B_dB[12];
    B_dB[0]       = p->B_r;
    B_dB[1]       = p->B_r_dr;
    B_dB[2]       = p->B_r_dphi;
    B_dB[3]       = p->B_r_dz;
    B_dB[4]       = p->B_phi;
    B_dB[5]       = p->B_phi_dr;
    B_dB[6]       = p->B_phi_dphi;
    B_dB[7]       = p->B_phi_dz;
    B_dB[8]       = p->B_z;
    B_dB[9]       = p->B_z_dr;
    B_dB[10]      = p->B_z_dphi;
    B_dB[11]      = p->B_z_dz;
    real prtpos[6];
    phys_gctoprt(p->mass, p->charge, p->r, p->phi, p->z, p->vpar, p->mu, p->theta, B_dB, prtpos);

    p->rprt       = prtpos[0];     
    p->phiprt     = prtpos[1];      
    p->zprt       = prtpos[2];   
    p->rdot       = prtpos[3]/p->mass; 
    p->phidot     = (prtpos[4]/p->mass)/p->rprt;    ;     
    p->zdot       = prtpos[5]/p->mass;
}

void particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml, int j, 
			  B_field_data* Bdata) {
    
    p_ml->r[j]          = p->r;
    p_ml->phi[j]        = p->phi;
    p_ml->z[j]          = p->z;

    p_ml->pitch[j]      = 2*(p->vpar >= 0) - 1.0;
    p_ml->time[j]       = p->time;
    p_ml->weight[j]     = p->weight;
    p_ml->id[j]         = p->id;
    p_ml->cputime[j]    = p->cputime;
    p_ml->rho[j]        = p->rho;
    p_ml->pol[j]        = p->pol;
    p_ml->endcond[j]    = p->endcond; 
    p_ml->walltile[j]   = p->walltile;

    p_ml->B_r[j]        = p->B_r;
    p_ml->B_r_dr[j]     = p->B_r_dr;
    p_ml->B_r_dphi[j]   = p->B_r_dphi;
    p_ml->B_r_dz[j]     = p->B_r_dz;

    p_ml->B_phi[j]      = p->B_phi;
    p_ml->B_phi_dr[j]   = p->B_phi_dr;
    p_ml->B_phi_dphi[j] = p->B_phi_dphi;
    p_ml->B_phi_dz[j]   = p->B_phi_dz;

    p_ml->B_z[j]        = p->B_z;
    p_ml->B_z_dr[j]     = p->B_z_dr;
    p_ml->B_z_dphi[j]   = p->B_z_dphi;
    p_ml->B_z_dz[j]     = p->B_z_dz;

    p_ml->running[j] = 1;
    if(p->endcond) {
	p_ml->running[j] = 0;
    }
    p_ml->index[j] = i;
    
}

void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p, 
			  B_field_data* Bdata) {

    p->rprt       = p_ml->r[j];
    p->phiprt     = p_ml->phi[j];
    p->zprt       = p_ml->z[j];
    p->rdot       = 0;
    p->phidot     = 0;
    p->zdot       = 0;

    p->r          = p_ml->r[j];
    p->phi        = p_ml->phi[j];
    p->z          = p_ml->z[j];
    p->vpar       = p_ml->pitch[j];
    p->mu         = 0;
    p->theta      = 0;
    p->mass       = 0;
    p->charge     = 0;
    p->time       = p_ml->time[j];
    p->weight     = p_ml->weight[j];
    p->id         = p_ml->id[j];
    p->cputime    = p_ml->cputime[j]; 
    p->rho        = p_ml->rho[j]; 
    p->pol        = p_ml->pol[j];
    p->endcond    = p_ml->endcond[j]; 
    p->walltile   = p_ml->walltile[j];

    p->B_r        = p_ml->B_r[j];
    p->B_r_dr     = p_ml->B_r_dr[j];
    p->B_r_dphi   = p_ml->B_r_dphi[j];
    p->B_r_dz     = p_ml->B_r_dz[j];

    p->B_phi      = p_ml->B_phi[j];
    p->B_phi_dr   = p_ml->B_phi_dr[j];
    p->B_phi_dphi = p_ml->B_phi_dphi[j];
    p->B_phi_dz   = p_ml->B_phi_dz[j];

    p->B_z        = p_ml->B_z[j];
    p->B_z_dr     = p_ml->B_z_dr[j];
    p->B_z_dphi   = p_ml->B_z_dphi[j];
    p->B_z_dz     = p_ml->B_z_dz[j];
}
