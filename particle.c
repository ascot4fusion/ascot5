/**
 * @file particle.c
 * @brief Particle representations and helper functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "error.h"
#include "consts.h"
#include "math.h"
#include "physlib.h"
#include "particle.h"
#include "B_field.h"
#include "E_field.h"


/**
 * @brief Makes a dummy full-orbit simulation struct
 *
 * @param p_fo  pointer to particle_simd_fo array
 * @param j     index of the new fo struct in the SIMD array
 */
void particle_to_fo_dummy(particle_simd_fo* p_fo, int j){
    p_fo->r[j]        = 1;
    p_fo->phi[j]      = 1;
    p_fo->z[j]        = 1;
    p_fo->rdot[j]     = 1;
    p_fo->phidot[j]   = 1;
    p_fo->zdot[j]     = 1;
    p_fo->mass[j]     = 1;
    p_fo->charge[j]   = 1;
    p_fo->weight[j]   = 0;
    p_fo->time[j]     = 0;
    p_fo->id[j]       = -1; 
    p_fo->running[j]  = 0;
    p_fo->endcond[j]  = 0; 
    p_fo->walltile[j] = 0;
    p_fo->B_r[j]      = 1;					  
    p_fo->B_phi[j]    = 1;				
    p_fo->B_z[j]      = 1;	        
    p_fo->index[j]    = -1;
    p_fo->err[j]      = 0;
}

/**
 * @brief Makes a dummy guiding center simulation struct
 *
 * @param p_gc  pointer to particle_simd_gc array
 * @param j     index of the new fo struct in the SIMD array
 */
void particle_to_gc_dummy(particle_simd_gc* p_gc, int j) {
    p_gc->r[j]          = 1;
    p_gc->phi[j]        = 1;
    p_gc->z[j]          = 1;
    p_gc->vpar[j]       = 1;
    p_gc->mu[j]         = 1;
    p_gc->theta[j]      = 1;
    p_gc->mass[j]       = 1;
    p_gc->charge[j]     = 1;
    p_gc->time[j]       = 0;
    p_gc->weight[j]     = 0;
    p_gc->id[j]         = -1;
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
    p_gc->index[j]      = -1;
    p_gc->err[j]        = 0;
}

/**
 * @brief Makes a dummy magnetic field line simulation struct
 *
 * @param p_fo  pointer to particle_simd_ml array
 * @param j     index of the new ml struct in the SIMD array
 */
void particle_to_ml_dummy(particle_simd_ml* p_ml, int j){
    p_ml->r[j]          = 1;
    p_ml->phi[j]        = 1;
    p_ml->z[j]          = 1;
    p_ml->time[j]       = 0;
    p_ml->id[j]         = -1; 
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
    p_ml->index[j]      = -1;
    p_ml->err[j]        = 0;
}

/**
 * @brief Clean finished markers from the SIMD struct and 
 * fetch a fresh one from the simulation queue
 *
 * @param q  particle queue storing all simulated particles
 * @param p  SIMD structure of particles
 * @param Bdata pointer to magnetic field data
 * @param cycle NSIMD integer array indicating what was done for each marker:
 *        0 : Nothing 
 *       -1 : Finished marker replaced with a dummy (queue is empty)
 *        1 : Finished marker replaced with a fresh one
 */
int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle) {    
    int f1 = q->finished;
    int f2 = 0;

    for(int i = 0; i < NSIMD; i++) {
	int i_prt;
	int newmarker = 0;
	cycle[i] = 0;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    newmarker = 1;
	}

	/* This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_fo_to_state(p, i, q->p[p->index[i]], Bdata);
	    newmarker = 1;
	    #pragma omp critical
	    f2 = q->finished++;
	}

	while(newmarker) {
            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		a5err err = particle_state_to_fo(q->p[i_prt], i_prt, p, i, Bdata);
		if(!err) {
		    cycle[i] = 1;
		    newmarker = 0;
		}
		else {
		    #pragma omp critical
 		    f2 = q->finished++;
		}
            }
            else {
                p->running[i] = 0;
                p->id[i] = -1;
		cycle[i] = -1;
		newmarker = 0;
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
        n_running += p->running[i];
    }
    
    if(f2 > f1) {
        printf("Progress: %d/%d %le\n", f2, q->n, ((real) f2)/q->n);
    }

    return n_running;
}

/**
 * @brief Clean finished markers from the SIMD struct and 
 * fetch a fresh one from the simulation queue
 *
 * @param q  guiding center queue storing all simulated guiding centers
 * @param p  SIMD structure of guiding centers
 * @param Bdata pointer to magnetic field data
 * @param cycle NSIMD integer array indicating what was done for each marker:
 *        0 : Nothing 
 *       -1 : Finished marker replaced with a dummy (queue is empty)
 *        1 : Finished marker replaced with a fresh one
 */
int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle) {    
    int f1 = q->finished;
    int f2 = 0;

    for(int i = 0; i < NSIMD; i++) {
	int i_prt;
	int newmarker = 0;
	cycle[i] = 0;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    newmarker = 1;
	}

	/* This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_gc_to_state(p, i, q->p[p->index[i]], Bdata);
	    newmarker = 1;
	    #pragma omp critical
	    f2 = q->finished++;
	}

	while(newmarker) {
            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		a5err err = particle_state_to_gc(q->p[i_prt], i_prt, p, i, Bdata);
		if(!err) {
		    cycle[i] = 1;
		    newmarker = 0;
		}
		else {
		    #pragma omp critical
		    f2 = q->finished++;
		}
            }
            else {
                p->running[i] = 0;
                p->id[i] = -1;
		cycle[i] = -1;
		newmarker = 0;
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
        n_running += p->running[i];
    }

    if(f2 > f1) {
      printf("Progress: %d/%d %le\n", f2, q->n, ((real) f2)/q->n);
    }
    
    return n_running;
}

/**
 * @brief Clean finished markers from the SIMD struct and 
 * fetch a fresh one from the simulation queue
 *
 * @param q  field line queue storing all simulated field line tracers
 * @param p  SIMD structure of field line tracers
 * @param Bdata pointer to magnetic field data
 * @param cycle NSIMD integer array indicating what was done for each marker:
 *        0 : Nothing 
 *       -1 : Finished marker replaced with a dummy (queue is empty)
 *        1 : Finished marker replaced with a fresh one
 */
int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      B_field_data* Bdata, int* cycle) {    
    int f1 = q->finished;
    int f2 = 0;

    for(int i = 0; i < NSIMD; i++) {
	int i_prt;
	int newmarker = 0;
	cycle[i] = 0;

	/* If there are markers in queue and this position is dummy,
	 * init a marker here */
	if(p->id[i] < 0 && q->next < q->n) {
	    newmarker = 1;
	}

	/* This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
	    particle_ml_to_state(p, i, q->p[p->index[i]], Bdata);
	    newmarker = 1;
	    #pragma omp critical
	    f2 = q->finished++;
	}

	while(newmarker) {
            #pragma omp critical
            i_prt = q->next++;
            if(i_prt < q->n) {
		a5err err = particle_state_to_ml(q->p[i_prt], i_prt, p, i, Bdata);
		if(!err) {
		    cycle[i] = 1;
		    newmarker = 0;
		}
		else {
		    #pragma omp critical
		    f2 = q->finished++;
		}
            }
            else {
                p->running[i] = 0;
                p->id[i] = -1;
		cycle[i] = -1;
		newmarker = 0;
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
        n_running += p->running[i];
    }

    if(f2 > f1) {
      printf("Progress: %d/%d %le\n", f2, q->n, ((real) f2)/q->n);
    }
    
    return n_running;
}

/**
 * @brief Transforms input markers to simulation state
 *
 * @param p  marker input
 * @param ps corresponding state that was filled
 * @param Bdata pointer to magnetic field data
 */
void particle_input_to_state(input_particle* p, particle_state* ps, B_field_data* Bdata) {
    a5err err = 0;
    integer id;

    if(p->type == input_particle_type_p) {
	/* Check that input is valid */
	if(!err && ( isnan(p->p.r) || p->p.r <= 0 ))                         {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.phi))                                          {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.z))                                            {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.v_r))                                          {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.v_phi))                                        {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.v_z))                                          {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && ( math_normc(p->p.v_r,p->p.v_phi,p->p.v_z) >= CONST_C2 )) {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.time))                                         {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(p->p.charge))                                       {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && ( isnan(p->p.mass) || p->p.mass <= 0 ))                   {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && ( isnan(p->p.weight) || p->p.weight <= 0 ))               {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && p->p.id <= 0)                                             {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}

	/* Particle to state */
	p->type = input_particle_type_s;
	id = p->p.id;

	if(!err) {
	    ps->rprt   = p->p.r;     
	    ps->phiprt = p->p.phi;      
	    ps->zprt   = p->p.z;   
	    ps->rdot   = p->p.v_r; 
	    ps->phidot = p->p.v_phi/ps->rprt;     
	    ps->zdot   = p->p.v_z;
	    ps->mass   = p->p.mass;      
	    ps->charge = p->p.charge;  
	    ps->weight = p->p.weight;    
	    ps->time   = p->p.time; 
	    ps->pol    = atan2(ps->zprt-B_field_get_axis_r(Bdata),
			       ps->rprt-B_field_get_axis_z(Bdata));   
	    ps->id       = id;    
	    ps->endcond  = 0; 
	    ps->walltile = 0;
	    ps->cputime  = 0;
	}

	/* Guiding center transformation */
	real B_dB[12], r, phi, z, vpar, mu, theta, gamma, ppar, psi[1], rho[1];
	if(!err) {err = B_field_eval_B_dB(B_dB, ps->rprt, ps->phiprt, ps->zprt, Bdata);}

	if(!err) {
	    gamma = physlib_relfactorv_fo(math_normc(p->p.v_r, p->p.v_phi, p->p.v_z));
	    physlib_fo2gc(ps->mass, ps->charge, B_dB, ps->rprt, ps->phiprt, ps->zprt, 
			  gamma*ps->mass*p->p.v_r, gamma*ps->mass*p->p.v_phi, gamma*ps->mass*p->p.v_z,
			  &r, &phi, &z, &mu, &ppar, &theta);
	}
	if(!err && ( isnan(r) || r <= 0 ))  {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(phi))              {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(z))                {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(mu) || mu < 0 )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(theta))            {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

	/* Update magnetic field at gc position */
	if(!err) {err = B_field_eval_B_dB(B_dB, r, phi, z, Bdata);}
	if(!err) {err = B_field_eval_psi(psi, r, phi, z, Bdata);}
	if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}

	if(!err) {
	    gamma = physlib_relfactorp_gc(ps->mass, mu, ppar, math_normc(B_dB[0], B_dB[4], B_dB[8]));
	    vpar = ppar/(ps->mass*gamma);
	}
	if(!err && ( isnan(vpar) || vpar >= CONST_C )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

	if(!err) {
	    ps->r     = r; 
	    ps->phi   = phi;
	    ps->z     = z;
	    ps->vpar  = vpar;
	    ps->mu    = mu;
	    ps->theta = theta;

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

	    ps->err = 0;
	}
	
    }
    else if(p->type == input_particle_type_gc) {
	/* Check that input is valid */
	if(!err && ( isnan(p->p_gc.r) || p->p_gc.r <= 0 ))              {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(p->p_gc.phi))                                  {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(p->p_gc.z))                                    {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(p->p_gc.pitch) || fabs(p->p_gc.pitch) > 1 )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(p->p_gc.energy) || p->p_gc.energy <= 0 ))    {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(p->p_gc.theta))                                {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(p->p_gc.time))                                 {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && isnan(p->p_gc.charge))                               {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(p->p_gc.mass) || p->p_gc.mass <= 0 ))        {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(p->p_gc.weight) || p->p_gc.weight <= 0 ))    {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && p->p_gc.id <= 0)                                     {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

        /* Guiding center to state */
	p->type = input_particle_type_s;
	id = p->p_gc.id;

	real B_dB[12], psi[1], rho[1];
	if(!err) {err = B_field_eval_B_dB(B_dB, p->p_gc.r, p->p_gc.phi, p->p_gc.z, Bdata);}
	if(!err) {err = B_field_eval_psi(psi, p->p_gc.r, p->p_gc.phi, p->p_gc.z, Bdata);}
	if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}

	if(!err) {
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

	/* Input is in (Ekin,xi) coordinates but state needs (mu,vpar) so we need to do that
	 * transformation first. */
	real gamma, mu, vpar;
	if(!err) {
	    /* From kinetic energy we get Lorentz factor as gamma = 1 + Ekin/mc^2 */
	    gamma = 1 + p->p_gc.energy / (p->p_gc.mass * CONST_C2);

	    /* And then we can use the usual formula for Lorentz factor to get total velocity */
	    real v = sqrt(1 - 1.0 / (gamma * gamma)) * CONST_C;
       
	    /* Now we can use library functions for transformation */
	    real B_norm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
	    physlib_gc_vxi2muvpar(p->p_gc.mass, B_norm, v, p->p_gc.pitch, &mu, &vpar);
	}
	if(!err && ( isnan(mu) || mu < 0 ))            {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
	if(!err && ( isnan(vpar) || vpar >= CONST_C )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

	if(!err) {
	    ps->r          = p->p_gc.r;     
	    ps->phi        = p->p_gc.phi;      
	    ps->z          = p->p_gc.z;   
	    ps->mu         = mu; 
	    ps->vpar       = vpar;     
	    ps->theta      = p->p_gc.theta;
	    ps->mass       = p->p_gc.mass;      
	    ps->charge     = p->p_gc.charge;  
	    ps->weight     = p->p_gc.weight;    
	    ps->time       = p->p_gc.time;
	    ps->pol        = atan2(ps->z-B_field_get_axis_z(Bdata),
				   ps->r-B_field_get_axis_r(Bdata));   
	    ps->id         = id;    
	    ps->endcond    = 0; 
	    ps->walltile   = 0;
	    ps->cputime    = 0;
	}
	  
	/* Guiding center transformation to get particle coordinates */
	real rprt, phiprt, zprt, pR, pphi, pz;
	if(!err) {
	    physlib_gc2fo(ps->mass, ps->charge, B_dB,
			  ps->r, ps->phi, ps->z, ps->mu, gamma*ps->mass*ps->vpar, ps->theta,
			  &rprt, &phiprt, &zprt, &pR, &pphi, &pz);
	    gamma = physlib_relfactorp_fo(ps->mass, math_normc(pR, pphi, pz));
	}
	if(!err && ( isnan(rprt) || rprt <= 0 ))  {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(phiprt))                 {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && isnan(zprt))                   {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	if(!err && ( isnan(gamma) || gamma < 1 )) {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
	
	if(!err) {
	    ps->rprt   = rprt;
	    ps->phiprt = phiprt;
	    ps->zprt   = zprt;
	    ps->rdot   = pR/(gamma*ps->mass); 
	    ps->phidot = pphi/(gamma*ps->mass*ps->rprt);     
	    ps->zdot   = pz/(gamma*ps->mass);
           
	    ps->err = 0;
	}
    }
    else if(p->type == input_particle_type_ml) {
	/* Check that input is valid */
	if(!err && ( isnan(p->p_ml.r) || p->p_ml.r <= 0 ))           {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && isnan(p->p_ml.phi))                               {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && isnan(p->p_ml.z))                                 {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && isnan(p->p_ml.pitch))                             {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && isnan(p->p_ml.time))                              {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && ( isnan(p->p_ml.weight) || p->p_ml.weight <= 0 )) {err = error_raise(ERR_UNPHYSICAL_ML, __LINE__);}
	if(!err && p->p_ml.id <= 0)  

	/* Magnetic field line to state */
	p->type = input_particle_type_s;
	id = p->p_ml.id;

	real B_dB[12], psi[1], rho[1];
	if(!err) {err = B_field_eval_B_dB(B_dB, p->p_ml.r, p->p_ml.phi, p->p_ml.z, Bdata);}
	if(!err) {err = B_field_eval_psi(psi, p->p_ml.r, p->p_ml.phi, p->p_ml.z, Bdata);}
	if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}

	if(!err) {
	    ps->rprt       = p->p_ml.r;     
	    ps->phiprt     = p->p_ml.phi;      
	    ps->zprt       = p->p_ml.z;   
	    ps->rdot       = 0; 
	    ps->phidot     = 0;     
	    ps->zdot       = 0;
     
	    ps->mass       = 0;      
	    ps->charge     = 0;  
	    ps->weight     = p->p_ml.weight;    
	    ps->time       = p->p_ml.time;     
	    ps->id         = id; 
	    ps->pol        = atan2(p->p_ml.z-B_field_get_axis_z(Bdata), 
				   p->p_ml.r-B_field_get_axis_r(Bdata)); 
	    ps->endcond    = 0; 
	    ps->walltile   = 0;
	    ps->cputime    = 0;

	    ps->r          = p->p_ml.r;
	    ps->phi        = p->p_ml.phi;
	    ps->z          = p->p_ml.z;
	    ps->vpar       = p->p_ml.pitch >= 0;
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

	    ps->err = 0;
	}
    }

    /* If particle was rejected, ensure that these fields are filled */
    if(err) {
	ps->id      = id;    
	ps->endcond = 0;
	ps->err     = error_module(err, ERRMOD_REJECTED);
    }
}

/**
 * @brief Transform state into a fo SIMD struct
 * 
 * @param p  pointer to state being transformed
 * @param i  index of this state in the state array
 * @param p_fo SIMD structure where marker is transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo, int j, 
			  B_field_data* Bdata) {
    a5err err = p->err;

    if(!err) {
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
    }

    /* Magnetic field stored in state is for the gc position */
    real B_dB[12], psi[1], rho[1];
    if(!err) {err = B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, Bdata);}
    if(!err) {err = B_field_eval_psi(psi, p->rprt, p->phiprt, p->zprt, Bdata);}
    if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}

    if(!err) {
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

	p_fo->err[j] = 0;
    }

    return err;
}

/**
 * @brief Transform fo struct into a state on its original location 
 * in the array of states
 * 
 * @param p_fo SIMD structure being transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param p  pointer state array where marker is transformed
 * @param Bdata pointer to magnetic field data
 */
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p, 
			  B_field_data* Bdata) {
    a5err err = p_fo->err[j];
    int simerr = 0; /* Error occurred during simulation */
    if(err) {simerr = 1;}

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
    real B_dB[12], psi[1], rho[1];
    rho[0]        = p_fo->rho[j];
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

    /* Guiding center transformation */
    real ppar, gamma;
    if(!err) {
	real vR   = p->rdot;
	real vphi = p->phidot * p->rprt;
	real vz   = p->zdot;
	
	gamma = physlib_relfactorv_fo(math_normc(vR, vphi, vz));
	physlib_fo2gc(p->mass, p->charge, B_dB, p->rprt, p->phiprt, p->zprt, 
		      gamma*p->mass*vR , gamma*p->mass*vphi, gamma*p->mass*vz,
		      &p->r, &p->phi, &p->z, &p->mu, &ppar, &p->theta);
    }
    if(!err && ( isnan(p->r) || p->r <= 0 ))  {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(p->phi))                 {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(p->z))                   {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && ( isnan(p->mu) || p->mu < 0 )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(p->theta))               {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

    if(!err) {err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, Bdata);}
    if(!err) {err = B_field_eval_psi(psi, p->r, p->phi, p->z, Bdata);}
    if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}

    if(!err) {
	gamma = physlib_relfactorp_gc(p->mass, p->mu, ppar, math_normc(B_dB[0], B_dB[4], B_dB[8]));
	p->vpar = ppar/(p->mass*gamma);
    }
    if(!err && ( isnan(p->vpar) || p->vpar >= CONST_C )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}

    /* Normally magnetic field data at gc position is stored here
     * but, if gc transformation fails, field at particle position is
     * stored instead. */
    p->rho        = rho[0];

    p->B_r        = B_dB[0];
    p->B_r_dr     = B_dB[1];
    p->B_r_dphi   = B_dB[2];
    p->B_r_dz     = B_dB[3];

    p->B_phi      = B_dB[4];
    p->B_phi_dr   = B_dB[5];
    p->B_phi_dphi = B_dB[6];
    p->B_phi_dz   = B_dB[7];

    p->B_z        = B_dB[8];
    p->B_z_dr     = B_dB[9];
    p->B_z_dphi   = B_dB[10];
    p->B_z_dz     = B_dB[11];

    if(!simerr && err) {err = error_module(err, ERRMOD_STATE);}
    p->err = err;
}

/**
 * @brief Transform state into a gc SIMD struct
 * 
 * @param p  pointer to state being transformed
 * @param i  index of this state in the state array
 * @param p_gc SIMD structure where marker is transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc, int j, 
			   B_field_data* Bdata) {
    a5err err = p->err;
    
    if(!err) {
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
	p_gc->err[j] = 0;
    }

    return err;
}

/**
 * @brief Transform gc struct into a state on its original location 
 * in the array of states
 * 
 * @param p_gc SIMD structure being transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param p  pointer state array where marker is transformed
 * @param Bdata pointer to magnetic field data
 */
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p, 
			  B_field_data* Bdata) {
    a5err err  = p_gc->err[j];
    int simerr = 0; /* Error occurred during simulation */
    if(err) {simerr = 1;}

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

    real pR, pphi, pz, gamma;
    if(!err) {
	gamma = physlib_relfactorv_gc(p->mass, p->mu, p->vpar, 
					   math_normc(B_dB[0], B_dB[4], B_dB[8]));
	physlib_gc2fo(p->mass, p->charge, B_dB,
		      p->r, p->phi, p->z, p->mu, gamma*p->mass*p->vpar, p->theta,
		      &p->rprt, &p->phiprt, &p->zprt, &pR, &pphi, &pz);
	gamma = physlib_relfactorp_fo(p->mass, math_normc(pR, pphi, pz));
    }
    if(!err && ( isnan(p->rprt) || p->rprt <= 0 )) {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
    if(!err && isnan(p->phiprt))                   {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
    if(!err && isnan(p->zprt))                     {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
    if(!err && ( isnan(gamma) || gamma < 1 ))      {err = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}

    if(!err) {
	p->rdot       = pR/(gamma*p->mass); 
	p->phidot     = pphi/(gamma*p->mass*p->rprt);     
	p->zdot       = pz/(gamma*p->mass);
    }
    

    if(!simerr && err) {err = error_module(err, ERRMOD_STATE);}
    p->err = err;
}

/**
 * @brief Transform state into a ml SIMD struct
 * 
 * @param p  pointer to state being transformed
 * @param i  index of this state in the state array
 * @param p_ml SIMD structure where marker is transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml, int j, 
			   B_field_data* Bdata) {
    a5err err = p->err;
    
    if(!err) {
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

	p_ml->err[j] = 0;
    }
    
    return err;
}

/**
 * @brief Transform ml struct into a state on its original location 
 * in the array of states
 * 
 * @param p_ml SIMD structure being transformed
 * @param j  index where in the SIMD structure marker is stored
 * @param p  pointer state array where marker is transformed
 * @param Bdata pointer to magnetic field data
 */
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p, 
			  B_field_data* Bdata) {
    a5err err = p_ml->err[j];
    int simerr = 0; /* Error occurred during simulation */
    if(err) {simerr = 1;}

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
    p->err        = p_ml->err[j];

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

    if(!simerr && err) {err = error_module(err, ERRMOD_STATE);}
    p->err = err;
}

/**
 * @brief Transform fo struct into a gc struct
 * 
 * @param p_fo  fo SIMD structure being transformed
 * @param j     index where in the SIMD structure marker is stored
 * @param p_gc  gc SIMD structure where marker is transformed
 * @param Bdata pointer to magnetic field data
 */
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc, 
		      B_field_data* Bdata) {
    a5err err = p_fo->err[j];
    int simerr = 0; /* Error has already occurred */
    if(err) {simerr = 1;}
    p_gc->id[j]      = p_fo->id[j]; 
    p_gc->index[j]   = p_fo->index[j];

    real r, phi, z, gamma, ppar, mu, theta, B_dB[12];
    if(!err) {
	real Rprt   = p_fo->r[j];
	real phiprt = p_fo->phi[j];
	real zprt   = p_fo->z[j];
	real vR     = p_fo->rdot[j];
	real vphi   = p_fo->phidot[j] * p_fo->r[j];
	real vz     = p_fo->zdot[j];
	real mass   = p_fo->mass[j];
	real charge = p_fo->charge[j];

	p_gc->mass[j]     = p_fo->mass[j];
	p_gc->charge[j]   = p_fo->charge[j];
	p_gc->weight[j]   = p_fo->weight[j];
	p_gc->time[j]     = p_fo->time[j];
	p_gc->pol[j]      = p_fo->pol[j]; // This is not accurate
	p_gc->endcond[j]  = p_fo->endcond[j];
	p_gc->walltile[j] = p_fo->walltile[j];
	p_gc->cputime[j]  = p_fo->cputime[j];

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

	/* Guiding center transformation */
	gamma = physlib_relfactorv_fo(math_normc(vR, vphi, vz));
	physlib_fo2gc(mass, charge, B_dB, Rprt, phiprt, zprt, 
		      gamma*mass*vR , gamma*mass*vphi, gamma*mass*vz,
		      &r, &phi, &z, &mu, &ppar, &theta);
    }
    if(!err && ( isnan(r) || r <= 0 ))  {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(phi))              {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(z))                {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && ( isnan(mu) || mu < 0 )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    if(!err && isnan(theta))            {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    
    real psi[1], rho[1];
    if(!err) {err = B_field_eval_B_dB(B_dB, p_gc->r[j], p_gc->phi[j], p_gc->z[j], Bdata);}
    if(!err) {err = B_field_eval_psi(psi, p_gc->r[j], p_gc->phi[j], p_gc->z[j], Bdata);}
    if(!err) {err = B_field_eval_rho(rho, psi[0], Bdata);}
    if(!err) {
	gamma = physlib_relfactorp_gc(p_gc->mass[j], mu, ppar, math_normc(B_dB[0], B_dB[4], B_dB[8]));
    }
    if(!err && ( isnan(gamma) || gamma < 1 )) {err = error_raise(ERR_UNPHYSICAL_GC, __LINE__);}
    
    if(!err) {
	p_gc->r[j]          = r;
	p_gc->phi[j]        = phi;
	p_gc->z[j]          = z;
	p_gc->mu[j]         = mu;
	p_gc->theta[j]      = theta;
	p_gc->vpar[j]       = ppar/(p_gc->mass[j]*gamma);
	p_gc->rho[j]        = rho[0];

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
    }
    if(!simerr) {err = error_module(err, ERRMOD_STATE);}
    p_gc->err[j] = err;
    if(p_gc->err[j]) {
	p_gc->running[j] = 0;
	p_gc->endcond[j] = 0;
    }

    return err > 0;
}
