/**
 * @file diag_orb.c
 * @brief Functions to write particle and guiding center information. 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ascot5.h"
#include "consts.h"
#include "phys_orbit.h"
#include "particle.h"
#include "B_field.h"
#include "diag_orb.h"

void diag_orb_intervalTrigger(diag_orb_data* data, integer* particleId, real* prevWriteTime,
			      integer* id, real* time, int* store);

void diag_orb_poincareTrigger(diag_orb_data* data, int* pol, int* tor, real* kpol, real* ktor, integer* id,
			      real* ftime, real* fpol, real* fphi,
			      real* itime, real* ipol, real* iphi);

void diag_orb_lastTrigger(diag_orb_data* data, integer* id, real* ftime, real* itime, int* store);

void diag_orb_init_offload(diag_orb_offload_data* data) {


}

void diag_orb_init(diag_orb_data* data, diag_orb_offload_data* offload_data) {
    
    int i;

    /* Set the mode (interval, poincare, N-last) */
    data->mode = offload_data->mode;

    /* Interval specific input */
    data->writeInterval = offload_data->writeInterval;

    /* Poincare specific input */
    data->ntoroidalplots = offload_data->ntoroidalplots;
    for(i=0; i<data->ntoroidalplots; i++) {
	data->toroidalangles[i] = offload_data->toroidalangles[i];
    }
    data->npoloidalplots = offload_data->npoloidalplots;
    for(i=0; i<data->ntoroidalplots; i++) {
	data->poloidalangles[i] = offload_data->poloidalangles[i];
    }
    data->writeNlast = offload_data->writeNlast;
    data->writelist = NULL;
    data->size = 0;
}

void diag_orb_update_fo(integer* particleId, real* prevWriteTime, int* nextN, diag_orb_dat** Nlist,
			diag_orb_data* data, particle_simd_fo* p_f, particle_simd_fo* p_i) {

    if(data->mode == DIAG_ORB_ORBIT) {
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data,  particleId, prevWriteTime, 
				 p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->fo.id     = p_f->id[i];
		new->fo.time   = p_f->time[i];
		new->fo.r      = p_f->r[i];
		new->fo.phi    = p_f->phi[i];
		new->fo.z      = p_f->z[i];
		new->fo.rdot   = p_f->rdot[i];
		new->fo.phidot = p_f->phidot[i];
		new->fo.zdot   = p_f->zdot[i];
		new->fo.rho    = p_f->rho[i];

		new->fo.mass   = p_f->mass[i];
		new->fo.charge = p_f->charge[i];
		new->fo.weight = p_f->weight[i];
		new->fo.B_r    = p_f->B_r[i];
		new->fo.B_phi  = p_f->B_phi[i];
		new->fo.B_z    = p_f->B_z[i];
	    }
	}
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
	/* Check if marker has crossed any Poincare planes.
	 * For interpolation, kpol and ktor ( in interval [0,1]) 
	 * indicate where between initial and final state the crossing 
	 * approximately occurred. */
	int pol[NSIMD];
	int tor[NSIMD];
	real kpol[NSIMD];
	real ktor[NSIMD];
        diag_orb_poincareTrigger(data, pol, tor, kpol, ktor, p_f->id,
				 p_f->time, p_f->pol, p_f->phi,
				 p_i->time, p_i->pol, p_i->phi);

	int i;
	for(i=0; i<NSIMD; i++) {
	    if(pol[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}
	    
		new->fo.id     = p_f->id[i];
		new->fo.time   = kpol[i] * p_f->time[i]   + (1 - kpol[i]) * p_i->time[i];
		new->fo.r      = kpol[i] * p_f->r[i]      + (1 - kpol[i]) * p_i->r[i];
		new->fo.phi    = kpol[i] * p_f->phi[i]    + (1 - kpol[i]) * p_i->phi[i];
		new->fo.z      = kpol[i] * p_f->z[i]      + (1 - kpol[i]) * p_i->z[i];
		new->fo.rdot   = kpol[i] * p_f->rdot[i]   + (1 - kpol[i]) * p_i->rdot[i];
		new->fo.phidot = kpol[i] * p_f->phidot[i] + (1 - kpol[i]) * p_i->phidot[i];
		new->fo.zdot   = kpol[i] * p_f->zdot[i]   + (1 - kpol[i]) * p_i->zdot[i];
		new->fo.rho    = kpol[i] * p_f->rho[i]    + (1 - kpol[i]) * p_i->rho[i];
		new->fo.B_r    = kpol[i] * p_f->B_r[i]    + (1 - kpol[i]) * p_i->B_r[i];
		new->fo.B_phi  = kpol[i] * p_f->B_phi[i]  + (1 - kpol[i]) * p_i->B_phi[i];
		new->fo.B_z    = kpol[i] * p_f->B_z[i]    + (1 - kpol[i]) * p_i->B_z[i];

		new->fo.mass   = p_f->mass[i];
		new->fo.charge = p_f->charge[i];
		new->fo.weight = p_f->weight[i];

		new->poincareId = pol[i];
	    }

	}
	for(i=0; i<NSIMD; i++) {
	    if(tor[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->fo.id     = p_f->id[i];
		new->fo.time   = ktor[i] * p_f->time[i]   + (1 - ktor[i]) * p_i->time[i];
		new->fo.r      = ktor[i] * p_f->r[i]      + (1 - ktor[i]) * p_i->r[i];
		new->fo.phi    = ktor[i] * p_f->phi[i]    + (1 - ktor[i]) * p_i->phi[i];
		new->fo.z      = ktor[i] * p_f->z[i]      + (1 - ktor[i]) * p_i->z[i];
		new->fo.rdot   = ktor[i] * p_f->rdot[i]   + (1 - ktor[i]) * p_i->rdot[i];
		new->fo.phidot = ktor[i] * p_f->phidot[i] + (1 - ktor[i]) * p_i->phidot[i];
		new->fo.zdot   = ktor[i] * p_f->zdot[i]   + (1 - ktor[i]) * p_i->zdot[i];
		new->fo.rho    = ktor[i] * p_f->rho[i]    + (1 - ktor[i]) * p_i->rho[i];
		new->fo.B_r    = ktor[i] * p_f->B_r[i]    + (1 - ktor[i]) * p_i->B_r[i];
		new->fo.B_phi  = ktor[i] * p_f->B_phi[i]  + (1 - ktor[i]) * p_i->B_phi[i];
		new->fo.B_z    = ktor[i] * p_f->B_z[i]    + (1 - ktor[i]) * p_i->B_z[i];

		new->fo.mass   = p_f->mass[i];
		new->fo.charge = p_f->charge[i];
		new->fo.weight = p_f->weight[i];

		new->poincareId = tor[i] + DIAG_ORB_MAXPOINCARES;
	    }
	}
	
    }
    else if(data->mode == DIAG_ORB_WRITELAST) {
	
	int store[NSIMD];
	diag_orb_lastTrigger(data, p_f->id, p_f->time, p_i->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		/* See if we should start a new Nlist */
		if(p_f->id[i] != particleId[i]) {
		    nextN[i] = -1;
		    particleId[i] = p_f->id[i];
		}

		/* See whether we need a new list or can we 
		   use some that exists */
		diag_orb_dat* new;
		if(nextN[i] < (data->writeNlast-1)) {
		    nextN[i]++;
		    new = malloc(sizeof(diag_orb_dat));
		    Nlist[i*data->writeNlast + nextN[i]] = new;

		    #pragma omp critical 
		    {
			new->prev = data->writelist;
			if(data->writelist != NULL) {
			    data->writelist->next = new;
			}
			data->writelist = new;
			data->size = data->size + 1;
		    }
		}
		else {
		    nextN[i]++;
		    if(nextN[i] == 2*data->writeNlast) {
			nextN[i] = data->writeNlast;
		    }
		    
		    new = Nlist[i*data->writeNlast + (nextN[i]-data->writeNlast)];
		}

		new->fo.id     = p_f->id[i];
		new->fo.time   = p_f->time[i];
		new->fo.r      = p_f->r[i];
		new->fo.phi    = p_f->phi[i];
		new->fo.z      = p_f->z[i];
		new->fo.rdot  = p_f->rdot[i];
		new->fo.phidot = p_f->phidot[i];
		new->fo.zdot   = p_f->zdot[i];
		new->fo.rho    = p_f->rho[i];

		new->fo.mass   = p_f->mass[i];
		new->fo.charge = p_f->charge[i];
		new->fo.weight = p_f->weight[i];
		new->fo.B_r    = p_f->B_r[i];
		new->fo.B_phi  = p_f->B_phi[i];
		new->fo.B_z    = p_f->B_z[i];
	    }
	}
    }
}

void diag_orb_update_gc(integer* particleId, real* prevWriteTime, int* nextN, diag_orb_dat** Nlist,
			diag_orb_data* data, particle_simd_gc* p_f, particle_simd_gc* p_i) {

    if(data->mode == DIAG_ORB_ORBIT) {
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data, particleId, prevWriteTime,
				 p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->gc.id     = p_f->id[i];
		new->gc.time   = p_f->time[i];
		new->gc.r      = p_f->r[i];
		new->gc.phi    = p_f->phi[i];
		new->gc.z      = p_f->z[i];
		new->gc.mu     = p_f->mu[i];
		new->gc.vpar   = p_f->vpar[i];
		new->gc.theta  = p_f->theta[i];
		new->gc.rho    = p_f->rho[i];

		new->gc.mass   = p_f->mass[i];
		new->gc.charge = p_f->charge[i];
		new->gc.weight = p_f->weight[i];
		new->gc.B_r    = p_f->B_r[i];
		new->gc.B_phi  = p_f->B_phi[i];
		new->gc.B_z    = p_f->B_z[i];
	    }
	}
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
	/* Check if marker has crossed any Poincare planes.
	 * For interpolation, kpol and ktor ( in interval [0,1]) 
	 * indicate where between initial and final state the crossing 
	 * approximately occurred. */
	int pol[NSIMD];
	int tor[NSIMD];
	real kpol[NSIMD];
	real ktor[NSIMD];
        diag_orb_poincareTrigger(data, pol, tor, kpol, ktor, p_f->id,
				 p_f->time, p_f->pol, p_f->phi,
				 p_i->time, p_i->pol, p_i->phi);

	int i;
	for(i=0; i<NSIMD; i++) {
	    if(pol[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));
		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}
	    
		new->gc.id     = p_f->id[i];
		new->gc.time   = kpol[i] * p_f->time[i]  + (1 - kpol[i]) * p_i->time[i];
		new->gc.r      = kpol[i] * p_f->r[i]     + (1 - kpol[i]) * p_i->r[i];
		new->gc.phi    = kpol[i] * p_f->phi[i]   + (1 - kpol[i]) * p_i->phi[i];
		new->gc.z      = kpol[i] * p_f->z[i]     + (1 - kpol[i]) * p_i->z[i];
		new->gc.mu     = kpol[i] * p_f->mu[i]    + (1 - kpol[i]) * p_i->mu[i];
		new->gc.vpar   = kpol[i] * p_f->vpar[i]  + (1 - kpol[i]) * p_i->vpar[i];
		new->gc.rho    = kpol[i] * p_f->rho[i]   + (1 - kpol[i]) * p_i->rho[i];
		new->gc.B_r    = kpol[i] * p_f->B_r[i]   + (1 - kpol[i]) * p_i->B_r[i];
		new->gc.B_phi  = kpol[i] * p_f->B_phi[i] + (1 - kpol[i]) * p_i->B_phi[i];
		new->gc.B_z    = kpol[i] * p_f->B_z[i]   + (1 - kpol[i]) * p_i->B_z[i];

		new->gc.theta  = p_f->theta[i];
		new->gc.mass   = p_f->mass[i];
		new->gc.charge = p_f->charge[i];
		new->gc.weight = p_f->weight[i];

		new->poincareId = pol[i];
	    }

	}
	for(i=0; i<NSIMD; i++) {
	    if(tor[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));
		#pragma omp critical 
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->gc.id     = p_f->id[i];
		new->gc.time   = ktor[i] * p_f->time[i]  + (1 - ktor[i]) * p_i->time[i];
		new->gc.r      = ktor[i] * p_f->r[i]     + (1 - ktor[i]) * p_i->r[i];
		new->gc.phi    = ktor[i] * p_f->phi[i]   + (1 - ktor[i]) * p_i->phi[i];
		new->gc.z      = ktor[i] * p_f->z[i]     + (1 - ktor[i]) * p_i->z[i];
		new->gc.mu     = ktor[i] * p_f->mu[i]    + (1 - ktor[i]) * p_i->mu[i];
		new->gc.vpar   = ktor[i] * p_f->vpar[i]  + (1 - ktor[i]) * p_i->vpar[i];
		new->gc.rho    = ktor[i] * p_f->rho[i]   + (1 - ktor[i]) * p_i->rho[i];
		new->gc.B_r    = ktor[i] * p_f->B_r[i]   + (1 - ktor[i]) * p_i->B_r[i];
		new->gc.B_phi  = ktor[i] * p_f->B_phi[i] + (1 - ktor[i]) * p_i->B_phi[i];
		new->gc.B_z    = ktor[i] * p_f->B_z[i]   + (1 - ktor[i]) * p_i->B_z[i];

		new->gc.theta  = p_f->theta[i];
		new->gc.mass   = p_f->mass[i];
		new->gc.charge = p_f->charge[i];
		new->gc.weight = p_f->weight[i];

		new->poincareId = tor[i] + DIAG_ORB_MAXPOINCARES;
	    }
	}
	
    }
    else if(data->mode == DIAG_ORB_WRITELAST) {
	
	int store[NSIMD];
	diag_orb_lastTrigger(data, p_f->id, p_f->time, p_i->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		/* See if we should start a new Nlist */
		if(p_f->id[i] != particleId[i]) {
		    nextN[i] = -1;
		    particleId[i] = p_f->id[i];
		}

		/* See whether we need a new list or can we 
		   use some that exists */
		diag_orb_dat* new;
		if(nextN[i] < (data->writeNlast-1)) {
		    nextN[i]++;

		    new = malloc(sizeof(diag_orb_dat));
		    #pragma omp critical 
		    {
			Nlist[i*data->writeNlast + nextN[i]] = new;

			new->prev = data->writelist;
			if(data->writelist != NULL) {
			    data->writelist->next = new;
			}
			data->writelist = new;
			data->size = data->size + 1;
		    }
		}
		else {
		    nextN[i]++;
		    if(nextN[i] == 2*data->writeNlast) {
			nextN[i] = data->writeNlast;
		    }
		    
		    new = Nlist[i*data->writeNlast + (nextN[i]-data->writeNlast)];
		}

		new->gc.id     = p_f->id[i];
		new->gc.time   = p_f->time[i];
		new->gc.r      = p_f->r[i];
		new->gc.phi    = p_f->phi[i];
		new->gc.z      = p_f->z[i];
		new->gc.mu     = p_f->mu[i];
		new->gc.vpar   = p_f->vpar[i];
		new->gc.theta  = p_f->theta[i];
		new->gc.rho    = p_f->rho[i];

		new->gc.mass   = p_f->mass[i];
		new->gc.charge = p_f->charge[i];
		new->gc.weight = p_f->weight[i];
		new->gc.B_r    = p_f->B_r[i];
		new->gc.B_phi  = p_f->B_phi[i];
		new->gc.B_z    = p_f->B_z[i];
	    }
	}
    }
}

void diag_orb_update_ml(integer* particleId, real* prevWriteTime, int* nextN, diag_orb_dat** Nlist,
			diag_orb_data* data, particle_simd_ml* p_f, particle_simd_ml* p_i) {
    if(data->mode == DIAG_ORB_ORBIT) {
	
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data, particleId, prevWriteTime,
				 p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		#pragma omp critical
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->ml.id     = p_f->id[i];
		new->ml.time   = p_f->time[i];
		new->ml.r      = p_f->r[i];
		new->ml.phi    = p_f->phi[i];
		new->ml.z      = p_f->z[i];
		new->ml.rho    = p_f->rho[i];
		new->ml.weight = p_f->weight[i];
		new->ml.B_r    = p_f->B_r[i];
		new->ml.B_phi  = p_f->B_phi[i];
		new->ml.B_z    = p_f->B_z[i];
	    }
	}
    }
    else if(data->mode == DIAG_ORB_POINCARE) {
	/* Check if marker has crossed any Poincare planes.
	 * For interpolation, kpol and ktor ( in interval [0,1]) 
	 * indicate where between initial and final state the crossing 
	 * approximately occurred. */
	int pol[NSIMD];
	int tor[NSIMD];
	real kpol[NSIMD];
	real ktor[NSIMD];
        diag_orb_poincareTrigger(data, pol, tor, kpol, ktor, p_f->id,
				 p_f->time, p_f->pol, p_f->phi,
				 p_i->time, p_i->pol, p_i->phi);

	int i;
	
	for(i=0; i<NSIMD; i++) {
	    if(pol[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));
		
		#pragma omp critical
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;

		    data->size = data->size + 1;
		}

		new->ml.id     = p_f->id[i];
		new->ml.time   = kpol[i] * p_f->time[i]   + (1 - kpol[i]) * p_i->time[i];
		new->ml.r      = kpol[i] * p_f->r[i]      + (1 - kpol[i]) * p_i->r[i];
		new->ml.phi    = kpol[i] * p_f->phi[i]    + (1 - kpol[i]) * p_i->phi[i];
		new->ml.z      = kpol[i] * p_f->z[i]      + (1 - kpol[i]) * p_i->z[i];
		new->ml.rho    = kpol[i] * p_f->rho[i]    + (1 - kpol[i]) * p_i->rho[i];
		new->ml.weight = kpol[i] * p_f->weight[i] + (1 - kpol[i]) * p_i->weight[i];
		new->ml.B_r    = kpol[i] * p_f->B_r[i]    + (1 - kpol[i]) * p_i->B_r[i];
		new->ml.B_phi  = kpol[i] * p_f->B_phi[i]  + (1 - kpol[i]) * p_i->B_phi[i];
		new->ml.B_z    = kpol[i] * p_f->B_z[i]    + (1 - kpol[i]) * p_i->B_z[i];
		new->ml.id     = p_f->id[i];
		new->ml.time   = p_f->time[i];
		new->ml.r      = p_f->r[i];
		new->ml.phi    = p_f->phi[i];
		new->ml.z      = p_f->z[i];
		new->ml.rho    = p_f->rho[i];
		new->ml.weight = p_f->weight[i];
		new->ml.B_r    = p_f->B_r[i];
		new->ml.B_phi  = p_f->B_phi[i];
		new->ml.B_z    = p_f->B_z[i];
	    
		new->poincareId = pol[i];		
	    }

	}
	for(i=0; i<NSIMD; i++) {
	    if(tor[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));
		#pragma omp critical
		{
		    new->prev = data->writelist;
		    if(data->writelist != NULL) {
			data->writelist->next = new;
		    }
		    data->writelist = new;
		    data->size = data->size + 1;
		}

		new->ml.id     = p_f->id[i];
		new->ml.time   = ktor[i] * p_f->time[i]   + (1 - ktor[i]) * p_i->time[i];
		new->ml.r      = ktor[i] * p_f->r[i]      + (1 - ktor[i]) * p_i->r[i];
		new->ml.phi    = ktor[i] * p_f->phi[i]    + (1 - ktor[i]) * p_i->phi[i];
		new->ml.z      = ktor[i] * p_f->z[i]      + (1 - ktor[i]) * p_i->z[i];
		new->ml.rho    = ktor[i] * p_f->rho[i]    + (1 - ktor[i]) * p_i->rho[i];
		new->ml.weight = ktor[i] * p_f->weight[i] + (1 - ktor[i]) * p_i->weight[i];
		new->ml.B_r    = ktor[i] * p_f->B_r[i]    + (1 - ktor[i]) * p_i->B_r[i];
		new->ml.B_phi  = ktor[i] * p_f->B_phi[i]  + (1 - ktor[i]) * p_i->B_phi[i];
		new->ml.B_z    = ktor[i] * p_f->B_z[i]    + (1 - ktor[i]) * p_i->B_z[i];
		new->ml.id     = p_f->id[i];
		new->ml.time   = p_f->time[i];
		new->ml.r      = p_f->r[i];
		new->ml.phi    = p_f->phi[i];
		new->ml.z      = p_f->z[i];
		new->ml.rho    = p_f->rho[i];
		new->ml.weight = p_f->weight[i];
		new->ml.B_r    = p_f->B_r[i];
		new->ml.B_phi  = p_f->B_phi[i];
		new->ml.B_z    = p_f->B_z[i];

		new->poincareId = tor[i] + DIAG_ORB_MAXPOINCARES;
	    }
	}
    }
    else if(data->mode == DIAG_ORB_WRITELAST) {
	
	int store[NSIMD];
	diag_orb_lastTrigger(data, p_f->id, p_f->time, p_i->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		/* See if we should start a new Nlist */
		if(p_f->id[i] != particleId[i]) {
		    nextN[i] = -1;
		    particleId[i] = p_f->id[i];
		}

		/* See whether we need a new list or can we 
		   use some that exists */
		diag_orb_dat* new;
		if(nextN[i] < (data->writeNlast-1)) {
		    nextN[i]++;
		    new = malloc(sizeof(diag_orb_dat));


		    Nlist[i*data->writeNlast + nextN[i]] = new;
		    #pragma omp critical
		    {
			new->prev = data->writelist;
			if(data->writelist != NULL) {
			    data->writelist->next = new;
			}
			data->writelist = new;
			data->size = data->size + 1;
		    }
		}
		else {
		    nextN[i]++;
		    if(nextN[i] == 2*data->writeNlast) {
			nextN[i] = data->writeNlast;
		    }
		    
		    new = Nlist[i*data->writeNlast + (nextN[i]-data->writeNlast)];
		}

		new->ml.id     = p_f->id[i];
		new->ml.time   = p_f->time[i];
		new->ml.r      = p_f->r[i];
		new->ml.phi    = p_f->phi[i];
		new->ml.z      = p_f->z[i];
		new->ml.rho    = p_f->rho[i];
		new->ml.weight = p_f->weight[i];
		new->ml.B_r    = p_f->B_r[i];
		new->ml.B_phi  = p_f->B_phi[i];
		new->ml.B_z    = p_f->B_z[i];

		
	    }
	}
    }
}


void diag_orb_clean(diag_orb_data* data) {
    while(data->writelist != NULL) {
	diag_orb_dat* list = data->writelist->prev;
	free(data->writelist);
	data->writelist = list;
    }
    data->writelist = NULL;
    data->size = 0;
}

/**
 * @brief Check whether a marker qualifies for Interval-mode writing
 *
 * The marker qualifies if following conditions are met:
 * - It is not a dummy marker
 * - The marker has not been written before OR enough time has passed from the last write.
 *
 * Checks are done for NSIMD markers simultaneously.
 */
void diag_orb_intervalTrigger(diag_orb_data* data, integer* particleId, real* prevWriteTime,
			      integer* id, real* time, int* store) {
    
    #pragma omp simd
    for(int i= 0; i < NSIMD; i++) {
	store[i] = 0;
	if( (id[i] != -1) && (                                           // Check if dummy
		(id[i] != particleId[i]) ||                        // Check if first time this particle
		(time[i] - prevWriteTime[i] > data->writeInterval) // Check if enough time has passed from previous write
		)) {
	    
	    store[i] = 1;
	    particleId[i] = id[i];
	    prevWriteTime[i] = time[i];
	}
    }
}

/**
 * @brief Check whether a marker qualifies for Poincare-mode writing
 *
 * The marker qualifies if following conditions are met:
 * - It is not a dummy marker
 * - The time step was accepted
 * - Marker has crossed one of the specified poloidal or toroidal planes
 *
 * Marker is assumed to cross only maximum of one poloidal and one toroidal plane. Other crosses are ignored.
 *
 * Checks are done for NSIMD markers simultaneously.
 */
void diag_orb_poincareTrigger(diag_orb_data* data, int* pol, int* tor, real* kpol, real* ktor, integer* id,
			      real* ftime, real* fpol, real* fphi,
			      real* itime, real* ipol, real* iphi){
    #pragma omp simd
    for(int i= 0; i < NSIMD; i++) {
	pol[i] = -1;
	tor[i] = -1;
        if( (id[i] != -1) && (ftime[i] != itime[i]) ) { // Check marker is not dummy and time step was accepted
            // Check if the particle has crossed one of the poloidal planes
	    
	    for(int ip = 0; ip < data->npoloidalplots; ip++) {
		
		/* The phi coordinate we use is "unmodulated", i.e., it is not limited to interval [0,2pi).
		 * We can then find whether this poloidal plane was crossed by adding that plane's toroidal
		 * coordinate on marker initial and final position, and see if the division with 2pi gives the
		 * same value (no crossing) or not (marker has crossed the plane) */

		if( floor( (fphi[i] + data->poloidalangles[ip])/CONST_2PI ) !=
		    floor( (iphi[i] + data->poloidalangles[ip])/CONST_2PI ) 
		    ) {
		    pol[i] = ip;
		    kpol[i] = ( data->poloidalangles[ip] - ( iphi[i] - CONST_2PI*floor(iphi[i] / CONST_2PI) ) ) / (fphi[i] - iphi[i]);
		    break;
		}
	    }

	    // Check if the particle has crossed one of the toroidal planes
	    for(int ip = 0; ip < data->ntoroidalplots; ip++) {

		/* The pol coordinate we use is "unmodulated", i.e., it is not limited to interval [0,2pi).
		 * We can then find whether this toroidal plane was crossed by adding that plane's poloidal
		 * coordinate on marker initial and final position, and see if the division with 2pi gives the
		 * same value (no crossing) or not (marker has crossed the plane) */

		if( floor( (fpol[i] + data->toroidalangles[ip])/CONST_2PI ) !=
		    floor( (ipol[i] + data->toroidalangles[ip])/CONST_2PI ) 
		    ) {
		    tor[i] = ip;
		    ktor[i] = ( data->toroidalangles[ip] - ( ipol[i] - CONST_2PI*floor(ipol[i] / CONST_2PI) ) ) / (fpol[i] - ipol[i]);
		    break;
		}
	    }
	}
    }


}

/**
 * @brief Check whether a marker qualifies for N-last-mode writing
 *
 * The marker qualifies if following conditions are met:
 * - It is not a dummy marker
 * - The time step was accepted
 *
 * Checks are done for NSIMD markers simultaneously.
 */
void diag_orb_lastTrigger(diag_orb_data* data, integer* id, real* ftime, real* itime, int* store) {
    #pragma omp simd
    for(int i= 0; i < NSIMD; i++) {
	store[i] = 0;
	if( (id[i] != -1) &&               // Check if dummy
		(ftime[i] != itime[i]) ) { // Check if time step was accepted
		store[i] = 1;
	}
    }
}
