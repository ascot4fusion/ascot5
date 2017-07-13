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

void diag_orb_intervalTrigger(diag_orb_data* data, integer* id, real* time, int* ok);

void diag_orb_poincareTrigger(diag_orb_data* data, int* pol, int* tor, real* kpol, real* ktor, integer* id,
			      real* ftime, real* fr, real* fphi, real* fz,
			      real* itime, real* ir, real* iphi, real* iz);

void diag_orb_lastTrigger(diag_orb_data* data, integer* id, real* time, int* ok);

void diag_orb_init_offload(diag_orb_offload_data* data) {


}

void diag_orb_init(diag_orb_data* data, diag_orb_offload_data* offload_data) {
    data->writeInterval = offload_data->writeInterval;
    data->mode = offload_data->mode;

    int i;
    data->ntoroidalplots = offload_data->ntoroidalplots;
    for(i=0; i<data->ntoroidalplots; i++) {
	data->toroidalangles[i] = offload_data->toroidalangles[i];
    }
    data->npoloidalplots = offload_data->npoloidalplots;
    for(i=0; i<data->ntoroidalplots; i++) {
	data->poloidalangles[i] = offload_data->poloidalangles[i];
    }

    for(i=0; i < NSIMD; i++) {
	data->particleId[i] = -2;
    }

    if(data->mode == DIAG_ORB_ORBIT) {
	data->writelist = NULL;
    }
    if(data->mode == DIAG_ORB_POINCARE) {
	data->poincarelist = NULL;
    }
    if(data->mode == DIAG_ORB_WRITELAST) {
	
    }
    data->size = 0;
}

void diag_orb_update_gc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data) {


    if(data->mode == DIAG_ORB_ORBIT) {
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data, p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->gc.id     = p_f->id[i];
		new->gc.time   = p_f->time[i];
		new->gc.r      = p_f->r[i];
		new->gc.phi    = p_f->phi[i];
		new->gc.z      = p_f->z[i];
		new->gc.mu     = p_f->mu[i];
		new->gc.vpar   = p_f->vpar[i];
		new->gc.theta  = p_f->theta[i];

		new->gc.mass   = p_f->mass[i];
		new->gc.charge = p_f->charge[i];
		new->gc.weight = p_f->weight[i];
		new->gc.B_r    = p_f->B_r[i];
		new->gc.B_phi  = p_f->B_phi[i];
		new->gc.B_z    = p_f->B_z[i];

		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
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
				 p_f->time, p_f->r, p_f->phi, p_f->z,
				 p_i->time, p_i->r, p_i->phi, p_i->z);

	int i;
	for(i=0; i<NSIMD; i++) {
	    if(pol[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->gc.time = kpol[i] * p_f->time[i] + (1 - kpol[i]) * p_i->time[i];
		new->gc.r    = kpol[i] * p_f->r[i]    + (1 - kpol[i]) * p_i->r[i];
		new->gc.phi  = kpol[i] * p_f->phi[i]  + (1 - kpol[i]) * p_i->z[i];
		new->gc.z    = kpol[i] * p_f->z[i]    + (1 - kpol[i]) * p_i->phi[i];
		new->gc.mu   = kpol[i] * p_f->mu[i]   + (1 - kpol[i]) * p_i->mu[i];
		new->gc.vpar = kpol[i] * p_f->vpar[i] + (1 - kpol[i]) * p_i->vpar[i];
	    
		new->poincareId = pol[i];
		
		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
	    }

	}
	for(i=0; i<NSIMD; i++) {
	    if(tor[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->gc.time = ktor[i] * p_f->time[i] + (1 - ktor[i]) * p_i->time[i];
		new->gc.r    = ktor[i] * p_f->r[i]    + (1 - ktor[i]) * p_i->r[i];
		new->gc.phi  = ktor[i] * p_f->phi[i]  + (1 - ktor[i]) * p_i->z[i];
		new->gc.z    = ktor[i] * p_f->z[i]    + (1 - ktor[i]) * p_i->phi[i];
		new->gc.mu   = ktor[i] * p_f->mu[i]   + (1 - ktor[i]) * p_i->mu[i];
		new->gc.vpar = ktor[i] * p_f->vpar[i] + (1 - ktor[i]) * p_i->vpar[i];
	    
		new->poincareId = tor[i];

		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
	    }
	}
	
    }
    else if(data->mode == DIAG_ORB_WRITELAST) {
	//diag_orb_lastTrigger(data, p_f->id, p_f->time, p_i->time, ok);
    }
    

}

void diag_orb_update_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data) {
    if(data->mode == DIAG_ORB_ORBIT) {
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data, p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->fo.id     = p_f->id[i];
		new->fo.time   = p_f->time[i];
		new->fo.r      = p_f->r[i];
		new->fo.phi    = p_f->phi[i];
		new->fo.z      = p_f->z[i];
		new->fo.rdot   = p_f->rdot[i];
		new->fo.phidot = p_f->phidot[i];
		new->fo.zdot   = p_f->zdot[i];

		new->fo.mass   = p_f->mass[i];
		new->fo.charge = p_f->charge[i];
		new->fo.weight = p_f->weight[i];
		new->fo.B_r    = p_f->B_r[i];
		new->fo.B_phi  = p_f->B_phi[i];
		new->fo.B_z    = p_f->B_z[i];

		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
	    }
	}
    }
}

void diag_orb_update_ml(particle_simd_ml* p_f, particle_simd_ml* p_i, diag_orb_data* data) {
    if(data->mode == DIAG_ORB_ORBIT) {
	
	/* Check first whether a marker should be stored */
	int store[NSIMD]; 
        diag_orb_intervalTrigger(data, p_f->id, p_f->time, store);

	/* Store marker data */
	for(int i= 0; i < NSIMD; i++) {
	    if(store[i]) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->ml.id     = p_f->id[i];
		new->ml.time   = p_f->time[i];
		new->ml.r      = p_f->r[i];
		new->ml.phi    = p_f->phi[i];
		new->ml.z      = p_f->z[i];

		new->ml.weight = p_f->weight[i];
		new->ml.B_r    = p_f->B_r[i];
		new->ml.B_phi  = p_f->B_phi[i];
		new->ml.B_z    = p_f->B_z[i];

		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
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
				 p_f->time, p_f->r, p_f->phi, p_f->z,
				 p_i->time, p_i->r, p_i->phi, p_i->z);

	int i;
	
	for(i=0; i<NSIMD; i++) {
	    //printf("%d\n",pol[i]);
	    if(pol[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->ml.time = kpol[i] * p_f->time[i] + (1 - kpol[i]) * p_i->time[i];
		new->ml.r    = kpol[i] * p_f->r[i] + (1 - kpol[i]) * p_i->r[i];
		new->ml.phi  = kpol[i] * p_f->phi[i]  + (1 - kpol[i]) * p_i->z[i];
		new->ml.z    = kpol[i] * p_f->z[i] + (1 - kpol[i]) * p_i->z[i];
	    
		new->poincareId = pol[i];
		
		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
	    }

	}
	for(i=0; i<NSIMD; i++) {
	    if(tor[i] > -1) {
		diag_orb_dat* new = malloc(sizeof(diag_orb_dat));

		new->ml.time = ktor[i] * p_f->time[i] + (1 - ktor[i]) * p_i->time[i];
		new->ml.r    = ktor[i] * p_f->r[i] + (1 - ktor[i]) * p_i->r[i];
		new->ml.phi  = ktor[i] * p_f->phi[i]  + (1 - ktor[i]) * p_i->z[i];
		new->ml.z    = ktor[i] * p_f->z[i] + (1 - ktor[i]) * p_i->z[i];
	    
		new->poincareId = tor[i];

		new->prev = data->writelist;
		if(data->writelist != NULL) {
		    data->writelist->next = new;
		}
		data->writelist = new;
		data->size = data->size + 1;
	    }
	}
    }
}


void diag_orb_clean(diag_orb_data* data) {

}

void diag_orb_intervalTrigger(diag_orb_data* data, integer* id, real* time, int* store) {
    
    #pragma omp simd
    for(int i= 0; i < NSIMD; i++) {
	store[i] = 0;
	if( (id[i] != -1) && (                                           // Check if dummy
		(id[i] != data->particleId[i]) ||                        // Check if first time this particle
		(time[i] - data->prevWriteTime[i] > data->writeInterval) // Check if enough time has passed from previous write
		)) {
	    
	    store[i] = 1;
	    data->particleId[i] = id[i];
	    data->prevWriteTime[i] = time[i];
	}
    }
}

void diag_orb_poincareTrigger(diag_orb_data* data, int* pol, int* tor, real* kpol, real* ktor, integer* id,
			      real* ftime, real* fr, real* fphi, real* fz,
			      real* itime, real* ir, real* iphi, real* iz){
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
	    real axisRz[2];
	    axisRz[0] = 6.2;
	    axisRz[1] = 0.6; // TODO remove hard coded values
	    real fpol = atan2(fz[i] - axisRz[1], fr[i] - axisRz[0]); // Marker final poloidal angle
	    real ipol = atan2(iz[i] - axisRz[1], ir[i] - axisRz[0]); // Marker initial poloidal angle
	    for(int ip = 0; ip < data->ntoroidalplots; ip++) {
		
		/* The poloidal coordinate is "modulated", i.e., it lies within interval [0,2pi). 
		 * We don't know direction marker is travelling but we assume it has taken the 
		 * shortest route. */

		if( ( (fpol - data->toroidalangles[ip] < 0) && (ipol - data->toroidalangles[ip] > 0) && 
		      (ipol - fpol < CONST_PI ) ) ||
		    ( (fpol - data->toroidalangles[ip] > 0) && (ipol - data->toroidalangles[ip] < 0) &&
		      (fpol - ipol < CONST_PI ) ) ) {
		    tor[i] = ip + DIAG_ORB_MAXPOINCARES;
		    ktor[i] = ( data->toroidalangles[ip] - ipol ) / (fpol - ipol);
		    break;
		}
	    }
	}
    }


}

void diag_orb_lastTrigger(diag_orb_data* data, integer* id, real* time, int* ok) {

}
