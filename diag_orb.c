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


void diag_orb_writegc(FILE* out, particle_simd_gc* p){
    fprintf(out, "%d, %le, %le, %le, %le, %le, %le\n",
	    (int) p->id[0], p->time[0], p->r[0], p->phi[0], p->z[0],
	    p->vpar[0], p->mu[0]);
}

void diag_orb_writegcarray(FILE* out, diag_orb_gcarray* w){
    for(int i = 0; i < w->size; i++){
	fprintf(out, "%d, %le, %le, %le, %le, %le, %le\n",
		w->id[i], w->time[i], w->r[i], w->phi[i], w->z[i],
		w->vpar[i], w->mu[i]);
    }
}

void diag_orb_writegclist(FILE* out, diag_orb_gclist** w){
    diag_orb_gclist* temp;
    while(*w != NULL){
	fprintf(out, "%d, %le, %le, %le, %le, %le, %le\n",
		(*w)->id, (*w)->time, (*w)->r, (*w)->phi, (*w)->z,
		(*w)->vpar, (*w)->mu);
	temp = *w;
	(*w) = (*w)->prev;
	free(temp);
    }
}

void diag_orb_storegcarray(particle_simd_gc* p, diag_orb_gcarray* w, int* write){
    int size = w->size;
    int i;
    for(i=0; i < NSIMD; i=i+1){
	if(write[i]){
	    w->id[size]   = p->id[i];
	    w->time[size] = p->time[i];
	    w->r[size]    = p->r[i];
	    w->phi[size]  = p->phi[i];
	    w->z[size]    = p->z[i];
	    w->vpar[size] = p->vpar[i];
	    w->mu[size]   = p->mu[i];
	    size = size + 1;
	}
    }
    w->size = size;
}


void diag_orb_storegclist(particle_simd_gc* p, diag_orb_gclist** winout, int* write){
    diag_orb_gclist* w;
    int i;
    for(i=0; i < NSIMD; i=i+1){
	if(write[i]){
	    w = malloc(sizeof(diag_orb_gclist));
	    w->id   = p->id[i];
	    w->time = p->time[i];
	    w->r    = p->r[i];
	    w->phi  = p->phi[i];
	    w->z    = p->z[i];
	    w->vpar = p->vpar[i];
	    w->mu   = p->mu[i];
	    w->prev = *winout;
	    *winout = w;
	}
    }
}

diag_orb_gcarray* diag_orb_initgcarray(){
    diag_orb_gcarray* w = malloc(sizeof(diag_orb_gcarray));
    w->size = 0;
    return w;
}

diag_orb_gclist* diag_orb_initgclist(){
    return NULL;
}

void diag_orb_poincaregc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data, int* write){


}

void diag_orb_orbitgc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data, int* write){
    int i;
    #pragma omp simd
    for(i=0; i < NSIMD; i=i+1){
	write[i] = (p_f->id[i] > 0) && (p_f->time[i] != p_i->time[i]) 
	    && ( (fabs(p_f->time[i] - data->lastWriteTime[i]) > data->writeInterval) 
		 || data->particleId[i] != p_f->id[i]);
	if(write[i]) {
	    data->lastWriteTime[i] = p_f->time[i];
	    data->particleId[i] = p_f->id[i];
	}
    }
}

void diag_orb_updategc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data){
   
    int write[NSIMD];
    switch(data->storewhat){
	
    case DIAG_ORB_ORBIT:
	diag_orb_orbitgc(p_f, p_i, data, write);
	break;

    case(DIAG_ORB_POINCARE):
	diag_orb_poincaregc(p_f, p_i, data, write);
	break;
    }


    if(data->islist){
	diag_orb_storegclist(p_f, &data->gclist, write);
    }
    else{
	diag_orb_storegcarray(p_f, data->gcarray, write);
    }
}

void diag_orb_writefo(FILE* out, particle_simd_fo* p){
    fprintf(out, "%d, %le, %le, %le, %le, %le, %le, %le\n",
	    (int) p->id[0], p->time[0], p->r[0], p->phi[0], p->z[0],
	    p->rdot[0], p->phidot[0]*p->r[0], p->zdot[0]);
}

void diag_orb_writefoarray(FILE* out, diag_orb_foarray* w){
    for(int i = 0; i < w->size; i++){
	fprintf(out, "%d, %le, %le, %le, %le, %le, %le, %le\n",
		w->id[i], w->time[i], w->r[i], w->phi[i], w->z[i],
	        w->rdot[i], w->phidot[i] * w->r[i], w->zdot[i]);
    }
}

void diag_orb_writefolist(FILE* out, diag_orb_folist** w){
    diag_orb_folist* temp;
    while(*w != NULL){
	fprintf(out, "%d, %le, %le, %le, %le, %le, %le, %le\n",
		(*w)->id, (*w)->time, (*w)->r, (*w)->phi, (*w)->z,
		(*w)->rdot, (*w)->phidot*(*w)->r, (*w)->zdot);
	temp = *w;
	(*w) = (*w)->prev;
	free(temp);
    }
}

void diag_orb_storefoarray(particle_simd_fo* p, diag_orb_foarray* w, int* write){
    int size = w->size;
    int i;
    for(i=0; i < NSIMD; i=i+1){
	if(write[i]){
	    w->id[size]     = p->id[i];
	    w->time[size]   = p->time[i];
	    w->r[size]      = p->r[i];
	    w->phi[size]    = p->phi[i];
	    w->z[size]      = p->z[i];
	    w->rdot[size]   = p->rdot[i];
	    w->phidot[size] = p->phidot[i];
	    w->zdot[size]   = p->zdot[i];
	    size = size + 1;
	}
    }
    w->size = size;
}


void diag_orb_storefolist(particle_simd_fo* p, diag_orb_folist** winout, int* write){
    diag_orb_folist* w;
    int i;
    for(i=0; i < NSIMD; i=i+1){
	if(write[i]){
	    w = malloc(sizeof(diag_orb_folist));
	    w->id     = p->id[i];
	    w->time   = p->time[i];
	    w->r      = p->r[i];
	    w->phi    = p->phi[i];
	    w->z      = p->z[i];
	    w->rdot   = p->rdot[i];
	    w->phidot = p->phidot[i];
	    w->zdot   = p->zdot[i];
	    w->prev   = *winout;
	    *winout   = w;
	}
    }
}

diag_orb_foarray* diag_orb_initfoarray(){
    diag_orb_foarray* w = malloc(sizeof(diag_orb_foarray));
    w->size = 0;
    return w;
}

diag_orb_folist* diag_orb_initfolist(){
    return NULL;
}

void diag_orb_poincarefo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data, int* write){


}

void diag_orb_orbitfo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data, int* write){
    int i;
    #pragma omp simd
    for(i=0; i < NSIMD; i=i+1){
	write[i] = (p_f->id[i] > 0) && (p_f->time[i] != p_i->time[i]) 
	    && ( (fabs(p_f->time[i] - data->lastWriteTime[i]) > data->writeInterval) 
		 || data->particleId[i] != p_f->id[i]);
	if(write[i]) {
	    data->lastWriteTime[i] = p_f->time[i];
	    data->particleId[i] = p_f->id[i];
	}
    }
}

void diag_orb_updatefo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data){
   
    int write[NSIMD];
    switch(data->storewhat){
	
    case DIAG_ORB_ORBIT:
	diag_orb_orbitfo(p_f, p_i, data, write);
	break;

    case(DIAG_ORB_POINCARE):
	diag_orb_poincarefo(p_f, p_i, data, write);
	break;
    }


    if(data->islist){
	diag_orb_storefolist(p_f, &data->folist, write);
    }
    else{
	diag_orb_storefoarray(p_f, data->foarray, write);
    }
}

void diag_orb_write(diag_orb_data* data){
    if(data->markertype == 0) {
	FILE* out = fopen("foorbits.test","w");
	if(data->islist){
	    diag_orb_writefolist(out, &data->folist);
	}
	else{
	    diag_orb_writefoarray(out, data->foarray);
	}
    }

    if(data->markertype == 1) {
	FILE* out = fopen("gcorbits.test","w");
	if(data->islist){
	    diag_orb_writegclist(out, &data->gclist);
	}
	else{
	    diag_orb_writegcarray(out, data->gcarray);
	}
    }
}

diag_orb_data* diag_orb_init(){
    diag_orb_data* data = malloc(sizeof(diag_orb_data));

    data->storewhat = DIAG_ORB_ORBIT;
    data->islist = 0;
    data->markertype = 0;
    data->writeInterval = 1.0e-8;
    data->gclist = NULL;
    data->gcarray = NULL;
    data->folist = NULL;
    data->foarray = NULL;

    if(data->markertype == 0) {
	if(data->islist) {
	    data->folist = diag_orb_initfolist();
	}
	else{
	    data->foarray = diag_orb_initfoarray();
	}
    }
    else {
	if(data->islist) {
	    data->gclist = diag_orb_initgclist();
	}
	else{
	    data->gcarray = diag_orb_initgcarray();
	}
    }

    for(int i=0; i < NSIMD; i++) {
	data->particleId[i] = -1;
    }

    return data;
}

void diag_orb_clean(diag_orb_data* data){
    if(data->gclist != NULL){
	diag_orb_gclist* temp;
	while(data->gclist != NULL){
	    temp = data->gclist;
	    data->gclist = data->gclist->prev;
	    free(temp);
	}
    }
    
    if(data->gcarray != NULL){
	free(data->gcarray);
    }

    if(data->folist != NULL){
	diag_orb_folist* temp;
	while(data->folist != NULL){
	    temp = data->folist;
	    data->folist = data->folist->prev;
	    free(temp);
	}
    }
    
    if(data->foarray != NULL){
	free(data->foarray);
    }

    free(data);

}
